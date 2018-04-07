"""
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
"""

__author__    = 'Roman Geisthövel'
__date__      = '2018-03-05'
__copyright__ = '(C) 2018 by Roman Geisthövel'

import os
import ctypes as ct
from ctypes.util import find_library

__all__ = """
    App
    Common
    FFI
    Raster
""".split()

def running_qgis():
    try:
        from qgis.core import Qgis
        return Qgis.QGIS_VERSION_INT
    except ImportError:
        return 0

class Common:
    "Functionality common to all platforms"

    @staticmethod
    def folder():
        return os.path.dirname(os.path.realpath(__file__))

    @staticmethod
    def mkpath(*args):
        return os.path.join(*args)

    @staticmethod
    def file_size(path, unit="byte"):
        if os.path.isfile(path):
            x = os.path.getsize(path)
            unit = unit.lower()
            if unit in ("mb", "mega", "megabyte"):
                return x >> 20
            elif unit in ("gb", "giga", "gigabyte"):
                return x >> 30
            else:
                return x
        else:
            return 0


class FFI:
    "C function interface utilities"

    _type_map = {}

    @staticmethod
    def load(name):
        try:
            if not os.path.isfile(name):
                name = find_library(name)
            if name:
                dll = ct.windll if os.name == "nt" else ct.cdll
                return dll.LoadLibrary(name)
        except OSError:
            pass
        return None

    @staticmethod
    def c_t(type_name):
        if not FFI._type_map:
            m = FFI._type_map

            # Integer types
            for aliases in ("char", "int i", "short", "long l", "longlong ll",
                            "bool ?", "float f", "double d", "longdouble ld"):
                aliases = aliases.split()
                T = aliases[0]
                for _ in aliases:
                    m[_]          = getattr(ct, "c_%s" % T)
                    if T == "char":
                        # ctypes has no c_uchar
                        m["uchar"] = ct.c_ubyte
                    elif T in "char int short long longlong":
                        # unsigned T
                        m["u%s" % _]  = getattr(ct, "c_u%s" % T)

                    # Pointer to T
                    m["%s*" % _]  = ct.POINTER(m[T])

                    if T not in "bool float double longdouble":
                        # Pointer to unsigned T
                        m["u%s*" % _] = ct.POINTER(m["u%s" % T])

            # Exact size integer (pointer) types
            for n in (1,2,4,8):
                s = "i%i" % n
                u = "u%i" % n
                m[s] = getattr(ct, "c_int%i" % (n*8))
                m[u] = getattr(ct, "c_uint%i" % (n*8))
                m["%s*" % s] = ct.POINTER(m[s])
                m["%s*" % u] = ct.POINTER(m[u])

            # Void pointer
            m["void*"] = m["*"] = ct.c_void_p

        return FFI._type_map[type_name]

    @staticmethod
    def array_t(T, size):
        return T * size

    @staticmethod
    def func_t(res_type, *arg_types):
        return ct.CFUNCTYPE(res_type, *arg_types)

    @staticmethod
    def sizeof(obj):
        return ct.sizeof(obj)


if running_qgis():

    from qgis.core import QgsMessageLog, QgsRasterDataProvider, QgsRasterLayer
    from qgis.core import Qgis
    from osgeo import gdal


    #-------------------------------------------------------------------------
    class App:

        @staticmethod
        def info(*msg, **opt):
            App.log(*msg, level=Qgis.Info, **opt)

        @staticmethod
        def warn(*msg, **opt):
            App.log(*msg, level=Qgis.Warning, **opt)

        @staticmethod
        def critical(*msg, **opt):
            App.log(*msg, level=Qgis.Critical, **opt)

        @staticmethod
        def log(*msg, **opt):
            src = opt.get("src", "")
            lvl = opt.get("level", Qgis.Info)
            sep = opt.get("sep", " ")
            QgsMessageLog.logMessage(sep.join(str(_) for _ in msg), tag=src, level=lvl)

    #-------------------------------------------------------------------------
    class Raster:

        gdal2numpy_type = dict(Byte="u1", Int16="i2", UInt16="u2", Int32="i4",
                                UInt32="u4", Float32="f4", Float64="f8")

        @staticmethod
        def geo_transform(x):
            dy,dx = Raster.cellsize(x)
            p = x.extent()
            return (p.xMinimum(),dx,0,p.yMinimum(),0,-dy)

        @staticmethod
        def crs(x):
            return x.crs().toWkt()

        # Returns tuple (dy, dx)
        @staticmethod
        def cellsize(x):
            return (x.rasterUnitsPerPixelY(), x.rasterUnitsPerPixelX())

        @staticmethod
        def num_bands(x):
            return x.bandCount()

        @staticmethod
        def from_numpy(x, **opt):
            ras = QgsRasterLayer()
            dp = QgsRasterDataProvider()
            ras.setDataProvider(dp)

        @staticmethod
        def numpy_to_file(x, file_name, **opt):
            h,w = x.shape
            src = opt.get("src", None)
            if src is not None:
                srcds = gdal.Open(src, gdal.GA_ReadOnly)
                crs   = srcds.GetProjectionRef()
                gt    = srcds.GetGeoTransform()
            else:
                crs = opt.get("crs", None)
                gt  = opt.get("geo_transform", (0,1,0,0,0,1))

            drv = gdal.GetDriverByName("GTiff")
            ds  = drv.Create(file_name, w, h, 1, gdal.GDT_Float32,
                    """COMPRESS=DEFLATE
                        ZLEVEL=4
                        BIGTIFF=IF_SAFER
                        PREDICTOR=3
                        NUM_THREADS=ALL_CPUS""".split())
            if crs:
                ds.SetProjection(crs)
            ds.SetGeoTransform(gt)
            ds.GetRasterBand(1).WriteArray(x)
            ds.FlushCache()

        @staticmethod
        def to_numpy(x, **opt):
            """
            Options     Description             Default
            band        Band number             1
            dtype       Data type of result     None, i.e. dtype of raster layer
            """
            ds    = gdal.Open(str(x.source()), gdal.GA_ReadOnly)
            band  = ds.GetRasterBand(opt.get("band", 1))
            out   = band.ReadAsArray()
            dtype = opt.get("dtype", out.dtype)
            if dtype != out.dtype:
                out = out.astype(dtype, copy=False)
            return out

        @staticmethod
        def shape(x):
            return (x.height(), x.width())
