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

__all__ = """
    Common
    App
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
