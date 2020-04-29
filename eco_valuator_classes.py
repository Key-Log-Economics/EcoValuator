

from osgeo import gdal, gdalnumeric, gdal_array
import numpy as np
from numpy import copy
import os
from os.path import splitext
import processing
import re
import datetime
import sqlite3

from qgis.core import (QgsUnitTypes,  
                       QgsFields,
                       QgsField,
                       QgsFeature,
                       QgsFeatureSink,
                       QgsColorRampShader,
                       QgsRasterBandStats,
                       QgsRasterFileWriter)

from PyQt5.QtCore import QCoreApplication

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
__esv_data_location__ = os.path.join(__location__, "esv_data")

class LULC_dataset:
    """Custom class to handle summarizing LULC rasters with ecosystem service values"""

    def __init__(self, source_name, raster):
        self.source_name = source_name
        self.raster = raster
        self.ras_array = self.get_raster_as_array()

    def is_valid(self):
        """Check that the raster is valid for the selected source type
           returns True if the raster is valid if raster is not valid returns a string with the reason it is invalid"""

        with ESV_dataset() as esv:
            expected_values = esv.get_pixel_options(self.source_name)

        # Find any unexpected raster values
        raster_vals = np.unique(self.ras_array)
        difference = set(raster_vals).difference(expected_values)
        if len(difference) > 0:
            error_message = f'The following raster values are invalid for the selected LULC source ({self.source_name}): {difference}. Please check that you have selected the correct data source.'
            return(error_message)
        else:
            return(True)

    def get_raster_as_array(self):
        """Returns the raster values as an array 
        by creating a GDAL Dataset object from the input LULC raster"""

        raster_source_path = self.raster.source()
        file_prepend = 'file://'
        if raster_source_path.find(file_prepend)==0:
            # Strip out 'file://' from begining of file path
            #  to pass to gdal
            raster_source_path = raster_source_path[len(file_prepend):]

        LULC_ds = gdal.Open(raster_source_path)
        LULC_ds_band = LULC_ds.GetRasterBand(1)

        nodata_val = LULC_ds_band.GetNoDataValue()
        raster_as_array = LULC_ds_band.ReadAsArray()

        if nodata_val is not None:
            # If the raster has a no_data value, filter that value
            #  from the array
            raster_as_array = raster_as_array[raster_as_array != nodata_val]
        
        # Remove 0 values
        raster_as_array = raster_as_array[raster_as_array != 0]
        return(raster_as_array)

    def cell_size(self):
        """Returns the size of the raster cell in hectares"""
        # Get the cell dimensions and units of the LULC raster
        cell_dim_x = self.raster.rasterUnitsPerPixelX()
        cell_dim_y = self.raster.rasterUnitsPerPixelY()
        raster_units = self.raster.crs().mapUnits()

        # Figure out the conversion factor from map units to hectares
        raster_units_square = QgsUnitTypes.distanceToAreaUnit(raster_units)
        hectare_unit_code = QgsUnitTypes.stringToAreaUnit("Hectares")[0]
        conversion_factor = QgsUnitTypes.fromUnitToUnitFactor(raster_units_square, hectare_unit_code)

        # Return raster cell size in hectares
        # Look out for floating point issues!
        return(cell_dim_x * cell_dim_y * conversion_factor)

    def summarize_raster_values(self):
        """Returns an array summarizing the raster values. Pixel counts are converted to area (hectares)"""

        #Check for NaN in the array        
        if not isinstance(self.ras_array.dtype.type, np.float):
            #If dtype is object, convert to float and filter Nan
            test = self.ras_array.astype(float)
            test = test[~np.isnan(test)]
        else:
            test = self.ras_array[~np.isnan(self.ras_array)]
        
        raster_values, value_pixel_count = np.unique(test, return_counts=True)
        value_area = value_pixel_count * self.cell_size()

        col_defs = np.dtype({'names':('raster_values',
                                    'value_pixel_count',
                                    'value_area'),
                            'formats':(raster_values.dtype,
                                    value_pixel_count.dtype,
                                    value_area.dtype)
                            }
                           )
        data_output = list(zip(raster_values,
                           value_pixel_count,
                           value_area)
                       )
        return(data_output)


    def convert_LULC_pixels_to_ESV_valuation(self,
                                            service,
                                            aggregation_method,
                                            output_raster_path):

        """Converts the raster of LULC land use class values to a per-pixel ecosystem service
        valueation based on the input value type and aggragation method using the ESV values from the 
        ESV csv file.
        @arg service: Name of the ecosystem service of intereste (from the tool arguments)
        @arg aggregation_method: Which ESV aggregation will be used (Min/Max/Avg from tool arguments)
        @arg output_raster_path: path to the output ESV summary raster (from the tool arguments)

        """

        """1. Take ESV csv file and transform it into a dictionary with LULC code as the key and the per-pixel
        ecosystem service valuation as a value. The dictionary reflects values for the ecosystem service and 
        aggregation method specified. The dictionary is used to map LULC pixel values to the valuation in the 
        output raster"""

        ESV_array = self.ESV_to_np_array()
        col_names = ESV_array.dtype.names

        value_field = self.get_field_name(aggregation_method, col_names)
        if not service in ESV_array['Service']:
            raise Exception(f'Selected service {service} is not valid for this LULC dataset type!')
        filtered_array = ESV_array[ESV_array['Service']==service][['LULC Value',value_field]]
        cell_size = self.cell_size()
        ESV_dict = {int(code):self.make_num(value) * cell_size for code, value in filtered_array}

        """2. Use the ESV_dict dictionary to map the LULC values numpy array to a per-pixel valuation
        numpy array"""


        def map_lulc_val(cell_val):
            """Helper function. Caluclates an LULC value for each individual element in a Numpy array
            @arg cell_val: Value of the raster cell (array element)
            @arg ESV_dict: per-pixel ecosystem service value (value) for each LULC type (key)"""
            ESV_value = ESV_dict.get(cell_val)
            if ESV_value:
                return(ESV_value)
            else:
                return(None)
        # Vecorize the mapping function to use on the numpy array
        vfunc = np.vectorize(map_lulc_val)

        ESV_valued_pixels = vfunc(self.ras_array)
        import sqlite3
        ESV_valued_pixels = np.vstack(ESV_valued_pixels[:, :]).astype(np.float)

        driver = gdal.GetDriverByName("GTiff") # will need to ensure raster file passed is a .tiff
        #dataType = gdal_array.NumericTypeCodeToGDALTypeCode(ESV_valued_pixels.dtype.type) #Figure out data type of output raster from array
        dataType = gdal.GDT_Float64
        dsOut = driver.Create(output_raster_path, self.LULC_ds.RasterXSize,  self.LULC_ds.RasterYSize, 1, dataType)
        gdalnumeric.CopyDatasetInfo(self.LULC_ds, dsOut)
        bandOut=dsOut.GetRasterBand(1)
        gdalnumeric.BandWriteArray(bandOut, ESV_valued_pixels)

        # Raster writes when the band and dataset variables are removed
        bandOut = None
        dsOut = None

        # Get range of output raster values for map styling
        self.output_min_val = int(round(np.nanmin(ESV_valued_pixels)))
        self.output_max_val = int(round(np.nanmax(ESV_valued_pixels)))

class ESV_dataset:
    def __init__(self):

        sqlite_file = os.path.join(__esv_data_location__, 'ESV_data.sqlite')
        if os.path.isfile(sqlite_file):
            self._conn = sqlite3.connect(sqlite_file)
            self._cursor = self._conn.cursor()
        else:
            raise OSError(sqlite_file)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.commit()
        self.connection.close()

    @property
    def connection(self):
        return self._conn

    @property
    def cursor(self):
        return self._cursor

    def commit(self):
        self.connection.commit()

    def execute(self, sql, params=None):
        self.cursor.execute(sql, params or ())

    def executemany(self, sql, params=None):
        self.cursor.executemany(sql, params or ())

    def fetchall(self):
        return self.cursor.fetchall()

    def query(self, sql, params=None):
        self.cursor.execute(sql, params or ())
        return self.fetchall()

    def test_query(self):
        result = self.query("""SELECT name FROM sqlite_master WHERE type='table';""")
        return(result)

    def create_temp_table_raster_area_summary(self, data):
        """Takes an input of data that summarizes LULC raster values in the AOI and converts into 
        a temporary table in the sqlite database. This allows the values to be summarized and used
        for later queries.
        
        Arguments:
            data {list of tuples} -- Data to be passed in to the temporary table:
                   (LULC Value, Pixel Count, Area) where
                      LULC Value is the integer LULC type code
                      Pixel Count is the integer count of pixels in the LULC type
                      Area is the total area of that type in the AOI (hectares)
        """
        # Make sure incoming data is the correct type (i.e. not Numpy data types)
        cleaned_data = [(d[0].item(),d[1].item(),d[2].item()) for d in data]

        self.execute("""DROP TABLE IF EXISTS raster_area_summary;""")
        self.execute("""CREATE TEMPORARY TABLE raster_area_summary(
            lulc_value INTEGER NOT NULL,
            pixel_count INTEGER NOT NULL,
            area REAL NOT NULL
        ); """)
        self.executemany("""INSERT INTO raster_area_summary(lulc_value,pixel_count,area) VALUES (?,?,?);""", cleaned_data)

    def get_LULC_evaluation_data(self, data, source_name):
        """ Final output for step 1. Exports data and column names used to build the
        output table for step 1.
        Arguments:
            data {list} -- output from LULC_dataset.summarize_raster_values() method
            source_name {str} -- name of the LULC source (e.g. NLCD, NALCMS...)

        Returns:
            dictionary of column names and data to populate the output table
        """
        self.create_temp_table_raster_area_summary(data)
        #self.create_area_summmay_pivot_view(source_name)

        #result = self.query("""SELECT * FROM esv_val_pivot""")

        query = """ SELECT esv_estimates.lulc_value AS LULC_code,
                           lulc_legend.description LULC_Description,
                           service_names.service_name AS Ecosystem_Service_Name,
                           ROUND(estimate_min * area) AS Minimum_Value_Estimate,
                           ROUND(estimate_max * area) AS Maximum_Value_Estimate,
                           ROUND(estimate_avg * area) AS Average_Value_Estimate
                    FROM esv_estimates JOIN raster_area_summary ON esv_estimates.lulc_value = raster_area_summary.lulc_value
                          JOIN service_names ON esv_estimates.service_id = service_names.service_id
                          JOIN (SELECT * FROM lulc_legend WHERE source = ?) AS lulc_legend ON esv_estimates.lulc_value = lulc_legend.value
                    WHERE esv_estimates.lulc_source = ?"""

        result = self.query(query, (source_name,source_name))

        for r in result:
            print(r)
        column_names = [d[0] for d in self.cursor.description]
        print(column_names)
        return({'data':result,
                'column_names':column_names})

    def make_reclassify_table(self,
                              cell_size,
                              source_name,
                              summary_type,
                              service):
        """Creates a table to be passed to the QGIS reclassify by table algorithm
        for step 2. The table is a list with the repeated pattern: 
            target LULC value-1, target LULC_value, ecosystem service valuation * cell_size.
            see documentation of the reclassify by table algorithm and the main script
            for step 2.
        
        Arguments:
            cell_size {float}  -- cell size in hectares of the raster to be reclassifed. Obtained from the LULC_dataset.cellsize() method
            source_name {str}  -- name of the LULC source (e.g. NLCD, NALCMS...)
            summary_type {str} -- aggregation of ecosystem service values to use (min, max, avg). text appended to "estimate_" to select a column
            service {str}      -- name of the ecosystem service to be used
        """

        query = f"""SELECT lulc_value, estimate_{summary_type}
        FROM esv_estimates JOIN service_names ON esv_estimates.service_id = service_names.service_id
        WHERE lulc_source = (?) AND service_name = (?)"""

        result = self.query(query, (source_name, service))

        LULC_max = [r[0] for r in result]
        LULC_min = [r-1 for r in LULC_max]
        ESV_val = [r[1]*cell_size for r in result]
        reclass_table = [row for sublist in zip(LULC_min, LULC_max, ESV_val) for row in sublist]

        return(reclass_table)

    def get_lulc_sources(self):
        """Returns a list with the LULC sources (e.g. NLCD, NLCMS) available in the dataset"""
        result = self.query("""SELECT DISTINCT source FROM lulc_legend""")
        sources = [r[0] for r in result]
        return(sources)
    
    def get_pixel_options(self, source):
        """Given a LULC source (e.g. NLCD, NLCMS), returns a list of the potential
        pixel values for that source. This is used to validate input rasters 
        against the selected source."""

        result = self.query("""SELECT value FROM lulc_legend WHERE source = (?)""",(source,))
        pixel_values = [r[0] for r in result]
        return(pixel_values)

    def get_ecosystem_service_names(self):
        """Returns a list of ecosystem services available in the ESV data"""
        result = self.query("""SELECT service_name FROM service_names""")
        services = [r[0] for r in result]
        return(services)
    
