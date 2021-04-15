

from osgeo import gdal, gdalnumeric, gdal_array
import numpy as np
from numpy import copy
import os
from os.path import splitext
import processing
import re
import datetime
import sqlite3
import csv
import tempfile

from qgis.core import (QgsUnitTypes,  
                       QgsFields,
                       QgsField,
                       QgsFeature,
                       QgsFeatureSink,
                       QgsColorRampShader,
                       QgsRasterBandStats,
                       QgsRasterFileWriter)

from qgis.core import *

from PyQt5.QtCore import QCoreApplication
from qgis.PyQt.QtGui import QColor


__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
__esv_data_location__ = os.path.join(__location__, "esv_data")

class LULC_dataset:
    """Custom class to handle summarizing LULC rasters with ecosystem service values"""

    def __init__(self, source_name, raster):
        self.source_name = source_name
        self.raster = raster
        self.raster_file_path = self.get_input_raster_path()
        self.raster_summary = self.summarize_raster_values()


    def get_input_raster_path(self):
        """Returns the file path of the input raster removing 'file://' from the front"""

        raster_source_path = self.raster.source()
        file_prepend = 'file://'
        if raster_source_path.find(file_prepend)==0:
            # Strip out 'file://' from begining of file path
            #  to pass to gdal
            raster_source_path = raster_source_path[len(file_prepend):]
        return(raster_source_path)

    def is_valid(self):
        """Check that the raster is valid for the selected source type
           returns True if the raster is valid if raster is not valid returns a string with the reason it is invalid"""

        with ESV_dataset() as esv:
            expected_values = esv.get_pixel_options(self.source_name)

        # Find any unexpected raster values (remove zero values)
        raster_vals = [r[0] for r in self.raster_summary if r[0] != 0]
        difference = set(raster_vals).difference(expected_values)
        if difference:
            error_message = f'The following raster values are invalid for the selected LULC source ({self.source_name}): {difference}. Please check that you have selected the correct data source.'
            return(error_message)
        else:
            return(True)

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


        with tempfile.TemporaryDirectory() as tmpdirname:
            csv_name = 'TEMP_SUMMARY.csv'
            summary_file = os.path.join(tmpdirname,csv_name)

            summary_params = {'BAND':1,
                              'INPUT':self.raster_file_path,
                              'OUTPUT_TABLE':summary_file }

            summarize = processing.run("native:rasterlayeruniquevaluesreport", summary_params)

            val = lambda x: int(float(x))

            with open(summary_file, 'r', encoding = 'LATIN') as csvfile:
                reader = csv.DictReader(csvfile)
                data = []
                for row in reader:
                    data.append([val(row['value']),
                                val(row['count']),
                                val(row['count'])*self.cell_size()])

        return(data)



class ESV_dataset:
    """Custom class which mostly handles rocessing of SQLite file to produce ESV_data file  """
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
        self.execute("""DROP TABLE IF EXISTS raster_area_summary;""")
        self.execute("""CREATE TEMPORARY TABLE raster_area_summary(
            lulc_value INTEGER NOT NULL,
            pixel_count INTEGER NOT NULL,
            area REAL NOT NULL
        ); """)
        self.executemany("""INSERT INTO raster_area_summary(lulc_value,pixel_count,area) VALUES (?,?,?);""", data)

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
                           ROUND(area, 2) AS land_cover_area,
                           esv_estimates.service_name AS Ecosystem_Service_Name,
                           ROUND(estimate_min * area) AS Minimum_Value_Estimate,
                           ROUND(estimate_max * area) AS Maximum_Value_Estimate,
                           ROUND(estimate_avg * area) AS Average_Value_Estimate
                    FROM esv_estimates JOIN raster_area_summary ON esv_estimates.lulc_value = raster_area_summary.lulc_value
                          JOIN (SELECT * FROM lulc_legend WHERE source = ?) AS lulc_legend ON esv_estimates.lulc_value = lulc_legend.value
                    WHERE esv_estimates.lulc_source = ?"""

        result = self.query(query, (source_name,source_name))

        column_names = [d[0] for d in self.cursor.description]
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
        FROM esv_estimates
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
        result = self.query("""SELECT DISTINCT service_name FROM esv_estimates""")
        services = [r[0] for r in result]
        return(services)
    



class Symbology:
    """ Custom class to handle symbology of input dataset 
    Returns color ramp for symbology of each respective input_esv_field"""

    def __init__(self, raster_stats, input_esv_field):
        self.raster_stats = raster_stats
        self.input_esv_field = input_esv_field

    def symbolize_input_raster(self):
        """Calculates statistics for input raster and returns color ramp corresponding to each Ecosystem service choice"""
        
        min_val = self.raster_stats.minimumValue            #minimum pixel value in layer
        max_val = self.raster_stats.maximumValue            #maximum pixel value in layer

        #first need to get range of values in raster layer
        value_range = list(range(int(min_val), int(max_val+1)))           #Range of values in raster layer. Without +1 doesn't capture highest value
        value_range.sort()
        for value in value_range:                   #deletes 0 value from value range so as not to skew shading in results
            if value < min_val:
                del value

        #we will categorize pixel values into 5 quintiles, based on value_range of raster layer
        #defining min and max values for each quintile. 
        #Also, values are rounded to 2 decimal places
        first_quintile_max = round(np.percentile(value_range, 20), 2)
        first_quintile_min = round(min_val, 2)
        second_quintile_max = round(np.percentile(value_range, 40), 2)
        second_quintile_min = round((first_quintile_max + .01), 2)
        third_quintile_max = round(np.percentile(value_range, 60), 2)
        third_quintile_min = round((second_quintile_max + .01), 2)
        fourth_quintile_max = round(np.percentile(value_range, 80), 2)
        fourth_quintile_min = round((third_quintile_max + .01), 2)
        fifth_quintile_max = round(np.percentile(value_range, 100), 2)
        fifth_quintile_min = round((fourth_quintile_max + .01), 2)

        #green color ramp
        if self.input_esv_field == 'Aesthetic':
            colors_list = [ QgsColorRampShader.ColorRampItem(0, QColor(255, 255, 255, .5), 'No Value'), \
                       QgsColorRampShader.ColorRampItem(first_quintile_max, QColor(204, 255, 204), f"${first_quintile_min}0 - ${first_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(second_quintile_max, QColor(153, 255, 153), f"${second_quintile_min} - ${second_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(third_quintile_max, QColor(51, 255, 51), f"${third_quintile_min} - ${third_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(fourth_quintile_max, QColor(0, 204, 0), f"${fourth_quintile_min} - ${fourth_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(fifth_quintile_max+1, QColor(0, 102, 0), f"${fifth_quintile_min} - ${fifth_quintile_max}0")]

        #green color ramp
        elif self.input_esv_field == 'Biodiversity':
            colors_list = [ QgsColorRampShader.ColorRampItem(0, QColor(255, 255, 255), 'No Value'), \
                       QgsColorRampShader.ColorRampItem(first_quintile_max, QColor(204, 255, 204), f"${first_quintile_min}0 - ${first_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(second_quintile_max, QColor(153, 255, 153), f"${second_quintile_min} - ${second_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(third_quintile_max, QColor(51, 255, 51), f"${third_quintile_min} - ${third_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(fourth_quintile_max, QColor(0, 204, 0), f"${fourth_quintile_min} - ${fourth_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(fifth_quintile_max+1, QColor(0, 102, 0), f"${fifth_quintile_min} - ${fifth_quintile_max}0")]
        
        #orange color ramp
        elif self.input_esv_field == 'Climate Regulation':
            colors_list = [ QgsColorRampShader.ColorRampItem(0, QColor(255, 255, 255, .5), 'No Value'), \
                       QgsColorRampShader.ColorRampItem(first_quintile_max, QColor(255, 229, 204), f"${first_quintile_min}0 - ${first_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(second_quintile_max, QColor(255, 204, 153), f"${second_quintile_min} - ${second_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(third_quintile_max, QColor(255, 153, 51), f"${third_quintile_min} - ${third_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(fourth_quintile_max, QColor(204, 102, 0), f"${fourth_quintile_min} - ${fourth_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(fifth_quintile_max+1, QColor(102, 51, 0), f"${fifth_quintile_min} - ${fifth_quintile_max}0")]
         
        #brown color ramp
        elif self.input_esv_field == 'Erosion Control':
            colors_list = [ QgsColorRampShader.ColorRampItem(0, QColor(255, 255, 255, .5), 'No Value'), \
                       QgsColorRampShader.ColorRampItem(first_quintile_max, QColor(220,187,148), f"${first_quintile_min}0 - ${first_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(second_quintile_max, QColor(198,168,134), f"${second_quintile_min} - ${second_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(third_quintile_max, QColor(169,144,115), f"${third_quintile_min} - ${third_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(fourth_quintile_max, QColor(138,117,93), f"${fourth_quintile_min} - ${fourth_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(fifth_quintile_max+1, QColor(100,85,67), f"${fifth_quintile_min} - ${fifth_quintile_max}0")]
        
        #pink color ramp
        elif self.input_esv_field == 'Food/Nutrition':
            colors_list = [ QgsColorRampShader.ColorRampItem(0, QColor(255, 255, 255, .5), 'No Value'), \
                       QgsColorRampShader.ColorRampItem(first_quintile_max, QColor(255,204,229), f"${first_quintile_min}0 - ${first_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(second_quintile_max, QColor(255,153,204), f"${second_quintile_min} - ${second_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(third_quintile_max, QColor(255,51,153), f"${third_quintile_min} - ${third_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(fourth_quintile_max, QColor(204,0,102), f"${fourth_quintile_min} - ${fourth_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(fifth_quintile_max+1, QColor(102,0,51), f"${fifth_quintile_min} - ${fifth_quintile_max}0")]        

        #yellow color ramp
        elif self.input_esv_field == 'Pollination':
            colors_list = [ QgsColorRampShader.ColorRampItem(0, QColor(255, 255, 255, .5), 'No Value'), \
                       QgsColorRampShader.ColorRampItem(first_quintile_max, QColor(255,245,145), f"${first_quintile_min}0 - ${first_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(second_quintile_max, QColor(255,245,6), f"${second_quintile_min} - ${second_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(third_quintile_max, QColor(205,195,4), f"${third_quintile_min} - ${third_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(fourth_quintile_max, QColor(168,160,3), f"${fourth_quintile_min} - ${fourth_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(fifth_quintile_max+1, QColor(92,87,2), f"${fifth_quintile_min} - ${fifth_quintile_max}0")]        

        #gray/black color ramp
        elif self.input_esv_field in ('Extreme Events', 'Protection from Extreme Events'):
            colors_list = [ QgsColorRampShader.ColorRampItem(0, QColor(255, 255, 255, .5), 'No Value'), \
                       QgsColorRampShader.ColorRampItem(first_quintile_max, QColor(224,224,224), f"${first_quintile_min}0 - ${first_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(second_quintile_max, QColor(192,192,192), f"${second_quintile_min} - ${second_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(third_quintile_max, QColor(128,128,128), f"${third_quintile_min} - ${third_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(fourth_quintile_max, QColor(64,64,64), f"${fourth_quintile_min} - ${fourth_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(fifth_quintile_max+1, QColor(0,0,0), f"${fifth_quintile_min} - ${fifth_quintile_max}0")]          

        #purple color ramp
        elif self.input_esv_field == 'Raw materials':
            colors_list = [ QgsColorRampShader.ColorRampItem(0, QColor(255, 255, 255, .5), 'No Value'), \
                       QgsColorRampShader.ColorRampItem(first_quintile_max, QColor(229,204,255), f"${first_quintile_min}0 - ${first_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(second_quintile_max, QColor(204,153,255), f"${second_quintile_min} - ${second_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(third_quintile_max, QColor(153,51,255), f"${third_quintile_min} - ${third_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(fourth_quintile_max, QColor(102,0,204), f"${fourth_quintile_min} - ${fourth_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(fifth_quintile_max+1, QColor(51,0,102), f"${fifth_quintile_min} - ${fifth_quintile_max}0")]
        
        #red color ramp
        elif self.input_esv_field == 'Recreation':
            colors_list = [ QgsColorRampShader.ColorRampItem(0, QColor(255, 255, 255, .5), 'No Value'), \
                       QgsColorRampShader.ColorRampItem(first_quintile_max, QColor(253,204,184), f"${first_quintile_min}0 - ${first_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(second_quintile_max, QColor(252,143,111), f"${second_quintile_min} - ${second_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(third_quintile_max, QColor(244,77,55), f"${third_quintile_min} - ${third_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(fourth_quintile_max, QColor(197,22,27), f"${fourth_quintile_min} - ${fourth_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(fifth_quintile_max+1, QColor(103,0,13), f"${fifth_quintile_min} - ${fifth_quintile_max}0")]
        
        #brown color ramp
        elif self.input_esv_field == 'Soil Formation':
            colors_list = [ QgsColorRampShader.ColorRampItem(0, QColor(255, 255, 255, .5), 'No Value'), \
                       QgsColorRampShader.ColorRampItem(first_quintile_max, QColor(220,187,148), f"${first_quintile_min}0 - ${first_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(second_quintile_max, QColor(198,168,134), f"${second_quintile_min} - ${second_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(third_quintile_max, QColor(169,144,115), f"${third_quintile_min} - ${third_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(fourth_quintile_max, QColor(138,117,93), f"${fourth_quintile_min} - ${fourth_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(fifth_quintile_max+1, QColor(100,85,67), f"${fifth_quintile_min} - ${fifth_quintile_max}0")]       
        
        #blue/purple color ramp
        elif self.input_esv_field == 'Waste Assimilation':
            colors_list = [ QgsColorRampShader.ColorRampItem(0, QColor(255, 255, 255, .5), 'No Value'), \
                       QgsColorRampShader.ColorRampItem(first_quintile_max, QColor(204,204,255), f"${first_quintile_min}0 - ${first_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(second_quintile_max, QColor(153,153,255), f"${second_quintile_min} - ${second_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(third_quintile_max, QColor(51,51,255), f"${third_quintile_min} - ${third_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(fourth_quintile_max, QColor(0,0,204), f"${fourth_quintile_min} - ${fourth_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(fifth_quintile_max+1, QColor(0,0,102), f"${fifth_quintile_min} - ${fifth_quintile_max}0")]
        
        #medium blue color ramp
        elif self.input_esv_field in ('Water Supply', 'Air Quality','Air quality'):
            colors_list = [ QgsColorRampShader.ColorRampItem(0, QColor(255, 255, 255, .5), 'No Value'), \
                       QgsColorRampShader.ColorRampItem(first_quintile_max, QColor(204,229,255), f"${first_quintile_min}0 - ${first_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(second_quintile_max, QColor(153,204,255), f"${second_quintile_min} - ${second_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(third_quintile_max, QColor(51,153,205), f"${third_quintile_min} - ${third_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(fourth_quintile_max, QColor(0,102,204), f"${fourth_quintile_min} - ${fourth_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(fifth_quintile_max+1, QColor(0,51,102), f"${fifth_quintile_min} - ${fifth_quintile_max}0")]       
        
        else:
            # Catchall if for some reason the ESV type does match any of the other values
            colors_list = [ QgsColorRampShader.ColorRampItem(0, QColor(255, 255, 255, .5), 'No Value'), \
                       QgsColorRampShader.ColorRampItem(first_quintile_max, QColor(204,229,255), f"${first_quintile_min}0 - ${first_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(second_quintile_max, QColor(153,204,255), f"${second_quintile_min} - ${second_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(third_quintile_max, QColor(51,153,205), f"${third_quintile_min} - ${third_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(fourth_quintile_max, QColor(0,102,204), f"${fourth_quintile_min} - ${fourth_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(fifth_quintile_max+1, QColor(0,51,102), f"${fifth_quintile_min} - ${fifth_quintile_max}0")]       


        return colors_list








