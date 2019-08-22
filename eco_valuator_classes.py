import csv
from osgeo import gdal
import numpy as np
from numpy import copy
import os
from os.path import splitext
import processing


from collections import defaultdict
from qgis.core import (QgsUnitTypes,  
                       QgsFields,
                       QgsField,
                       QgsFeature,
                       QgsFeatureSink,
                       QgsColorRampShader,
                       QgsRasterShader,
                       QgsSingleBandPseudoColorRenderer,
                       QgsRasterBandStats,
                       QgsRasterFileWriter)

from PyQt5.QtCore import QVariant, QCoreApplication
from PyQt5.QtGui import *

from .appinter import (Raster, App)



class LULC_dataset:
    """Custom class to handle summarizing LULC rasters with ecosystem service values"""

    def __init__(self, source_name, raster, esv_data_directory):
        self.source_name = source_name
        self.raster = raster
        self.esv_data_directory = esv_data_directory

        print('converting to array')
        self.ras_array = self.get_raster_as_array()
        print(np.unique(self.ras_array))

        print('summarizing raster')
        ras_summary = self.summarize_raster_values()
        self.ras_summary_array = ras_summary['array']
        self.ras_summary_dict = ras_summary['dict']

    def is_valid(self):
        """Check that the raster is valid for the selected source type
             returns True if the raster is valid
             if raster is not valid returns a string with the reason it is invalid"""

        code_library = {
            'NLCD':[11, 21, 22, 23, 24, 31, 41, 42, 43, 52, 71, 81, 82, 90, 95],
            'NALCMS':[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
                       }
        
        # Find any unexpected raster values
        expected_values = code_library[self.source_name]
        raster_vals = self.ras_summary_dict.keys()
        difference = set(raster_vals).difference(expected_values)
        if len(difference) > 0:
            error_message = f'The following raster values are invalid for the selected LULC source ({self.source_name}): {difference}. Please check that you have selected the correct data source.'
            return(error_message)
        else:
            return(True)

    def get_raster_as_array(self):
        """returns the raster values as an array"""

        raster_source_path = self.raster.source()
        file_prepend = 'file://'
        if raster_source_path.find(file_prepend)==0:
            # Strip out 'file://' from begining of file path
            #  to pass to gdal
            raster_source_path = raster_source_path[len(file_prepend):]

        ds = gdal.Open(raster_source_path)
        ds_band = ds.GetRasterBand(1)
        nodata_val = ds_band.GetNoDataValue()
        raster_as_array = ds_band.ReadAsArray()

        #  Mask 0 pixel values from the array. These are either NoData vals or are introduced by GDAL
        raster_as_array = raster_as_array[raster_as_array > 0]

        if nodata_val:
            # If the raster has a no_data value, filter that value
            #  from the array
            raster_as_array = raster_as_array[raster_as_array != nodata_val]

        return(raster_as_array)

    def ESV_to_np_array(self):
        """Selects the most up to date ESV csv file based on the source type
            and returns it as a structured numpy array"""

        CSV_options = [csv_file for csv_file in os.listdir(self.esv_data_directory) if csv_file.endswith(".csv")]
        source_match_CSVs = [csv_file for csv_file in CSV_options if csv_file.find(self.source_name)>-1]
        latest_CSV = source_match_CSVs[0] #TODO: add latest date logic!
        csv_file = os.path.join(self.esv_data_directory, latest_CSV)

        with open(csv_file) as f:
            reader = csv.reader(f, delimiter=',')
            text = [tuple(r) for r in reader]

        first_line = 0 # If adding metadata header to the CSV file, change this value to the row with attribute names
        col_line = text[first_line]
        data = text[first_line+1:]

        max_string_len = max([max([len(entry) for entry in r]) for r in text])
        str_dtype = f'<U{max_string_len}'

        col_names = tuple(col_name for col_name in col_line)


        col_defs = np.dtype({'names':(col_names),
                            'formats':([str_dtype]*len(col_names)) })
        summary = np.array(data, dtype=col_defs)

        return(summary)

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

        raster_values, value_pixel_count = np.unique(self.ras_array, return_counts=True)
        value_area = value_pixel_count * self.cell_size()

        col_defs = np.dtype({'names':('raster_values',
                                    'value_pixel_count',
                                    'value_area'),
                            'formats':(raster_values.dtype,
                                    value_pixel_count.dtype,
                                    value_area.dtype)
                            }
                           )

        col_vals = list(zip(raster_values,
                           value_pixel_count,
                           value_area)
                       )
        summary = np.array(col_vals, dtype=col_defs)

        raster_area_dict = {k:v for k,v in list(summary[['raster_values','value_area']])} 

        ras_summary_output = {'array':summary,
                              'dict':raster_area_dict}
        return(ras_summary_output)

    def get_service_value_dict(self):
        """Returns a dictionary with the value summaries for each ecosystem service in the raster"""

        def get_field_name(matcher):
            # Finds the name of the column of the ESV file
            #  that contains the passed string: i.e. minimum, average, maximum
            fields = col_names
            field = [field for field in fields if field.lower().find(matcher)>-1][0]
            return(field)

        def make_num(x):
            # Converts an input string for ESV CSV into an integer
            #  removing thousand separator commas
            return(float(x.replace(',','')))

        ESV = self.ESV_to_np_array()
        col_names = ESV.dtype.names

        # remove reserved characters from service name so that it can be
        #  used as an attribute name in the output table
        clean = lambda x: x.lower().replace(" ", "-").replace(",", "")

        # Create a nested dictionary from the ESV and raster summary arrays
        #  first level keys are the output columns of *service*_*summary* (e.g. aesthetic_min)
        #  second level keys are the LULC type with the total value of the service for that LULC type
        output_data = defaultdict(dict)

        print(self.ras_summary_dict)
        for row in list(ESV):
            LULC_val = row[col_names.index('LULC Value')]
            if int(LULC_val) in self.ras_summary_dict:
                # Check if LULC type exists in raster. If it does add the values to the output dictionary
                service = row[col_names.index('Service')]

                for col_append, agg_esv_col in {'min':'minimum value',
                                                'avg':'average value',
                                                'max':'maximum value'}.items():
                    

                    out_col_name = f'{clean(service)}_{col_append}'
                    ESV_multiplier = row[col_names.index(get_field_name(agg_esv_col))]
                    value_out = make_num(ESV_multiplier) * self.ras_summary_dict[int(LULC_val)]
                    output_data[out_col_name][LULC_val] = value_out

        return(output_data)

    def get_output_QgsFields(self):
        """Return fields for output table as a list of QgsField objects"""

        # Create list of output fieldter
        prepend_field_names = ['LULC_code','LULC_description','pixel_count','area_hectares']
        output_table_field_names = prepend_field_names + list(self.get_service_value_dict().keys())

        # Create QgsFields object from column names
        output_table_QgsFields = QgsFields()
        for field in output_table_field_names:
            output_table_QgsFields.append(QgsField(field, QVariant.String, len=50))
        return(output_table_QgsFields)


    def create_output_table(self, sink):
        """ Reshape and write data to feature sink"""

        output_data = self.get_service_value_dict()
        ESV = self.ESV_to_np_array()
        output_table_QgsFields = self.get_output_QgsFields()


        LULC_descriptions = {i:k for i, k in ESV[['LULC Value', 'LULC Description']]}
        # Reshape the data with a row for each LULC type
        for row in self.ras_summary_dict.keys():
            # Add LULC code and descriptions to feature
            new_feature = QgsFeature(output_table_QgsFields)
            new_feature['LULC_code'] = int(row)
            new_feature['LULC_description'] = str(LULC_descriptions.get(str(row))) # Returns None if the LULC code does not exist in the ESV CSV file

            if(new_feature['LULC_description']):
                # Check if the current LULC value in the raster exists in the ESV dataset
                #  If it exists, add the total pixel count and area for each value
                raster_summary_for_LULC_val = self.ras_summary_array[self.ras_summary_array['raster_values'] == int(row)]
                pixel_count = raster_summary_for_LULC_val['value_pixel_count']
                area_hectares = raster_summary_for_LULC_val['value_area']
                new_feature['pixel_count'] = str(pixel_count[0])
                new_feature['area_hectares'] = str(area_hectares[0])

            for colname, colvalue in output_data.items():
                # Add service values to the feature
                service_value = colvalue.get(str(row))
                new_feature[colname] = str(service_value)

            sink.addFeature(new_feature, QgsFeatureSink.FastInsert)

        return(True)
       

    def check_output_format(output_raster_destination):
        """Check output file format to make sure it is a geotiff"""

        output_format = QgsRasterFileWriter.driverForExtension(splitext(output_raster_destination)[1])
        if not output_format or output_format.lower() != "gtiff":
            error_message = "CRITICAL: Currently only GeoTIFF output format allowed, exiting!"
            return error_message
        else:
            message = "Output file is GeoTIFF. Check"
            return message


    def check_for_nlcd_raster(grid, nlcd_codes, input_nodata_value, raster_value_mapping_dict):
        """Check to make sure the input raster is an NLCD raster, ie: has the right kinds of pixel values"""
        
        unique_pixel_values_of_input_raster = np.unique(grid)
        nlcd_codes.append(str(input_nodata_value))

        if all(str(i) in nlcd_codes for i in unique_pixel_values_of_input_raster):
            message = "The input raster has the correct NLCD codes for pixel values. Check"
        else:
            error_message = "The input raster's pixels aren't all legitimate NLCD codes. They must all be one of these values: " + str(nlcd_codes) + ". The raster you input had these values: " + str(unique_pixel_values_of_input_raster)
            return error_message
        output_array = copy(grid)
        for key, value in raster_value_mapping_dict.items():
            output_array[grid == key] = value
        
        return message, output_array


    def check_for_nalcms_raster(grid, nalcms_codes, input_nodata_value, raster_value_mapping_dict):
        """Check to make sure the input raster is an NALCMS raster, ie: has the right kind of pixel values"""
        
        unique_pixel_values_of_input_raster = np.unique(grid)
        nalcms_codes.append(str(input_nodata_value))

        if all(str(i) in nalcms_codes for i in unique_pixel_values_of_input_raster):
            message = "The input raster has the correct NALCMS codes for pixel values. Check"
        else:
            error_message = "The input raster's pixels aren't all legitimate NALCMS codes. They must all be one of these values: " + str(nalcms_codes) + ". The raster you input had these values: " + str(unique_pixel_values_of_input_raster)
            return error_message
        output_array = copy(grid)
        for key, value in raster_value_mapping_dict.items():
            output_array[grid == key] = value
        
        return message, output_array
        

    def check_nlcd_codes(input_esv_field, input_esv_table, input_esv_stat, input_nodata_value):
        """Checks to make sure all NLCD land use codes are valid"""
        
        raster_value_mapping_dict = {}

        input_esv_table_features = input_esv_table.getFeatures()
        nlcd_codes = ['11', '21', '22', '23', '24', '31', '41', '42', '43', '52', '71', '81', '82', '90', '95']

        for input_esv_table_feature in input_esv_table_features:
            nlcd_code = str(input_esv_table_feature.attributes()[0])
#            # Check to make sure this is a legit nlcd code. If it's not throw and error and abort the algorithm
            if nlcd_code not in nlcd_codes:
                error_message = "Found a value in the first column of the input ESV table that isn't a legitimate NLCD code: " + str(nlcd_code) + ". All the values in the first column of the input ESV table must be one of these: " + str(nlcd_codes)
                return error_message  
            try:
                selected_esv = input_esv_table_feature.attribute(input_esv_field.lower().replace(" ", "-").replace(",", "") + "_" + input_esv_stat)
            except KeyError:
                error_message = ("The Input ESV field you specified (" + input_esv_field + "_" + input_esv_stat + ") doesn't exist in this dataset. Please enter one of the fields that does exist: ")
#                error_message = (f'The input ESV field you specified: {input_esv_field}_{input_esv_stat} does not exist in this dataset. Please enter one of the fields that does exist: \n Input table fields: {str(input_esv_table.fields().names()[4:]}')
                error_data = ("Input table fields: ", str(input_esv_table.fields().names()[4:]))
                return error_message, error_data
            
#            # If there is no ESV for tis particular NLCD-ES combo Then
#            # the cell will be Null (i.e. None) and so we're dealing with
#            # that below by setting the value to 255, which is the value
#            # of the other cells that don't have values (at least for this
#            # data)

            if selected_esv == 'None':
                selected_esv = input_nodata_value

#            # If it's not null then we need to convert the total ESV for
#            # the whole area covered by that land cover (which is in USD/hectare)
#            # to the per pixel ESV (USD/pixel)
            else:
                num_pixels = int(input_esv_table_feature.attributes()[2])
                selected_esv = float(selected_esv) / num_pixels
            raster_value_mapping_dict.update({int(nlcd_code): selected_esv})
            
        return raster_value_mapping_dict, nlcd_codes


    def check_nalcms_codes(input_esv_field, input_esv_table, input_esv_stat, input_nodata_value):  #TODO: debug this function
        """Checks to make sure all NALCMS land use codes are valid"""
        
        raster_value_mapping_dict = {}

        input_esv_table_features = input_esv_table.getFeatures()
        nalcms_codes = ['1', '2', '3', '4' , '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19']

        for input_esv_table_feature in input_esv_table_features:  
            nalcms_code = str(input_esv_table_feature.attributes()[0])
#            # Check to make sure this is a legit nalcms code. If it's not throw and error and abort the algorithm
            if nalcms_code not in nalcms_codes:
                error_message = "Found a value in the first column of the input ESV table that isn't a legitimate NALCMS code: " + str(nalcms_code) + ". All the values in the first column of the input ESV table must be one of these: " + str(nalcms_codes)
                return error_message  
            try:
                selected_esv = input_esv_table_feature.attribute(input_esv_field.lower().replace(" ", "-").replace(",", "") + "_" + input_esv_stat)
            except KeyError:
                error_message = ("The Input ESV field you specified (" + input_esv_field + "_" + input_esv_stat + ") doesn't exist in this dataset. Please enter one of the fields that does exist: ")
#                error_message = (f'The input ESV field you specified: {input_esv_field}_{input_esv_stat} does not exist in this dataset. Please enter one of the fields that does exist: \n Input table fields: {str(input_esv_table.fields().names()[4:]}')
                error_data = ("Input table fields: ", str(input_esv_table.fields().names()[4:]))
                return error_message, error_data
            
#            # If there is no ESV for tis particular NLCD-ES combo Then
#            # the cell will be Null (i.e. None) and so we're dealing with
#            # that below by setting the value to 255, which is the value
#            # of the other cells that don't have values (at least for this
#            # data)

            if selected_esv == 'None':
                selected_esv = input_nodata_value

#            # If it's not null then we need to convert the total ESV for
#            # the whole area covered by that land cover (which is in USD/hectare)
#            # to the per pixel ESV (USD/pixel)
            else:
                num_pixels = int(input_esv_table_feature.attributes()[2])
                selected_esv = float(selected_esv) / num_pixels
            raster_value_mapping_dict.update({int(nalcms_code): selected_esv})
            
        return raster_value_mapping_dict, nalcms_codes

        
    def check_esv_table_length(input_esv_table):
        """Check to make sure the input ESV table has at least 4 columns"""
        
        input_esv_table_col_names = input_esv_table.fields().names()
        if len(input_esv_table_col_names) <= 4:
            error_message = "The Input ESV table should have at least 5 columns, the one you input only has " + str(len(input_esv_table_col_names))
            return error_message
        else:
            message = "Input ESV table has at least 5 columns. Check"
            return message, input_esv_table_col_names
           
        
    def check_for_esv_stats(input_esv_table_col_names):
        """Check to make sure the input ESV table has columns with valid ESV stats"""
        
        stats = ['min', 'avg', 'max']                   #changed from 'mean' to 'avg' as 'avg' was already in names of columns
        input_esv_table_esv_stat_col_names = input_esv_table_col_names[4:]
        input_esv_table_name_stats = []
        for name in input_esv_table_esv_stat_col_names:
            if len(name.split('_', 1)) > 1:
                input_esv_table_name_stats.append(name.split('_', 1)[1])
            else:
                error_message = "One or more of the columns in your Input ESV table doesn't appear to be an ESV stat. Columns 5 through the last column should all have an underscore between the ecosystem service and the statistic, e.g. aesthetic_min."
                return error_message
            
        if all(str(i) in stats for i in input_esv_table_name_stats):
            message = "The table appears to include ESV stats columns. Check"
            return message
        else:
            error_message = "One or more of the columns in your Input ESV table doesn't appear to be an ESV stat. Columns 5 through the last column should all end with \"_min\", \"_avg\", or \"_max\"."
            return error_message


    def compute_range_of_values(layer, provider, extent):
        """Takes information from active raster layer and uses that to build range of values
        for that layer. Then breaks values into 5 evenly spaced quintiles (minimum and maximum
        value for each quintile) and returns those as a tuple."""
        
        raster_stats = provider.bandStatistics(1, QgsRasterBandStats.All) 
        min_val = raster_stats.minimumValue            #minimum pixel value in layer
        max_val = raster_stats.maximumValue            #maximum pixel value in layer

        value_range = list(range(int(min_val), int(max_val+1)))           #Range of values in raster layer. Without +1 doesn't capture highest value
        value_range.sort()
        for value in value_range:                   #deletes 0 value from value range so as not to skew shading in results
            if value < raster_stats.minimumValue:
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

        return first_quintile_max, first_quintile_min, second_quintile_max, second_quintile_min, third_quintile_max, third_quintile_min, fourth_quintile_max, fourth_quintile_min, fifth_quintile_max, fifth_quintile_min
        
        
        
    def create_color_ramp_and_shade_output(layer, input_esv_field, first_quintile_max, first_quintile_min, second_quintile_max, second_quintile_min,
                          third_quintile_max, third_quintile_min, fourth_quintile_max, fourth_quintile_min, 
                          fifth_quintile_max, fifth_quintile_min):
        """Takes values for each quintile and builds raster shader with discrete color for each quintile.
        Unique color ramp chosen for each ESV value. I tried to be intuitive with the colors.
        Lastly, shades output in QGIS."""

        #green color ramp
        if input_esv_field == 'aesthetic':
            raster_shader = QgsColorRampShader()
            raster_shader.setColorRampType(QgsColorRampShader.Discrete)           #Shading raster layer with QgsColorRampShader.Discrete
            colors_list = [ QgsColorRampShader.ColorRampItem(0, QColor(255, 255, 255, .5), 'No Value'), \
                       QgsColorRampShader.ColorRampItem(first_quintile_max, QColor(204, 255, 204), f"${first_quintile_min}0 - ${first_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(second_quintile_max, QColor(153, 255, 153), f"${second_quintile_min} - ${second_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(third_quintile_max, QColor(51, 255, 51), f"${third_quintile_min} - ${third_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(fourth_quintile_max, QColor(0, 204, 0), f"${fourth_quintile_min} - ${fourth_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(fifth_quintile_max, QColor(0, 102, 0), f"${fifth_quintile_min} - ${fifth_quintile_max}0")]

        #light blue color ramp
        elif input_esv_field == 'air quality':
            raster_shader = QgsColorRampShader()
            raster_shader.setColorRampType(QgsColorRampShader.Discrete)           #Shading raster layer with QgsColorRampShader.Discrete
            colors_list = [ QgsColorRampShader.ColorRampItem(0, QColor(255, 255, 255, .5), 'No Value'), \
                       QgsColorRampShader.ColorRampItem(first_quintile_max, QColor(204, 255, 255), f"${first_quintile_min}0 - ${first_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(second_quintile_max, QColor(153, 255, 255), f"${second_quintile_min} - {second_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(third_quintile_max, QColor(51, 255, 255), f"${third_quintile_min} - {third_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(fourth_quintile_max, QColor(0, 204, 204), f"${fourth_quintile_min} - {fourth_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(fifth_quintile_max, QColor(0,102,102), f"${fifth_quintile_min} - {fifth_quintile_max}0")]
            
        #green color ramp
        elif input_esv_field == 'biodiversity':
            raster_shader = QgsColorRampShader()
            raster_shader.setColorRampType(QgsColorRampShader.Discrete)           #Shading raster layer with QgsColorRampShader.Discrete
            colors_list = [ QgsColorRampShader.ColorRampItem(0, QColor(255, 255, 255, .5), 'No Value'), \
                       QgsColorRampShader.ColorRampItem(first_quintile_max, QColor(204, 255, 229), f"${first_quintile_min}0 - ${first_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(second_quintile_max, QColor(153, 255, 204), f"${second_quintile_min} - ${second_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(third_quintile_max, QColor(51, 255, 153), f"${third_quintile_min} - ${third_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(fourth_quintile_max, QColor(0, 204, 102), f"${fourth_quintile_min} - ${fourth_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(fifth_quintile_max, QColor(0, 102, 51), f"${fifth_quintile_min} - ${fifth_quintile_max}0")]
        
        #orange color ramp
        elif input_esv_field == 'climate regulation':
            raster_shader = QgsColorRampShader()
            raster_shader.setColorRampType(QgsColorRampShader.Discrete)           #Shading raster layer with QgsColorRampShader.Discrete
            colors_list = [ QgsColorRampShader.ColorRampItem(0, QColor(255, 255, 255, .5), 'No Value'), \
                       QgsColorRampShader.ColorRampItem(first_quintile_max, QColor(255, 229, 204), f"${first_quintile_min}0 - ${first_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(second_quintile_max, QColor(255, 204, 153), f"${second_quintile_min} - ${second_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(third_quintile_max, QColor(255, 153, 51), f"${third_quintile_min} - ${third_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(fourth_quintile_max, QColor(204, 102, 0), f"${fourth_quintile_min} - ${fourth_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(fifth_quintile_max, QColor(102, 51, 0), f"${fifth_quintile_min} - ${fifth_quintile_max}0")]
        
        #orange color ramp
        elif input_esv_field == 'cultural, Other':
            raster_shader = QgsColorRampShader()
            raster_shader.setColorRampType(QgsColorRampShader.Discrete)           #Shading raster layer with QgsColorRampShader.Discrete
            colors_list = [ QgsColorRampShader.ColorRampItem(0, QColor(255, 255, 255, .5), 'No Value'), \
                       QgsColorRampShader.ColorRampItem(first_quintile_max, QColor(255, 229, 204), f"${first_quintile_min}0 - ${first_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(second_quintile_max, QColor(255, 204, 153), f"${second_quintile_min} - ${second_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(third_quintile_max, QColor(255, 153, 51), f"${third_quintile_min} - ${third_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(fourth_quintile_max, QColor(204, 102, 0), f"${fourth_quintile_min} - ${fourth_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(fifth_quintile_max, QColor(102, 51, 0), f"${fifth_quintile_min} - ${fifth_quintile_max}0")]
        
        #brown color ramp
        elif input_esv_field == 'erosion control':
            raster_shader = QgsColorRampShader()
            raster_shader.setColorRampType(QgsColorRampShader.Discrete)           #Shading raster layer with QgsColorRampShader.Discrete
            colors_list = [ QgsColorRampShader.ColorRampItem(0, QColor(255, 255, 255, .5), 'No Value'), \
                       QgsColorRampShader.ColorRampItem(first_quintile_max, QColor(220,187,148), f"${first_quintile_min}0 - ${first_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(second_quintile_max, QColor(198,168,134), f"${second_quintile_min} - ${second_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(third_quintile_max, QColor(169,144,115), f"${third_quintile_min} - ${third_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(fourth_quintile_max, QColor(138,117,93), f"${fourth_quintile_min} - ${fourth_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(fifth_quintile_max, QColor(100,85,67), f"${fifth_quintile_min} - ${fifth_quintile_max}0")]
        
        #pink color ramp
        elif input_esv_field == 'food/nutrition':
            raster_shader = QgsColorRampShader()
            raster_shader.setColorRampType(QgsColorRampShader.Discrete)           #Shading raster layer with QgsColorRampShader.Discrete
            colors_list = [ QgsColorRampShader.ColorRampItem(0, QColor(255, 255, 255, .5), 'No Value'), \
                       QgsColorRampShader.ColorRampItem(first_quintile_max, QColor(255,204,229), f"${first_quintile_min}0 - ${first_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(second_quintile_max, QColor(255,153,204), f"${second_quintile_min} - ${second_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(third_quintile_max, QColor(255,51,153), f"${third_quintile_min} - ${third_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(fourth_quintile_max, QColor(204,0,102), f"${fourth_quintile_min} - ${fourth_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(fifth_quintile_max, QColor(102,0,51), f"${fifth_quintile_min} - ${fifth_quintile_max}0")]        

        #pink color ramp
        elif input_esv_field == 'medicinal':
            raster_shader = QgsColorRampShader()
            raster_shader.setColorRampType(QgsColorRampShader.Discrete)           #Shading raster layer with QgsColorRampShader.Discrete
            colors_list = [ QgsColorRampShader.ColorRampItem(0, QColor(255, 255, 255, .5), 'No Value'), \
                       QgsColorRampShader.ColorRampItem(first_quintile_max, QColor(255,204,229), f"${first_quintile_min}0 - ${first_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(second_quintile_max, QColor(255,153,204), f"${second_quintile_min} - ${second_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(third_quintile_max, QColor(255,51,153), f"${third_quintile_min} - ${third_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(fourth_quintile_max, QColor(204,0,102), f"${fourth_quintile_min} - ${fourth_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(fifth_quintile_max, QColor(102,0,51), f"${fifth_quintile_min} - ${fifth_quintile_max}0")]


        #yellow color ramp
        elif input_esv_field == 'pollination':
            raster_shader = QgsColorRampShader()
            raster_shader.setColorRampType(QgsColorRampShader.Discrete)           #Shading raster layer with QgsColorRampShader.Discrete
            colors_list = [ QgsColorRampShader.ColorRampItem(0, QColor(255, 255, 255, .5), 'No Value'), \
                       QgsColorRampShader.ColorRampItem(first_quintile_max, QColor(255,255,204), f"${first_quintile_min}0 - ${first_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(second_quintile_max, QColor(255,255,153), f"${second_quintile_min} - ${second_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(third_quintile_max, QColor(255,255,51), f"${third_quintile_min} - ${third_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(fourth_quintile_max, QColor(204,204,0), f"${fourth_quintile_min} - ${fourth_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(fifth_quintile_max, QColor(102,102,0), f"${fifth_quintile_min} - ${fifth_quintile_max}0")]        

        #gray/black color ramp
        elif input_esv_field == 'protection from extreme events':
            raster_shader = QgsColorRampShader()
            raster_shader.setColorRampType(QgsColorRampShader.Discrete)           #Shading raster layer with QgsColorRampShader.Discrete
            colors_list = [ QgsColorRampShader.ColorRampItem(0, QColor(255, 255, 255, .5), 'No Value'), \
                       QgsColorRampShader.ColorRampItem(first_quintile_max, QColor(224,224,224), f"${first_quintile_min}0 - ${first_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(second_quintile_max, QColor(192,192,192), f"${second_quintile_min} - ${second_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(third_quintile_max, QColor(128,128,128), f"${third_quintile_min} - ${third_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(fourth_quintile_max, QColor(64,64,64), f"${fourth_quintile_min} - ${fourth_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(fifth_quintile_max, QColor(0,0,0), f"${fifth_quintile_min} - ${fifth_quintile_max}0")]       
        
        #purple color ramp
        elif input_esv_field == 'raw materials':
            raster_shader = QgsColorRampShader()
            raster_shader.setColorRampType(QgsColorRampShader.Discrete)           #Shading raster layer with QgsColorRampShader.Discrete
            colors_list = [ QgsColorRampShader.ColorRampItem(0, QColor(255, 255, 255, .5), 'No Value'), \
                       QgsColorRampShader.ColorRampItem(first_quintile_max, QColor(229,204,255), f"${first_quintile_min}0 - ${first_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(second_quintile_max, QColor(204,153,255), f"${second_quintile_min} - ${second_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(third_quintile_max, QColor(153,51,255), f"${third_quintile_min} - ${third_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(fourth_quintile_max, QColor(102,0,204), f"${fourth_quintile_min} - ${fourth_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(fifth_quintile_max, QColor(51,0,102), f"${fifth_quintile_min} - ${fifth_quintile_max}0")]
        
        #red color ramp
        elif input_esv_field == 'recreation':
            raster_shader = QgsColorRampShader()
            raster_shader.setColorRampType(QgsColorRampShader.Discrete)           #Shading raster layer with QgsColorRampShader.Discrete
            colors_list = [ QgsColorRampShader.ColorRampItem(0, QColor(255, 255, 255, .5), 'No Value'), \
                       QgsColorRampShader.ColorRampItem(first_quintile_max, QColor(255,102,102), f"${first_quintile_min}0 - ${first_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(second_quintile_max, QColor(255,51,51), f"${second_quintile_min} - ${second_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(third_quintile_max, QColor(255,0,0), f"${third_quintile_min} - ${third_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(fourth_quintile_max, QColor(204,0,0), f"${fourth_quintile_min} - ${fourth_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(fifth_quintile_max, QColor(153,0,0), f"${fifth_quintile_min} - ${fifth_quintile_max}0")]

        #red color ramp
        elif input_esv_field == 'renewable energy':
            raster_shader = QgsColorRampShader()
            raster_shader.setColorRampType(QgsColorRampShader.Discrete)           #Shading raster layer with QgsColorRampShader.Discrete
            colors_list = [ QgsColorRampShader.ColorRampItem(0, QColor(255, 255, 255, .5), 'No Value'), \
                       QgsColorRampShader.ColorRampItem(first_quintile_max, QColor(255,102,102), f"${first_quintile_min}0 - ${first_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(second_quintile_max, QColor(255,51,51), f"${second_quintile_min} - ${second_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(third_quintile_max, QColor(255,0,0), f"${third_quintile_min} - ${third_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(fourth_quintile_max, QColor(204,0,0), f"${fourth_quintile_min} - ${fourth_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(fifth_quintile_max, QColor(153,0,0), f"${fifth_quintile_min} - ${fifth_quintile_max}0")]

        
        #brown color ramp
        elif input_esv_field == 'soil formation':
            raster_shader = QgsColorRampShader()
            raster_shader.setColorRampType(QgsColorRampShader.Discrete)           #Shading raster layer with QgsColorRampShader.Discrete
            colors_list = [ QgsColorRampShader.ColorRampItem(0, QColor(255, 255, 255, .5), 'No Value'), \
                       QgsColorRampShader.ColorRampItem(first_quintile_max, QColor(220,187,148), f"${first_quintile_min}0 - ${first_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(second_quintile_max, QColor(198,168,134), f"${second_quintile_min} - ${second_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(third_quintile_max, QColor(169,144,115), f"${third_quintile_min} - ${third_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(fourth_quintile_max, QColor(138,117,93), f"${fourth_quintile_min} - ${fourth_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(fifth_quintile_max, QColor(100,85,67), f"${fifth_quintile_min} - ${fifth_quintile_max}0")]       
        
        #blue/purple color ramp
        elif input_esv_field == 'waste assimilation':
            raster_shader = QgsColorRampShader()
            raster_shader.setColorRampType(QgsColorRampShader.Discrete)           #Shading raster layer with QgsColorRampShader.Discrete
            colors_list = [ QgsColorRampShader.ColorRampItem(0, QColor(255, 255, 255, .5), 'No Value'), \
                       QgsColorRampShader.ColorRampItem(first_quintile_max, QColor(204,204,255), f"${first_quintile_min}0 - ${first_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(second_quintile_max, QColor(153,153,255), f"${second_quintile_min} - ${second_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(third_quintile_max, QColor(51,51,255), f"${third_quintile_min} - ${third_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(fourth_quintile_max, QColor(0,0,204), f"${fourth_quintile_min} - ${fourth_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(fifth_quintile_max, QColor(0,0,102), f"${fifth_quintile_min} - ${fifth_quintile_max}0")]
        
        #medium blue color ramp
        elif input_esv_field == 'water supply':
            raster_shader = QgsColorRampShader()
            raster_shader.setColorRampType(QgsColorRampShader.Discrete)           #Shading raster layer with QgsColorRampShader.Discrete
            colors_list = [ QgsColorRampShader.ColorRampItem(0, QColor(255, 255, 255, .5), 'No Value'), \
                       QgsColorRampShader.ColorRampItem(first_quintile_max, QColor(204,229,255), f"${first_quintile_min}0 - ${first_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(second_quintile_max, QColor(153,204,255), f"${second_quintile_min} - ${second_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(third_quintile_max, QColor(51,153,205), f"${third_quintile_min} - ${third_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(fourth_quintile_max, QColor(0,102,204), f"${fourth_quintile_min} - ${fourth_quintile_max}0"), \
                       QgsColorRampShader.ColorRampItem(fifth_quintile_max, QColor(0,51,102), f"${fifth_quintile_min} - ${fifth_quintile_max}0")]       

        raster_shader.setColorRampItemList(colors_list)         #applies colors_list to raster_shader
        shader = QgsRasterShader()
        shader.setRasterShaderFunction(raster_shader)       

        renderer = QgsSingleBandPseudoColorRenderer(layer.dataProvider(), 1, shader)    #renders selected raster layer
        layer.setRenderer(renderer)
        layer.triggerRepaint()
        
        

                    
                    
        
    
    
    
    
    
    
    
    
    
    
    