import csv
from osgeo import gdal
import numpy as np
import os
from collections import defaultdict

from qgis.core import (QgsUnitTypes,
                       QgsRasterLayer,  
                       QgsFields,
                       QgsField,
                       QgsFeature,
                       QgsFeatureSink)
                       

from PyQt5.QtCore import (QVariant,
                          QCoreApplication,
                          QFileInfo)

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

    def check_pixel_size(self, input_raster):
        """Checks to make sure pixel size of input raster are valid 30x30m pixels"""
        units_per_pixel_x = input_raster.rasterUnitsPerPixelX()
        units_per_pixel_y = input_raster.rasterUnitsPerPixelY()
        if units_per_pixel_x != 30 or units_per_pixel_y != 30:
            if round(units_per_pixel_x) == 30 and round(units_per_pixel_y) == 30:
                feedback = "Your input raster pixels weren't exactly 30x30 meters, but were close enough that the program will continue to run. Your input raster pixels were " + str(units_per_pixel_x) + "x" + str(units_per_pixel_y) + "."
                return feedback
            else:
                error_message = "The input raster should have 30x30 meter pixels. The one you input has " + str(units_per_pixel_x) + "x" + str(units_per_pixel_y) + "."
                return {'error': error_message}
        else:
            feedback = "The input raster's pixel size is correct: 30x30. Check"
            return feedback
        
    def esv_research_data_list(self, input_esv_data_location, input_esv_file):
        input_esv_table = []
        with open(os.path.join(input_esv_data_location, input_esv_file), newline='') as f:
            reader = csv.reader(f)
            for row in reader:
                input_esv_table.append(row)
                
        return input_esv_table
    
    def check_number_of_columns_in_esv_table(self, esv_table):
        if len(esv_table[0]) != 6:
            error_message = "The Input table of ESV research data should have 6 columns, the one you input has " + str(len(esv_table[0]))
            return {'error': error_message}
        else:
            message = "Input table of ESV research data has 6 columns. Check"
            return message  
        

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
        
        
    def build_color_ramp(self, input_esv_field):
    	"""Gets input_esv_field string and returns custom color ramp based on that input. Each value range is broken into 5 quintiles and each
    	quintile has a discrete color along a spectrum/color ramp (light green to dark green for example). Hex codes are used to for each color. Some 
    	ecosystem services naturally lent themselves to a color choice (ex: biodiversity = green color ramp). Others were just random colors that 
    	had not been used yet."""
    	pass