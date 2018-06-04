# -*- coding: utf-8 -*-

"""
/***************************************************************************
 EcosystemServiceValuator
                                 A QGIS plugin
 Calculate ecosystem service values for a given area
 Generated by Plugin Builder: http://g-sherman.github.io/Qgis-Plugin-Builder/
                              -------------------
        begin                : 2018-04-02
        copyright            : (C) 2018 by Phil Ribbens/Key-Log Economics
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
"""

__author__ = 'Phil Ribbens/Key-Log Economics'
__date__ = '2018-04-02'
__copyright__ = '(C) 2018 by Phil Ribbens/Key-Log Economics'

# This will get replaced with a git SHA1 when you do a git archive

__revision__ = '$Format:%H$'

import numpy as np
from numpy import copy
import processing

from os.path import splitext

from PyQt5.QtCore import QCoreApplication
from qgis.core import (QgsProcessing,
                       QgsProcessingAlgorithm,
                       QgsProcessingParameterRasterLayer,
                       QgsProcessingParameterFeatureSource,
                       QgsProcessingParameterRasterDestination,
                       QgsRasterFileWriter,
                       QgsRasterLayer,
                       QgsProcessingParameterString,
                       QgsProcessingParameterNumber,
                       QgsProcessingOutputLayerDefinition,
                       QgsProcessingParameterEnum,
                       QgsMapLayerStyle
                       )

from .appinter import (Raster, App)

class CreateEcosystemServiceValueRasterAlgorithm(QgsProcessingAlgorithm):
    # Constants used to refer to parameters and outputs. They will be
    # used when calling the algorithm from another algorithm, or when
    # calling from the QGIS console.
    INPUT_RASTER = 'INPUT_RASTER'
    INPUT_NODATA_VALUE = 'INPUT_NODATA_VALUE'
    INPUT_ESV_TABLE = 'INPUT_ESV_TABLE'
    INPUT_ESV_FIELD = 'INPUT_ESV_FIELD'
    INPUT_ESV_FIELD_OPTIONS = ['aesthetic','air quality', 'biocontrol', 'climate', 'cognitive', 'cultural (other)', 'cultural service [general]', 'energy', 'erosion', 'extreme events', 'food', 'genepool', 'genetic', 'medical', 'nursery', 'pollination', 'raw materials', 'recreation', 'science and education', 'soil fertility', 'soil formation', 'tev', 'total', 'various', 'waste', 'water', 'water flows']
    INPUT_ESV_STAT = 'INPUT_ESV_STAT'
    STATS = ['min','mean','max']
    OUTPUT_RASTER = 'OUTPUT_RASTER'
    OUTPUT_RASTER_FILENAME_DEFAULT = 'Output esv raster'


    def initAlgorithm(self, config):
        """
        Here we define the inputs and output of the algorithm
        """
        # Add a parameter for the clipped raster layer
        self.addParameter(
            QgsProcessingParameterRasterLayer(
                self.INPUT_RASTER,
                self.tr('Input NLCD raster')
            )
        )

        #Add a parameter where the user can specify what the nodata value of
        # the raster they're inputting is.
        # Must be an integer
        # Default is 255
        #This will be used later to make sure that any pixels in the incoming rasterany
        # that have this value will continue to have this value in the output rasterself.
        #It's also used to give this value to any pixels that would otherwise be Null
        # in the output raster
        self.addParameter(
            QgsProcessingParameterNumber(
                self.INPUT_NODATA_VALUE,
                self.tr('Nodata value of input raster'),
                QgsProcessingParameterNumber.Integer,
                255
            )
        )

        self.addParameter(
            QgsProcessingParameterFeatureSource(
                self.INPUT_ESV_TABLE,
                self.tr('Input ESV table'),
                [QgsProcessing.TypeFile]
            )
        )

        self.addParameter(
            QgsProcessingParameterEnum(
                self.INPUT_ESV_FIELD,
                self.tr('Input ecosystem service to create raster for'),
                self.INPUT_ESV_FIELD_OPTIONS
            )
        )

        self.addParameter(
            QgsProcessingParameterEnum(
                self.INPUT_ESV_STAT,
                self.tr('Statistic to create ESV raster for'),
                self.STATS
            )
        )

        # Add a parameter for the output raster layer
        self.addParameter(
            QgsProcessingParameterRasterDestination(
                self.OUTPUT_RASTER,
                self.tr(self.OUTPUT_RASTER_FILENAME_DEFAULT),
                ".tif"
            )
        )

    def processAlgorithm(self, parameters, context, feedback):
        """
        Here is where the processing itself takes place.
        """

        log = feedback.setProgressText
        input_raster = self.parameterAsRasterLayer(parameters, self.INPUT_RASTER, context)
        input_nodata_value = self.parameterAsInt(parameters, self.INPUT_NODATA_VALUE, context)
        input_esv_field_index = self.parameterAsEnum(parameters, self.INPUT_ESV_FIELD, context)
        input_esv_field = self.INPUT_ESV_FIELD_OPTIONS[input_esv_field_index]
        input_esv_stat_index = self.parameterAsEnum(parameters, self.INPUT_ESV_STAT, context)
        input_esv_stat = self.STATS[input_esv_stat_index]

        #Append input raster filename to end of output raster filename
        if isinstance(parameters['OUTPUT_RASTER'], QgsProcessingOutputLayerDefinition):
            dest_name = self.OUTPUT_RASTER_FILENAME_DEFAULT.replace(" ", "_") + "-" + input_esv_field.replace(" ", "_") + "_" + input_esv_stat + "-" + input_raster.name()
            setattr(parameters['OUTPUT_RASTER'], 'destinationName', dest_name)

        output_raster_destination = self.parameterAsOutputLayer(parameters, self.OUTPUT_RASTER, context)
        result = { self.OUTPUT_RASTER : output_raster_destination }

        #Check that the input raster is in the right CRS
        input_raster_crs = input_raster.crs().authid()
        if input_raster_crs == "EPSG:102003":
            log("The input raster is in the right CRS: EPSG:102003. Check")
        else:
            error_message = "The input raster isn't in the right CRS. It must be in EPSG:102003. The one you input was in " + str(input_raster_crs) + "."
            feedback.reportError(error_message)
            log("")
            return {'error': error_message}

        #Check that the input raster has the right pixel size
        units_per_pixel_x = input_raster.rasterUnitsPerPixelX()
        units_per_pixel_y = input_raster.rasterUnitsPerPixelY()
        if units_per_pixel_x != 30 or units_per_pixel_y != 30:
            if round(units_per_pixel_x) == 30 and round(units_per_pixel_y) == 30:
                feedback.pushDebugInfo("Your input raster pixels weren't exactly 30x30 meters, but were close enough that the program will continue to run. Your input raster pixels were " + str(units_per_pixel_x) + "x" + str(units_per_pixel_y) + ".")
            else:
                error_message = "The input raster should have 30x30 meter pixels. The one you input has " + str(units_per_pixel_x) + "x" + str(units_per_pixel_y) + "."
                feedback.reportError(error_message)
                log("")
                return {'error': error_message}
        else:
            log("The input raster's pixel size is correct: 30x30. Check")

        input_esv_table = self.parameterAsSource(parameters, self.INPUT_ESV_TABLE, context)

        #Check to make sure the input ESV table has at least 4 columns
        input_esv_table_col_names = input_esv_table.fields().names()
        if len(input_esv_table_col_names) <= 4:
            feedback.reportError("The Input ESV table should have at least 5 columns, the one you input only has " + str(len(input_esv_table_col_names)))
            log("")
            return result
        else:
            log("Input ESV table has at least 5 columns. Check")

        #Check to make sure the input ESV table appears to have columns with ESV stats
        stats = ['min','mean','max']
        input_esv_table_esv_stat_col_names = input_esv_table_col_names[4:]
        input_esv_table_name_stats = [i.split('_', 1)[1] for i in input_esv_table_esv_stat_col_names]
        if all(str(i) in stats for i in input_esv_table_name_stats):
            log("The table appears to include ESV stats columns. Check")
        else:
            feedback.reportError("One or more of the columns in your Input ESV table doesn't appear to be an ESV stat. Columns 4 through the last column should all end with \"_min\", \"_mean\", or \"_max\".")
            log("")
            return result

        # Check output format
        output_format = QgsRasterFileWriter.driverForExtension(splitext(output_raster_destination)[1])
        if not output_format or output_format.lower() != "gtiff":
            log("CRITICAL: Currently only GeoTIFF output format allowed, exiting!")
            return result

        raster_value_mapping_dict = {}

        input_esv_table_features = input_esv_table.getFeatures()

        nlcd_codes = ['11', '21', '22', '23', '24', '31', '41', '42', '43', '52', '71', '81', '82', '90', '95']

        for input_esv_table_feature in input_esv_table_features:
            nlcd_code = input_esv_table_feature.attributes()[0]
            #Check to make sure this is a legit nlcd code. If it's not throw and error and abort the alg
            if nlcd_code not in nlcd_codes:
                error_message = "Found a value in the first column of the input ESV table that isn't a legitimate NLCD code: " + str(nlcd_code) + ". All the values in the first column of the input ESV table must be one of these: " + str(nlcd_codes)
                feedback.reportError(error_message)
                log("")
                return {'error': error_message}
            try:
                selected_esv = input_esv_table_feature.attribute(input_esv_field.replace(" ","-") + "_" + input_esv_stat)
            except KeyError:
                feedback.reportError("The Input ESV field you specified (" + input_esv_field + "_" + input_esv_stat + ") doesn't exist in this dataset. Please enter one of the fields that does exist: ")
                feedback.pushDebugInfo(str(input_esv_table.fields().names()[4:]))
                log("")
                return result
            #If there is no ESV for tis particular NLCD-ES combo Then
            # the cell will be Null (i.e. None) and so we're dealing with
            # that below by setting the value to 255, which is the value
            # of the other cells that don't have values (at least for this
            # data)
            if selected_esv is None:
                selected_esv = input_nodata_value
            #If it's not null then we need to convert the total ESV for
            # the whole area covered by that land cover (which is in USD/hectare)
            # to the per pixel ESV (USD/pixel)
            else:
                num_pixels = input_esv_table_feature.attributes()[2]
                selected_esv = int(selected_esv[1:].replace(',','')) / 0.0001 / int(num_pixels)
            raster_value_mapping_dict.update({int(nlcd_code): selected_esv})

        #Create a new raster whose pixel values are, instead of being NLCD code values, the per-pixel ecosystem service values corresponding to the NLCD codes
        log(self.tr("Reading input raster into numpy array ..."))
        grid = Raster.to_numpy(input_raster, band=1, dtype=int)
        #Check to make sure the input raster is an NLCD raster, i.e. has the right kinds of pixel values
        unique_pixel_values_of_input_raster = np.unique(grid)
        nlcd_codes.append(str(input_nodata_value))
        if all(str(i) in nlcd_codes for i in unique_pixel_values_of_input_raster):
            log("The input raster has the correct NLCD codes for pixel values. Check")
        else:
            error_message = "The input raster's pixels aren't all legitimate NLCD codes. They must all be one of these values: " + str(nlcd_codes) + ". The raster you input had these values: " + str(unique_pixel_values_of_input_raster)
            feedback.reportError(error_message)
            log("")
            return {'error': error_message}
        log(self.tr("Array read"))
        log(self.tr("Mapping values"))
        output_array = copy(grid)
        for key, value in raster_value_mapping_dict.items():
            if feedback.isCanceled():
                return result
            output_array[grid==key] = value
        log(self.tr("Values mapped"))
        log(self.tr("Saving output raster ..."))
        Raster.numpy_to_file(output_array, output_raster_destination, src=str(input_raster.source()))
        log(self.tr("Done!\n"))

        # Return the results of the algorithm. In this case our only result is
        # the feature sink which contains the processed features, but some
        # algorithms may return multiple feature sinks, calculated numeric
        # statistics, etc. These should all be included in the returned
        # dictionary, with keys matching the feature corresponding parameter
        # or output names.
        return result

    def name(self):
        """
        Returns the algorithm name, used for identifying the algorithm. This
        string should be fixed for the algorithm, and must not be localised.
        The name should be unique within each provider. Names should contain
        lowercase alphanumeric characters only and no spaces or other
        formatting characters.
        """
        return 'Step 2: Create ESV raster'

    def displayName(self):
        """
        Returns the translated algorithm name, which should be used for any
        user-visible display of the algorithm name.
        """
        return self.tr(self.name())

    def group(self):
        """
        Returns the name of the group this algorithm belongs to. This string
        should be localised.
        """
        return self.tr(self.groupId())

    def groupId(self):
        """
        Returns the unique ID of the group this algorithm belongs to. This
        string should be fixed for the algorithm, and must not be localised.
        The group id should be unique within each provider. Group id should
        contain lowercase alphanumeric characters only and no spaces or other
        formatting characters.
        """
        return 'Ecosystem service valuator'

    def shortHelpString(self):
        """
        Returns a localised short helper string for the algorithm. This string
        should provide a basic description about what the algorithm does and the
        parameters and outputs associated with it..
        """
        return self.tr("Short description of " + self.name() + " algorithm")

    def tr(self, string):
        return QCoreApplication.translate('Processing', string)

    def createInstance(self):
        return CreateEcosystemServiceValueRasterAlgorithm()
