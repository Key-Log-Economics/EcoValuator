﻿# -*- coding: utf-8 -*-

"""
/***************************************************************************
 EcoValuator
                                 A QGIS plugin
 Calculate ecosystem service values for a given area
 Generated by Plugin Builder: http://g-sherman.github.io/Qgis-Plugin-Builder/
                              -------------------
        begin                : 2018-04-02
        copyright            : (C) 2018 by Key-Log Economics
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

__author__ = 'Key-Log Economics'
__date__ = '2018-04-02'
__copyright__ = '(C) 2018 by Key-Log Econnomics'

# This will get replaced with a git SHA1 when you do a git archive

__revision__ = '$Format:%H$'

import os
import processing

from PyQt5.QtCore import (QCoreApplication,
                          QFileInfo,
                          QVariant
                          )

from qgis.core import (QgsProcessing,
                       QgsProcessingAlgorithm,
                       QgsProcessingParameterFile,
                       QgsProcessingParameterFeatureSink,
                       QgsFields,
                       QgsField,
                       QgsFeature,
                       QgsFeatureSink,
                       QgsProcessingParameterRasterLayer,
                       QgsProcessingParameterFileDestination,
                       QgsProcessingOutputLayerDefinition,
                       QgsRasterLayer,
                       QgsProcessingParameterEnum,
                       QgsProcessingParameterVectorLayer,
                       QgsProcessingParameterRasterDestination
                       )

from .eco_valuator_classes import LULC_dataset, ESV_dataset

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
__esv_data_location__ = os.path.join(__location__, "esv_data")


class EstimateEcosystemServiceValuesForStudyRegion(QgsProcessingAlgorithm):
    """ Constants used to refer to parameters and outputs. They will be
     used when calling the algorithm from another algorithm, or when
     calling from the QGIS console."""
    INPUT_RASTER = 'INPUT_RASTER'
    MASK_LAYER = 'MASK_LAYER'
    CLIPPED_RASTER = 'CLIPPED_RASTER'
    CLIPPED_RASTER_FILENAME_DEFAULT = 'Output clipped raster'
    HTML_OUTPUT_PATH = 'HTML_OUTPUT_PATH'
    INPUT_ESV = 'INPUT_ESV'
    INPUT_LULC_SOURCE = 'INPUT_LULC_SOURCE'

    # The LULC_SOURCES dictionary contains data settings that depend
    #  on the data source such as projection data and ESV validation
    # map unit of 0 = meters
    with ESV_dataset() as esv:
        LULC_SOURCES = esv.get_lulc_sources()

    OUTPUT_RASTER_SUMMARY_TABLE = 'OUTPUT_RASTER_SUMMARY_TABLE'
    OUTPUT_RASTER_SUMMARY_TABLE_FILENAME_DEFAULT = 'Output table of raster unique values'
    OUTPUT_ESV_TABLE = 'OUTPUT_ESV_TABLE'
    OUTPUT_ESV_TABLE_FILENAME_DEFAULT = 'Output ESV table'

    def initAlgorithm(self, config):
        """
        Here we define the inputs and output of the algorithm
        """
        self.addParameter(
            QgsProcessingParameterEnum(
                self.INPUT_LULC_SOURCE,
                self.tr('Select land use/land cover data source'),
                self.LULC_SOURCES
            )
        )
        self.addParameter(
            QgsProcessingParameterRasterLayer(
                self.INPUT_RASTER,
                self.tr('Input land cover raster')
            )
        )
        # Input vector to be mask for raster
        self.addParameter(
            QgsProcessingParameterVectorLayer(
                self.MASK_LAYER,
                self.tr('Input mask layer'),
                [QgsProcessing.TypeVectorAnyGeometry]
            )
        )

        self.addParameter(
            QgsProcessingParameterRasterDestination(
                self.CLIPPED_RASTER,
                self.tr('Clipped raster layer')
            )
        )

        self.addParameter(
            QgsProcessingParameterFeatureSink(
                self.OUTPUT_ESV_TABLE,
                self.tr('Output ESV table')
            )
        )

    def processAlgorithm(self, parameters, context, feedback):
        """
        Here is where the processing itself takes place.
        """
        log = feedback.setProgressText

        input_raster = self.parameterAsRasterLayer(parameters, self.INPUT_RASTER, context)

        # Check that the input raster has been loaded correctlyvalidate
        if not input_raster.isValid():
            error_message = "Layer failed to load."
            feedback.reportError(error_message)
            return {'error': error_message}

        # Set the LULC data source and input vector layer
        input_lulc_source_index = self.parameterAsEnum(parameters, self.INPUT_LULC_SOURCE, context)
        input_lulc_source = self.LULC_SOURCES[input_lulc_source_index]
        input_vector = self.parameterAsVectorLayer(parameters, self.MASK_LAYER, context)
        
        feedback.pushDebugInfo(f'Processing for {input_lulc_source} data')

        # Append input raster and input vector filenames to end of output clipped raster filename
        if isinstance(parameters['CLIPPED_RASTER'], QgsProcessingOutputLayerDefinition):
            dest_name = input_raster.name() + "-CLIPPED_BY-" + input_vector.name()
            setattr(parameters['CLIPPED_RASTER'], 'destinationName', dest_name)
        elif isinstance(parameters['CLIPPED_RASTER'], str):  # for some reason when running this as part of a model parameters['OUTPUT_ESV_TABLE'] isn't a QgsProcessingOutputLayerDefinition object, but instead is just a string
            if parameters['CLIPPED_RASTER'][0:7] == "memory:":
                parameters['CLIPPED_RASTER'] = input_raster.name() + "-CLIPPED_BY-" + input_vector.name()

        clipped_raster_destination = self.parameterAsOutputLayer(parameters, self.CLIPPED_RASTER, context)

        # Clip the input raster by the input mask layer (vector)
        log("Clipping raster...")
        processing.run("gdal:cliprasterbymasklayer", {'INPUT': input_raster, 'MASK': input_vector.source(), 'ALPHA_BAND': False, 'CROP_TO_CUTLINE': True, 'KEEP_RESOLUTION': False, 'DATA_TYPE': 0, 'OUTPUT': clipped_raster_destination}, context=context, feedback=feedback)
        log("Done clipping raster.")

        # Summarize the raster, i.e. calculate the pixel counts and total area for each NLCD value
        log("Summarizing raster...")
        clipped_raster = QgsRasterLayer(clipped_raster_destination)

        #Make instance of LULC dataset from clipped layer here
        LULC_clipped_raster = LULC_dataset(input_lulc_source, clipped_raster)

        # Check that raster is valid for given data source
        valid = LULC_clipped_raster.is_valid()
        if isinstance(valid, str):
            #If is instance returns a string it is not valid. The string contains the error message
            error_message = valid
            feedback.reportError(error_message)
            return {'error': error_message}

        with ESV_dataset() as ESV_data:
            # pass the LULC area summary data from LULC raster object to ESV object
            LULC_area_summary = LULC_clipped_raster.raster_summary
            data_for_table = ESV_data.get_LULC_evaluation_data(LULC_area_summary, input_lulc_source)

        # Make output table
        output_table_QgsFields = QgsFields()
        for field in data_for_table['column_names']:
            output_table_QgsFields.append(QgsField(field, QVariant.String, len=50))

        (sink, dest_id) = self.parameterAsSink(parameters, self.OUTPUT_ESV_TABLE, context, output_table_QgsFields)

        result = {self.CLIPPED_RASTER: clipped_raster_destination,
                  self.OUTPUT_ESV_TABLE: dest_id}

        for record in data_for_table['data']:
            # convert from tuple to list  for setAttributes method
            record = list(record)
            new_feature = QgsFeature(output_table_QgsFields)
            new_feature.setAttributes(record)
            sink.addFeature(new_feature, QgsFeatureSink.FastInsert)

        # Return the results of the algorithm, which includes the clipped raster
        # and the output esv table
        return(result)

    def name(self):
        """
        Returns the algorithm name, used for identifying the algorithm. This
        string should be fixed for the algorithm, and must not be localised.
        The name should be unique within each provider. Names should contain
        lowercase alphanumeric characters only and no spaces or other
        formatting characters.
        """
        return 'Step 1: Estimate ecosystem service values for study region'

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
        return 'EcoValuator'

    def shortHelpString(self):
        """
        Returns a localised short helper string for the algorithm. This string
        should provide a basic description about what the algorithm does and the
        parameters and outputs associated with it..
        """
        return self.tr("This algorithm does several things:\n 1. Clips the input land use/land cover (LULC) raster to the extent of the input mask layer.\n 2. Calculates how much area each type of land cover accounts for in the now-clipped LULC raster layer.\n 3. Multiplies those areas by each of the associated per-hectare ecosystem service values (ESV), which is provided on the back end.\n Outputs: Include the clipped raster layer and a table of aggregate ESV values for each land cover type in the study region.\n ~~~~~~~~~~~~~~~~ \n Inputs:\n Select land use/land cover data source: Choose either NLCD (National Land Cover Dataset) or NALCMS (North American Land Change Monitoring System).\n Input land cover raster: Supply your NLCD or NALCMS raster layer (Info Below)\n Input mask layer: Supply mask layer for your area of interest. This should be a vector polygon.\n Clipped raster layer: This output is the result of the clipping process. You can specify an output location or save to temporary file.\n Output ESV table:\n (LEAVE THIS FIELD BLANK) \nThis table contains land cover vales and descriptions as rows and associated ecosystem service values based on the minimum, mean, and maximum values per hectare from the ESV research data. Not that many NULL values will appear in th etable due to a lack of existing research on certain ecosystem services in each land cover type. NULL does not correspond to a dollar value of 0. It simply indicated a current lack of primary studies to determine the dollar amount.\n ~~~~~~~~~~~~~~~~ \n Note on data sources:\n NLCD data available for download here:\n https://www.usgs.gov/centers/eros/science/national-land-cover-database?qt-science_center_objects=0#qt-science_center_objects\n NALCMS data available for download here:\n http://www.cec.org/tools-and-resources/map-files/land-cover-30m-2015-landsat-and-rapideye")
        


    def tr(self, string):
        return QCoreApplication.translate('Processing', string)

    def helpUrl(self):
        """
        Returns the location of the help file for this algorithm. This is the
        location that will be followed when the user clicks the Help button
        in the algorithm's UI.
        """
        return "http://www.keylogeconomics.com/ecovaluator.html"

    def createInstance(self):
        return EstimateEcosystemServiceValuesForStudyRegion()
