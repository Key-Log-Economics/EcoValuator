# -*- coding: utf-8 -*-

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
__copyright__ = '(C) 2018 by Key-Log Economics'

# This will get replaced with a git SHA1 when you do a git archive

__revision__ = '$Format:%H$'

import os
import csv
import processing
import numpy



from PyQt5.QtGui import *


from qgis.core import (QgsProcessing,
                       QgsProcessingAlgorithm,
                       QgsProcessingParameterString,
                       QgsProcessingParameterFileDestination,
                       QgsProcessingOutputLayerDefinition,
                       QgsRasterLayer,
                       QgsProject,
                       QgsPrintLayout,
                       QgsLayoutItemMap,
                       QgsUnitTypes,
                       QgsLayoutPoint,
                       QgsLayoutSize,
                       QgsLayoutItemLegend,
                       QgsLayoutItemLabel,
                       QgsLayerTree,
                       QgsRasterBandStats,
                       QgsLayoutExporter,
                       QgsColorRampShader,
                       QgsRasterShader,
                       QgsSingleBandPseudoColorRenderer
                       )


from qgis.utils import *


__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))



class CreatePrintLayoutAndExportMap(QgsProcessingAlgorithm):
    # Constants used to refer to parameters and outputs. They will be
    # used when calling the algorithm from another algorithm, or when
    # calling from the QGIS console.
    INPUT_VECTOR = 'INPUT_VECTOR'
    INPUT_TITLE = 'INPUT_TITLE'
    INPUT_SUBTITLE = 'INPUT_SUBTITLE'
    INPUT_CREDIT_TEXT = 'INPUT_CREDIT_TEXT'
    OUTPUT_PDF_PATH = 'OUTPUT_PDF_PATH'
    OUTPUT_PDF_FILENAME_DEFAULT = 'Choose file path for pdf output'
    
    def initAlgorithm(self, config):
        """
        Here we define the inputs and output of the algorithm
        """
        #Add String as input
        self.addParameter(
            QgsProcessingParameterString(
                self.INPUT_TITLE,
                self.tr('Input title string (Optional)'),
                " "                                         #this is in place of making the dialog box "optional". Instead just gives default value as blank string
            )
        )

        #Add String as input
        self.addParameter(
            QgsProcessingParameterString(
                self.INPUT_SUBTITLE,
                self.tr('Input Subtitle (this should be returned from the ESV choice in step 2)(Optional)'),
                " "
            )
        )

        #Add String as input
        self.addParameter(
            QgsProcessingParameterString(
                self.INPUT_CREDIT_TEXT,
                self.tr('Input Credit Text (Optional)'),
                " "
            )
        )

        #Add file path as input
        self.addParameter(
            QgsProcessingParameterFileDestination(
                self.OUTPUT_PDF_PATH,
                self.tr(self.OUTPUT_PDF_FILENAME_DEFAULT),
                ".pdf"
            )
        )

    
    def processAlgorithm(self, parameters, context, feedback):
        """This actually does the processing for creating the print layout and exporting as .pdf document"""
        #needs all the arguments (self, parameters, context, feedback)
        
        log = feedback.setProgressText
        
        input_title = self.parameterAsString(parameters, self.INPUT_TITLE, context)
        input_subtitle = self.parameterAsString(parameters, self.INPUT_SUBTITLE, context)
        input_credit_text = self.parameterAsString(parameters, self.INPUT_CREDIT_TEXT, context)
        output_pdf_path = self.parameterAsString(parameters, self.OUTPUT_PDF_PATH, context)
        
        log(f"Title: {input_title}")                       


        #This creates a new print layout
        project = context.project()             
        manager = project.layoutManager()           
        layout = QgsPrintLayout(project)            
        layoutName = 'PrintLayout'                 #layoutName is going to be name of Title. Change this later

        layouts_list = manager.printLayouts()
        for layout in layouts_list:
            if layout.name() == layoutName:
                manager.removeLayout(layout)
        
        layout = QgsPrintLayout(project)
        layout.initializeDefaults()                 #create default map canvas
        layout.setName(layoutName)
        manager.addLayout(layout)



        #This adds a map item to the Print Layout
        map = QgsLayoutItemMap(layout)
        map.setRect(20, 20, 20, 20)  
        
        #Set Extent
        canvas = iface.mapCanvas()
        map.setExtent(canvas.extent())                  #sets map extent to current map canvas
        layout.addLayoutItem(map)

        #Move & Resize
        map.attemptMove(QgsLayoutPoint(5, 27, QgsUnitTypes.LayoutMillimeters))
        map.attemptResize(QgsLayoutSize(234, 178, QgsUnitTypes.LayoutMillimeters))
        
        
        #Gather visible layers in project layer tree and create a list of the map layer objects
        #Those which are not active (layers_to_remove) will subsequently remove from the legend model
        tree_layers = project.layerTreeRoot().children()
        active_layers = [layer.name() for layer in tree_layers if layer.isVisible()]
        layers_to_remove = [layer for layer in project.mapLayers().values() if layer.name() not in active_layers]


        #This adds a legend item to the Print Layout
        legend = QgsLayoutItemLegend(layout)
        layout.addLayoutItem(legend)
        legend.attemptMove(QgsLayoutPoint(235, 5, QgsUnitTypes.LayoutMillimeters))        
        #Get reference to existing legend model and root group then remove the unchecked layers
        legend.setAutoUpdateModel(False) #not sure if this line is required
        model = legend.model()
        group = model.rootGroup()
        for layer in layers_to_remove:
            group.removeLayer(layer)
        legend.adjustBoxSize()
        
        
        #this symbolizes the raster layer in the map
        #defining raster layer to work with (active layer in layer panel)
        layer = iface.activeLayer()
        print("Active Layer: ", layer.name())
        provider = layer.dataProvider()
        extent = layer.extent()
        #Using RasterBandStats to find range of values in raster layer
        stats = provider.bandStatistics(1, QgsRasterBandStats.All) 
        min_val = stats.minimumValue            #minimum pixel value in layer
        max_val = stats.maximumValue            #maximum pixel value in layer
        print("min value =", min_val)
        print("max value =", max_val)

        value_range = list(range(int(min_val), int(max_val+1)))           #Range of values in raster layer. Without +1 doesn't capture highest value
        value_range.sort()
        for value in value_range:                   #deletes 0 value from value range so as not to skew shading in results
            if value < stats.minimumValue:
                del value

        #we will categorize pixel values into 5 quintiles, based on value_range of raster layer
        #defining min and max values for each quintile. 
        #Also, values are rounded to 2 decimal places
        first_quintile_max = round(numpy.percentile(value_range, 20), 2)
        first_quintile_min = round(min_val, 2)
        second_quintile_max = round(numpy.percentile(value_range, 40), 2)
        second_quintile_min = round((first_quintile_max + .01), 2)
        third_quintile_max = round(numpy.percentile(value_range, 60), 2)
        third_quintile_min = round((second_quintile_max + .01), 2)
        fourth_quintile_max = round(numpy.percentile(value_range, 80), 2)
        fourth_quintile_min = round((third_quintile_max + .01), 2)
        fifth_quintile_max = round(numpy.percentile(value_range, 100), 2)
        fifth_quintile_min = round((fourth_quintile_max + .01), 2)


        #builds raster shader with colors_list. 
        raster_shader = QgsColorRampShader()
        raster_shader.setColorRampType(QgsColorRampShader.Discrete)           #Shading raster layer with QgsColorRampShader.Discrete
        colors_list = [ QgsColorRampShader.ColorRampItem(0, QColor(255, 255, 255), 'No Value'), \
                       QgsColorRampShader.ColorRampItem(first_quintile_max, QColor(204, 219, 255), f"{first_quintile_min} - {first_quintile_max}"), \
                       QgsColorRampShader.ColorRampItem(second_quintile_max, QColor(153, 184, 255), f"{second_quintile_min} - {second_quintile_max}"), \
                       QgsColorRampShader.ColorRampItem(third_quintile_max, QColor(102, 148, 255), f"{third_quintile_min} - {third_quintile_max}"), \
                       QgsColorRampShader.ColorRampItem(fourth_quintile_max, QColor(51, 113, 255), f"{fourth_quintile_min} - {fourth_quintile_max}"), \
                       QgsColorRampShader.ColorRampItem(fifth_quintile_max, QColor(0, 77, 255), f"{fifth_quintile_min} - {fifth_quintile_max}")]

        raster_shader.setColorRampItemList(colors_list)         #applies colors_list to raster_shader
        shader = QgsRasterShader()
        shader.setRasterShaderFunction(raster_shader)       

        renderer = QgsSingleBandPseudoColorRenderer(layer.dataProvider(), 1, shader)    #renders selected raster layer
        layer.setRenderer(renderer)
        layer.triggerRepaint()
        
        
        #This adds labels to the map
        title = QgsLayoutItemLabel(layout)
        title.setText(input_title)
        title.setFont(QFont("Arial", 28))
        title.adjustSizeToText()
        layout.addLayoutItem(title)
        title.attemptMove(QgsLayoutPoint(10, 4, QgsUnitTypes.LayoutMillimeters))

        subtitle = QgsLayoutItemLabel(layout)
        subtitle.setText(input_subtitle)
        subtitle.setFont(QFont("Arial", 17))
        subtitle.adjustSizeToText()
        layout.addLayoutItem(subtitle)
        subtitle.attemptMove(QgsLayoutPoint(11, 20, QgsUnitTypes.LayoutMillimeters))

        credit_text = QgsLayoutItemLabel(layout)
        credit_text.setText(input_credit_text)
        credit_text.setFont(QFont("Arial", 10))
        credit_text.adjustSizeToText()
        layout.addLayoutItem(credit_text)
        credit_text.attemptMove(QgsLayoutPoint(246, 190, QgsUnitTypes.LayoutMillimeters))


        #This exports a Print Layout as an image
        manager = QgsProject.instance().layoutManager()     #this is a reference to the layout Manager, which contains a list of print layouts

        layout = manager.layoutByName(layoutName)         #this accesses a specific layout, by name (which is a string)

        exporter = QgsLayoutExporter(layout)                #this creates a QgsLayoutExporter object
        exporter.exportToPdf(output_pdf_path, QgsLayoutExporter.PdfExportSettings())


        results = {}                    #All I know is processAlgorithm wants to return a dictionary
        return results

    def flags(self):
        """
        From documentation: Algorithm is not thread safe and cannot be run in a
        background thread, e.g. algorithms which manipulate the current project,
        layer selections, or with external dependencies which are not thread safe.
        """
        return super().flags() | QgsProcessingAlgorithm.FlagNoThreading

    def name(self):
        """
        Returns the algorithm name, used for identifying the algorithm. This
        string should be fixed for the algorithm, and must not be localised.
        The name should be unique within each provider. Names should contain
        lowercase alphanumeric characters only and no spaces or other
        formatting characters.
        """
        return 'Step 3: Create Print Layout and Export as .pdf'

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
        return self.tr("This step takes an output raster layer from step 2 as input and automatically produces a finished map output as a .pdf. The output will contain the map (zoomed to the extent of your current screen) and a legend which contains the active layers in the project (***NEEDS WORK***)")

    def tr(self, string):
        return QCoreApplication.translate('Processing', string)

    def helpUrl(self):
        """
        Returns the location of the help file for this algorithm. This is the
        location that will be followed when the user clicks the Help button
        in the algorithm's UI.
        """
        return "http://keylogeconomics.com/ecovaluator-help/"

    def createInstance(self):
        return CreatePrintLayoutAndExportMap()
