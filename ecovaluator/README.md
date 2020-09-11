# Overview

The QGIS EcoValuator provides a simple means of estimating the dollar value of ecosystem services, such as for recreation, water supply, food, and others.  

In the first step (estimate_ecosystem_service_values_for_study_region.py), the EcoValuator clips and attaches raster data from the National Land Cover Database (NLCD) to a user-supplied polygon or polygon layer and attaches a set of per-unit-area values based on the land cover type of each pixel. Output from this step includes a new raster including just the clipped data and a table of aggregate ecosystem service values (in dollars per year) for each land cover-and-ecosystem service combination.

In the second step (map_the_value_of_individual_ecosystem_services.py), the EcoValuator generates a raster for a chosen ecosystem service so that users can view the spatial variation in that service’s value across the study region. By repeating this step for other ecosystem services, the user can create a series of maps of different values that can be layered or stacked to show combined value for multiple services. The raster is then shaded according to the ecosystem service chosen and its values are symbolized in 5 even quintiles.

In the third step (create_print_layout_and_export_as_pdf.py), the EcoValuator creates a finished map output as a .pdf file. A print layout is created called "EcoValuator Layout" and includes a map, legend, title, subtitle, and credit text. The map is created based on the extent of the active map window in the project. The legend is created including all active layers in the layers toolbar. The title, subtitle, legend, and output are specified by the user in the Step 3 dialog box.

# More Information

See the [Key-Log Economics website](http://keylogeconomics.com/ecovaluator) or the plugin's sidebar help for more information.
