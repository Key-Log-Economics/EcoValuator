# Overview

The QGIS EcoValuator provides a simple means of estimating the dollar value of ecosystem services, such as for recreation, water supply, food, and others.  

In the first step (estimate_ecosystem_service_values_for_study_region.py), the EcoValuator clips and attaches raster data from the National Land Cover Database (NLCD) to a user-supplied polygon or polygon layer and attaches a set of per-unit-area values based on the land cover type of each pixel. Output from this step includes a new raster including just the clipped data and a table of aggregate ecosystem service values (in dollars per year) for each land cover-and-ecosystem service combination.

In the second step (map_the_value_of_individual_ecosystem_services.py), the EcoValuator generates a raster for a chosen ecosystem service so that users can view the spatial variation in that service’s value across the study region. By repeating this step for other ecosystem services, the user can create a series of maps of different values that can be layered or stacked to show combined value for multiple services.

In the third step (create_print_layout.py), which is optional, users can choose to use raster output from step 2 and automatically generate a print layout as a .pdf document. This step takes one raster layer as input. The .pdf will show the current extent at which users are viewing the map and a legend with all active map layers.



# More Information

See the [Key-Log Economics website](http://keylogeconomics.com/ecovaluator-help/) or the plugin's sidebar help for more information.
