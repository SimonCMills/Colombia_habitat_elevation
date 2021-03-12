# Export elevation and forest cover rasters for Eastern Cordillera for plotting
# maps. 

## Packages ----
library(rgee); library(sf); library(dplyr)

## remove prior versions
file.remove(list.files("D:/Drive/rgee_exports/", 
                       pattern="forest cover E Cord_1km|elevation E Cord_500m", 
                       full.names = TRUE))

## Set up rgee session ----
ee_Initialize()

## EE datasets ----
ALOS <- ee$Image("JAXA/ALOS/AW3D30/V2_2")$select("AVE_DSM")
tc <- ee$Image("UMD/hansen/global_forest_change_2019_v1_7")$select("treecover2000")

## Export ele and fc from earth engine ----
df_ptInfo <- readRDS("data/point_ele_Eastern Cordillera.rds") %>%
    filter(forest==1) %>%
    st_transform(crs = flat)

bbox <- st_bbox(df_ptInfo) %>%
    # st_geometry() %>%
    st_as_sfc %>%
    st_buffer(50000) %>% # note: buffer by 50 km, but in plot buffer by 30 km
    st_bbox() %>%
    st_as_sfc()

# convert to earth engine geometry
bbox_rgee <- sf_as_ee(bbox)

# clip rasters
# ALOS_flattened <- ALOS$clip(bbox_rgee)
ALOS_clipped <- ALOS$clip(bbox_rgee)
tc_clipped <- tc$clip(bbox_rgee)

# calculate forest cover at a 1000 m resolution (~corresponds to the 500 m radius)
forest <- tc$gte(50)

forestCover <- forest$reduceResolution(
    reducer = ee$Reducer$mean(),
    maxPixels = 2000
)$reproject(crs = forest$projection(), scale=1000)

# note dropping scale to 5k for purposes of plot
export_fc <- ee_image_to_drive(image = forestCover, 
                               description = "forest cover E Cord_1km", 
                               folder = "rgee_exports",
                               region = bbox_rgee, 
                               crs = "epsg:3117", # projected CRS for Colombia
                               maxPixels = 1e9)
export_fc$start()
ee_monitoring(export_fc)

## Elevation 
export_ele <- ee_image_to_drive(image = ALOS_clipped, 
                                description = "elevation E Cord_500m", 
                                folder = "rgee_exports",
                                scale = 500, # drop resolution
                                region = bbox_rgee, 
                                crs = "epsg:3117", # projected CRS for Colombia
                                maxPixels = 1e9)
export_ele$start()
ee_monitoring(export_ele)
