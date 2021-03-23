# Export elevation and forest cover rasters for Eastern Cordillera for plotting
# maps. 

## Packages ----
library(rgee); library(sf); library(dplyr)

## Set up rgee session ----
ee_Initialize()

## EE datasets ----
ALOS <- ee$Image("JAXA/ALOS/AW3D30/V2_2")$select("AVE_DSM")
tc <- ee$Image("UMD/hansen/global_forest_change_2019_v1_7")$select("treecover2000")

## Export ele and fc from earth engine ----
df_ptInfo <- readRDS("data/point_ele_Eastern Cordillera.rds") %>%
    filter(forest==1) %>%
    st_transform(crs = "epsg:3117")

bbox <- st_bbox(df_ptInfo) %>%
    # st_geometry() %>%
    st_as_sfc %>%
    st_buffer(50000) %>% 
    st_bbox() %>%
    st_as_sfc()

# convert to earth engine geometry
bbox_rgee <- sf_as_ee(bbox)

# clip rasters
ALOS_clipped <- ALOS$reproject(crs=tc$projection())$clip(bbox_rgee)
tc_clipped <- tc$clip(bbox_rgee)

# calculate forest cover at a 500 m resolution 
forest <- tc_clipped$gte(50)

forestCover <- tc_clipped$reduceResolution(
    reducer = ee$Reducer$mean(),
    maxPixels = 2000
)$reproject(crs = "epsg:3117", scale=500)

ele_500m <- ALOS_clipped$reduceResolution(
    reducer = ee$Reducer$mean(),
    maxPixels = 2000
)$reproject(crs = "epsg:3117", scale=500)
    

fc_stars <- ee_as_stars(forestCover)
ele_stars <- ee_as_stars(ele_500m)

write_stars(ele_stars, "data/ele_fc_rasters/elevation_500m.tif")
write_stars(fc_stars, "data/ele_fc_rasters/fc_500m.tif")
