# Extract treecover metrics for buffer around points 

# IMPORTANT: requires GEE account, Python installation, and setting up 
# authentification keys (which was somewhat of a faff on windows..)

# Code adapted from https://philippgaertner.github.io/2019/12/earth-engine-rstudio-reticulate/
# & heavily borrows from Jacob & Jorgen's code 

## Packages ----
library(reticulate); library(dplyr)

## Set up GEE session ----
# point reticulate to the conda environment created in GEE_setup.sh
use_condaenv('gee_interface', conda = "auto", required = TRUE)
ee <- import("ee")          # Import the Earth Engine library
ee$Initialize()             # Trigger the authentication

# Get E Cord pts ----
library(dplyr)
df_bird <- readRDS("../Colombia/data/bird data_Jan&Jun2019/analysisDataset_EasternCordillera.rds")
unique_pt <- df_bird$Point %>% unique

pts <- read.csv("data/CO_sampling_points_metafile_ver2.csv") %>%
    filter(point_id %in% unique_pt)
pts <- pts[!pts$long > 0,]   # Fairly sure this is redundant now

# EE datasets ----
# Load raster files (not using most currently, but have retained for future use)
# ALOS <- ee$Image('JAXA/ALOS/AW3D30/V2_2')
# SRTM30 <- ee$Image("USGS/SRTMGL1_003")
tc <- ee$Image("UMD/hansen/global_forest_change_2019_v1_7")
epsg <- ee$Projection('EPSG:32719')


## Extract Hansen treecover ----
# get 2000 treecover and loss metrics
tc2000 <- tc$select("treecover2000")
loss <- tc$select("lossyear")$unmask()

# functions for buffering and getting metrics
do_buffer <- function(x) {
    ee$Geometry$Point(c(pts$long[x], pts$lat[x]))$buffer(buffer.width, max.error)
}

get_tc <- function(geometry_i) {
    reduced <- latlng$reduceRegion(reducer = ee$Reducer$toList(), 
                                   geometry = geometry_i,  
                                   maxPixels = 1e8, 
                                   scale = 30)
    as.data.frame(reduced$getInfo())
}

#* get buffers ----
buffer.width <- 500      # radius, in meters
max.error <- 1             # maximum error (controls number of vertices), in meters
geomcircs <- lapply(1:nrow(pts), do_buffer)

#* create object ----
# note: for lat long coords, replace pixelCoordinates with pixelLonLat
latlng = ee$Image$pixelCoordinates(projection=epsg)$addBands(tc2000)$addBands(loss)

#* extract tc values ----
# note: this is fairly slow (takes 2:2.5 minutes for 450 pts with 500m buffer)
# You can much more quickly extract all values outside of a apply loop, but it 
# wasn't immediately obvious how to then append the point identifier back to this
# total vector of values. 
tcs <- lapply(geomcircs, get_tc)
tcs_df <- tcs %>%
    bind_rows(., .id = "id") %>%
    mutate(id = as.numeric(id), 
           point_id = pts$point_id[id])

# get example plots
tcs_df %>%
    # as.data.frame %>% 
    filter(id %in% 1:20) %>% 
    mutate(tc = cut(treecover2000, c(-.1, 25, 50, 75, 100.1), ordered=T)) %>%
    ggplot(aes(x, y, fill=tc)) + geom_tile() + 
    facet_wrap(~point_id, scales="free") +
    theme_bw() +
    theme(aspect.ratio=1, 
          strip.text=element_text(hjust=0, face="bold"), 
          strip.background = element_blank(), 
          panel.border = element_blank(), 
          panel.grid=element_blank(), 
          axis.text = element_blank(), 
          axis.ticks=element_blank()) +
    labs(x="", y="")
ggsave("figures/example_tc2000.png")

# save df
saveRDS(tcs_df, "data/hansen_500m_buffer.rds")
