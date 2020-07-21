# Extract treecover metrics for buffer around points 

# IMPORTANT: requires GEE account, Python installation, and setting up 
# authentification keys (which was somewhat of a faff on windows..)

# Code adapted from https://philippgaertner.github.io/2019/12/earth-engine-rstudio-reticulate/
# & heavily borrows from Jacob & Jorgen's code 

## Packages ----
library(reticulate); library(dplyr); library(ggplot2); library(sf)

## Set up GEE session ----
# point reticulate to the conda environment created in GEE_setup.sh
use_condaenv('gee_interface', conda = "auto", required = TRUE)
ee <- import("ee")          # Import the Earth Engine library
ee$Initialize()             # Trigger the authentication

# Get E Cord pts ----
pts <- readRDS("../Colombia/data/points/point_ele_Eastern Cordillera.rds")

# EE datasets ----
# Load raster files (not using most currently, but have retained for future use)
# ALOS <- ee$Image('JAXA/ALOS/AW3D30/V2_2')
# SRTM30 <- ee$Image("USGS/SRTMGL1_003")
tc <- ee$Image("UMD/hansen/global_forest_change_2019_v1_7")
epsg <- ee$Projection('EPSG:4326')


## Extract Hansen treecover ----
# get 2000 treecover and loss metrics
tc2000 <- tc$select("treecover2000")
loss <- tc$select("lossyear")$unmask()

# functions for buffering and getting metrics
do_buffer <- function(x) {
    ee$Geometry$Point(c(pts$lon[x], pts$lat[x]))$buffer(buffer.width, max.error)
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
           point = pts$point[id])

# get example plots
tcs_df %>%
    # as.data.frame %>% 
    filter(id %in% 1:42) %>% 
    mutate(tc = cut(treecover2000, c(-.1, 25, 50, 75, 100.1), ordered=T)) %>%
    ggplot(aes(x, y, fill=tc)) + geom_tile() + 
    facet_wrap(~point, scales="free") +
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

# save df ----
# saveRDS(tcs_df, "data/hansen_100m_buffer.rds")
# saveRDS(tcs_df, "data/hansen_200m_buffer.rds")
saveRDS(tcs_df, "data/hansen_500m_buffer.rds")

## forest metrics (apply specifically to forest points)
library(landscapemetrics); library(raster)

forest_points <- pts %>% filter(forest == 1) %>% pull(point)

tc_metrics <- tcs_df %>%
    filter(point %in% forest_points) %>%
    mutate(tc_2000 = treecover2000 > 50, 
           tc_2018 = as.numeric(tc_2000 - (lossyear != 0) == 1)) %>%
    group_by(point) %>%
    summarise(amt_tc = sum(tc_2018), pct_tc = amt_tc/n())

for(i in 1:nrow(tc_metrics)) {
    point_i <- tc_metrics$point[i]
    x <- tcs_df %>% filter(point == point_i)
    y <- x %>%
        mutate(tc_2000 = treecover2000 > 50, 
               tc_2018 = tc_2000 - (lossyear != 0) == 1, 
               not_tc = ifelse(!tc_2018, 1, NA)) %>%
        dplyr::select(x, y, not_tc, tc_2018)
    r <- rasterFromXYZ(y %>% dplyr::select(-tc_2018))
    cl <- clump(r)
    # get clusters > 10 pixels, classify these as edge
    cl_tab <- table(cl[]) > 10
    edge <- as.numeric(names(cl_tab)[cl_tab])
    r[!(cl[] %in% edge)] <- NA 
    # get distance from edge
    d <- try(distance(r), silent = T)
    if(class(d) == "try-error") {
        tc_metrics$dist_from_edge[i] = 500   
    } else {
        tc_metrics$dist_from_edge[i] = extract(d, pts %>% filter(point == point_i))
    }
    # now get core area metric
    r2 <- rasterFromXYZ(y %>% dplyr::select(-not_tc))
    tc_metrics$pct_core_area[i] <- lsm_c_cai_mn(r2, directions=8) %>%
        filter(class==1) %>%
        pull(value)
}

saveRDS(tcs_metrics, "data/forest_cover_metrics_500m.rds")
