##
library(mapview); library(leaflet); library(dplyr); library(sf); library(htmltools)

library(ggspatial)
# create colours for points
col_lookup <- data_frame(type = c("F", "P", "A"),
                         cols = c("#32827B", "#EEDFB8", "#9F7032"))

df_bird <- readRDS("../Colombia/data/bird data_Jan&Jun2019/analysisDataset_EasternCordillera.rds")
unique_pt <- df_bird$Point %>% unique

# weird- corrupted or something
df_ptInfo <- read.csv("data/CO_sampling_points_metafile_ver2.csv") %>%
    as_tibble() %>%
    rename(Point = point_id) %>%
    filter(Point %in% unique_pt) %>%
    select(Point, lat, long, ele=ALOSelev) %>%
    mutate(type = substr(Point, 3, 3))

df_ptInfo %>% filter(Point %in% c("VIF4", "VIF5", "VIF6"))
pts %>% filter(name %in% c("VIF4", "VIF5", "VIF6"))

WGS84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
flat <- "epsg:3117"

pts <- st_as_sf(df_ptInfo, coords=c("long", "lat"), crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") %>%
    st_transform(crs = flat)
pts <- bind_cols(pts, as_tibble(st_coordinates(pts)))

library(raster)
ele <- raster("data/elevation_Colombia_3.tif")
ele2 <- flip(ele, direction="y")
ele_df <- rasterToPoints(ele) %>% 
    as_tibble() #%>%
    #st_as_sf(., coords=c("x", "y"), crs=WGS84) #%>%
    # st_transform(crs = flat)
ele_df <- bind_cols(ele_df, as_tibble(st_coordinates(ele_df)))

ggplot(ele_df) + 
    geom_tile(aes(x, y, fill=elevation_Colombia_3)) + 
    scale_fill_gradientn(colours = c("#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", "#FFFFBF", "#FEE090", "#FDAE61", "#F46D43", "#D73027")) +
    geom_point(data=pts, aes(X, Y), pch=21) +
    scale_x_continuous(expand=c(0,0)) + #limits = c(5.9e+05, NA)) + 
    scale_y_continuous(expand=c(0,0)) +
    labs(x="Lat", y="Long", fill="Elevation") +
    coord_equal() +
    geom_line(data=scale, aes(x, y)) +
    geom_label(data=data.frame(x=875000, y=910000, label="50 km"), aes(x=x, y=y, label=label), label.size=0, fill=NA) +
    
    geom_point(data=Bogota, aes(X, Y), size=3) + 
    geom_label(data=Bogota, aes(X,Y+20000), label="Bogotá", label.size=0, fill=NA) +
    theme(axis.ticks = element_blank(), axis.title = element_blank(), axis.text=element_blank()) 

#scalebar
scale <- tibble(x=c(850000, 900000), y=c(900000,900000))
900000 - 850000

Bogota <- data_frame(y=4.7110, x=-74.0721) %>%
    st_as_sf(., coords=c("x", "y"), crs=WGS84) %>%
    st_transform(., crs=flat) %>%
    st_coordinates(.) %>%
    as_tibble


ggsave("figures/map_firstdraft.png")

df_ptInfo %>%
    mutate(type_full = case_when(type=="F" ~ "Forest", 
                                 type == "A" ~ "Paramo", 
                                 type == "P" ~ "Pasture"), 
           type_full = factor(type_full, levels=c("Forest", "Pasture", "Paramo"))) %>%
    ggplot(., aes(ele)) + 
    geom_histogram(binwidth=100, fill="white", colour="black") + 
    facet_wrap(~type_full, ncol=1) +
    theme(strip.text=element_text(hjust=0, face="bold"), 
          strip.background = element_rect(fill=NA, colour=NA)) +
    labs(x="Elevation", y="Frequency")

ggsave("figures/elevations.png")

