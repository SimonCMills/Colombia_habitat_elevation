## plotting map of points. 
# note: needs post-processing to move inset figure to overlap plot (not 
# straightforward to do this without leaving plot margins of first plot running
# through inset figure; solution was to plot adjacent to figure and move afterwards,
# outside of R)

# housekeeping
library(dplyr); library(sf); library(stars); library(ggplot2); 
library(rnaturalearthhires); library(rnaturalearth)
WGS84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
flat <- "epsg:3117"
recalc_rasters <- FALSE

# only run if you want to regenerate rasters (takes 5 minutes or so)
if(recalc_rasters) {
    source("code/get_elevation_forestcover_raster.R")
}

# read rasters ----
ele <- read_stars("data/ele_fc_rasters/elevation_500m.tif")
fc <- read_stars("data/ele_fc_rasters/fc_500m.tif")

# read points ----
df_ptInfo <- readRDS("data/point_ele_Eastern Cordillera.rds") %>%
    filter(forest==1, !grepl("^CC", point), ele_jaxa >= 875) %>%
    st_transform(crs = "epsg:3117")

df_ptInfo %>% 
  st_union %>% 
  st_convex_hull %>% 
  st_area %>% 
  units::set_units(., "km^2")

pts <- st_coordinates(df_ptInfo) %>% as_tibble

# convert rasters and clip to bbox ----
bbox <- st_bbox(df_ptInfo) %>%
    st_as_sfc %>%
    st_buffer(45000) %>% # note: buffer by 50 km, but in plot buffer by 30 km
    st_bbox() %>%
    st_as_sfc()

ele_df <- ele %>% 
  st_crop(bbox) %>%
  st_as_sf(., as_points = TRUE, merge = FALSE)

fc_df <- fc %>% 
  st_crop(bbox) %>%
  st_as_sf(., as_points = TRUE, merge = FALSE) 

both <- st_join(fc_df, ele_df) %>% 
    rename(fc = 1, ele = 2) %>%
    bind_cols(., as_tibble(st_coordinates(.))) %>%
    mutate(ele_f = cut(ele, seq(0, 5000, 1000)))

# colour scales ----
greens <- RColorBrewer::brewer.pal(9, "Greens")
n_levels <- length(levels(both$ele_f))
fc_cols <- colorRampPalette(greens[3:9])(n_levels)
ele_cols <- colorRampPalette(c("grey90", "grey10"))(n_levels)

# create separate rasters for plot ----
fc_df_for_plot <- both %>%
    filter(fc > 30) %>%
  as_tibble

ele_clip <- both %>%
    filter(ele >= 4000)

# other features ----
V_border <- ne_countries(country="Colombia", scale=10) %>% 
  st_as_sf() %>%
  st_transform(., flat) %>%
  st_crop(., bbox) %>%
  as_Spatial() %>%
  fortify(.)

Bogota_df <- data_frame(y=4.7110, x=-74.0721) %>%
  st_as_sf(., coords=c("x", "y"), crs=WGS84) %>%
  st_transform(., crs=flat) %>%
  mutate(label = "Bogotá")%>%
  st_coordinates(.) %>%
  as_tibble

# get coordinates for margins ----
dimensions_WGS84 <- st_bbox(bbox %>% st_transform(., WGS84))

x_coords <- expand.grid(long_WGS84 = seq(-74.5, -72.5, 1), lat_WGS84 = dimensions_WGS84["ymin"]) %>%
  as.matrix
y_coords <- expand.grid(long_WGS84 = dimensions_WGS84["xmin"], lat_WGS84 = seq(4, 8, 2)) %>%
  as.matrix

x_sp <- x_coords %>% 
  st_multipoint() %>%
  st_geometry
st_crs(x_sp) <- WGS84
x_coords_flattened <- st_transform(x_sp, flat) %>% 
  st_coordinates() %>% 
  as_tibble

## now y
y_sp <- y_coords %>% 
  st_multipoint() %>%
  st_geometry
st_crs(y_sp) <- WGS84

y_coords_flattened <- st_transform(y_sp, flat) %>% 
  st_coordinates() %>% 
  as_tibble

y_trans <- bind_cols(as_tibble(y_coords), y_coords_flattened) %>%
  mutate(lat_WGS84 = paste0(lat_WGS84, "\u00B0"))
x_trans <- bind_cols(as_tibble(x_coords), x_coords_flattened) %>%
  mutate(long_WGS84 = paste0(long_WGS84, "\u00B0"))

# plot 1 ----
map_fc <- ggplot(fc_df_for_plot) + 
  geom_tile(data=ele_clip, aes(X, Y), fill="grey80") +
  geom_tile(aes(X, Y, fill=ele_f)) +
  geom_polygon(data = V_border, aes(long, lat), fill=NA, col="black") + #geom_tile(aes(X, Y, fill = ele_f)) +
  scale_fill_manual(values=fc_cols) +
  geom_point(data=pts, aes(X, Y), pch=21, fill="grey90", size=2) +
  scale_x_continuous(expand=c(0,0), 
                     breaks=pull(x_trans, X), 
                     labels = pull(x_trans, long_WGS84)) + 
  scale_y_continuous(expand=c(0,0),
                     breaks=pull(y_trans, Y), 
                     labels = pull(y_trans, lat_WGS84)) +
  labs(x="Lat", y="Long", fill="Elevation") +
  geom_point(data=Bogota_df, aes(X, Y), size=3, shape=17) + 
  # geom_label(data=Bogota_df, aes(X,Y+20000), label="Bogotá", label.size=0, fill=NA) +
  theme_bw() +
  theme(axis.title = element_blank(), 
        axis.text=element_text(colour="black"), 
        panel.grid = element_blank()) +
  guides(fill="none") +
  coord_equal() 

# plot 2 ----
hist_fc <- ggplot(df_ptInfo, aes(ele_jaxa, fill=cut(ele_jaxa, seq(0, 5000, 1000)))) + 
  geom_histogram(binwidth=200, col="black", boundary=0) +
  scale_fill_manual(values = fc_cols) +
  scale_y_continuous(expand=c(0,0), breaks=c(0, 20, 40)) +
  scale_x_continuous(expand = c(0,0)) +
  labs(x="", y="") +
  guides(fill="none") +
  theme(axis.text = element_text(colour="black"), 
        panel.background = element_blank(), 
        plot.background = element_blank(), 
        panel.grid = element_blank(), 
        axis.text.x = element_text(margin=margin(2,0,0,0))) 

# plot 3 ----
# create bounding box for inset figure (all of CO plus outlying areas)
# just removes far-away areas to prevent self intersections when projected for 
# Colombia
left <- -90
right <- -60
top <- 13
bottom <- -8
bbox_inset <- matrix(c(left, top, 
                       left, bottom, 
                       right, bottom, 
                       right, top, 
                       left, top), 
                     byrow=T, ncol=2) %>%
    list() %>%
    st_polygon() %>%
    st_sfc(., crs="WGS84")

# get Colombia bounding box (this will be the extent of the inset plot)
CO_buffer <- ne_countries(country="Colombia") %>%
    st_as_sf %>%
    st_transform(., crs=flat) %>%
    st_buffer(50000) %>%
    st_bbox()

# get country outlines
country_borders <- rnaturalearth::ne_countries(scale=10) %>%
  st_as_sf %>%
  st_crop(., bbox_inset) %>%
  st_transform(., crs=flat) %>% 
  st_crop(., CO_buffer)

Bogota <- data_frame(y=4.7110, x=-74.0721) %>%
    st_as_sf(., coords=c("x", "y"), crs=WGS84) %>%
    st_transform(., crs=flat) %>%
    mutate(label = "Bogotá")
    
# get mountain ranges
mountains <- readRDS("../range_shifts/data/mountain_polygons.rds") %>% 
  st_union() %>%
  st_transform(., crs=flat) %>%
  st_buffer(1) %>%
  st_crop(., country_borders) %>%
  st_transform(., crs=flat) %>%
  st_crop(CO_buffer)
  
mountains_ftfy <- mountains %>%
  as_Spatial %>%
  fortify()

# Plot inset figure (3) ---
mapplot <- ggplot(Bogota) + 
    geom_sf(data=mountains, fill="grey90", col=NA) +
    geom_sf(data=test, fill=NA, col="black", alpha=1) +
    geom_sf(data=bbox, fill=NA, col="black", alpha=1) +
    theme_void() +
    geom_sf(size=2, shape=17, col="black", alpha=1)

ggsave("figures/map_inset.png", mapplot, height=200*.2, width=150*.2, units="mm")

# Join all 3 plots ----
# because the plot is coord_equal() scaled, you have to maintain the dimensions
# when inserting the grob. 
# note: inset is plotted adjacent to main figure- this is moved to overlap figure
# in post-processing. 
dims <- bbox %>% st_bbox()

plot_all3 <- map_fc +
  theme(plot.margin = unit(c(2, 0, 1, 0), "cm")) +
  annotation_custom(
    ggplotGrob(hist_fc),
    xmin = 909666 - 200000, xmax = 900000,
    ymax = y_trans$Y[1] + 8e4, ymin= y_trans$Y[1] -4e4
  ) +
  annotation_custom(
    ggplotGrob(mapplot), 
    xmin = 7e+05 - 3e5 , xmax = 7e+05 - 1e5, 
    ymax = 15e5 + 2e5/(dims["ymax"] - dims["ymin"]), ymin= 13e5 #bbox2["ymin"] 
  ) 

ggsave("figures/map_fullplot_v3.png", plot_all3, height=200*.8, width=250*.8, 
       units="mm", dpi=400)
