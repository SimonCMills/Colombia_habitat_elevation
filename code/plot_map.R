##
library(mapview); library(leaflet); library(dplyr); library(sf); library(htmltools)
library(ggspatial)
WGS84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
flat <- "epsg:3117"
# create colours for points
col_lookup <- data_frame(type = c("F", "P", "A"),
                         cols = c("#32827B", "#EEDFB8", "#9F7032"))

df_bird <- readRDS("../Colombia/data/bird data_Jan&Jun2019/analysisDataset_EasternCordillera.rds") %>%
    st_transform(crs = flat)
unique_pt <- df_bird$point %>% unique

# weird- corrupted or something
df_ptInfo <- readRDS("data/point_ele_Eastern Cordillera.rds") %>%
    filter(forest==1) %>%
    st_transform(crs = flat)
#df_ptInfo[287,]

library(raster)
ele <- raster("data/elevation_Colombia_3.tif")
ele2 <- flip(ele, direction="y")
ele_df <- rasterToPoints(ele) %>% 
    as_tibble() #%>%
    #st_as_sf(., coords=c("x", "y"), crs=WGS84) #%>%
    # st_transform(crs = flat)
# ele_df <- bind_cols(ele_df, as_tibble(st_coordinates(ele_df)))

Bogota <- data_frame(y=4.7110, x=-74.0721) %>%
    st_as_sf(., coords=c("x", "y"), crs=WGS84) %>%
    st_transform(., crs=flat) %>%
    st_coordinates(.) %>%
    as_tibble
#scalebar
scale <- tibble(x=c(945000 - 105000, 945000 - 55000), y=c(826000 + 90000,826000 + 90000))

rgdal::ogrListLayers("../Colombia/col_admbnda_adm0-1-2.gdb/col_admbnda_adm0-1-2.gdb")
polys <- st_read("../Colombia/col_admbnda_adm0-1-2.gdb/col_admbnda_adm0-1-2.gdb",
                 layer="col_admbnda_adm1_unodc_ocha") #%>%
    # filter(admin1Name_es %in% c("Boyacá", "Cundinamarca"))

clipped <- st_crop(polys %>% st_transform(crs(ele)), ele) %>%
    st_transform(flat) %>% 
    st_union()

pts <- st_coordinates(df_ptInfo) %>% as_tibble
p1 <- ggplot(ele_df) + 
    geom_tile(aes(x, y, fill=elevation_Colombia_3)) + 
    scale_fill_gradientn(colours = c("#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", "#FFFFBF", "#FEE090", "#FDAE61", "#F46D43", "#D73027")) +
    geom_sf(data=clipped, fill=NA, col="grey20") +
    geom_point(data=pts, aes(X, Y), pch=21, fill="#A5D6A7", size=2) +
    scale_x_continuous(expand=c(0,0), breaks=-72.5:-74.5) + #limits = c(5.9e+05, NA)) + 
    scale_y_continuous(expand=c(0,0)) +
    labs(x="Lat", y="Long", fill="Elevation") +
    #geom_line(data=scale, aes(x, y)) +
    #geom_label(data=data.frame(x=865000, y=916000, label="50 km"), aes(x=x, y=y, label=label), vjust=0, label.size=0, fill=NA) +
    geom_point(data=Bogota, aes(X, Y), size=2, shape=17) + 
    geom_label(data=Bogota, aes(X,Y+20000), label="Bogotá", label.size=0, fill=NA) +
    theme(axis.title = element_blank(), axis.text=element_text(colour="black")) +
    coord_sf(xlim = c(540000 + 60000, 945000 - 50000), 
             ylim = c(826000 + 80000, 1426000 - 1000)) +
    guides(fill=F)
p1
# ggsave("figures/map_firstdraft.png")


ncols <- length(levels(cut(df_ptInfo$ele_jaxa, seq(0, 4000, 100), right=F)))
cols <- colorRampPalette(c("#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", "#FFFFBF", 
                           "#FEE090", "#FDAE61", "#F46D43", "#D73027"))(ncols)

library(ggthemes)
p2 <- ggplot(df_ptInfo, aes(ele_jaxa, fill=cut(ele_jaxa, seq(0, 4000, 100), right=F))) + 
    geom_histogram(binwidth=100, col="black", boundary=200) +
    scale_fill_manual(values = cols) +
    scale_y_continuous(breaks=seq(0, 100, 10), expand = expansion(mult = 0.01)) +
    scale_x_continuous(expand=expansion(mult = 0.01)) +
    labs(x="", y="") +
    guides(fill=F) +
    theme(axis.text = element_text(colour="black"), 
          panel.background = element_blank(), 
          plot.background = element_blank(), 
          panel.grid = element_blank()) 

dims <- c(x= c(540000 + 60000, 945000 - 50000), y = c(826000 + 80000, 1426000 - 1000))
x_inset <- c(dims["x1"] + diff(dims[c("x1", "x2")])/3, dims["x2"] - 8000)
y_inset <- c(dims["y1"] - 20000, dims["y1"] + 100000 - 20000)

p_a_both <- p1 + 
    annotation_custom(
        ggplotGrob(p2), 
        xmin = x_inset[1] + 5000, xmax = x_inset[2], ymax = y_inset[1], ymin= y_inset[2]
    )

ggsave("figures/map_and_elevations.png", p_a_both, height=200*.8, width=150*.8, units="mm")
plot
