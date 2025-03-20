library(here)
library(sf)
library(tidyverse)
library(rnaturalearth)
library(smoothr)
library(raster)

# doing this with a high-resolution map of the coastline gives too much resolution around the Chesapeake, Long Island, etc. I want a really coarse shape (similar to the 10m isobath approach)

# xmin=-78
# xmax=-66
ymin=35
ymax=45

xmin=-77
xmax=-66

usamap <- rnaturalearth::ne_countries(scale = "small", country = "united states of america", returnclass = "sf")[1] %>% 
  st_cast("MULTILINESTRING") %>%
  dplyr::select(geometry)

bbox1 <- st_set_crs(st_as_sf(as(raster::extent(xmin, xmax, ymin, ymax), "SpatialPolygons")), st_crs(usamap))
bbox2 <- st_set_crs(st_as_sf(as(raster::extent(-78, -74, 42, 45), "SpatialPolygons")), st_crs(usamap))
bbox3 <- st_set_crs(st_as_sf(as(raster::extent(-78,-71, 43, 48), "SpatialPolygons")), st_crs(usamap))

neusmap <- usamap %>% 
  st_intersection(bbox1)  %>%  # can replace with st_crop when CRAN version of sf updates 
  st_difference(bbox3) # get rid of extra non coastal line 
  
plot(neusmap)

smoothmap <- neusmap %>% 
  smoothr::smooth(method="ksmooth", smoothness = 8)
# smoother was applied incrementally more until the Chesapeake went away 
# https://cran.r-project.org/web/packages/smoothr/vignettes/smoothr.html

st_length(smoothmap)

plot(smoothmap)

smoothgeom <- smoothmap %>% 
  as("Spatial") %>% 
  geom()

geomdists <- pointDistance(smoothgeom[-nrow(smoothgeom), c("x", "y")], smoothgeom[-1, c("x", "y")], lonlat=TRUE) 
coastdistdat <- data.frame(smoothgeom[, c('x','y')], seglength=c(0, geomdists))

coastdistdat$lengthfromhere <- rev(cumsum(rev(coastdistdat[,"seglength"])))
# first row should match st_length(smoothmap) 
## csl: it's not

st_write(smoothmap, here("processed-data","coastline.shp"), append=FALSE)
write_rds(coastdistdat, here("Data","coastdistdat.rds"))

rm(list=ls())

####
# Snap to R
####
get_length <- function(lon, lat, distdf) {
  tmp <- distdf %>% 
    mutate(abs.diff.x2 = abs(x-lon)^2,
           abs.diff.y2 = abs(y-lat)^2,
           abs.diff.xy = sqrt(abs.diff.x2 + abs.diff.y2
    )) %>% 
    filter(abs.diff.xy == min(abs.diff.xy)) %>% 
    dplyr::select(lengthfromhere) %>% 
    pull()
  return(tmp)
}

# With COG dataframe, you would then do...
# coastdist_km = (get_length(lon=cog$lon, lat=cog$lat, distdf = coastdistdat))/1000
