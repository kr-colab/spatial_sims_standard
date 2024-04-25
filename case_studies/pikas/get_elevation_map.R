library(elevatr)
library(raster)
library(png)


# coordinates of target region
center_point = c(40.471094, -105.650122)   # crystal lake, RMNP
lat_range = 0.135
long_range = 0.165

# create df
examp_df <- data.frame(rbind(cbind(center_point[2]-(long_range/2), center_point[1]-(lat_range/2)),
                             cbind(center_point[2]-(long_range/2), center_point[1]+(lat_range/2)),
                             cbind(center_point[2]+(long_range/2), center_point[1]+(lat_range/2)),
                             cbind(center_point[2]+(long_range/2), center_point[1]-(lat_range/2))))
colnames(examp_df) = c("x","y")

# get elevation
rast = get_elev_raster(examp_df, prj = 4326, z=14)  # WGS84 projection; max zoom=14
#image(rast)  # note extra pixels along periphery
rast = rasterToPoints(rast) 

# crop to focal region
cropped = rast[(rast[,1] > center_point[2]-(long_range/2)) & (rast[,1] < center_point[2]+(long_range/2)),] 
cropped = cropped[cropped[,2] > center_point[1]-lat_range/2 & cropped[,2] < center_point[1]+lat_range/2,]  

# wrap into 2d 
map = matrix(cropped[,3], nrow = length(unique(cropped[,1])), byrow = FALSE)

# crop some more to make sqare
map = map[0:dim(map)[2],]

# convert vals to (0,1)
map = (map - min(map)) / (max(map)-min(map))

# shrink
map = (map - min(map)) / (max(map)-min(map))

# rotate
map = apply(map, 1, rev)
map = apply(map, 2, rev)

# write PNG
writePNG(map, "elevation_map.png")

