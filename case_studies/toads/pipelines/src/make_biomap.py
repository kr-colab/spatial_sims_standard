import sys
import argparse
import pandas as pd
import numpy as np
import geopandas as gpd
import rasterio
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import matplotlib.cm as cm
import matplotlib.mlab as mlab
from shapely.geometry import Polygon, LineString, Point, MultiPoint
import rasterio


def parse_args(args):
    parser = argparse.ArgumentParser(
        description="Little script to plot bioclim variable onto a shapefile."
    )

    parser.add_argument(
        "--bioclim_tiff", 
        type=str, 
        help="tiff file with desired bioclim variable"
    )

    parser.add_argument(
        "--shapefile", 
        type=str, 
        help="shapefile of region to plot. location_string should match"
    )

    parser.add_argument(
        "-r", "--resolution", 
        type=int,
        default = 200,
        required= False,
        help=
        """
        How finely mapped the climatic variable should be.
        Default is 200, results in 200 X 200 grid over desired area.
        """
    )

    parser.add_argument(
        "-l", "--location_string", 
        type = str,
        default = "110,155,-45,-10",
        required= False,
        help = 
        """
        comma separated values specifying locations to map bioclim values over
        i.e. '110,155,-45-10' for australia
        """
    )

    parser.add_argument(
        "--outfile", 
        type=str, 
        help="Name of the output file containing the map"
    )

    return parser.parse_args()


#converting to and from lat/long to locations on a map image

def map_to_scale(a, c, x, y, z):
    """
    solve for b with relationships between two intervals
    a---b---c 
    x---------y---------z
    preserving relative distances across points according to:
    (z-y)/(z-x) = (c-b)/(c-a)
    solving or b:
    b = c - (z-y)/(z-x)*(c-a) 
    """
    return c - (z-y)/(z-x)*(c-a) 

def scaled_location(map_location, map_dim, image_dim):
    """
    given dimensions on a map and image,
    convert location on a map to location on an image
    """
    map_x, map_y = map_location
    map_x0, map_x1, map_y0, map_y1 = map_dim
    image_x0, image_x1, image_y0, image_y1 = image_dim
    image_x = map_to_scale(image_x0, image_x1, map_x0, map_x, map_x1)
    image_y = map_to_scale(image_y0, image_y1, map_y0, map_y, map_y1)
    return image_x, image_y

def test_scaled_location():
    location = (147.0, -19.0)
    map_dims = (110.0, 155.0, -45.0, -10.0)
    image_dims = (0.0, 1999.0, 0.0, 1999.0)
    #from latlong to image coords
    to_image = scaled_location(
        location,
        map_dims,
        image_dims
    )
    #from output image coord back to latlong
    to_map = scaled_location(
        to_image,
        image_dims,
        map_dims,
    )

    assert(location == to_map)

def make_envmap(latitude_range, longitude_range, raster, shape, resolution = [200, 200]):
    """
    Construct map of environment at specfic coordinates

    #:param str config: Configuration file listing input parameters shared
    #    across commands.
    :param latitude_range: The min and max latitude for the location.
    :param longitude_range: The min and max longitude for the location.
    :param raster: Raster object containing environment values by location 
        generated using rasterio.open(). 
        Has only been tested with bioclim tif files as input.
    :param shape: geopandas shape file to clip map to. 
        Must fall with the ranges specified.
    :param List[int] resolution: width (longitude) and length (latitude) of rectangle 
        where environmental values should be extracted from. 
    """
    x = np.linspace(longitude_range[0], longitude_range[1], resolution[0])
    y = np.linspace(latitude_range[0], latitude_range[1], resolution[1])
    pairs = np.meshgrid(x, y)
    dims = len(x) * len(y)
    coord_list = [(x,y) for x,y in zip(pairs[0].reshape(dims,), pairs[1].reshape(dims,))]
    out = [x for x in raster.sample(coord_list)]
    env_df = pd.DataFrame([(coord_list[i][1],coord_list[i][0],out[i][0]) for (i, _) in enumerate(out)])
    env_df.columns = ['lat', 'long', 'variable']

    geo_df = gpd.GeoDataFrame(
            env_df, geometry=gpd.points_from_xy(env_df.long, env_df.lat)) \
            .clip(shape)

    return geo_df


if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    aus = gpd.read_file(args.shapefile)
    fp = args.bioclim_tiff
    img = rasterio.open(fp)
    
    matplotlib.rcParams['figure.figsize'] = 20,20

    rez = args.resolution
    x1, x2, y1, y2 = [int(i) for i in args.location_string.split(",")]

    gdf = make_envmap(
        latitude_range = [y1,y2], 
        longitude_range = [x1,x2], 
        raster = img, 
        shape = aus, 
        resolution = [rez, rez]
    )

    fig, ax = plt.subplots(1, 1)
    ax.scatter(gdf.long, gdf.lat, c = gdf.variable, marker = 's', cmap="gray")
    ax.axis('off')
    ax.set_xlim(x1, x2)
    ax.set_ylim(y1, y2)
    ax.margins(0)
    plt.savefig(args.outfile, pad_inches=0, bbox_inches='tight')
    plt.close()

    land_fig_name = args.outfile.split(".")[0].split("_")[0] + "_land.png"
    fig, ax = plt.subplots(1, 1)
    ax.scatter(gdf.long, gdf.lat, c = [1] * len(gdf.variable), marker = 's', cmap="gray")
    ax.axis('off')
    ax.set_xlim(x1, x2)
    ax.set_ylim(y1, y2)
    ax.margins(0)
    plt.savefig(land_fig_name, pad_inches=0, bbox_inches='tight')
    plt.close()


