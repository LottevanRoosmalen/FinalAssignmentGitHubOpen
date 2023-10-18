def extract_garden (address):
    import os
    import pandas as pd
    import geopandas as gpd
    from owslib.wfs import WebFeatureService
    from geopy.geocoders import Nominatim                                         #source: https://www.askpython.com/python/python-geopy-to-find-geocode-of-an-address//www.askpython.com/python/python-geopy-to-find-geocode-of-an-address
    from shapely.geometry import Point
    import pyproj as proj
    
    ##Geocoding address
    #address we need to geocode
    loc = address
    #making an instance of Nominatim class
    geolocator = Nominatim(user_agent="my_request") 
    #applying geocode method to get the location
    location = geolocator.geocode(loc)
     
    ##convert lat lon to x y epsg:28992                                           #source: https://gis.stackexchange.com/questions/212723/how-can-i-convert-lon-lat-coordinates-to-x-y
    # setup projections
    crs_wgs = proj.Proj(init='epsg:4326') 
    crs_bng = proj.Proj(init='epsg:28992') # use locally appropriate projected CRS
    # then cast geographic coordinate pair to the projected system
    x, y = proj.transform(crs_wgs, crs_bng, location.longitude, location.latitude)
    
    ## create bounding box
    xmin, xmax, ymin, ymax = x - 50, x + 50, y - 50, y + 50
    
    ##convert geocoded location to point geometry
    addressDF = pd.DataFrame({'x': [x], 'y': [y]})                                #source: https://gis.stackexchange.com/questions/174159/converting-pandas-dataframe-to-geodataframe
    geometry = [Point(xy) for xy in zip(addressDF['x'], addressDF['y'])]
    pointGDF = gpd.GeoDataFrame(addressDF, geometry=geometry, crs="EPSG:28992")
    pointGDF.to_file(filename='data/point_address.geojson', driver = 'GeoJSON')
    
    ## load and save WFS data with bbox
    # load WFS from url
    wfsUrl = 'https://service.pdok.nl/kadaster/kadastralekaart/wfs/v5_0?request=GetCapabilities' 
    wfs = WebFeatureService(url=wfsUrl, version='2.0.0')
    # get data from WFS
    response = wfs.getfeature(typename=list(wfs.contents)[3], bbox=(xmin, ymin, xmax, ymax))
    with open('data/KadasterGebouwen.xml', 'wb') as file:
        file.write(response.read())
    response = wfs.getfeature(typename=list(wfs.contents)[0], bbox=(xmin, ymin, xmax, ymax))
    with open('data/KadasterPercelen.xml', 'wb') as file:
        file.write(response.read())
    
    ## create geodataframes of buildings and their parcels
    Kadaster_ParcelsGDF = gpd.read_file('data/KadasterPercelen.xml')
    Kadaster_BuildingsGDF = gpd.read_file('data/KadasterGebouwen.xml')
    
    ## select house and building of address 
    Address_Building = gpd.sjoin(Kadaster_BuildingsGDF, pointGDF, how="inner")
    Address_Parcel = gpd.sjoin(Kadaster_ParcelsGDF, pointGDF, how="inner")
    
    # save building to file for later use
    Address_Building.to_file(filename='data/AddressBuilding.geojson', driver='GeoJSON')
    
    ## extract garden shape from geodataframes
    Garden = Address_Parcel.difference(Address_Building, align=False)
    #write garden shape to data file
    Garden.to_file(filename='data/Garden.geojson', driver='GeoJSON')
    
    return location.longitude, location.latitude
    
    ## Remove unecessary data 
    os.remove('data/KadasterGebouwen.gfs')
    os.remove('data/KadasterGebouwen.xml')
    os.remove('data/KadasterPercelen.gfs')
    os.remove('data/KadasterPercelen.xml')
    os.remove('data/point_address.geojson')
    
def extract_BAG ():
    from owslib.wfs import WebFeatureService
    import geopandas as gpd

    ## Buffer garden (200m)
    # Take vector file from data folder
    garden = gpd.read_file("data/Garden.geojson")
    # Buffer 200m around garden                                                               
    buffered_garden_multi = garden.buffer(15)                                                   
    
    ## Transform multipolygon to single
    buffered_garden = buffered_garden_multi.unary_union
    
    ## Create bounding box based on garden with buffer
    # Obtain bounding coordinates from buffered_garden
    min_x, min_y, max_x, max_y = buffered_garden.bounds                           #source: https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoSeries.bounds.html
    # Create bounding box
    bbox = (min_x, min_y, max_x, max_y)
    
    ## Download data BAG using bounding box
    # put wfs url in variable
    wfsUrl = 'https://data.3dbag.nl/api/BAG3D/wfs?request=getcapabilities'
    # Create a WFS object
    wfs = WebFeatureService(url=wfsUrl, version='2.0.0')
    # get data for garden with buffer
    BAG = wfs.getfeature(typename=list(wfs.contents)[0], bbox=bbox, outputFormat='json')
    
    # save to disk
    with open('data/BAG.geojson', 'wb') as file:
        file.write(BAG.read())
    
    

def clip_SoilRaster():
    import requests
    import geopandas as gpd
    import fiona
    import rasterio
    import rasterio.mask
    import numpy as np
    
    ## Load garden variable
    Garden = gpd.read_file('data/Garden.geojson')
    
    ## Create bounds based on garden
    min_x, min_y, max_x, max_y = list(Garden.total_bounds)                       #source: https://geopandas.org/en/stable/docs/reference/api/geopandas.GeoSeries.total_bounds.html
    
    ## Retrieve wcs info from url
    url = f'https://data.rivm.nl/geo/ank/wcs?SERVICE=WCS&VERSION=1.0.0&REQUEST=GetCoverage&FORMAT=GeoTIFF&COVERAGE=ank:rivm_r201d_20170608_gm_fysische_geschiktheid_groenten_10m&BBOX={min_x},{min_y},{max_x},{max_y}&CRS=EPSG:28992&RESPONSE_CRS=EPSG:28992&WIDTH=1000&HEIGHT=1000'
    
    ## Load url and write to file
    img_data = requests.get(url).content
    with open('data/soilRasterGarden.tif', 'wb') as handler:
        handler.write(img_data)
        
    ## Clip soil raster to garden
    # Load garden contour with fiona 
    with fiona.open('data/Garden.geojson', "r") as shapefile:
        shapes = [feature["geometry"] for feature in shapefile]
       
    # Read soil raster, clip with garden
    with rasterio.open("data/soilRasterGarden.tif") as src:
        out_image, out_transform = rasterio.mask.mask(src, shapes, crop=True)
        out_meta = src.meta   
    
    # Save clipped raster
    out_meta.update({"driver": "GTiff",
                 "height": out_image.shape[1],
                 "width": out_image.shape[2],
                 "transform": out_transform})

    with rasterio.open('data/SoilRasterClipped.tif', "w", **out_meta) as dest:
        dest.write(out_image) 

    ## Reclassify soil Raster
    # reclassification    
    with rasterio.open('data/SoilRasterClipped.tif') as src:    
    # Read as numpy array
        array = src.read()
        profile = src.profile

    # Reclassify                                                                  #source: https://www.datacamp.com/tutorial/python-numpy-tutorial
    array[np.where((0.0 <= array) & (array < 0.2))] = 1     # Unsuitable
    array[np.where((0.2 <= array) & (array < 0.4))] = 2     # Low
    array[np.where((0.4 <= array) & (array < 0.6))] = 3     # Fair
    array[np.where((0.6 <= array) & (array < 0.8))] = 4     # Good
    array[np.where((0.8 <= array) & (array <= 1.0))] = 5    # Excellent
    array[np.where(array == -9999)] = np.nan
    
    with rasterio.open('data/SoilRasterReclassified.tif', 'w', **profile) as dst:
        # Write to disk
        dst.write(array)  
    
    from osgeo import gdal 
    import matplotlib.pyplot as plt 
      
    
    dataset = gdal.Open(r'data/SoilRasterReclassified.tif') 
    raster = dataset.ReadAsArray() 
    
    ## create image and save png
    from matplotlib.patches import Patch
    from matplotlib.colors import ListedColormap
    import matplotlib.colors as colors
    import numpy as np
    
    fig, ax = plt.subplots(figsize=(6, 6))
    
    # Define the colors you want
    cmap = ListedColormap(["Orange", "Yellow", "lightgreen", "green", "darkgreen"])
    
    # Define a normalization from values -> colors
    norm = colors.BoundaryNorm([1, 2, 3, 4, 5], 5)
    
    chm_plot = ax.imshow(raster, cmap=cmap,norm=norm)
    
    # Add a legend for labels
    legend_labels = {"Orange":      "1", 
                     "Yellow":      "2", 
                     "lightgreen":  "3",
                     "green":       "4",
                     "darkgreen":   "5"
                     }
    
    patches = [Patch(color=color, label=label)
               for color, label in legend_labels.items()]
    
    ax.legend(handles=patches,
              bbox_to_anchor=(1.35, 1),
              facecolor="white")
    
    ax.set_axis_off()
    
    plt.savefig('data/soilraster.png',bbox_inches='tight') 
    plt.close(fig)

    ## Remove starting data
    


def sunshine_one_day (start_date):
    import pybdshadow                                                                #source: https://pybdshadow.readthedocs.io/en/latest/analysis.html
    sunshine = pybdshadow.cal_sunshine(buildings,
                                       day=start_date,
                                       roof=False,
                                       accuracy=1,
                                       precision=3600)
    return (sunshine)



def calculate_sun_exposure(start_date, run_time):
    import pandas as pd
    import geopandas as gpd
    import pybdshadow
    import fiona
    from dateutil.relativedelta import relativedelta
    from geocube.api.core import make_geocube

    
    ##Prepare buidling geodataframe and garden for final clip
    #Read building GeoJSON data
    buildings = gpd.read_file(r'data/BAG.geojson') 
    
    #change crs and column names so the geodataframe is compatible with the package pybdshadow
    buildings = buildings.to_crs("EPSG:4326")
    buildings.rename(columns={"b3_h_max": "height"}, inplace=True)
    buildings.rename(columns={"fid": "building_id"}, inplace=True)
    #step required from pybdshadow documentation
    buildings = pybdshadow.bd_preprocess(buildings)                                #source: https://pybdshadow.readthedocs.io/en/latest/preprocess.html

    
    
    ##perform sunshine analysis with pybdshadow for three options: a day, the commented code after the code for one day is a work in progress 
    if run_time == 'day':
        sunshine = pybdshadow.cal_sunshine(buildings,
                                           day=start_date,
                                           roof=False,
                                           accuracy=1,
                                           precision=900)
        #transform to RD new
        sunshine.crs = 'EPSG:4326'
        sunshine_transformed = sunshine.to_crs('epsg:28992')


        #clip to garden bounds
        Garden = gpd.read_file('data/Garden.geojson')
        sunshine_vis_clip = sunshine_transformed.clip(Garden)
        
        #load buildings
        building = gpd.read_file('data/AddressBuilding.geojson')
        
        #Visualize buildings and sunshine time using matplotlib
        import matplotlib.pyplot as plt
        fig = plt.figure(1,(10,5))
        ax = plt.subplot(111)
        #define colorbar
        cax = plt.axes([0.15, 0.33, 0.02, 0.3])
        plt.title('Hour')
        #plot the sunshine time
        sunshine_vis_clip.plot(ax = ax,cmap = 'plasma',column ='Hour',alpha = 1,legend = True,cax = cax,)
        #Buildings
        building.plot(ax = ax,edgecolor='k',color=('grey'))
        plt.sca(ax)
        plt.title('Sunshine time')
        plt.savefig("data/sunshine.png")
        plt.close(fig)
        

        #write to geojson file because we need it for the final visualisation
        sunshine_vis_clip.to_file('data/Sunshine_clip.geojson', driver="GeoJSON")
        
        ## convert the sunshine_vis_clip geodataframe to a raster for later analysis: namely combining it with the soil raster
        # Using GeoCube to rasterize the Vector
        sunshine_raster = make_geocube(
            vector_data=sunshine_vis_clip,
            measurements=["Hour"],
            resolution=(-1, 1),
            fill = 0
        )

        # Save raster census raster
        sunshine_raster.rio.to_raster('data/sunshine_raster.tiff')
    ##This code below aims to enable to calculate the sun-exposure for a longer time period but was too slow: so enhancements or changes need to be made to make the code faster
'''    
    elif run_time == 'Three months':

        # Define the start date
        start_date = pd.to_datetime(start_date)

        # Calculate the end date by adding 3 months
        end_date = start_date + relativedelta(days=3)                                   #source: https://www.geeksforgeeks.org/python-pandas-tseries-offsets-dateoffset/

        # Create a date range for the three-month period
        date_range = pd.date_range(start_date, end_date)
        date_range = pd.Series(date_range.format())

        # Use the apply function to calculate the sun exposure map for each date in the range
        sun_exposure_maps = date_range.apply(sunshine_one_day)
        
        #make a list out of the created dataframes in order to stack them
        df_list = []
        for df in sun_exposure_maps:
            df_list.append(df)
        #stack the dataframes in order to calculate the mean sun hours per grid cell    
        stacked_df = pd.concat(df_list, ignore_index=True)                              #source: https://pandas.pydata.org/docs/reference/api/pandas.concat.html

        #calculate mean hour value for each grid cell
        mean_hours = stacked_df.groupby(['LONCOL', 'LATCOL'], as_index=False)['Hour'].mean()
        mean_hours = pd.merge(df_list[0], mean_hours, on=['LONCOL', 'LATCOL'], how='left', suffixes=('','_mean'))
        sunshine = mean_hours
        
        #transform to RD new
        sunshine.crs = 'EPSG:4326'
        sunshine_transformed = sunshine.to_crs('epsg:28992')


        #clip to garden bounds
        Garden = gpd.read_file('data/Garden.geojson')
        sunshine_vis_clip = sunshine_transformed.clip(Garden)
        sunshine_vis_clip.to_file('data/Sunshine_clip.geojson', driver="GeoJSON")
 
'''
def combinedata(preference):

    import requests
    import geopandas as gpd
    import fiona
    import rasterio
    import rasterio.mask
    import numpy as np
    import rioxarray as riox
    from rioxarray.merge import merge_arrays
    from rasterio.warp import reproject, Resampling, calculate_default_transform
    
    
    ## load rasters
    # Open Sun Raster
    with rasterio.open('data/sunshine_raster.tiff', 'r') as sunras:     
        ## Reclassify sunRaster   
        # Read as numpy array
        Sunarray = sunras.read()
        profile = sunras.profile
        
        # Reclassify
        Sunarray[np.where((0.0 < Sunarray) & (Sunarray < 3))] = 1  # Unsuitable
        Sunarray[np.where((3 <= Sunarray) & (Sunarray < 5))] = 2    # Low
        Sunarray[np.where((5 <= Sunarray) & (Sunarray < 7))] = 3    # Fair
        Sunarray[np.where((7 <= Sunarray) & (Sunarray < 10))] = 4   # Good
        Sunarray[np.where((10 <= Sunarray))] = 5                    # Excellent
        
    with rasterio.open('data/SunRasterReclassified.tif', 'w', **profile) as dst:
        # Write to disk
        dst.write(Sunarray)  
        
    
    ## Resample soil raster                                 source: https://pygis.io/docs/e_raster_resample.html
    def reproj_match(infile, match, outfile):
        """Reproject a file to match the shape and projection of existing raster. 
        
        Parameters
        ----------
        infile : (string) path to input file to reproject
        match : (string) path to raster with desired shape and projection 
        outfile : (string) path to output file tif
        """
        # open input
        with rasterio.open(infile) as src:
            src_transform = src.transform
            
            # open input to match
            with rasterio.open(match) as match:
                dst_crs = match.crs
                
                # calculate the output transform matrix
                dst_transform, dst_width, dst_height = calculate_default_transform(
                    src.crs,     # input CRS
                    dst_crs,     # output CRS
                    match.width,   # input width
                    match.height,  # input height 
                    *match.bounds,  # unpacks input outer boundaries (left, bottom, right, top)
                )
    
            # set properties for output
            dst_kwargs = src.meta.copy()
            dst_kwargs.update({"crs": dst_crs,
                               "transform": dst_transform,
                               "width": dst_width,
                               "height": dst_height,
                               "nodata": 0})
            print("Coregistered to shape:", dst_height,dst_width,'\n Affine',dst_transform)
            # open output
            with rasterio.open(outfile, "w", **dst_kwargs) as dst:
                # iterate through bands and write using reproject function
                for i in range(1, src.count + 1):
                    reproject(
                        source=rasterio.band(src, i),
                        destination=rasterio.band(dst, i),
                        src_transform=src.transform,
                        src_crs=src.crs,
                        dst_transform=dst_transform,
                        dst_crs=dst_crs,
                        resampling=Resampling.nearest)
    
    reproj_match('data/SoilRasterReclassified.tif', 'data/SunRasterReclassified.tif', 'data/soilResampled.tif')
    with rasterio.open('data/soilResampled.tif', 'r') as soilResampled:
        # read sunraster
        with rasterio.open('data/SunRasterReclassified.tif') as sunreclas:
        
   
            
            input_prioritization = preference
            
            # give soil and sun value
            if input_prioritization == "focus on sunhours":
                sunvalue = 2
            elif input_prioritization == "prioritize sunhours":
                sunvalue = 1.5
            elif input_prioritization == "focus on sunhours and soil quality equally":
                sunvalue = 1
            elif input_prioritization == "prioritize soil quality":
                sunvalue = 0.5
            elif input_prioritization == "focus on soil quality":
                sunvalue = 0
            else:
                sunvalue = 1  # Set a default value
            
            soilvalue = 2 - sunvalue
            
            
            
            #oilwithvalue = soilResampled.read(1) * Soilvale
            
            soil = soilResampled.read(1) * soilvalue
            sun = sunreclas.read(1) * sunvalue
            
            merged = soil + sun
    
            #write data to file
            with rasterio.open('data/merged.tif', 'w', **profile) as dst:
                # Write to disk
                dst.write(merged, 1) 

       
    from osgeo import gdal 
    import matplotlib.pyplot as plt 
    
    dataset = gdal.Open(r'data/merged.tif')
    rastermerged = dataset.ReadAsArray() 
 
    ## create image and save png
    from matplotlib.patches import Patch
    from matplotlib.colors import ListedColormap
    import matplotlib.colors as colors
    import numpy as np
     
    fig, ax = plt.subplots(figsize=(6, 6))
     
    # Define the colors you want
    cmap = ListedColormap(["Red", "orangered","Yellow", "lightgreen", "green", "darkgreen"])
     
    # Define a normalization from values -> colors
    norm = colors.BoundaryNorm([1, 2, 3, 4, 5,6,7,8,9,10], 10)
     
    chm_plot = ax.imshow(rastermerged, cmap=cmap,norm=norm)
     
    # Add a legend for labels
    legend_labels = {"Red":      "1", 
                     "Orangered": "2",
                     "darkorange": "3", 
                     "orange":  "4",
                     "gold":       "5",
                     "yellow":   "6",
                     "greenyellow": "7",
                     "limegreen": "8",
                     "green": "9",
                     "darkgreen": "10"}
     
    patches = [Patch(color=color, label=label)
               for color, label in legend_labels.items()]
     
    ax.legend(handles=patches,
              bbox_to_anchor=(1.35, 1),
              facecolor="white")
     
    ax.set_axis_off()
     
    plt.savefig('data/combinedraster.png',bbox_inches='tight') 
    plt.close(fig)
     

