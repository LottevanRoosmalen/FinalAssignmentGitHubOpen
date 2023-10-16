def extract_garden (address):
    import os
    import pandas as pd
    import geopandas as gpd
    from owslib.wfs import WebFeatureService
    from geopy.geocoders import Nominatim      #source: https://www.askpython.com/python/python-geopy-to-find-geocode-of-an-address//www.askpython.com/python/python-geopy-to-find-geocode-of-an-address
    from shapely.geometry import Point
    import pyproj as proj
    
    ##Geocoding address
    #address we need to geocode
    loc = address
    #making an instance of Nominatim class
    geolocator = Nominatim(user_agent="my_request") 
    #applying geocode method to get the location
    location = geolocator.geocode(loc)
     
    ##convert lat lon to x y epsg:28992                  source: https://gis.stackexchange.com/questions/212723/how-can-i-convert-lon-lat-coordinates-to-x-y
    # setup projections
    crs_wgs = proj.Proj(init='epsg:4326') 
    crs_bng = proj.Proj(init='epsg:28992') # use locally appropriate projected CRS
    # then cast geographic coordinate pair to the projected system
    x, y = proj.transform(crs_wgs, crs_bng, location.longitude, location.latitude)
    
    ## create bounding box
    xmin, xmax, ymin, ymax = x - 50, x + 50, y - 50, y + 50
    
    ##convert geocoded location to point geometry
    addressDF = pd.DataFrame({'x': [x], 'y': [y]})      #source: https://gis.stackexchange.com/questions/174159/converting-pandas-dataframe-to-geodataframe
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
    buffered_garden_multi = garden.buffer(30)                                                   
    
    ## Transform multipolygon to single
    buffered_garden = buffered_garden_multi.unary_union
    
    ## Create bounding box based on garden with buffer
    # Obtain bounding coordinates from buffered_garden
    min_x, min_y, max_x, max_y = buffered_garden.bounds
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
    
    ## Load garden variable
    Garden = gpd.read_file('data/Garden.geojson')
    
    ## Create bounds based on garden
    min_x, min_y, max_x, max_y = list(Garden.total_bounds)
    
    ## Retrieve wcs info from url
    url = f'https://data.rivm.nl/geo/ank/wcs?SERVICE=WCS&VERSION=1.0.0&REQUEST=GetCoverage&FORMAT=GeoTIFF&COVERAGE=ank:rivm_r201d_20170608_gm_fysische_geschiktheid_groenten_10m&BBOX={min_x},{min_y},{max_x},{max_y}&CRS=EPSG:28992&RESPONSE_CRS=EPSG:28992&WIDTH=2&HEIGHT=2'
    
    ## Load url and write to file
    img_data = requests.get(url).content
    with open('data/soilRasterGarden.tif', 'wb') as handler:
        handler.write(img_data)
            
    ## Remove starting data