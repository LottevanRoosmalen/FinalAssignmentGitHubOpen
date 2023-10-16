#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 15:49:02 2023

@author: osboxes
"""

import streamlit as st                             # Source: https://docs.streamlit.io/
from streamlit_folium import st_folium             # Source: https://github.com/randyzwitch/streamlit-folium
import folium
import geopandas as gp
import Functions as funcs 
import datetime
import time

#Set with of the app
st.set_page_config(layout="wide")                 # Source: https://discuss.streamlit.io/t/how-to-increase-the-width-of-web-page/7697              

st.title('Tomatoes in the sun')
st.divider()
st.subheader ('This app will provide you....')
st.caption('Made by....')

# input box address buttom
address = st.text_input('What is your adress? (streetname and housenumber, City, Zipcode)', 'Lawickse Allee 26, Wageningen, 6707AG')
st.write('your adress is', address)

# starting date buttom
d = st.date_input("Startingdate", datetime.date(2019, 7, 6))
st.write('Your startingdate is:', d)

# period of sunexposure buttom
option = st.selectbox(
    'For what period do you want to know the sunhours?',
    ('Day', 'Week', 'Month'))

st.write('You selected:', option)

# Status
progress_text = "Operation in progress. Please wait."
my_bar = st.progress(0, text=progress_text)

for percent_complete in range(100):
    time.sleep(0.01)
    my_bar.progress(percent_complete + 1, text=progress_text)
time.sleep(1)
my_bar.empty()

st.button("Rerun")

## Create map garden
# functions, extract location
lon, lat = funcs.extract_garden(address)

#load data
garden = 'data/Garden.geojson'
gardenGDF = gp.read_file(garden)
gardenGDF = gardenGDF.to_crs(4326)

#create map, zoom to address
gardenmap = folium.Map(location=[lat, lon], zoom_start=20)

# Add garden geojson to map
folium.GeoJson(gardenGDF).add_to(gardenmap)


# Make 2 colums for sun exposure and soil data maps
col1, col2 = st.columns(2)

with col1:
   st.header("Sun exposure per hour")
   # add map to colum
   st_data = st_folium(gardenmap, width = 750)

with col2:
   st.header("Soil quality for vegetables")
   # st_data = st_folium(....., width = 750)
   
# Add combined map
st.header("Ideal places to make your vegetable garden")
# st_data = st_folium(....., width = 750)
