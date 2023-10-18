#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 14:42:14 2023

@author: osboxes
"""
import os
import Functions as funcs 
import streamlit as st                             # Source: https://docs.streamlit.io/
import geopandas as gp
import datetime
import time

### Set layout app
## Set width of the app
st.set_page_config(layout="wide")                 # Source: https://discuss.streamlit.io/t/how-to-increase-the-width-of-web-page/7697              
    
st.title('Your garden plan')
st.divider()
st.subheader ('Sun exposure and soil quality of your garden')
st.caption('Made by....')
    
## input box address buttom
address = st.text_input('WHAT IS YOUR ADDRESS? (streetname and housenumber, City, Zipcode)', 'Nassauweg 16, 6703 CH Wageningen')
st.write('your adress is', address)
    
## starting date buttom
date = st.date_input("STARTING DATE", datetime.date(2019, 7, 6))
dateString = date.strftime('%m/%d/%Y')
st.write('Your startingdate is:', dateString) 
st.caption('take into account the growing period of the vegetables (may till august)')

## Status running time
progress_text = "Operation in progress. Please wait."
my_bar = st.progress(0, text=progress_text)

for percent_complete in range(100):
    time.sleep(0.01)
    my_bar.progress(percent_complete + 1, text=progress_text)
time.sleep(1)
my_bar.empty()
    
st.button("Rerun")   

## announcement goal when we would have more time
st.divider()
st.subheader('UNDER CONSTRUCTION')
st.caption('We are busy expanding our sun hours analysis. With the next update of the app, you will be able to choose periods during which the sun hours are calculated for you. Like a one-month period or the growing period of your favourite vegetables.')
st.divider()

   
### Run main script part
# Create data and output folders
if not os.path.exists('data'):
    os.makedirs('data')
if not os.path.exists('output'):
    os.makedirs('output')
  
## Call functions from main With input from streamlit    
# call function Extract Garden
funcs.extract_garden(address)
# call function Clip BAG
funcs.extract_BAG()
# call function CLip soil raster
funcs.clip_SoilRaster()
# call function sunexposure
funcs.calculate_sun_exposure(dateString, 'day')
    
## Make visualization
import matplotlib.image as mpimg

img = mpimg.imread('data/sunshine.png')
img2 = mpimg.imread('data/soilraster.png')


### Part two app layout
## Make 2 colums for sun exposure and soil data maps
col1, col2 = st.columns(2)

with col1:
    st.header("Sun exposure per hour")
    st.image(img, caption=None, width=750, use_column_width=None, clamp=False, channels="RGB", output_format="auto")
    # add map to colum
    #st_data = st_folium(basemap1, width = 750)

with col2:
    st.header("Soil quality for vegetables")
    st.image(img2, caption=None, width=750, use_column_width=None, clamp=False, channels="RGB", output_format="auto")
    # add map to colum
    #st_data = st_folium(basemap2, width = 750)
    # add info
    st.write("There is offcourse a possibility that the soil in your garden does not have a good soil quality for vegetables. If this is the case, do not let it scare you!! You can still make your own vegetable garden, but maybe consider to make a planter.")
    st.caption("if your soil quality is low, you can consider to not pay attention to it looking at suitable places vor your garden (set preference to only sunhours)")

## Add header
st.divider()
st.header("Ideal places to make your vegetable garden")

## Add slider for prioritization of the variables
Choice = st.select_slider(
    'Select your preferences',
    options=['focus on sunhours', 'prioritize sunhours', 'focus on sunhours and soil quality equally', 'prioritize soil quality', 'focus on soil quality'])
st.write('you decided to', Choice, 'for the location of your vegetable garden')

## Add text
st.subheader("based on your preferences we provide you with the most suitable locations of your garden to make a vegetable garden. We wish you all the best in making your garden and enjoy your freshly grown veggies!!")


### run code and make map for optimal place vegetable garden
## run function combine data to location vegetable garden
#funs.   (Choice)
#img3 = mpimg.imread('data/.........png')


### part three app layout
## add map
#st.image(img3, caption=None, width=750, use_column_width=None, clamp=False, channels="RGB", output_format="auto")
