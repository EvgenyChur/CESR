# -*- coding: utf-8 -*-
"""
Script for reading and calculating information about LAI from MODIS or 
COPERNICUS satellite data

Autors of project: Evgenii Churiulin
                                                   
Current Code Owner: CESR, Evgenii Churiulin
phone:  +49  561 804-6142
fax:    +49  561 804-6116
email:  evgenychur@uni-kassel.de


History:
Version    Date       Name
---------- ---------- ----                                                   
    1.1    2021-12-10 Evgenii Churiulin, Center for Enviromental System Research (CESR)
           Initial release
"""

import os
import pandas as pd
os.environ['PROJ_LIB'] = 'C:/Users/Churiulin/Anaconda3/Library/share/proj'


from netCDF4 import Dataset
import numpy as np
import PIL as pil
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from netCDF4 import num2date, date2num 



# Get general information about NetCDF
#------------------------------------------------------------------------------
# Input  parameteres: ncfile - netcdf file
# Output parameters: none
#------------------------------------------------------------------------------

def get_start_nc_info(ncfile):
    #Get general information about NetCF:
    print('Main parameteres of current NetCDF: ', ncfile.variables.keys(), '\n' )
    print(f'NetCDF attributes for each parameter: {ncfile.dimensions.keys()} \n')
    print(f'NetCDF Title: ---> {ncfile.title} \n ')                                   # Get NetCDF title

    # Parameters list 
    keys = ncfile.variables.keys()
    for key in keys:
        if key in ['time', 'lon', 'lat', 'LAI']:
            print (f'{key} ---> {ncfile.variables[key]} \n')
        else:
            print(f'Parameter is available, but not interesting: {key} \n')
#------------------------------------------------------------------------------



# get_lai_data - the function for downloading, opening and creating the 
#                dataframes with LAI information
#------------------------------------------------------------------------------
# Input parameters : data_type      - MODIS or COPERNICUS data
#                    reseactch_area - domain
# Output parameters: df_lai         - dataframe with information about lai
#------------------------------------------------------------------------------
def get_lai_data(data_type, research_area):

    lai_list = []
    for dom_id in range(len(research_area)):
        
        # Get information about domains
        if research_area[dom_id] == 'parc':
            domain = 'PARC'     
        elif research_area[dom_id] == 'linden':
            domain = 'LINDEN'
        else:
            domain = 'LINDENBERG'
        
        # Get initial information
        if data_type == 'COPERNICUS':       
            # Name of NetCDF file
            filename = f'LAI_{domain}_2010_2015'
            outname  = f'COP_{filename}'
            # Get actual keys from NetCDF
            lon_var   = 'lon'
            lat_var   = 'lat'
            param_var = 'LAI'
            time_var  = 'time'  
        
        else:                                                                           # need settings      
            # Name of NetCDF file
            filename = f'LAI_{domain}_2010_2015'
            outname  = f'MOD_{filename}'
            # Get actual keys from NetCDF
            lon_var   = 'lon'
            lat_var   = 'lat'
            param_var = 'lai'
            time_var  = 'time'    
        
        # Get actual path to data
        main = 'C:/Users/Churiulin/LAI_DATA/LAI_SATELLITE/'
        # The absolute path to file + file
        file = main + f'{data_type}/NETCDF_DATA/{filename}.nc'
        # the path for results
        
        path_exit = main + f'RESULT/{domain}/'
        #----------------------------------------------------------------------
        
        # Main program
        #----------------------------------------------------------------------
        
        # Open NetCDF file
        fh = Dataset(file, mode='r')
        
        # Get information about NetCDF
        netcdf_info = get_start_nc_info(fh)
        
        # Get data from NetCDF
        time   = fh.variables[time_var][:]                                             # Get timesteps
        tunits = fh.variables[time_var].units                                          # Get timeunits
        # Get actual parameter --> LAI
           
        # Copernicus data
        if data_type == 'COPERNICUS':
            pname = fh.variables[param_var].long_name 
        # Modis data
        else: 
            pname = fh.variables[param_var].long_name                                  # Get LAI long name
            punits = fh.variables[param_var].units                                     # Get LAI units   
           
        # Get more information about time:
        #print(f'Get actual timesteps from NetCDF --->: {time} \n')                     
        #print(f'Get NetCDF timeunits: {tunits} \n')
        
        # Time numbers from netCDF to dates
        dates = num2date(time[:], units = tunits)
        
        # Get data from netcdf file for research parameter
        llist = []                                                                     # list for parameter data (timesteps) - lai
        ylist = []                                                                     # list for actual year                    
        mlist = []                                                                     # list for actual month
        dlist = []                                                                     # list for actual day
        for time_step in range(len(time)):
            # Get parameter values form NetCDF for research domain
            temp = fh.variables[param_var][time_step, :, :] 
            llist.append(pd.DataFrame(temp).mean(axis = 1).mean())
            # Get actual date
            ylist.append(dates[time_step].year )
            mlist.append(dates[time_step].month)
            dlist.append(dates[time_step].day  )
        
        # Close NetCDF file
        fh.close()
            
        # Get dataframes with actual data
        date_index = []
        for step_id in range(len(ylist)): 
            current_step = pd.to_datetime(f'{ylist[step_id]}-{mlist[step_id]}-{dlist[step_id]}', format='%Y-%m-%d')
            date_index.append(current_step)
            
        # Create LAI dataframe for actual domai
        df_lai = pd.DataFrame(llist, index = date_index)
        df_lai.columns = [f'lai_{research_area[dom_id]}']
        # Save data
        df_lai.to_csv(path_exit + f'{outname}.csv', float_format='%.3f', sep = ',',
                      header = ['lai'], index_label = 'date')
    
        lai_list.append(df_lai)
    df = pd.concat(lai_list, axis = 1)
    return df
        
