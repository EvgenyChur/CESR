# -*- coding: utf-8 -*-
"""
Script works with NetCDF data based on MODIS satellite information about LAI.

The program allows to open NetCDF data, get data based on key elements and
create plots with actual satellite values for Germany.

The program has:
    1. plot_netcdf_mod - the function for visialization of netcdf data - lai


Acknowledgements: Stefan Kern


Current Code Owner: CESR, Evgenii Churiulin
phone:  +49  561 804-6142
fax:    +49  561 804-6116
email:  evgenychur@uni-kassel.de


History:
Version    Date       Name
---------- ---------- ----                                                   
    1.1    2021-12-13 Evgenii Churiulin, Center for Enviromental System Research (CESR)
           Initial release
"""


import os
os.environ['PROJ_LIB'] = 'C:/Users/Churiulin/Anaconda3/Library/share/proj'


from netCDF4 import Dataset
import numpy as np
import PIL as pil
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from netCDF4 import num2date, date2num


def plot_netcdf_mod(main, lGIF = False):
    '''
    

    Parameters
    ----------
    main : string
        the main path for data project.
    lGIF : logic, optional
        logical paramters for GIF plot visualization. The default is False.

    Returns
    -------
    None.

    '''
    
    year_list = [2010, 2011, 2012, 2013, 2014, 2015]
        
    for year in year_list:
        file      = main + f'MODIS/NETCDF_DATA/MODIS_LAI_GERMANY_{year}.nc'
        path_exit = main + 'RESULT/FIGURES/MODIS/'
        
        # Open NetCDF file
        fh = Dataset(file, mode='r')
                
        # Get variables from NetCDF file based on keys
        lons         = fh.variables['lon'][:]                                          # Get longitude            
        lats         = fh.variables['lat'][:]                                          # Get latitude
        lai          = fh.variables['lai'][:]                                          # Get LAI
        lai_units    = fh.variables['lai'].units                                       # Get LAI units   
        lai_longname = fh.variables['lai'].long_name                  
        time         = fh.variables['time'][:]                                         # Get timesteps
        tunits       = fh.variables['time'].units                                      # Get timeunits
                                  
            
        # Преобразование непонятных чисел в даты
        dates = num2date(time[:], units = tunits)
        
        # Close NetCDF file
        fh.close()
        
        # Visualization of GLOBAL LAI data
        for time_index in range(len(time)):
            
            # Get current date
            current_year  = (dates[time_index].year)
            current_month = (dates[time_index].month)
            current_day   = (dates[time_index].day)
            cur_data      = f'{current_day}/{current_month}/{current_year}'
            out_name      = f'_{current_year}_{current_month}_{current_day}'
            
            # Additional elements for work controlling
            print('-' * 40)
            print(f'Actual time step is {cur_data}') 
            print('-' * 40)    
            
            # Start visualization
            fig = plt.figure(figsize=(8, 8))
            
            # Select projections:
            m = Basemap(width = 1200000   , height = 600000, resolution='l', 
                        projection='stere', lat_ts = 52.5, lat_0 = 51.5, lon_0 = 10)
            
            # Because our lon and lat variables are 1D, use meshgrid to create 2D arrays
            # Not necessary if coordinates are already in 2D arrays.
            lon, lat = np.meshgrid(lons, lats)
            xi, yi   = m(lon, lat)
            
            cs = m.contourf(xi, yi, np.squeeze(lai[time_index]))
        
            # Add Grid Lines
            m.drawparallels(np.arange(40.0, 60.1, 1.0), 
                            labels = [1,0,0,0]        , 
                            fontsize = 10             )
            m.drawmeridians(np.arange(-10., 81.0, 2.5), 
                            labels = [0,0,0,1]        ,
                            fontsize = 10             )
            
            # Add Coastlines, States, and Country Boundaries
            m.drawlsmask(land_color  = 'antiquewhite',
                         ocean_color = 'aqua'        ,
                         lakes = True)
            m.drawcoastlines()
            m.drawstates()
            m.drawcountries() 
            
            # Add Colorbar
            #cbar = m.colorbar(cs, location='bottom', pad="10%")  
            
            colormesh = m.pcolormesh(lon, lat, lai[time_index], vmin = 0, vmax = 7)
            cbar = m.colorbar(colormesh, location = 'bottom', pad="10%")    
            cbar.set_label(lai_units)    
            
            # Add Title
            plt.title(f'{lai_longname} {cur_data}')
        
            plt.savefig(path_exit + f'MODIS{out_name}.png', format = 'png', dpi = 300) 
            plt.close(fig)        
            plt.gcf().clear() 
           
        if lGIF == True:
            # Create GIF animation
            imege_list = []
            for time_index in range(len(time)):  
                # Create GIF
                current_year  = (dates[time_index].year)
                current_month = (dates[time_index].month)
                current_day   = (dates[time_index].day)
                out_name      = f'_{current_year}_{current_month}_{current_day}'
                
                imege_list.append(pil.Image.open(r'C:/Users/Churiulin/LAI_DATA/LAI/RESULT/FIGURES/MODIS/' + 'MODIS' + str(out_name) + '.png'))
                
            imege_list[0].save('C:/Users/Churiulin/LAI_DATA/LAI/RESULT/FIGURES/GIF/' + 'MODIS' + str(year) + '.gif', format = 'GIF',
                               append_images = imege_list[1:],
                               save_all = True, duration = 300,
                               loop = 0)
            print('GIF animation was created')
        
        else:
            print('GIF animation was not created')




