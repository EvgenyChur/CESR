# -*- coding: utf-8 -*-
"""
Created on Mon Dec  6 15:14:53 2021

@author: churiulin
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

# Get actual path to data
main = 'C:/Users/Churiulin/LAI_DATA'

file      = main + '/LAI_SATELLITE/MODIS/NETCDF_DATA/MODIS_LAI_GERMANY_2010.nc'
path_exit = main + '/YANDEX/'


lGIF = False

# Open NetCDF file
fh = Dataset(file, mode='r')

#Get general information about NetCF:
    
print('Main parameteres of current NetCDF: ', fh.variables.keys(), '\n')
print(f'NetCDF attributes for each parameter: {fh.dimensions.keys()} \n')
#print(f'References: ---> {fh.references} \n')                                  # Get info about NetCDF
#print(f'History: ---> {fh.history} \n')                                        # Get info about NetCDF history
#print(f'Description: ---> {fh.description} \n')                                # Get description of NetCDF file
print (f'NetCDF Title: ---> {fh.title} \n ')                                   # Get NetCDF title

# Список переменных в виде ключей
keys = fh.variables.keys()
for key in keys:
    if key in ['time', 'lon', 'lat', 'fpar', 'lai']:
        print (f'{key} ---> {fh.variables[key]} \n')
    else:
        print(f'Parameter is available, but not interesting: {key} \n')
    

# Get variables from NetCDF file based on keys
lons      = fh.variables['lon'][:]                                             # Get longitude            
lats      = fh.variables['lat'][:]                                             # Get latitude
lai       = fh.variables['lai'][:]                                             # Get LAI
lai_units = fh.variables['lai'].units                                          # Get LAI units   
#lai_1moment = fh.variables['lai'][1, :, :]
#lai_longname = fh.variables['lai'].long_name                  

# Get variables for time
time      = fh.variables['time'][:]                                            # Get timesteps
tunits    = fh.variables['time'].units                                         # Get timeunits
                      

# Get more information about time:
print(f'Get actual timesteps from NetCDF --->: {time} \n')                     # Это значения переменной "Время"
print(f'Get NetCDF timeunits: {tunits} \n')
#print(f'Время в файле netCDF выглядит как-то так: {time[0]} \n')

# Преобразование непонятных чисел в даты
dates = num2date(time[:], units = tunits)

print('Get dates type %s' % type(dates[0]), '\n')
print('Get actual date %s' % format(dates[0]), '\n')

for i in range(10):
    print ('{} {}'.format(i, dates[i]))


# Select area
sel_lon, sel_lat = [5.75, 7.25], [50.25, 51.25]

# latitude lower and upper index
latli = np.argmin( np.abs( lats - sel_lat[0] ) )
latui = np.argmin( np.abs( lats - sel_lat[1] ) ) 

# longitude lower and upper index
lonli = np.argmin( np.abs( lons - sel_lon[0] ) )
lonui = np.argmin( np.abs( lons - sel_lon[1] ) )

lai_list = []
year_list = []
month_list = []
day_list = []
for time_index in range(len(time)):
    # Get LAI for domain
    lai_2 = fh.variables['lai'][ time_index, latli:latui , lonli:lonui ] 
    lai_list.append(pd.DataFrame(lai_2).mean(axis = 1).mean())
    # Get actual date
    year_list.append(dates[time_index].year)
    month_list.append(dates[time_index].month)
    day_list.append(dates[time_index].day)
    
# Close NetCDF file
fh.close()

date_index = []
for i in range(len(year_list)): 
    current_day = pd.to_datetime(f'{year_list[i]}-{month_list[i]}-{day_list[i]}', format='%Y-%m-%d')
    date_index.append(current_day)
    
# Create LAI dataframe for actual domai
df_lai = pd.DataFrame(lai_list, index = date_index)

# Get plot
df_lai.plot(figsize = (14,10), grid = True)


# Visualization of GLOBAL LAI data
for time_index in range(len(time)):
        
    # Additional elements for work controlling
    print('Actual time step N %d:' % (time_index))
    print('----------------------')
    print('Год %d'  % dates[time_index].year  )
    print('Месяц %d'% dates[time_index].month )
    print('День %d' % dates[time_index].day   )
    print('Час %d'  % dates[time_index].hour  )
    print('Сек %d'  % dates[time_index].second)
    print('----------------------')
    
    
    # Start visualization
    fig = plt.figure(figsize=(8, 8))
    
    
    # Select projections:
    # Option 1
    m = Basemap(width = 5000000, height = 3500000, resolution='l', projection='stere', 
                lat_ts = 40, lat_0 = 50, lon_0 = 50)
            
    # Because our lon and lat variables are 1D, use meshgrid to create 2D arrays
    # Not necessary if coordinates are already in 2D arrays.
    lon, lat = np.meshgrid(lons, lats)
    xi, yi   = m(lon, lat)
    
    cs = m.contourf(xi, yi, np.squeeze(lai[time_index]))

    # Add Grid Lines
    m.drawparallels(np.arange(-80., 81., 5.),   labels = [1,0,0,0], fontsize = 10)
    m.drawmeridians(np.arange(-180., 181., 5.), labels = [0,0,0,1], fontsize = 10)
    
    # Add Coastlines, States, and Country Boundaries
    m.drawcoastlines()
    m.drawstates()
    m.drawcountries()
    
    # Add Colorbar
    cbar = m.colorbar(cs, location='bottom', pad="10%")  
    cbar.set_label(lai_units)

    current_year  = (dates[time_index].year)
    current_month = (dates[time_index].month)
    current_day   = (dates[time_index].day)
    cur_data      = f'MODIS_LAI: {current_day} {current_month} {current_year}'    
    
    # Add Title
    plt.title(f'{cur_data}')
    plt.clim(0,7)
    
     
    plt.savefig(path_exit + f'MODIS{time_index}.png', format = 'png', dpi = 300) 
    plt.close(fig)        
    plt.gcf().clear()
    
if lGIF == True:
    # Create GIF animation
    imege_list = []
    for time_index in range(len(time)):  
        # Create GIF
        imege_list.append(pil.Image.open(r'C:/Users/Churiulin/LAI_DATA/LAI/RESULT/' + 'MODIS' + str(time_index) + '.png'))
    
    imege_list[0].save('C:/Users/Churiulin/LAI_DATA/LAI/' + 'LAI.gif', format = 'GIF',
                       append_images = imege_list[1:],
                       save_all = True, duration = 300,
                       loop = 0)
else:
    print('GIF animation was not create')




