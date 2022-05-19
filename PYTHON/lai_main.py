# -*- coding: utf-8 -*-
"""
Created on Fri Dec 10 17:17:01 2021

@author: churiulin
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Import personal libraries
import sys
sys.path.append('C:/Users/Churiulin/Python/scripts/CESR_project')
import vis_module     as vis                                                   # module for visualisation    
import get_lai_data   as gld                                                   # 
import lai_copernicus as lcop                                                  # module for visualisation of LAI data from COPERNICUS in NetCDF
import lai_modis      as lmod                                                  # module for visualisation of LAI data from MODIS      in NetCDF
import lai_cosmo      as lac                                                   # module for processing of LAI data from COSMO


#------------------------------------------------------------------------------
# Get initial information about data
#------------------------------------------------------------------------------

# Get actual path to data
mf_com = 'C:/Users/Churiulin/LAI_DATA'

# Data path for actual folders
sat_path = f'{mf_com}/LAI_SATELLITE/'                                          # Lai data based on satellite data
cos_path = f'{mf_com}/LAI_COSMO/'                                              # Lai data based on COSMO data
c_exit   = f'{cos_path}/RESULTS/'




# logical parameters:
lget_data     = True                                                           # get LAI data  based on satellite data
lcosmo_data   = True                                                           # get LAI data  based on COSMO data    
lplot_basemap = False                                                          # get LAI plots based on NetCDF data
lplot_cosmo   = True                                                           # get LAI plots based on COSMO data
lplot4line    = True

# COSMO-CLM parameters for analysis
#------------------------------------------------------------------------------
clm_param = ['SUR_LAI' , 'LAI_SHA'    , 
             'LAI_SUN' , 'LAI'        , 
             'PLCOV'   , 'RUNOFF_G'   , 
             'RUNOFF_S', 'SUR_BIOMASS']    

# Names of COSMO-CLM plot labels
#------------------------------------------------------------------------------
cplot_label = ['Leaf area index - SUR_LAI'       ,
               'vegetation_area_fraction - PLCOV',
               'Runoff amount - RUNOFF'          ,
               'Biomass amount - SUR_BIOMASS'    ]

# Names of y axis
#------------------------------------------------------------------------------ 
caxis_label   = ['LAI, m2 / m2'  ,
                 'PLCOV'         ,            
                 'RUNOFF, kg m-2',
                 'BIOMASS'       ]

# Common part of COSMO
fn_cosmo      = '_ts_mean_1999_2015.csv'  


# Settings parameteres
#------------------------------------------------------------------------------


research_area = ['parc']#, 'linden', 'lindenberg']                               

# Get LAI data - COPERNICUS and MODIS
if lget_data == True:
    lai_coper = gld.get_lai_data('COPERNICUS', research_area)
    lai_mod   = gld.get_lai_data('MODIS'     , research_area)
else:
    print('The data from NetCDF are not created')


# Plot LAI data - COPERNICUS and MODIS
if lplot_basemap == True:
    plot_lai_coper = lcop.plot_netcdf_cop(sat_path, lGIF = False)        
    plot_lai_coper = lmod.plot_netcdf_mod(sat_path, lGIF = False)
else:
    print('The basemap plots are not created')    

# Get LAI data - COSMO model results
if lcosmo_data == True:
    lai_cosmo = lac.get_cosmo(cos_path       ,                                 # Cosmo path
                              fn_cosmo       ,                                 # Cosmo file name
                              clm_param      ,                                 # COSMO parameteres
                              '2011-01-01'   ,                                 # start calculations
                              '2011-12-31 23').apply(lambda x: x.rolling(7).mean().dropna())                               #   end calculations   
else:
    print('The COSMO data are not loaded')
    
# Plot LAI data - COSMO plot separate
if lplot_cosmo == True:
    plot_lai_cosmo = lac.plot_cosmo(lai_cosmo      , 
                                    clm_param      ,
                                    cplot_label    ,
                                    caxis_label    ,
                                    c_exit         ,
                                    '2011-01-01'   ,
                                    '2011-12-31 23')


if lplot4line == True:
    # Visualization of the results (4 LAI on the one plot)
    #==============================================================================
    # Create time period
    time_start    = pd.to_datetime(['2010-12-30'   ])
    time_stop     = pd.to_datetime(['2011-12-31 23'])
    daily_period  = pd.date_range(time_start[0], time_stop[0], freq = 'D')     # dayly timesteps 
    
    # Get data for plots:
    lai_coper = lai_coper['lai_parc'][daily_period]
    lai_mod   = lai_mod['lai_parc'][daily_period]
                
    # Modis data correction based on interpolation of previous and next value
    for i in range(len(lai_mod)):
        year = 2011
        if lai_mod.index[i]   == pd.to_datetime([f'{year}-01-25']):
            lai_mod[i] = 0.46      
        elif lai_mod.index[i] == pd.to_datetime([f'{year}-05-09']):
            lai_mod[i] = 2.53
        elif lai_mod.index[i] == pd.to_datetime([f'{year}-06-10']):
            lai_mod[i] = 2.660
        elif lai_mod.index[i] == pd.to_datetime([f'{year}-06-18']):        
            lai_mod[i] = 2.661
        elif lai_mod.index[i] == pd.to_datetime([f'{year}-07-04']):        
            lai_mod[i] = 2.56     
        elif lai_mod.index[i] == pd.to_datetime([f'{year}-07-12']):        
            lai_mod[i] = 2.50
        elif lai_mod.index[i] == pd.to_datetime([f'{year}-07-20']):        
            lai_mod[i] = 2.43   
        elif lai_mod.index[i] == pd.to_datetime([f'{year}-08-05']):        
            lai_mod[i] = 2.34              
    
    # COSMO data
    lai_cclmv45    = lai_cosmo['SUR_LAI'][daily_period]
    lai_ref        = lai_cosmo['LAI'][daily_period]
    
    # Get plot
    fig = plt.figure(figsize = (14,10))
    ax  = fig.add_subplot(111)
    # Get plots
    ax.scatter(lai_coper.index, lai_coper   , label = 'COPERNICUS', color = 'blue' )
    ax.scatter(lai_mod.index  , lai_mod     , label = 'MODIS'     , color = 'green')
    ax.plot(lai_cclmv45.index , lai_cclmv45 , label = 'CCLMv4.5'  , color = 'brown')
    ax.plot(lai_ref.index     , lai_ref     , label = 'CCLMref'   , color = 'red'  )
    # Get plot settings
    vis.lplots_stomata2(ax, 'Leaf Area Index, m2 m-2', 'upper right' , 0.0, 3.01, 0.25)  
    # Create output plot name    
    output_name = 'LAI.png'
    # Save plot            
    plt.savefig(c_exit + output_name, format = 'png', dpi = 300) 
    # Clean memory
    plt.close(fig)        
    plt.gcf().clear()
#==============================================================================
