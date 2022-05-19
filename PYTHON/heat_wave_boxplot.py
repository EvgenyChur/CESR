# -*- coding: utf-8 -*-
"""
Created on Thu Feb 24 09:08:07 2022

@author: churiulin
"""

# Import standart liblaries 
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Import personal libraries                                                                                                                                                       
import cosmo_data       as csm_data                                            




# Calculations of the standardized anomaly
#------------------------------------------------------------------------------
def extreme_data(df):
    ''' 
    Parameters
    ----------
    df : Series
        The series with informataion about temperatures.
    Returns
    -------
    ts_data : Series
        Series with standardized anomaly.
    '''
    t2m_mean = df.mean()
    t2m_std  = df.std()    
    # Create a nan timeseries
    ts_data    = pd.Series(np.nan, index = df.index)
    # Calculate values 
    for i in range(len(df)):
        ts_data[i] = (df[i] - t2m_mean) / t2m_std     
    return ts_data
#------------------------------------------------------------------------------

# Additional parameters for plots
#------------------------------------------------------------------------------
def plot_set(cx, rot, size):
    '''
    Parameters
    ----------
    cx : figure
        Figure for plot.
    rot : int
        Rotation of x-axis.
    size : int
        Size of text.
    Returns
    -------
    None.
    '''
    # Settings for axis
    for label in cx.xaxis.get_ticklabels():
        label.set_color('black')
        label.set_rotation(rot)
        label.set_fontsize(size)
    for label in cx.yaxis.get_ticklabels():
        label.set_color('black')
        label.set_fontsize(size)
    # Add grid and save plot    
    cx.grid(True)   



#==============================================================================
# Get COSMO-CLM parameteres
#==============================================================================
mf_com    = 'C:/Users/Churiulin/Desktop/COSMO_RESULTS/'                        # Main path for project  
cosmo_old = mf_com + 'COSMO/OLD_data/PARC/CTR/'                                # Path
data_exit = mf_com + '/ANALYSIS/PLOTS/PARC/'                                   # Main path for project results

# Names of parameters for COSMO data
clm_name = ['T_2M', 'T_S', 'TMAX_2M', 'TMIN_2M', 'TOT_PREC', 'W_SO']

# Name of COSMO data
fn_cosmo = '_ts_mean_1999_2015.csv'

# Define the time intervals
years = ['1999', '2000', '2001', '2002', '2003', 
         '2004', '2005', '2006', '2007', '2008',
         '2009', '2010', '2011', '2012', '2013', 
         '2014', '2015', '2016']

# Define step for filtering data
t_start = 0                                                                    # parameter --> from  
t_stop  = 23                                                                   # parameter --> to 

# Define the step for rolling mean   
days    = 7

#==============================================================================
# Start calculations
#==============================================================================

# Get COSMO-CLM data
df_cclm_ref  = csm_data.cosmo_data(cosmo_old , fn_cosmo, clm_name)

# Create filter only for summer moths
time_start = []
time_stop  = []
for year in years:
    time_start.append(f'{year}-07-01'  )
    time_stop.append(f'{year}-08-31 23')
time_start = pd.to_datetime(pd.Series(time_start))
time_stop  = pd.to_datetime(pd.Series(time_stop))

# Create dataframe for filter data
year_list = []
for index in range(len(time_start)):
    # Create periods 
    hourly_period = pd.date_range(time_start[index], 
                                  time_stop[index], 
                                  freq = 'H')                                  # hourly timesteps      
    # Filter for hours
    res_period   = [x for x in hourly_period if x.hour >= t_start and x.hour <= t_stop]     
    df_cosmo     = csm_data.get_timeseries(df_cclm_ref, clm_name, res_period, 'D')
    year_list.append(df_cosmo)

# Create a DataFrame for analysis
cclm_ref = pd.concat(year_list, axis = 0)    
cclm_ref['extrem_t2m'] = extreme_data(cclm_ref['T_2M'])
cclm_ref['year']       = cclm_ref.index
cclm_ref['year']       = cclm_ref['year'].dt.year


#==============================================================================
# Create boxplots for the Standardized anomaly
#==============================================================================

fig = plt.figure(figsize = (14, 8))
ax  = fig.add_subplot(111)
# Creat boxplot 
ax = sns.boxplot(x    = 'year'      ,
                 y    = 'extrem_t2m', 
                 data = cclm_ref    ) 
      
# Settings for boxplot
ax.set_yticks(np.arange(-3.0, 4.01, 1.0))                                       
    
# Create labels
ax.set_ylabel('Standardized anomaly (s.d.)', 
              color    = 'black', 
              fontsize = 16     ,
              labelpad = 20     )
ax.set_xlabel('Years'           ,
              color    = 'black', 
              fontsize = 16     , 
              labelpad = 20     )
# Additional settings for axis (rotation and size)
plot_set(ax, 15, 16)
ax.figure.savefig(data_exit + '1. Heat_wave.png', dpi = 300)

#==============================================================================
# Create barplot for total precipitation
#==============================================================================
from numpy import sum

fig1 = plt.figure(figsize = (14, 8))
bx   = fig1.add_subplot(111)
# Create barplot
bx = sns.barplot(x         = 'year'    , 
                 y         = 'TOT_PREC',
                 data      = cclm_ref  ,
                 estimator = sum       ,
                 ci        = None      ) 
# Settings for boxplot
bx.set_yticks(np.arange(0.0, 275.01, 25.0))                                       
# Create labels
bx.set_ylabel('Total precipitation, mm', 
              color    = 'black',
              fontsize = 16     ,
              labelpad = 20     )
bx.set_xlabel('Years'           ,
              color    = 'black',
              fontsize = 16     ,
              labelpad = 20     )
# Additional settings for axis (rotation and size)
plot_set(bx, 15, 16)
bx.figure.savefig(data_exit + '2. Precipitations.png', dpi = 300)


#==============================================================================
# Create barplot for w_so
#==============================================================================
from numpy import median

fig2 = plt.figure(figsize = (14, 8))
cx   = fig2.add_subplot(111)
# Create barplot
bx = sns.barplot(x         = 'year'    , 
                 y         = 'W_SO',
                 data      = cclm_ref  ,
                 estimator = median       ) 
# Settings for boxplot
bx.set_yticks(np.arange(0.0, 0.06, 0.01))                                       
# Create labels
bx.set_ylabel('Soil moisture, m H2O', 
              color    = 'black',
              fontsize = 16     ,
              labelpad = 20     )
bx.set_xlabel('Years'           ,
              color    = 'black',
              fontsize = 16     ,
              labelpad = 20     )
# Additional settings for axis (rotation and size)
plot_set(cx, 15, 16)
cx.figure.savefig(data_exit + '3. Soil_moisture.png', dpi = 300)   
    