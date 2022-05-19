# -*- coding: utf-8 -*-
"""
The script for creation of boxplot based on the TRY database and T2m categories

Plant trait data in the TRY Database are public under the Creative Commons 
Attribution License ‘CC BY’ (https://creativecommons.org/licenses/by/4.0).

Data is available by reguests via the TRY website www.try‐db.org/TryWeb/Prop0.php

More information about data is available:
Kattge, J, Bönisch, G, Díaz, S, et al. TRY plant trait database – enhanced
coverage and open access. Glob Change Biol. 2020; 26: 119– 188. 
https://doi.org/10.1111/gcb.14904


The module contains several subroutines:
    1. annual_rstom    - Get annual mean data for stomatal resistance



Autors of scripts: Evgenii Churiulin - Center for Enviromental System Research (CESR) 

                                                   
Current Code Owner: CESR, Evgenii Churiulin
phone:  +49  561 804-6142
fax:    +49  561 804-6116
email:  evgenychur@uni-kassel.de


History:
Version    Date       Name
---------- ---------- ----                                                   
    1.1    2022-02-15 Evgenii Churiulin, Center for Enviromental System Research (CESR)
           Initial release
                 
"""


# Import standart liblaries 
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Import personal libraries
import sys
sys.path.append('C:/Users/Churiulin/Python/scripts/CESR_project')
import vis_module       as vis                                                                                 
import cosmo_data       as csm_data        
import get_new_data     as gnd



def annual_rstom(df, clm_name, time_start, time_stop, t_start, t_stop):
    df_list  = []
    for t_index in range(len(time_start)):
        hourly_period = pd.date_range(time_start[t_index],
                                      time_stop[t_index] , 
                                      freq = 'H'         )                    # hourly timesteps        
        # General time period for COSMO, FLUXNET and EURONET data
        res_period   = [x for x in hourly_period if x.hour >= t_start and x.hour <= t_stop]              
        # get data
        cclm_df  = csm_data.get_timeseries(df, clm_name, res_period, 'D')
        # Add to list
        df_list.append(cclm_df[['T_2M','RSTOM']].reset_index())
    # Get dataframe with all data
    df_ref  = pd.concat(df_list , axis = 1)
    # Rename columns
    df_ref.columns = ['Date' , 'T_2M_2010', 'RSTOM_2010', 
                      'Date2', 'T_2M_2011', 'RSTOM_2011',
                      'Date3', 'T_2M_2012', 'RSTOM_2012', 
                      'Date4', 'T_2M_2013', 'RSTOM_2013',
                      'Date5', 'T_2M_2014', 'RSTOM_2014', 
                      'Date6', 'T_2M_2015', 'RSTOM_2015']    
        
    # Delete dates and set new index
    df_ref = (df_ref.drop(['Date2', 'Date3',
                           'Date4', 'Date5', 'Date6'], axis = 1)
                    .set_index(df_ref['Date'])
                    .drop(['Date'], axis = 1)
             )
    # Get new annual mean values
    df_ref['RSTOM'] = (df_ref['RSTOM_2010'] + df_ref['RSTOM_2011'] + 
                       df_ref['RSTOM_2012'] + df_ref['RSTOM_2013'] +
                       df_ref['RSTOM_2014'] + df_ref['RSTOM_2015'] ) / 6.0
        
    df_ref['T_2M'] = (df_ref['T_2M_2010'] + df_ref['T_2M_2011'] + 
                      df_ref['T_2M_2012'] + df_ref['T_2M_2013'] +
                      df_ref['T_2M_2014'] + df_ref['T_2M_2015'] ) / 6.0
            
    return df_ref



# RSTOM for 1 year
def year_rstom(df, clm_name, time_start, time_stop, t_start, t_stop):
    df_list  = []
    for t_index in range(len(time_start)):
        hourly_period = pd.date_range(time_start[t_index],
                                      time_stop[t_index] , 
                                      freq = 'H'         )                    # hourly timesteps        
        # General time period for COSMO, FLUXNET and EURONET data
        res_period   = [x for x in hourly_period if x.hour >= t_start and x.hour <= t_stop]              
        # get data
        cclm_df  = csm_data.get_timeseries(df, clm_name, res_period, 'D')
        # Add to list
        df_list.append(cclm_df[['T_2M','RSTOM']].reset_index())
    # Get dataframe with all data
    df_ref  = pd.concat(df_list , axis = 1)
    # Rename columns
    df_ref.columns = ['Date', 'T_2M_2013', 'RSTOM_2013']    
        
    # Delete dates and set new index
    df_ref = df_ref.set_index(df_ref['Date']).drop(['Date'], axis = 1)
    # Get new annual mean values
    df_ref['RSTOM'] = df_ref['RSTOM_2013']    
    df_ref['T_2M'] = df_ref['T_2M_2013']
       
    return df_ref





# Choose your region for analysis (Possible options: 1 - PARC domain
#                                                    2 - LINDEN domain
#                                                    3 - LINDENBERG domain )

input_region = '1'


while True: 
    if input_region == '1':
        domain    = 'PARC'                                                     # Name of folder for COSMO data
        fn_region = 'parc'                                                     # Name of file for GLEAM data
        station_plot = 'Rollesbroich'
        break
    
    elif input_region == '2':
        domain         = 'LINDEN'
        fn_region      = 'linden'
        station_plot   = 'Linden'
        sf_in_situ     = 'IN-SITU/' + domain + '/'
        fn_in_situ     = 'EC4.csv'
        break
    
    elif input_region == '3':
        domain         = 'LINDENBERG'
        fn_region      = 'lindenberg'
        station_plot   = 'Lindenberg'
        sf_in_situ     = 'IN-SITU/' + domain + '/'
        fn_in_situ     = 'Lindenberg.csv'
        break
    else:
        print('Error: Incorrect format of region')

#------------------------------------------------------------------------------
# Setting paths for COSMO, EURONET, FLUXNET, GLEAM, HYRAS, E-OBS data
#------------------------------------------------------------------------------
mf_com    = 'C:/Users/Churiulin/Desktop/COSMO_RESULTS/'                        # Main path for project  
data_exit    = mf_com + '/ANALYSIS/PLOTS/' + domain + '/'                         # Main path for project results



exit_boxplot = mf_com + '/ANALYSIS/PLOTS/BOXPLOTS/'

sf_cclm_ref  = mf_com + 'COSMO/' + domain + '/CTR/'                            # Path for CCLM_ref  data                
sf_cclm_v35  = mf_com + 'COSMO/' + domain + '/v3.5/'                           # Path for CCLMv3.5  data                
sf_cclm_v45  = mf_com + 'COSMO/' + domain + '/v4.5/'                           # Path for CCLMv4.5  data               
sf_cclm_v45e = mf_com + 'COSMO/' + domain + '/v4.5e/'                          # Path for CCLMv4.5e data   


# Names of parameters for COSMO data
clm_name = ['AEVAP_S', 'ALHFL_S' , 'ASHFL_S', 'RSTOM'   ,
            'ZVERBO' , 'T_2M'    , 'T_S'    , 'TMAX_2M' ,
            'TMIN_2M', 'ZTRALEAV', 'W_SO'   , 'TOT_PREC']

# Name of COSMO data
fn_cosmo = '_ts_mean_1999_2015.csv'

# Input station
input_station = 'RuR'

# Timestep for data (need for pd.resample)
time_step     = '1D'


plant_species = ['Lolium perenne'       ,
                 'Poa pratensis'        ,
                 'Arrhenatherum elatius',
                 'Filipendula ulmaris'  ,
                 'Festuca rubra'        ]

#mode = 'annual'
mode = 'year'  
actual_step = 'day'
#actual_step = 'night'

#------------------------------------------------------------------------------
# Get initial data
#------------------------------------------------------------------------------
df_cclm_ref  = csm_data.cosmo_data(sf_cclm_ref , fn_cosmo, clm_name)           # Get COSMO_ref  data
df_cclm_v35  = csm_data.cosmo_data(sf_cclm_v35 , fn_cosmo, clm_name)           # Get COSMOv3.5  data
df_cclm_v45  = csm_data.cosmo_data(sf_cclm_v45 , fn_cosmo, clm_name)           # Get COSMOv4.5  data
df_cclm_v45e = csm_data.cosmo_data(sf_cclm_v45e, fn_cosmo, clm_name)           # Get COSMOv4.5e data

in_situ  = gnd.get_new_TRYdata(mf_com, plant_species, 1, ldf_save = False)     # stomatal resistance
#in_situ  = gnd.get_new_TRYdata(mf_com, plant_species, 2, ldf_save = False)     # Vcmax

# Add column with experiment name
in_situ['Name']  = 'in-situ' 
   

# Define summer months
if mode == 'annual':
    time_start = pd.to_datetime(['2010-05-01', '2011-05-01', 
                                 '2012-05-01', '2013-05-01',
                                 '2014-05-01', '2015-05-01'])
        
    time_stop  = pd.to_datetime(['2010-08-31 23', '2011-08-31 23', 
                                 '2012-08-31 23', '2013-08-31 23',
                                 '2014-08-31 23', '2015-08-31 23']) 
else:
    time_start = pd.to_datetime(['2013-07-01'])       
    time_stop  = pd.to_datetime(['2013-08-31 23']) 
   
    
    
if actual_step == 'day':
    # Day-time values
    t_start = 12                                                               
    t_stop  = 13     
else:
    # Night-time values
    t_start = 18                                                          
    t_stop  = 23 
    
# Get stomatal resistance data 
if mode == 'annual':
    cclm_ref  = annual_rstom(df_cclm_ref , clm_name, time_start, time_stop, t_start, t_stop)
    cclm_v35  = annual_rstom(df_cclm_v35 , clm_name, time_start, time_stop, t_start, t_stop)
    cclm_v45  = annual_rstom(df_cclm_v45 , clm_name, time_start, time_stop, t_start, t_stop)
    cclm_v45e = annual_rstom(df_cclm_v45e, clm_name, time_start, time_stop, t_start, t_stop)
else:
    cclm_ref  = year_rstom(df_cclm_ref , clm_name, time_start, time_stop, t_start, t_stop)
    cclm_v35  = year_rstom(df_cclm_v35 , clm_name, time_start, time_stop, t_start, t_stop)
    cclm_v45  = year_rstom(df_cclm_v45 , clm_name, time_start, time_stop, t_start, t_stop)
    cclm_v45e = year_rstom(df_cclm_v45e, clm_name, time_start, time_stop, t_start, t_stop)    
   
# Add column with experiment name
cclm_ref['Name']  = 'CCLMref'
cclm_v35['Name']  = 'CCLMv3.5'
cclm_v45['Name']  = 'CCLMv4.5'
cclm_v45e['Name'] = 'CCLMv4.5e'    
    
# Correct vulues
if actual_step == 'night':
    cclm_v35['RSTOM'] = (cclm_v35['RSTOM']
                                 .replace(cclm_v35['RSTOM'][cclm_v35['RSTOM_2014'] < 3000],
                                          5460.0)
                        )
# Correct vulues
#if actual_step == 'day':
#    cclm_v45['RSTOM'] = (cclm_v45['RSTOM']
#                                 .replace(cclm_v45['RSTOM'][cclm_v45['RSTOM'] > 1200],
#                                          1200.0)
#                        ) 
    
# Get data for boxplot
df_data = pd.concat([cclm_ref[['RSTOM', 'T_2M', 'Name']], 
                     cclm_v35[['RSTOM', 'T_2M', 'Name']],
                     cclm_v45[['RSTOM', 'T_2M', 'Name']], 
                    cclm_v45e[['RSTOM', 'T_2M', 'Name']],
                      in_situ[['RSTOM', 'T_2M', 'Name']]], axis = 0)        
    
df_data['t2m_cat'] = pd.cut(df_data['T_2M'], [10, 15, 20, 25, np.inf], 
                            labels = ['10-15',
                                      '15-20',
                                      '20-25',
                                      '>25'  ])   
# Create boxplots
fig = plt.figure(figsize = (12, 7))
ax  = fig.add_subplot(111)
    
# Create boxplot
    
#ax = sns.boxplot(x    = 'Name'  ,
#                 y    = 'RSTOM' , 
#                 data = rstom) 
  
ax = sns.boxplot(x    = 't2m_cat'  ,
                 y    = 'RSTOM' , 
                 hue  = 'Name',
                 data = df_data)      
        
if mode == 'annual':
    if actual_step == 'day':
        # Day-time values
        ax.set_yticks(np.arange(0, 1600.1, 200.0))#was 0.0, 1600.1, 200.0
    else:
        # Night-time values
        ax.set_yticks(np.arange(2000, 12000.1, 2000))
# get data for 1 year
else:
    if actual_step == 'day':    
        # Day-time values
        ax.set_yticks(np.arange(0, 4000.1, 500.0))
    else:
        # Night-time values
        ax.set_yticks(np.arange(2000, 12000.1, 2000))        
    
# Create labels
ax.set_ylabel('RSTOM, s m \u207b\u00B9'    , color = 'black', fontsize = 16, labelpad = 20)
ax.set_xlabel('Temperature ranges, \u2070C', color = 'black', fontsize = 16, labelpad = 20)
# Settings for axis
for label in ax.xaxis.get_ticklabels():
    label.set_color('black')
    label.set_rotation(0)
    label.set_fontsize(16)
for label in ax.yaxis.get_ticklabels():
    label.set_color('black')
    label.set_fontsize(16)
# Add grid and save plot    
ax.grid(True)      
ax.figure.savefig(exit_boxplot + f'{domain}_{actual_step}_T2M.png', dpi = 300)        
  




