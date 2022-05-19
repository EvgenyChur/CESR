# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
"""
The script for creation of heatmap based on COSMO-CLM data

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


def annual_rstom(df, clm_name, time_start, time_stop, t_start, t_stop):
    df_list  = []
    for t_index in range(len(time_start)):
        hourly_period = pd.date_range(time_start[t_index],
                                      time_stop[t_index] , 
                                      freq = 'H'      )                    # hourly timesteps
            
        # General time period for COSMO, FLUXNET and EURONET data
        res_period   = [x for x in hourly_period if x.hour >= t_start and x.hour <= t_stop]        
            
        # get data
        cclm_df  = csm_data.get_timeseries(df, clm_name, res_period, 'D')

        cclm_df['month'] = cclm_df.index
        cclm_df['month'] = cclm_df['month'].dt.month
        
        # Add to list
        df_list.append(cclm_df)
    # Get dataframe with all data
    df_ref  = pd.concat(df_list , axis = 0)
    return df_ref

#------------------------------------------------------------------------------
# Start program
#------------------------------------------------------------------------------          
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

if mode == 'annual':   
    # Define summer months
    time_start = pd.to_datetime(['2010-05-01', '2011-05-01', 
                                 '2012-05-01', '2013-05-01',
                                 '2014-05-01', '2015-05-01'])
        
    time_stop  = pd.to_datetime(['2010-08-31 23', '2011-08-31 23',
                                 '2012-08-31 23', '2013-08-31 23', 
                                 '2014-08-31 23', '2015-08-31 23'])  
else:
    # Define summer months
    time_start = pd.to_datetime(['2013-05-01'])       
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
cclm_ref  = annual_rstom(df_cclm_ref , clm_name, time_start, time_stop, t_start, t_stop)

mcclm_ref  = cclm_ref.groupby(['month']).agg({'AEVAP_S' : 'mean',
                                              'ALHFL_S' : 'mean',
                                              'ASHFL_S' : 'mean',
                                              'RSTOM'   : 'mean',
                                              'ZVERBO'  : 'mean',
                                              'T_2M'    : 'mean',
                                              #'TMAX_2M' : 'mean',
                                              #'TMIN_2M' : 'mean',
                                              'W_SO'    : 'mean',
                                              'TOT_PREC': 'mean'})
# Calculate correlation
df_corr = mcclm_ref.corr()
 
 # Create boxplots
fig = plt.figure(figsize = (14, 8))
ax  = fig.add_subplot(111)
ax  = sns.heatmap(df_corr, annot = True)
 
ax.figure.savefig(exit_boxplot + f'{domain}_{actual_step}_heatmap.png', dpi = 300)


#cclm_v35  = annual_rstom(df_cclm_v35 , clm_name, time_start, time_stop, t_start, t_stop)
#cclm_v45  = annual_rstom(df_cclm_v45 , clm_name, time_start, time_stop, t_start, t_stop)
#cclm_v45e = annual_rstom(df_cclm_v45e, clm_name, time_start, time_stop, t_start, t_stop)

#mcclm_v35  = cclm_v35.resample('1M').mean()
#mcclm_v45  = cclm_v45.resample('1M').mean()
#mcclm_v45e = cclm_v45e.resample('1M').mean()




