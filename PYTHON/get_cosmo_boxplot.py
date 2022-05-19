# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
"""
The script for creation of box and line plots based on COSMO-CLM data

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
import STOMATA          as stm 

# Annual
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
        # Add to list
        df_list.append(cclm_df['RSTOM'].reset_index())
    # Get dataframe with all data
    df_ref  = pd.concat(df_list , axis = 1)
    # Rename columns
    df_ref.columns = ['Date' , 'RSTOM_2010', 'Date2', 'RSTOM_2011',
                      'Date3', 'RSTOM_2012', 'Date4', 'RSTOM_2013',
                      'Date5', 'RSTOM_2014', 'Date6', 'RSTOM_2015']    
        
    # Delete dates and set new index
    df_ref = (df_ref.drop(['Date2', 'Date3',
                           'Date4', 'Date5', 'Date6'], axis = 1)
                    .set_index(df_ref['Date'])
                    .drop(['Date'], axis = 1)
             )
           
    # Get new annual mean values
    df_ref['RSTOM'] = df_ref.mean(axis = 1)
        
    return df_ref

# Year
def year_rstom(df, clm_name, time_start, time_stop, t_start, t_stop):
    df_list  = []
    for t_index in range(len(time_start)):
        hourly_period = pd.date_range(time_start[t_index],
                                      time_stop[t_index] , 
                                      freq = 'H'      )                    # hourly timesteps            
        # General time period for COSMO, FLUXNET and EURONET data
        res_period   = [x for x in hourly_period if x.hour >= t_start and x.hour <= t_stop]               
        # get data
        cclm_df  = csm_data.get_timeseries(df, clm_name, res_period, 'D')
        # Add to list
        df_list.append(cclm_df['RSTOM'].reset_index())
    # Get dataframe with all data
    df_ref  = pd.concat(df_list , axis = 1)
    # Rename columns
    df_ref.columns = ['Date', 'RSTOM_2013']    
    # Delete dates and set new index
    df_ref = df_ref.set_index(df_ref['Date']).drop(['Date'], axis = 1)      
    # Get new annual mean values
    df_ref['RSTOM'] = df_ref.mean(axis = 1)        
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
llplot = True                                                       
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
# Year 
else: 
    # Define summer months
    time_start = pd.to_datetime(['2013-07-01'])      
    time_stop  = pd.to_datetime(['2013-08-31 23'])     
    

    
if actual_step == 'day':
    # Day-time values
    t_start = 12                                                               # for control 12 - 13
    t_stop  = 13     
else:
    # Night-time values
    t_start = 18                                                          
    t_stop  = 23 

#==============================================================================
# Get COSMO-CLM data
#==============================================================================
if mode == 'annual':    
    # Get stomatal resistance data 
    cclm_ref  = annual_rstom(df_cclm_ref , clm_name , 
                             time_start  , time_stop, 
                             t_start     , t_stop   )
    cclm_v35  = annual_rstom(df_cclm_v35 , clm_name , 
                             time_start  , time_stop, 
                             t_start     , t_stop   )
    cclm_v45  = annual_rstom(df_cclm_v45 , clm_name ,
                             time_start  , time_stop,
                             t_start     , t_stop   )
    cclm_v45e = annual_rstom(df_cclm_v45e, clm_name , 
                             time_start  , time_stop,
                             t_start     , t_stop   )
else:
    # Get stomatal resistance data 
    cclm_ref  = year_rstom(df_cclm_ref , clm_name , 
                           time_start  , time_stop,
                           t_start     , t_stop   )
    cclm_v35  = year_rstom(df_cclm_v35 , clm_name ,
                           time_start  , time_stop,
                           t_start     , t_stop   ) 
    cclm_v45  = year_rstom(df_cclm_v45 , clm_name ,
                           time_start  , time_stop,
                           t_start     , t_stop   )
    cclm_v45e = year_rstom(df_cclm_v45e, clm_name , 
                           time_start  , time_stop,
                           t_start     , t_stop   )    
    
    
# Add column with experiment name
cclm_ref['Name']  = 'CCLMref'
cclm_v35['Name']  = 'CCLMv3.5'
cclm_v45['Name']  = 'CCLMv4.5'
cclm_v45e['Name'] = 'CCLMv4.5e'    
    

if mode == 'annual':  
    # Correct vulues
    if actual_step == 'night':
        cclm_v35['RSTOM'] = (cclm_v35['RSTOM']
                                     .replace(cclm_v35['RSTOM'][cclm_v35['RSTOM_2014'] < 3000],
                                              5460.0)
                            )  
    
    
    
    
# Get data for boxplot
rstom = pd.concat([cclm_ref[['RSTOM' , 'Name']], 
                   cclm_v35[['RSTOM' , 'Name']],
                   cclm_v45[['RSTOM' , 'Name']], 
                   cclm_v45e[['RSTOM', 'Name']]], axis = 0)
    
rstom['month'] = rstom.index
#rstom['month'] = rstom['month'].dt.month
rstom['month'] = rstom['month'].dt.strftime('%b')
      
#==============================================================================
# Create boxplots for stomatal resitance
#==============================================================================

fig = plt.figure(figsize = (12, 7))
ax  = fig.add_subplot(111)

#ax = sns.boxplot(x    = 'Name'  ,
#                 y    = 'RSTOM' , 
#                 data = rstom) 
   
ax = sns.boxplot(x    = 'month'  ,
                 y    = 'RSTOM' , 
                 hue  = 'Name',
                 data = rstom) 
 
if mode == 'annual':         
    if actual_step == 'day':
        ax.set_yticks(np.arange(0, 600.1, 200))                               # Day-time values
    else:
        ax.set_yticks(np.arange(2000, 12000.1, 2000))                          # Night-time values
# Year
else:         
    if actual_step == 'day':
        ax.set_yticks(np.arange(0, 4000.1, 500))                               # Day-time values
    else:
        ax.set_yticks(np.arange(2000, 12000.1, 2000))                          # Night-time values        
        
    
# Create labels
ax.set_ylabel('RSTOM, s m \u207b\u00B9', color = 'black', fontsize = 16, labelpad = 20)
ax.set_xlabel('Months'                 , color = 'black', fontsize = 16, labelpad = 20)
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
    
ax.figure.savefig(exit_boxplot + f'{domain}_{actual_step}.png', dpi = 300)





# Plot line figure
fig1 = plt.figure(figsize = (12, 7))
bx  = fig1.add_subplot(111)
    
plot4lines = vis.plots4(bx, cclm_ref['RSTOM'].rolling(1).mean() ,
                            cclm_v35['RSTOM'].rolling(1).mean() ,
                            cclm_v45['RSTOM'].rolling(1).mean() , 
                            cclm_v45e['RSTOM'].rolling(1).mean(),
                            'CCLMref', 'CCLMv3.5',
                            'CCLMv4.5', 'CCLMv4.5e')    

if mode == 'annual':        
    if actual_step == 'day':
        # Day-time values
        plot4lines = vis.lplots_stomata2(bx, 'RSTOM, s m \u207b\u00B9', 'upper left',
                                             0, 5000.1, 500)#0, 1600.1, 200.0)       
    else:
        # Night-time values    
        plot4lines = vis.lplots_stomata2(bx, 'RSTOM, s m \u207b\u00B9', 'upper left',
                                             2000, 8000.1, 1000.0)
else:
    if actual_step == 'day':
        plot4lines = vis.lplots_stomata2(bx, 'RSTOM, s m \u207b\u00B9', 'upper left',
                                                      0, 4000.1, 500.0)      
    else:
        # Night-time values    
        plot4lines = vis.lplots_stomata2(bx, 'RSTOM, s m \u207b\u00B9', 'upper left',
                                             2000, 12000.1, 2000.0)        
            
            
            
plt.savefig(exit_boxplot + f'{domain}_RSTOM_{actual_step}.png', format = 'png', dpi = 300)




#==============================================================================
# Create a line plot for comparison with in-situ data
#==============================================================================
if llplot == True:
    #--------------------------------------------------------------------------
    # Get in-situ data North America
    #--------------------------------------------------------------------------   
    if mode == 'annual':
        dates = [pd.Timestamp("2010-06-20"),
                 pd.Timestamp("2010-06-22"),
                 pd.Timestamp("2010-06-25"),         
                 pd.Timestamp("2010-07-04"),
                 pd.Timestamp("2010-07-13"),
                 pd.Timestamp("2010-07-20"),
                 pd.Timestamp("2010-07-29"),
                 pd.Timestamp("2010-08-09"),
                 pd.Timestamp("2010-08-18"),
                 pd.Timestamp("2010-08-27")]
    # Year
    else:
        dates = [pd.Timestamp("2013-06-20"),
                 pd.Timestamp("2013-06-22"),
                 pd.Timestamp("2013-06-25"),         
                 pd.Timestamp("2013-07-04"),
                 pd.Timestamp("2013-07-13"),
                 pd.Timestamp("2013-07-20"),
                 pd.Timestamp("2013-07-29"),
                 pd.Timestamp("2013-08-09"),
                 pd.Timestamp("2013-08-18"),
                 pd.Timestamp("2013-08-27")]
        
    # Create timeseries based on measurements
    in_situ = pd.Series([140 , 295, 172, 230, 196 ,
                         410 , 210, 280, 160, 221 ], dates)
       
    
    # Get data for statistic analisys    
    df_stom = pd.concat([cclm_ref['RSTOM'], 
                         cclm_v35['RSTOM'], 
                         cclm_v45['RSTOM'], 
                        cclm_v45e['RSTOM'], 
                          in_situ        ], axis = 1).dropna()
    
    #df_stom2 = df_stom.dropna()
    
    
    df_stom.columns = ['COSMO_CTR' ,'COSMO_v3.5' , 
                       'COSMO_v4.5','COSMO_v4.5e',
                       'OBS']
     
    stat_stom_ctr  = stm.statistic(df_stom['OBS'], df_stom['COSMO_CTR'  ], 'CCLM_ref')
    stat_stom_v35  = stm.statistic(df_stom['OBS'], df_stom['COSMO_v3.5' ], 'CCLMv3.5' )
    stat_stom_v45  = stm.statistic(df_stom['OBS'], df_stom['COSMO_v4.5' ], 'CCLMv4.5' )
    stat_stom_v45e = stm.statistic(df_stom['OBS'], df_stom['COSMO_v4.5e'], 'CCLMv4.5e')
    
    df_com_stat = pd.concat([stat_stom_ctr, stat_stom_v35 ,
                             stat_stom_v45, stat_stom_v45e], axis = 0)
    df_com_stat.to_excel(exit_boxplot + f'RSTOM_{fn_region}.xlsx', float_format='%.3f')   
    
    
    # Plot line figure
    fig1 = plt.figure(figsize = (12, 7))
    bx  = fig1.add_subplot(111)
    
    
    #plot5lines = vis.plots5_stomata(bx, cclm_ref['RSTOM'].rolling(7).mean(),
    #                                    cclm_v35['RSTOM'].rolling(7).mean(),
    #                                    cclm_v45['RSTOM'].rolling(7).mean(), 
    #                                    cclm_v45e['RSTOM'].rolling(7).mean(),
    #                                    in_situ, 'CCLMref', 'CCLMv3.5',
    #                                    'CCLMv4.5', 'CCLMv4.5e', 'OBS')    
    
    
    plot5lines = vis.plots4(bx, cclm_ref['RSTOM'].rolling(7).mean(),
                                cclm_v35['RSTOM'].rolling(7).mean(),
                                cclm_v45['RSTOM'].rolling(7).mean(), 
                                cclm_v45e['RSTOM'].rolling(7).mean(),
                                'CCLMref' , 'CCLMv3.5' ,
                                'CCLMv4.5', 'CCLMv4.5e')  
    
    if mode == 'annual':        
        if actual_step == 'day':
            # Day-time values
            plot5lines = vis.lplots_stomata2(bx, 'RSTOM, s m \u207b\u00B9', 'upper left',
                                             0, 5000.1, 500)#0, 1600.1, 200.0)       
        else:
            # Night-time values    
            plot5lines = vis.lplots_stomata2(bx, 'RSTOM, s m \u207b\u00B9', 'upper left',
                                             2000, 8000.1, 1000.0)
    else:
        if actual_step == 'day':
            plot5lines = vis.lplots_stomata2(bx, 'RSTOM, s m \u207b\u00B9', 'upper left',
                                                      0, 3500.1, 500.0)      
        else:
            # Night-time values    
            plot5lines = vis.lplots_stomata2(bx, 'RSTOM, s m \u207b\u00B9', 'upper left',
                                             2000, 8000.1, 1000.0)        
            
            
            
    plt.savefig(exit_boxplot + f'{domain}_RSTOM2_{actual_step}.png', format = 'png', dpi = 300)
    