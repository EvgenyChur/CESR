# -*- coding: utf-8 -*-
"""
Created on Mon Jul 26 14:58:58 2021

@author: churiulin
"""

# Import standart liblaries 

import pandas as pd
import matplotlib.pyplot as plt

# Import personal libraries
                                                            
import cosmo_data        as csm_data                                                                                    
import vis_module        as vsp                                                 

# Improt methods for statistical analysis
from sklearn.metrics import mean_absolute_error
from sklearn.metrics import mean_squared_error                                           


def stomata(data, years):
    # Create time period
    t_1 = pd.to_datetime([f'{years}-06-01 13'])
    t_2 = pd.to_datetime([f'{years}-08-31 13'])    
    
    x = []
    
    for tr in range(len(t_1)):
        period = pd.date_range(t_1[tr], t_2[tr], freq = '1D')     
        stom   = data.iloc[:,3][period]
        stom   = stom.reset_index()
        x.append(stom)
        
    df = pd.concat(x, axis = 1) 
              
    df.columns = [f'{years}',f'R{years}']
    #df = df.drop(['2012','2013', '2014', '2015'], axis = 1)  
    #df = df.set_index('2011')     
    df = df.set_index(f'{years}') 
    df_mean = df.mean(axis = 1)
    return df_mean



def statistic(x, y, column_name):
    mean_obs = x.mean()
    std_obs  = x.std()
            
    mean_mod = y.mean()
    std_mod  = y.std() 
             
    mae      = mean_absolute_error(x, y)
    rmse     = mean_squared_error(x , y)
    corr     = x.corr(y)
        
    df_stat = pd.DataFrame({'Parameter' : [column_name],
                              'Mean OBS': mean_obs,
                              'Mean MOD': mean_mod,
                              'STD OBS' : std_obs ,
                              'STD_MOD' : std_mod ,
                              'MAE'     : mae     ,
                              'RMSE'    : rmse    ,
                              'CORR'    : corr    })
        
    df_stat.set_index('Parameter', inplace = True)
    return df_stat
#------------------------------------------------------------------------------
# Main programm
#------------------------------------------------------------------------------

if __name__ == '__main__':
    # General path to the project folders                            

    # Choose your region for analysis (Possible options: 1 - PARC domain
    #                                                    2 - LINDEN domain
    #                                                    3 - LINDENBERG domain )
    
    input_region = '1'
    
    
    while True: 
        if input_region == '1':
            domain    = 'PARC'                                                     # Name of folder for COSMO data
            fn_region = 'parc'                                                     # Name of file for GLEAM data
            station_plot = 'Rollesbroich'   
            date_ind = 'Mean RSTOM from 2011 to 2015'         
            l_p      = 'upper left'
            nst      = 'Parc'  
    
            break
        
        elif input_region == '2':
            domain         = 'LINDEN'
            fn_region      = 'linden'
            station_plot   = 'Linden'
            sf_in_situ     = 'IN-SITU/' + domain + '/'
            fn_in_situ     = 'EC4.csv'        
            date_ind       = 'Mean RSTOM from 2011 to 2015'         
            l_p            = 'upper left'
            nst            = 'Linden'         
            break
        
        elif input_region == '3':
            domain         = 'LINDENBERG'
            fn_region      = 'lindenberg'
            station_plot   = 'Lindenberg'
            sf_in_situ     = 'IN-SITU/' + domain + '/'
            fn_in_situ     = 'Lindenberg.csv'       
            date_ind       = 'Mean RSTOM from 2011 to 2015'         
            l_p            = 'upper left'
            nst            = 'Lindenberg'  
            break
        else:
            print('Error: Incorrect format of region')
        
        
    # Path to the main project folder
    mf_com    = 'C:/Users/Churiulin/Desktop/COSMO_RESULTS/'
    data_exit = mf_com + '/ANALYSIS/PLOTS/' + domain + '/'
    
    #------------------------------------------------------------------------------
    # Setting for COSMO data
    #------------------------------------------------------------------------------
    
    # Path to COSMO subfolders
    sf_parc01_ctr     = mf_com + 'COSMO/' + domain + '/CTR/'                                         
    sf_parc_v35       = mf_com + 'COSMO/' + domain + '/v3.5/'                                         
    sf_parc_v45       = mf_com + 'COSMO/' + domain + '/v4.5/'                                        
    sf_parc_v45_evap  = mf_com + 'COSMO/' + domain + '/v4.5e/'                                    
    
    # Names of parameters for COSMO data
    clm_name = ['AEVAP_S', 'ALHFL_S' , 'ASHFL_S', 'RSTOM'   ,
                'ZVERBO' , 'T_2M'    , 'T_S'    , 'TMAX_2M' ,
                'TMIN_2M', 'ZTRALEAV']
    
    name_1 = ['Stomata resistance - RSTOM']    
    
    name_2 = ['RSTOM, s m  \u207b\u00B9']    
        
    
    # Name of COSMO data
    fn_cosmo = '_ts_mean_1999_2015.csv'
    
    exit_path = 'C:/Users/Churiulin/Desktop/COSMO_RESULTS/ANALYSIS/PLOTS/'
    
    
    # USER data with user timesteps
    #             RSTOM  
    y_min  = [     0.0 ]
    y_max  = [  3500.1 ]
    y_step = [   250.0 ]  
    
    year_list = [2010, 2011, 2012, 2013, 2014, 2015]  
    #------------------------------------------------------------------------------
    # Get initial COSMO data
    #------------------------------------------------------------------------------
    for year in year_list:
        
        dates = [pd.Timestamp(f"{year}-06-20 13"),
                 pd.Timestamp(f"{year}-06-22 13"),
                 pd.Timestamp(f"{year}-06-25 13"),         
                 pd.Timestamp(f"{year}-07-04 13"),
                 pd.Timestamp(f"{year}-07-13 13"),
                 pd.Timestamp(f"{year}-07-20 13"),
                 pd.Timestamp(f"{year}-07-29 13"),
                 pd.Timestamp(f"{year}-08-09 13"),
                 pd.Timestamp(f"{year}-08-18 13"),
                 pd.Timestamp(f"{year}-08-27 13")]#,
                 #pd.Timestamp(f"{year}-09-05 13")]
    
        #ts = pd.Series([290, 230, 275, 400, 280, 250, 265, 370, 420], dates)
        #ts = pd.Series([100 , 290, 160, 230, 86 ,
        #                410 ,  80,  55,  60, 180, 270], dates)
        
        #ts = pd.Series([140 , 295, 172, 230, 196 ,
        #                410 , 210, 280, 160, 221, 270], dates)
        
        ts = pd.Series([140 , 295, 172, 230, 196 ,
                        410 , 210, 280, 160, 221 ], dates)
    
        parc01_ctr = csm_data.cosmo_data(sf_parc01_ctr    , fn_cosmo, clm_name)
        parc_v35   = csm_data.cosmo_data(sf_parc_v35      , fn_cosmo, clm_name)
        parc_v45   = csm_data.cosmo_data(sf_parc_v45      , fn_cosmo, clm_name)
        parc_v45e  = csm_data.cosmo_data(sf_parc_v45_evap , fn_cosmo, clm_name)
        
        stom_ctr  = stomata(parc01_ctr, year)
        stom_v35  = stomata(parc_v35  , year)
        stom_v45  = stomata(parc_v45  , year)
        stom_v45e = stomata(parc_v45e , year )
    
    
        #------------------------------------------------------------------------------
        # Preparing data for stat analysis
        #------------------------------------------------------------------------------
        df_stom = pd.concat([stom_ctr, stom_v35, stom_v35, stom_v45e, ts], axis = 1)
        df_stom2 = df_stom.dropna()
        
        
        df_stom2.columns = ['COSMO_CTR' ,'COSMO_v3.5' , 
                            'COSMO_v4.5','COSMO_v4.5e',
                            'OBS']
           
    
    
    
    
    
        
        stat_stom_ctr  = statistic(df_stom2['OBS'], df_stom2['COSMO_CTR'  ])
        stat_stom_v35  = statistic(df_stom2['OBS'], df_stom2['COSMO_v3.5' ])
        stat_stom_v45  = statistic(df_stom2['OBS'], df_stom2['COSMO_v4.5' ])
        stat_stom_v45e = statistic(df_stom2['OBS'], df_stom2['COSMO_v4.5e'])
        
        df_com_stat = pd.concat([stat_stom_ctr, stat_stom_v35 ,
                                 stat_stom_v45, stat_stom_v45e], axis = 0)
        df_com_stat.to_excel(exit_path + f'RSTOM_{fn_region} {year}.xlsx', float_format='%.3f')
        
        
        fig1  = plt.figure(figsize = (14,10))
        ax1   = fig1.add_subplot(111)
                
        
        stom = vsp.plots5_stomata(ax1, stom_ctr ,   stom_v35,   stom_v45,   stom_v45e,    ts,
                                       'CCLMref', 'CCLMv3.5', 'CCLMv4.5', 'CCLMv4.5e', 'OBS')
        
        stom = vsp.lplots_stomata(ax1, name_2[0],  nst, l_p,  y_min[0]    ,  
                                        y_max[0], y_step[0],  station_plot)
        
        plt.savefig(exit_path + f'RSTOM_{fn_region} {year}.png', format = 'png', dpi = 300) 
