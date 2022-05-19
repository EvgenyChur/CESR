# -*- coding: utf-8 -*-
"""
The script for creation new dataframe with stomatal resistance and Vcmax data 
based on filtered datasets from the TRY database

Plant trait data in the TRY Database are public under the Creative Commons 
Attribution License ‘CC BY’ (https://creativecommons.org/licenses/by/4.0).

Data is available by reguests via the TRY website www.try‐db.org/TryWeb/Prop0.php

More information about data is available:
Kattge, J, Bönisch, G, Díaz, S, et al. TRY plant trait database – enhanced
coverage and open access. Glob Change Biol. 2020; 26: 119– 188. 
https://doi.org/10.1111/gcb.14904


The module contains several subroutines:
    1. get_data        - get data from the datasets with date, latitude and longitude
    2. get_new_TRYdata - Get new data for stomatal resistance or Vcmax based on TRY data 
                         (main function)
    3. get_timeseries - Get new data in spesial format (timeseries)


Autors of scripts: Evgenii Churiulin - Center for Enviromental System Research (CESR) 

                                                   
Current Code Owner: CESR, Evgenii Churiulin
phone:  +49  561 804-6142
fax:    +49  561 804-6116
email:  evgenychur@uni-kassel.de


History:
Version    Date       Name
---------- ---------- ----                                                   
    1.1    2022-02-16 Evgenii Churiulin, Center for Enviromental System Research (CESR)
           Initial release
                 
"""

import pandas as pd
import numpy as np
import datetime as dt
import get_TRY_data as gtd

# get data from the datasets with date, latitude and longitude
#------------------------------------------------------------------------------
def get_data(df, date_id, param_id, dt_base_name):
    '''
    Parameters
    ----------
    df             : dataframe
        table with data
    date_id        : integer 
        ID code for date.
    lat_id, lon_id : integer
        ID code for latitude and longitude
    param_id : integer
        ID code for parameters. The ID code can be found in DataID column of 
        input tables (st_data, vc_data). For example: ID code for:
                                                stomatal conductance - 50;
                                                Vcmax - 1275
    Returns
    -------
    dataset : dataframe
    '''
    # Date from database
    if dt_base_name == 'Ti Tree':
        date = pd.to_datetime(df['OrigValueStr'].loc[df['DataID'] == date_id],
                                      format='%b %Y').reset_index(drop = True)  
        
    elif dt_base_name in ('Leaf and Whole Plant Traits 2','Global Respiration'):
        date = df['OrigValueStr'].loc[df['DataID'] == date_id].reset_index(drop = True)
        
    elif dt_base_name in ('TROBIT West Africa',
                          'Leaf carbon exchange dataset for tropical, temperate, and boreal species of North and Central A'):
        date = pd.to_datetime(df['OrigValueStr'].loc[df['DataID'] == date_id],
                              format='%m/%d/%Y').reset_index(drop = True)
        
    elif dt_base_name in ('LBA ECO Tapajos'):
        date = pd.to_datetime(df['OrigValueStr'].loc[df['DataID'] == date_id],
                              format='%d/%m/%Y').reset_index(drop = True)
      
    # Actual parameters
    st_values  = (df[['AccSpeciesName',
                      'StdValue'      , 
                      'ErrorRisk'     ,
                      'UnitName'      ,
                      'Reference'     ,
                      'Dataset'       ]].loc[df['DataID'] == param_id]          
                                        .reset_index(drop = True))  
    # Create new datset
    dataset = pd.concat([date, st_values], axis = 1)
  
    return dataset
#------------------------------------------------------------------------------



# Get new data for stomatal resistance or Vcmax based on TRY data
#------------------------------------------------------------------------------
def get_new_TRYdata(main_path, plant_species, mode, ldf_save = True):
    '''
    Parameters
    ----------
    main_path : Objects
        The main path to the project.
    plant_species : List of objects
        The plant type species
    mode : Int
        The binary flag for data (1 - RSTOM,
                                  2 - VCMAX)
    ldf_save : Logical
        The parameter is responsible save or not the new dataframe
        (default - True)
    Returns
    -------
    df_new : DataFrame
        The new dataframe with data.

    '''
    # Define type of data
    if mode == 1:
        # The type of data: stomatal resistance
        data_request = 1
        data_type    = 'RSTOM'
        id_code      = 50                                                      # Stomatal resistance ID code in TRY database
    else:
        # The type of data: Vcmax 
        data_request = 2
        data_type = 'VCMAX'
        id_code      = 1275 
      
    # Path to data from TRY databases
    st_data   = main_path + f'/REANALYSIS/TRY/{data_type}/DATASETS/'  
    # Output data
    path_exit = main_path + '/REANALYSIS/TRY/DATASET/'

    # Format of the databases (IN, OUT) 
    dbf_in  = '.txt'
    dbf_out = '.csv'
    
    fnameTRY  = 'original_data.txt'                                            # Common name for TRY data (initial)

    # Get database names
    databases_list = (gtd.get_try_dataset(main_path   , 
                                         fnameTRY     , 
                                         plant_species, 
                                         data_request ,
                                         ldata = False)
                     ).tolist()

    # Create new DataFrame with stomatal resistance or Vcmax data     
    data_list = []
    
    for act_base in databases_list:
        # Define database 
        iPath_act  = st_data + act_base + dbf_in
                
        # Read data
        df = pd.read_csv(iPath_act, sep = '\t', header = 0, encoding = 'latin-1')
                
        # Filter data by species
        plant_filter = []
        for plant in plant_species:
            plant_filter.append(df.query('(AccSpeciesName == @plant)'))
        data_st = pd.concat(plant_filter, axis = 0)
            
        # Get date relating to stomatal resitance
        st_values   = (data_st[['AccSpeciesName', 'StdValue', 
                                'ErrorRisk'     , 'UnitName', 
                                'Reference'     , 'Dataset' ]].loc[data_st['DataID'] == id_code]
                                                              .reset_index(drop = True))   
        # get data related to temperature
        if act_base in ('Onoda 2017 leaf dataset', 
                        'The DIRECT Plant Trait Database'):
            st_temp = data_st['OrigValueStr'].loc[data_st['DataID'] == 322].reset_index(drop = True)
                
        elif act_base in ('Global Respiration Database'):
            st_temp = data_st['OrigValueStr'].loc[data_st['DataID'] ==  51].reset_index(drop = True)
                
        elif act_base in ('Photosynthesis Traits Worldwide', 
                              'GLOPNET - Global Plant Trait Network Database'):
            st_temp = data_st['OrigValueStr'].loc[data_st['DataID'] ==  62].reset_index(drop = True)
                
        else:
            st_temp = pd.Series(25.0, index = st_values.index)
            
        dataset_new = pd.concat([st_values, st_temp], axis = 1)
        dataset_new.columns = ['plant species', f'{data_type}', 
                               'error'        , 'units'       , 
                               'reference'    , 'dataset'     , 'T_2M' ]
        data_list.append(dataset_new)        
        
    # Prepare new dataset for furher work
    df_new = pd.concat(data_list, axis = 0).dropna()    
    df_new['T_2M'] = df_new['T_2M'].astype('float64')
    df_new[f'{data_type}'] = (df_new[f'{data_type}'].replace(df_new[f'{data_type}'][df_new['T_2M'] < 7.0], 347))  
    df_new['T_2M'] = (df_new['T_2M'].replace(df_new['T_2M'][df_new['T_2M'] < 6.0], 10.6)
                                    .replace(df_new['T_2M'][df_new['T_2M'] < 7.0], 11.6)
                     )  
    if ldf_save == True:
        # Save new dataframe
        df_new.to_csv(f'{path_exit}{data_type}_TRY_databases{dbf_out}',
                       sep=';', encoding = 'utf-8', 
                       float_format = '%9.3f')  
    return df_new
#------------------------------------------------------------------------------    

# Get timeseries based on TRY data 
#------------------------------------------------------------------------------
def get_timeseries(main_path, plant_species, database, mode):
    if mode == 1:
        data_type = 'RSTOM'
    else:
        data_type = 'VCMAX'
                  
    # # Path to data from TRY databases  
    st_data   = main_path + f'/REANALYSIS/TRY/{data_type}/DATASETS/'  
    # Output data
    path_exit = main_path + '/REANALYSIS/TRY/DATASET/TIMESERIES/' 
    
    # Get data
    df_TRY = pd.read_csv(st_data + f'{database}.txt', sep = '\t', header = 0, encoding = 'latin-1')
                 
    if mode == 1:
        if database == 'LBA ECO Tapajo':
            try_id  = [241, 2350]  
        elif database in ('Leaf and Whole Plant Traits 2',
                          'Global Respiration'           ):
            try_id  = [112,   50]       
        elif database in ('TROBIT West Africa',
                          'Leaf carbon exchange dataset for tropical, temperate, and boreal species of North and Central A',
                          'Ti Tree'):
            try_id  = [241,   50]
           
    elif mode == 2:
        if database == 'LBA ECO Tapajo':
            try_id  = [241, 550]    
        elif database == 'Leaf and Whole Plant Traits 2':
            try_id  = [241,  550]  
        elif database in ('Ti Tree',
                          'Global Respiration'):       
            try_id  = [112, 1275]         
        elif database in ('TROBIT West Africa',
                          'Leaf carbon exchange dataset for tropical, temperate, and boreal species of North and Central A'):  
            try_id  = [241, 1275]      
                             
    # Get data
    #                   df      dateID    ParamID   database_name
    df_new = get_data(df_TRY, try_id[0], try_id[1], database)        
    
    
     # Filter data by species
    plant_filter = []
    for plant in plant_species:
        plant_filter.append(df_new.query('(AccSpeciesName == @plant)'))
    data_st = pd.concat(plant_filter, axis = 0)
    
    data_st.columns = ['index', 'plant species', f'{data_type}', 'ERROR',
                       'units', 'reference'    , 'dataset']
    
    data_st = data_st.set_index(data_st['index']).drop(['index'], axis = 1)
    # Save the final dataset 
    data_st.to_csv(path_exit + database + f'_{data_type}.csv', 
                   sep = ';', encoding = 'utf-8', 
                   float_format = '%9.3f'       )
    #data_st = data_st.set_index(data_st['index'])

    return data_st
#------------------------------------------------------------------------------


if __name__ == '__main__':
    # General path to the project folders
    main_path     = 'C:/Users/Churiulin/Desktop/COSMO_RESULTS'
    plant_species = ['Lolium perenne'       ,
                     'Poa pratensis'        ,
                     'Arrhenatherum elatius',
                     'Filipendula ulmaris'  ,
                     'Festuca rubra'        ]
    '''
    mode = 1 ---> stomatal resistance data
    mode = 2 ---> Vcmax data
    '''
    df_stom  = get_new_TRYdata(main_path, plant_species, 1, ldf_save = True)
    df_vcmax = get_new_TRYdata(main_path, plant_species, 2, ldf_save = True)
    
    # Test 2
    plant_species = ['Miconia acinodendron',
                     'Brachiaria brizantha']
    
    # Name of the database
    database    = 'LBA ECO Tapajo' 
    
    # Get data
    rstom = get_timeseries(main_path, plant_species, database, 1)
    vcmax = get_timeseries(main_path, plant_species, database, 2)     
