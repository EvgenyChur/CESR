# -*- coding: utf-8 -*-
"""
The script for pre-processing of data from the TRY database

Plant trait data in the TRY Database are public under the Creative Commons 
Attribution License ‘CC BY’ (https://creativecommons.org/licenses/by/4.0).

Data is available by reguests via the TRY website www.try‐db.org/TryWeb/Prop0.php

More information about data is available:
Kattge, J, Bönisch, G, Díaz, S, et al. TRY plant trait database – enhanced
coverage and open access. Glob Change Biol. 2020; 26: 119– 188. 
https://doi.org/10.1111/gcb.14904


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


# Get the dataset names with information about species in the research
#------------------------------------------------------------------------------
def get_uniq_datasets(df, plant_species):
    '''
    Parameters
    ----------
    df : DataFrame
        The dataframe with all data from TRY database.
    Returns
    -------
    None.

    '''
    # Get uniq names of databases
    plants_databases = []
    for plants in plant_species:
        # Read data from request      
        data_grass = df.query('AccSpeciesName == @plants')
        plants_databases.append(data_grass['Dataset'].unique())
    # Names of dataset with information about species in the research
    uniq_datasets = pd.Series(np.concatenate(plants_databases, axis = 0)).unique()   
    return uniq_datasets
#------------------------------------------------------------------------------


# Get names of the datasets with actual plant species and create new dataframes
#------------------------------------------------------------------------------
def get_try_dataset(main_path, fnameTRY, plant_species, mode, ldata = True):
    '''
    Parameters
    ----------
    main_path     : Objects
        The main path for project folder.
    fnameTRY      : Objects
        Name of data TRY data.
    plant_species : List of objects
        Names of plant species
    mode          : Int
        The binary flag for data (1 - RSTOM,
                                  2 - VCMAX)
    ldata         : Logical
        The logical parameter for creation of output dataset
                                            (default - True)
    Returns
    -------
    uniq_datasets : Objects
        The unic names of datasets with actual plant species.
    Also the script save datasets in new DataFrames

    '''
   
    if mode == 1:
        # The type of data: stomatal resistance or Vcmax
        data_type = 'RSTOM'
    else:
        data_type = 'VCMAX'
    
    # Path to Vcmax data from TRY database
    st_data   = main_path + f'/REANALYSIS/TRY/{data_type}/{fnameTRY}'  
    # Output data
    path_exit = main_path + f'/REANALYSIS/TRY/{data_type}/DATASETS/'
        
    # Get original stomatal resistance data
    df = pd.read_csv(st_data, sep = '\t', header = 0, encoding='latin-1')
       
    uniq_datasets = get_uniq_datasets(df, plant_species)
    print('')
    print('The names of the datasets where there is information about actual plant species')
    print(uniq_datasets)
    
    
    # Get datasets only with the reseach species and save them in new dataframes
    #get_data = []
    if ldata == True:
        for i, dbase_name in enumerate(uniq_datasets):
            df_dataset = df.query('Dataset == @dbase_name')
            
            df_dataset.to_csv(path_exit + dbase_name + '.txt', 
                              sep      = '\t'                , 
                              encoding = 'latin-1'           )
        
        # Get data from a special datasets
        sp_dataset = "LBA ECO Tapajos: Leaf Characteristics and Photosynthesis"
        df_special = df.query('Dataset == @sp_dataset') 
        df_special.to_csv(path_exit + sp_dataset[0:14] + '.txt', 
                                      sep      = '\t'          , 
                                      encoding = 'latin-1'     ) 
    return uniq_datasets
#------------------------------------------------------------------------------

if __name__ == '__main__':
    # General path to the project folders
    main_path     = 'C:/Users/Churiulin/Desktop/COSMO_RESULTS'
    fnameTRY      = 'original_data.txt'  
    plant_species = ['Lolium perenne'       ,
                     'Poa pratensis'        ,
                     'Arrhenatherum elatius',
                     'Filipendula ulmaris'  ,
                     'Festuca rubra'        ]
    
    unic_name_rstom = get_try_dataset(main_path, fnameTRY, plant_species, 1)
    unic_name_vcmax = get_try_dataset(main_path, fnameTRY, plant_species, 2)
    