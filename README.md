# The results of my work in CESR (Kassel University) for VAINT project

The land surface processes significantly affect the conditions in the low-level atmosphere ([**Tölle and Churiulin, 2021**][1]). The main parameters, which determine the interactions between the land surface and atmosphere, are the soil water content ([**Koster et. al., 2002**][2]) and the surface roughness ([**de Noblet-Ducoudre and Pitman, 2021**][3]). The impact of surface processes is evident in the low-level temperature, humidity, the structure of the planetary boundary layer and precipitation ([**Arora, 2002**][4]). [**Tölle et al. 2014**][5] have shown in climate simulations at convection-permitting scale that vegetation type changes can have a significant impact on extreme temperatures. Associated changes in vegetation phenology influence the energy and water cycle. Therefore, atmospheric models have to represent the land surface processes in a realistic way. However, the evapotranspiration simulated by the multilayer land surface scheme *TERRA-ML* of the [**Consortium for Small-scale Modelling**][6] – *COSMO* was found to be systematically underestimated from April to October during the growing season. One of the possible reasons for the underestimation of evapotranspiration is connected with the fact that in *TERRA-ML* the vegetation is not sufficiently represented in the surface energy balance ([**Schulz et al., 2015**][7]).

**Project purpose:** Improving the current phenology of vegetation and photosynthesis in the COSMO model. Preparing the work for a future implementation in the ICON model. 

**Project tasks:**
1. Developing the new vegetation parameterisation scheme for [COSMO-CLM][8] v5.16 and v6.0 based on [CLM][9] and [SURFEX][10] models;
2. Adaptating and implemeting new vegetation parameterisation scheme into COSMO-CLM v5.16 and v6.0;
3. Developing the instruments for postprocessing model results;
4. Developing the instruments for data analysis and visualisation; 

**Project methods:**
Depending on the project task, the different methods and programing languages were used for achieving the final results:
1. **Task 1** - For realisation of the first task I have used the next methods:
    1. collecting information from literature review and COSMO-CLM, CLM and SURFEX model documentations;
    2. analysing and researching the model codes and other resources (e.g., data for verification and estimation of model results);
    3. creating the blockschemes of important for modernisation modules of different models (SURFEX, CLM and [TERRA-ML][11]).  I used a special [Diagram Editor][12] for creating blockschemes. 
2. **Task 2** - The COSMO-CLM model uses a *Fortran90* as a main programing language, because of that all new updates have been written in Fortran90.
3. **Task 3** - I have done my personal *Shell* scripts for initial postprocessing of COSMO-CLM output data;
4. **Task 4** - I have done my personal *Python* scripts for data analysis and visualisation.

## Project results:

1. **Task 1** - the blockshemes of different models are presented in Table 1.

Table 1. Blockschemes of different models used in the project:
<table>
<tr>
<td><b>Model:</b></td>
<td><b>Blockscheme:</b></td>  
<td><b>Model:</b></td>
<td><b>Blockscheme:</b></td>
<td><b>Model:</b></td>
<td><b>Blockscheme:</b></td> 
<tr>
    
<td>COSMO-CLM v5.16</td>
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/Blockschemes/COSMO-CLM/src_radiation.jpg" target="_blank"><b>src_radiation</b></a></td>      
<td>COSMO-CLM v5.16</td>
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/Blockschemes/COSMO-CLM/terra_ml.png" target="_blank"><b>terra_ml</b></a></td>      
<td>SURFEX</td>
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/Blockschemes/SURFEX/vegetation_evol.jpg" target="_blank"><b>vegetation_evol</b></a></td>      
<tr>    

<td>CLM3.5</td>
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/Blockschemes/CLM/1.%20CLM3.5%20model%20structure.jpg" target="_blank"><b>modules structure</b></a></td>      
<td>CLM3.5</td>
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/Blockschemes/CLM/2.%20Program_off%20(CLM).jpg" target="_blank"><b>program off</b></a></td>      
<td>CLM3.5</td>
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/Blockschemes/CLM/3.%20clm_comp%20(CLM).jpg" target="_blank"><b>clm comp</b></a></td>      
<tr>  

<td>CLM3.5</td>
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/Blockschemes/CLM/4.%20Initialize1%20(CLM).jpg" target="_blank"><b>initialize 1</b></a></td>      
<td>CLM3.5</td>
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/Blockschemes/CLM/5.%20Initialiaze2%20-%20important%20(CLM).jpg" target="_blank"><b>initialize 2</b></a></td> 
<td>CLM3.5</td>
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/Blockschemes/CLM/6.%20InitSurfalbMod%20-%20first%20step%20(CLM).jpg" target="_blank"><b>initSurfalbMod</b></a></td>      
<tr>
    
<td>CLM3.5</td>
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/Blockschemes/CLM/7.%20Driver1%20-%20calculations%20(CLM).jpg" target="_blank"><b>driver 1</b></a></td>      
<td>CLM3.5</td>
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/Blockschemes/CLM/8.%20Driver2%20-%20printing%20(CLM).jpg" target="_blank"><b>driver 2</b></a></td> 
<td>CLM3.5</td>
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/Blockschemes/CLM/9.%20CanopyFluxesMod.jpg" target="_blank"><b>CanopyFluxesMod</b></a></td>      
<tr>    
    
<td>CLM3.5</td>
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/Blockschemes/CLM/10.%20SurfaceRadiationMod.jpg" target="_blank"><b>SurfaceRadiationMod</b></a></td>      
<td>CLM3.5 STATICE</td>
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/Blockschemes/CLM/STATICE/STATICEcosysDynMod.jpg" target="_blank"><b>STATICEcosysdynMOD</b></a></td> 
<td>CLM3.5 DGVM</td>
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/Blockschemes/CLM/DGVM/DGVMMod.jpg" target="_blank"><b>modules structure</b></a></td>    
<tr>    
    
<td>CLM3.5 DGVM</td>
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/Blockschemes/CLM/DGVM/1.%20DGVMEcosystemDyniti.jpg" target="_blank"><b>DGVMEcosystemDyniti</b></a></td>      
<td>CLM3.5 DGVM</td>
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/Blockschemes/CLM/DGVM/2.%20DGVMEcosystemMod.jpg" target="_blank"><b>DGVMEcosystemDyn</b></a></td> 
<td>CLM3.5 DGVM</td>
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/Blockschemes/CLM/DGVM/3.%20Phenology%20%2B%20Fire.jpg" target="_blank"><b>phenology</b></a></td>     
<tr>     
    
<td>CLM3.5 DGVM</td>
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/Blockschemes/CLM/DGVM/4.%20DGVMRespiration.jpg" target="_blank"><b>respiration</b></a></td>    
<td>CLM3.5 DGVM</td>
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/Blockschemes/CLM/DGVM/5.%20LitterSOM.jpg" target="_blank"><b>litterSOM</b></a></td>     
<td>CLM3.5 CN</td>
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/Blockschemes/CLM/CN/1.%20CN_model_structure.jpg" target="_blank"><b>modules structure</b></a></td>  
<tr>
    
<td>CLM3.5 CN</td>
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/Blockschemes/CLM/CN/2.%20CNEcosystemDyn.jpg" target="_blank"><b>CNEcosystemDyn</b></a></td>    
<td>CLM3.5 CN</td>
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/Blockschemes/CLM/CN/3.%20CNSetValueMod.jpg" target="_blank"><b>CNSetValueMod</b></a></td>     
<td>CLM3.5 CN</td>
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/Blockschemes/CLM/CN/4.%20CNMResp.jpg" target="_blank"><b>CNMRespMod</b></a></td>  
<tr>

<td>CLM3.5 CN</td>
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/Blockschemes/CLM/CN/5.%20CNDecompMod.jpg" target="_blank"><b>CNDecompMod</b></a></td>    
<td>CLM3.5 CN</td>
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/Blockschemes/CLM/CN/6.%20CNAllocaionMod.jpg" target="_blank"><b>CNAllocationMod</b></a></td>     
<td>CLM3.5 CN</td>
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/Blockschemes/CLM/CN/7.%20CNStateUpdate1Mod.jpg" target="_blank"><b>CNStateUpdate1Mod</b></a></td>  
<tr>    
    
<td>CLM3.5 CN</td>
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/Blockschemes/CLM/CN/8.%20CNStateUpdate2Mod.jpg" target="_blank"><b>CNStateUpdate2Mod</b></a></td>    
<td>CLM3.5 CN</td>
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/Blockschemes/CLM/CN/9.%20CNStateUpdate3Mod.jpg" target="_blank"><b>CNStateUpdate3Mod</b></a></td>     
<td>CLM3.5 CN</td>
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/Blockschemes/CLM/CN/10.%20CNVegStructUpdate.jpg" target="_blank"><b>CNVegStructUpdate</b></a></td>    
</table>

2. **Task 2** - the new vegetation parameterisation scheme has been implemented and adapted for COSMO-CLM model versions v5.16 and v6.0. The new algorithms apply the physical [Bell-Berry approach][13] coupled with [Farquhar][14] and [Collatz][15] leaf photosynthesis models and [two big leaf approach][16]. In general, the 4 different versions have been developed: 
    1. COSMO-CLMv5.16 with vegetation algorithms adapted from the Community Land Model (CLM) version 3.5. This version has a name CCLMv3.5;
    2. COSMO-CLMv5.16 with vegetation algorithms adapted from the Community Land Model (CLM) version 4.5. This version has a name CCLMv4.5;
    3. COSMO-CLMv5.16 with vegetation algorithms adapted from the Community Land Model (CLM) version 4.5 and additional changes in algorithm for calculation of transpiration from the dry leaf surface. This version has a name CCLMv4.5e;
    4. [COSMO-CLMv6.0][18] with vegetation algorithms adapted from the Community Land Model (CLM) version 4.5. Also the additional LAI algorithm from SURFEX model have been implemented in the latest version of COSMO-CLM. 

The full code of COSMO-CLMv5.16 or COSMO-CLM 6.0 with updates is available by request (evgenychur@gmail.com). Also, updated code of COSMO-CLMv6.0 can be available after the registration for the official members of the [COSMO consortium][17]

Table 2. Updates implemented in COSMO-CLM v6.0:
<table>
<tr>
<td><b>Number:</b></td>
<td><b>Changes:</b></td>  
<td><b>Available code:</b></td>
<tr>

<td><b>1</b></td>
<td>Forthran code with updates. The updates have been implemented in 17 modules of COSMO-CLMv6.0 and use tiles structure as other COSMO-CLM soil parameters</td>
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/changes%20in%20COSMO-CLMv6.0.f90" target="_blank"><b>updates in COSMO-CLMv6.0</b></a></td>
<tr>

<td><b>2</b></td>
<td>Forthran code for the new vegetation parameterisation scheme. Module has 6 new subroutines for calculations stomatal resistance, leaf photosynthesis, two big leaf approach and other parameters required for correct work of the new scheme.</td>
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/sfc_phenology.f90"><b>sfc_phenology</b></a></td>
<tr>
    
<td><b>3</b></td>
<td>Forthran code with constant and PFT parameters requered for correct work of sfc_phenology module. In case, if you want to add new PFT the corresponding values should be added to table pft_CN_par and pft_SURFEX</td>
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/sfc_phenology_data.f90" target="_blank"><b>sfc_phenology_data</b></a></td>   
</table>    

3. **Task 3** - my personal shell scripts for postprocessing presented in folder ***script***. Folder SCRIPTS has subfolders:
    1. `COSMO_CLM_v6.0` – has scripts for postprocessing COSMO-CLMv6.0 output data;
    2. `COSMO_CLM_v5.16` – has scripts for postprocessing COSMO_v3.5, COSMO_v4.5, COSMO_v4.5_evap and COSMO-CLM-v5.0_clm16_orig output data;
    3. `DATASETS` – has scripts for postprocessing of EOBS, HYRAS and GLEAM data;
    4. `LAI` – has scripts for postprocessing of the simulations with the LAI algorithms implemented in COSMO-CLMv5.16 and satellite data;
    5. `LC_MAPS` – has scripts for statistical analysis of data. The results of these scripts were used for the manuscript in the Frontiers journal (https://doi.org/10.3389/feart.2021.722244);
    6. `SIMULATION_STAT` – has scripts for statistical analysis of data. The results of these scripts were used for the manuscript in the Biogeosciences journal (https://doi.org/10.5194/bg-2021-294).

Table 3. Shell scripts for postprocessing:
<table>
<tr>
<td><b>Number:</b></td>
<td><b>Folder:</b></td>  
<td><b>Purpose:</b></td> 
<td><b>Available code:</b></td>
<tr>

<td><b>1</b></td>
<td>COSMO_CLM_v6.0</td>
<td>Script for analysis of data from *COSMO_CLM_6.0_balance* version</td>    
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/script/COSMO_CLM_v6.0/balance_mode.sh" target="_blank"><b>balance_mode</b></a></td>
<tr>

<td><b>2</b></td>
<td>COSMO_CLM_v6.0</td>
<td>Script for analysis of the new vegetation algorithm implemented in *COSMO_CLM_6.0*</td>    
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/script/COSMO_CLM_v6.0/test_mode.sh" target="_blank"><b>test_mode </b></a></td>
<tr>

<td><b>3</b></td>
<td>COSMO_CLM_v5.16</td>
<td>Script works only for experimental simulations *CCLMv3.5*, *CCLMv4.5* and *CCLMv4.5_evap*. Versions can be corrected in the parameter `exp_version` and `ver_count` where the correct number of the versions should be written. COSMO-CLM output parameters in `params` are available only in the new simulations</td>    
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/script/COSMO_CLM_v5.16/add_params.sh" target="_blank"><b>add_params</b></a></td>
<tr>    

<td><b>4</b></td>
<td>COSMO_CLM_v5.16</td>
<td>Script works for all simulations due to the common output parameters of the simulations. The version can be corrected in the parameters `exp_version` and `ver_count` where the correct number of versions should be written</td>    
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/script/COSMO_CLM_v5.16/main_params.sh" target="_blank"><b>main_params</b></a></td>
<tr>  
    
<td><b>5</b></td>
<td>COSMO_CLM_v5.16</td>
<td>Script works only for the reference *v5.0_clm16_orig* simulation. The names of reference simulation (*ctr* or *cclmref*) can be corrected in the parameter `exp_version`. </td>    
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/script/COSMO_CLM_v5.16/ref_simul.sh" target="_blank"><b>ref_simul</b></a></td>
<tr>  
    
<td><b>6</b></td>
<td>DATASETS</td>
<td>Script for interpolation of EOBS dataset grid to COSMO-CLM rotated grid</td>    
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/script/DATASETS/EOBS_domain.sh" target="_blank"><b>EOBS_domain</b></a></td>
<tr>  
    
<td><b>7</b></td>
<td>DATASETS</td>
<td>Script for interpolation of HYRAS dataset grid to COSMO-CLM rotated grid</td>    
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/script/DATASETS/HYRAS_domain.sh" target="_blank"><b>HYRAS_domain</b></a></td>
<tr>  
    
<td><b>8</b></td>
<td>DATASETS</td>
<td>Script for interpolation of GLEAM dataset grid to COSMO-CLM rotated grid</td>    
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/script/DATASETS/GLEAM_domain.sh" target="_blank"><b>GLEAM_domain</b></a></td>
<tr>  
    
<td><b>9</b></td>
<td>DATASETS</td>
<td>Script for preparing NetCDF data of the EOBS dataset presented at COSMO-CLM rotated grid to python scripts and converting it to `.csv` format</td>    
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/script/DATASETS/EOBS_python.sh" target="_blank"><b>EOBS_python</b></a></td>
<tr>  
    
<td><b>10</b></td>
<td>DATASETS</td>
<td>Script for preparing NetCDF data of the HYRAS dataset presented at COSMO-CLM rotated grid to python scripts and converting it to `.csv` format</td>    
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/script/DATASETS/HYRAS_python.sh" target="_blank"><b>HYRAS_python</b></a></td>
<tr>  

<td><b>11</b></td>
<td>DATASETS</td>
<td>Script for preparing NetCDF data of the GLEAM dataset presented at COSMO-CLM rotated grid to python scripts and converting it to `.csv` format</td>    
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/script/DATASETS/GLEAM_python.sh" target="_blank"><b>GLEAM_python</b></a></td>
<tr>     

<td><b>12</b></td>
<td>LAI</td>
<td>Script for analysis of the results from 3 different simulations: *COSMO_LAI*, *_LAI2*, *_LAI3*. Script has a parameter exp_version where you can define the actual version of simulation</td>    
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/script/LAI/cosmo_lai.sh" target="_blank"><b>cosmo_lai</b></a></td>
<tr> 
    
<td><b>13</b></td>
<td>LAI</td>
<td>Script for postprocessing of MODIS satellite data with information about LAI. Script automatically work with *PARC*, *LINDEN* and *LINDENBERG* domains and require the **initial COSMO-CLM domains info**</td>    
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/script/LAI/modis_lai.sh" target="_blank"><b>modis_lai</b></a></td>
<tr>  
    
<td><b>14</b></td>
<td>LAI</td>
<td>Script for postprocessing of COPERNICUS satellite data with information about LAI. Script automatically work with *PARC*, *LINDEN* and *LINDENBERG* domains and require the **initial COSMO-CLM domains info**</td>    
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/script/LAI/copernicus_lai.sh" target="_blank"><b>copernicus_lai</b></a></td>
<tr>  
    
<td><b>15</b></td>
<td>LC_MAPS</td>
<td>Script for statistical analysis. The land cover map LU_G was set as a reference experiment and all statistical results were given from this assumption</td>    
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/script/LC_MAPS/G_refer.sh" target="_blank"><b>G_refer</b></a></td>
<tr>  
    
<td><b>16</b></td>
<td>LC_MAPS</td>
<td>Script for statistical analysis. The land cover map LU_GC was set as a reference experiment and all statistical results were given from this assumption<</td>    
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/script/LC_MAPS/GC_refer.sh" target="_blank"><b>GC_refer</b></a></td>
<tr>  
    
<td><b>17</b></td>
<td>LC_MAPS</td>
<td>Script for statistical analysis. HYRAS observational gridded dataset was set as a reference. Statistical results were given from this assumption;</td>    
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/script/LC_MAPS/HYRAS_refer.sh" target="_blank"><b>HYRAS_refer</b></a></td>
<tr>  
    
<td><b>18</b></td>
<td>SIMULATION_STAT</td>
<td>Comparison of <i>GLEAM</i> data with <i>AEVAP</i> and <i>ZVERBO</i> output parameters of <i>COSMO-CLMv5.16</i> simulations. Data for comparison is presented on COSMO-CLM rotated grid. Script automatically works only for one `domain` which is presented in domain parameter. Nevertheless, script is adapted for all domains</td>   
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/script/SIMULATION_STAT/gleam_data.sh" target="_blank"><b>gleam_data</b></a></td>
<tr>  
    
<td><b>19</b></td>
<td>SIMULATION_STAT</td>
<td>Comparison of <i>HYRAS</i>  and <i>E-OBS</i>  data with temperature (near surface – <i>T2m</i>, maximum – <i>T_MAX</i>, minimum – <i>T_MIN</i> and surface – <i>T_S</i>) output parameters of COSMO-CLMv5.16 simulations. Data for comparison is presented on COSMO-CLM rotated grid. Script automatically works only for one domain which is presented in domain parameter. Nevertheless, script is adapted for all domains.  The version of dataset can be chosen in refer parameter.</td>    
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/script/SIMULATION_STAT/eobs_hyras_data.sh" target="_blank"><b>eobs_hyras_data</b></a></td>  
</table>

4. **Task 4**

Table 4. Python scripts for postprocessing:

<table>
<tr>
<td><b>Number:</b></td>
<td><b>Purpose</b></td>  
<td><b>Available code:</b></td> 
<td><b>Connected with:</b></td>
<tr>
    
<td><b>1</b></td>
<td>Get data from TRY database. The data can be downloaded from the TRY Database by requests via  the TRY website www.try‐db.org/TryWeb/Prop0.php. After that this script can be used</td>
 <td><a href="https://github.com/EvgenyChur/CESR/blob/main/PYTHON/get_TRY_data.py" target="_blank"><b>get_TRY_data</b></a></td>   
<td> Script doesn't use other personal modules</td>    
<tr>     
 
<td><b>2</b></td>
<td>Script for creating new dataframes with stomatal resistance and Vcmax data based on filtered datasets from the TRY database</td>
 <td><a href="https://github.com/EvgenyChur/CESR/blob/main/PYTHON/get_new_data.py" target="_blank"><b>get_new_data</b></a></td>   
<td> Script uses the personal module: <b>get_TRY_data</b></td>    
<tr>     

<td><b>3</b></td>
<td>Script for creation of correlation heatmap based on the reference COSMO-CLM simulation</td>
 <td><a href="https://github.com/EvgenyChur/CESR/blob/main/PYTHON/get_heatmap.py" target="_blank"><b>get_heatmap</b></a></td>   
<td> Script uses the results of shell scripts for postprocessing of COSMO-CLM data. Also, the personal module: <b>cosmo_data</b> is used in it</td>    
<tr>

<td><b>4</b></td>
<td>Script for creation of stomatal resistance boxplot depending on data from TRY database. All TRY data have near surface air temperatures because of that the special T2m categories have created for data analysis</td>
 <td><a href="https://github.com/EvgenyChur/CESR/blob/main/PYTHON/get_boxplot_t2m.py" target="_blank"><b>get_boxplot_t2m</b></a></td>   
<td> Script uses the two personal modules: <b>cosmo_data</b> and <b>get_new_data</b></td>    
<tr>

<td><b>5</b></td>
<td>Script for creation of stomatal resistance boxplot and line plot depending on data from COSMO-CLM simulations. This script all to compare output COSMO-CLM stomatal resistance data with other simulations</td>
 <td><a href="https://github.com/EvgenyChur/CESR/blob/main/PYTHON/get_cosmo_boxplot.py" target="_blank"><b>get_cosmo_boxplot</b></a></td>   
<td> Script uses the three personal modules: <b>cosmo_data</b>, <b>vis_module</b> and <b>STOMATA</b></td>    

<td><b>6</b></td>
<td>Script for processing of NetCDF data. It is a special script for students of Yandex Practicum. </td>
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/PYTHON/netcdf_yandex_lai.py" target="_blank"><b>netcdf_yandex_lai</b></a></td>   
<td> Script doesn't use other personal modules</td>    

 <td><b>7</b></td>
<td>The main program for analysis of LAI satellite data from MODIS and COPERNICUS systems. Also, script has modules for analysis of COSMO-CLM output information about LAI</td>
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/PYTHON/lai_main.py" target="_blank"><b>lai_main</b></a></td>   
<td> Script uses next personal modules: <b>vis_module</b>, <b>get_lai_data</b>, <b>lai_copernicus</b>, <b>lai_modis</b> and <b>lai_cosmo</b></td>     
    
<td><b>8</b></td>
<td>Script for reading and calculating information about LAI from MODIS or COPERNICUS satellite data </td>
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/PYTHON/get_lai_data.py" target="_blank"><b>get_lai_data</b></a></td>   
<td>Script is under control of <b>lai_main</b></td>      
    
<td><b>9</b></td>
<td>Script works with NetCDF data based on COPERNICUS  satellite information about LAI</td>
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/PYTHON/lai_copernicus.py" target="_blank"><b>lai_copernicus</b></a></td>   
<td>Script is under control of <b>lai_main</b></td>  
    
<td><b>10</b></td>
<td>Script works with NetCDF data based on MODIS satellite information about LAI</td>
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/PYTHON/lai_modis.py" target="_blank"><b>lai_modis</b></a></td>   
<td>Script is under control of <b>lai_main</b></td>  
    
<td><b>11</b></td>
<td>Script works with COSMO-CLM simulations COSMO-CLM_LAI, _LAI2, _LAI3 and allows to create 4 verification linear plots for such output parameters as: SUR_LAI,
SUR_BIOMAS, RUNOFF_G, PLCOV</td>
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/PYTHON/lai_cosmo.py" target="_blank"><b>lai_cosmo</b></a></td>   
<td>Script uses the two personal modules: <b>cosmo_data</b> and <b>vis_module</b>. Script is under control of <b>lai_main</b></td>      
</table>

[1]: https://www.frontiersin.org/articles/10.3389/feart.2021.722244/full
[2]: https://journals.ametsoc.org/view/journals/hydr/3/3/1525-7541_2002_003_0363_ctdola_2_0_co_2.xml
[3]: https://oxfordre.com/climatescience/view/10.1093/acrefore/9780190228620.001.0001/acrefore-9780190228620-e-825
[4]: https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2001RG000103
[5]: https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2013JD020877
[6]: http://www.cosmo-model.org/
[7]: https://www.cosmo-model.org/content/model/documentation/newsLetters/newsLetter15/default.htm
[8]: https://www.cmcc.it/models/cosmo-clm-climate-limited-area-modelling-community
[9]: https://www.cesm.ucar.edu/models/cesm1.2/clm/CLM45_Tech_Note.pdf
[10]: https://www.umr-cnrm.fr/surfex/IMG/pdf/surfex_scidoc_v8.1.pdf
[11]: http://www.cosmo-model.org/content/model/documentation/core/default.htm
[12]: https://www.diagrameditor.com/
[13]: https://link.springer.com/chapter/10.1007%2F978-94-017-0519-6_48
[14]: https://link.springer.com/article/10.1007%2FBF00386231
[15]: https://www.sciencedirect.com/science/article/pii/0168192391900028?via%3Dihub
[16]: https://journals.ametsoc.org/view/journals/clim/20/15/jcli4222.1.xml
[17]: https://wiki.coast.hzg.de/clmcom
[18]: https://github.com/EvgenyChur/CESR/blob/main/changes%20in%20COSMO-CLMv6.0.f90
