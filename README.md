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

The full code of COSMO-CLMv5.16 or COSMO-CLM 6.0 with updates is available by requests. Also, updated code of COSMO-CLMv6.0 can be available after the registration for the official members of the [COSMO consortium][17]

Table 1. Blockschemes of different models used in the project:
<table>
<tr>
<td><b>Number:</b></td>
<td><b>Changes:</b></td>  
<td><b>Available code:</b></td>
<tr>

<td><b>1</b></td>
<td>Forthran code with updates. The updates have been implemented in 17 modules of COSMO-CLMv6.0 and use tiles structure as other COSMO-CLM soil parameters</td>
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/changes%20in%20COSMO-CLMv6.0.f90" target="_blank"><b>CNVegStructUpdate</b></a></td>
<tr>

<td><b>2</b></td>
<td>Forthran code for the new vegetation parameterisation scheme. Module has 6 new subroutines for calculations stomatal resistance, leaf photosynthesis, two big leaf approach and other parameters required for correct work of the new scheme.</td>
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/sfc_phenology.f90"><b>CNVegStructUpdate</b></a></td>
<tr>
    
<td><b>3</b></td>
<td>Forthran code with constant and PFT parameters requered for correct work of sfc_phenology module. In case, if you want to add new PFT the corresponding values should be added to table pft_CN_par and pft_SURFEX</td>
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/sfc_phenology_data.f90" target="_blank"><b>CNVegStructUpdate</b></a></td>   
</table>    

3. **Task 3**

4. **Task 4**



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
