# The results of my work in CESR (Kassel University) for VAINT project

**Project purpose:** Improving the current phenology of vegetation and photosynthesis in the COSMO model. Preparing the work for a future implementation in the ICON model. 

**Project tasks:**
1. Developing the new vegetation parameterisation scheme for [COSMO-CLM][1] v5.16 and v6.0 based on [CLM][2] and [SURFEX][3] models;
2. Adaptating and implemeting new vegetation parameterisation scheme into COSMO-CLM v5.16 and v6.0;
3. Developing the instruments for postprocessing model results;
4. Developing the instruments for data analysis and visualisation; 

**Project methods:**
Depending on the project task, the different methods and programing languages were used for achieving the final results:
1. **Task 1** - For realisation of the first task I have used the next methods:
    1. collecting information from literature review and COSMO-CLM, CLM and SURFEX model documentations;
    2. analysing and researching the model codes and other resources (e.g., data for verification and estimation of model results);
    3. creating the blockschemes of important for modernisation modules of different models (SURFEX, CLM and [TERRA-ML][4]).  I used a special [Diagram Editor][4] for creating blockschemes. 
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
</table>




[1]: https://www.cmcc.it/models/cosmo-clm-climate-limited-area-modelling-community
[2]: https://www.cesm.ucar.edu/models/cesm1.2/clm/CLM45_Tech_Note.pdf
[3]: https://www.umr-cnrm.fr/surfex/IMG/pdf/surfex_scidoc_v8.1.pdf
[4]: http://www.cosmo-model.org/content/model/documentation/core/default.htm
[5]: https://www.diagrameditor.com/
