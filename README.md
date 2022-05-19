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
<td><b>Number:</b></td>    
<td><b>Model:</b></td>
<td><b>Blockscheme:</b></td>
<tr>
    
<td><b>1</b></td>
<td>COSMO-CLM v5.16</td>
<td><a href="https://github.com/EvgenyChur/CESR/blob/main/src_radiation.jpg" target="_blank"><b>src_radiation.f90</b></a></td>      
<tr>    
</table>




[1]: https://www.cmcc.it/models/cosmo-clm-climate-limited-area-modelling-community
[2]: https://www.cesm.ucar.edu/models/cesm1.2/clm/CLM45_Tech_Note.pdf
[3]: https://www.umr-cnrm.fr/surfex/IMG/pdf/surfex_scidoc_v8.1.pdf
[4]: http://www.cosmo-model.org/content/model/documentation/core/default.htm
[5]: https://www.diagrameditor.com/
