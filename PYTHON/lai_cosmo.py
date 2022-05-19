# -*- coding: utf-8 -*-
"""
Script works with COSMO-CLM simulations COSMO-CLM_LAI, _LAI2, _LAI3 and allows
to create 4 verification linear plots for such output parameters as: SUR_LAI,
SUR_BIOMAS, RUNOFF_G, PLCOV


Script contains subroutines:
    get_cosmo - get COSMO-CLM data;
    plot_cosmo - plot 4 linear plots
       
    
Autors of project: Evgenii Churiulin, Center for Enviromental System
                                                   Research (CESR) 

Current Code Owner: CESR, Evgenii Churiulin
phone:  +49  561 804-6142
fax:    +49  561 804-6116
email:  evgenychur@uni-kassel.de


History:
Version    Date       Name
---------- ---------- ----                                                   
    1.1    2022-04-15 Evgenii Churiulin, Center for Enviromental System Research (CESR)
           Initial release
"""

# Import standart liblaries 
import pandas as pd
import matplotlib.pyplot as plt

# Import personal libraries
import sys
sys.path.append('C:/Users/Churiulin/Python/scripts/CESR_project')
import vis_module       as vis                                                                                 
import cosmo_data       as csm_data                                            
 

# The function for COSMO data preprocessing
def get_cosmo(data_path, fn_cosmo, clm_name, t_start, t_stop):
    '''
    Parameters
    ----------
    data_path : object
        The absolute path to the data folder.
    fn_cosmo : object
        The common part of the COSMO-CLM files.
    clm_name : list
        The names of COSMO-CLM parameters for analysis
    t_start  : object
        The date of start calculations
    t_stop   : object
        The date of end calculations

    Returns
    -------
    df - the dataframe with COSMO-CLM data.
    '''
    #--------------------------------------------------------------------------
    # Get initial data: COSMO
    #--------------------------------------------------------------------------
    # The COSMO data has a hourly timestep
    df_cclm  = csm_data.cosmo_data(data_path, fn_cosmo, clm_name)

    #--------------------------------------------------------------------------
    # Get time periods
    #--------------------------------------------------------------------------
    time_start = pd.to_datetime([t_start])
    time_stop  = pd.to_datetime([t_stop ])

    t_start = 0                                                                      
    t_stop  = 23   

    #--------------------------------------------------------------------------
    # Section: Work with data and data visualization
    #--------------------------------------------------------------------------  

    # Ged monthly date for time periods and plot maps
    for  index in range(len(time_start)):   
        #----------------------------------------------------------------------
        # Create time periods 
        #---------------------------------------------------------------------- 
        hourly_period = pd.date_range(time_start[index], 
                                      time_stop[index] ,  
                                      freq = 'H'       )                       # hourly timesteps
        
        daily_period  = pd.date_range(time_start[index],
                                      time_stop[index] ,
                                      freq = 'D'       )                       # dayly timesteps 
        
        # General time period for COSMO, FLUXNET and EURONET data
        res_period   = [x for x in hourly_period if x.hour >= t_start and x.hour <= t_stop]
           
        #----------------------------------------------------------------------
        # Subsection: Create data for plots
        #----------------------------------------------------------------------    
        # COSMO data --> daily mean
        df = csm_data.get_timeseries(df_cclm, clm_name, res_period, 'D') # CCLMref  --> original COSMO  
    
        df['COSMO_LAI'] = df['LAI_SHA'] + df['LAI_SUN'] 
    
    return df


def plot_cosmo(df, clm_name, plot_label, axis_label, path_exit, t_start, t_stop):
    ''' 
    Parameters
    ----------
    df : DateFrame
        The dataframe with COSMO-CLM data.
    clm_name   : list
        The names of COSMO-CLM parameters for analysis.
    plot_label : list
        The names of plots
    axis_label : list
        The names of y axis
    path_exit  : object
        The output path
    t_start    : object
        The date of start calculations
    t_stop     : object
        The date of end calculations

    Returns
    -------
    None.

    '''
    #--------------------------------------------------------------------------
    # Get time periods
    #--------------------------------------------------------------------------
    time_start = pd.to_datetime([t_start])
    time_stop  = pd.to_datetime([t_stop ])

    t_start = 0                                                                      
    t_stop  = 23 
    
    # Ged monthly date for time periods and plot maps
    for  index in range(len(time_start)):       
        # Additional information for plot title: General            
        time_int_1 = str(time_start[index])[0:10]                                  # The date of period start --> need only for print
        time_int_2 =  str(time_stop[index])[0:10]                                  # The date of period stop  --> need only for print
        date_ind   = f'Time step: {time_int_1} to {time_int_2}'                    # The full date of period  --> need for plot label          
        l_p        = 'upper left'                                                  # The position of legend
        nst        = 'Parc domain'  
        
        #----------------------------------------------------------------------
        # Simple plot for visualisation of LAI fields from COSMO-CLM simulation
        #----------------------------------------------------------------------     
        fig = plt.figure(figsize = (14,10))
        ax  = fig.add_subplot(111) 
        
        plot2par = vis.plots2(ax, df[clm_name[0]],  df[clm_name[3]]   ,
                                  'SURFEX_LAI'   ,  'COSMO_LAI extpar')
                    
        plot2par = vis.lplots(ax, plot_label[0], 
                                  axis_label[0],
                                  date_ind, nst, l_p,
                                  0.0, 3.01, 0.25   ,
                                  time_start[index] , 
                                  time_stop[index]  ) 
        
        output_name = f'{clm_name[0]}_{time_int_1}_{time_int_2}.png'
                
        plt.savefig(path_exit + output_name, format = 'png', dpi = 300) 

        #----------------------------------------------------------------------
        # Simple plot for visualisation of PLCOV field from COSMO-CLM simulation
        #---------------------------------------------------------------------- 
        fig2 = plt.figure(figsize = (14,10))
        bx  = fig2.add_subplot(111) 
        
        plot1par = vis.plots1(bx, df[clm_name[4]], 'PLCOV')
                    
        plot1par = vis.lplots(bx, plot_label[1],
                                  axis_label[1], 
                                  date_ind, nst, l_p,
                                  0.65, 0.851, 0.05 ,
                                  time_start[index] ,
                                  time_stop[index]  ) 
        
        output_name = f'{clm_name[4]}_{time_int_1}_{time_int_2}.png'
                
        plt.savefig(path_exit + output_name, format = 'png', dpi = 300) 
    
        #----------------------------------------------------------------------
        # Simple plot for visualisation of runoff fields from COSMO-CLM simulation
        #---------------------------------------------------------------------- 
        fig3 = plt.figure(figsize = (14,10))
        cx  = fig3.add_subplot(111) 
        # plot runioff field
        plot2par = vis.plots2(cx, df[clm_name[5]], df[clm_name[6]],
                                  'subsurface_runoff_amount - RUNOFF_G', 
                                  'surface_runoff_amount - RUNOFF_S'   )
        # settings for plot            
        plot2par = vis.lplots(cx, plot_label[2], 
                                  axis_label[2], 
                                  date_ind, nst, l_p,
                                  0.0, 0.61, 0.10   ,
                                  time_start[index] ,
                                  time_stop[index]  ) 
        # output name of plot
        output_name = f'{clm_name[5]}_{time_int_1}_{time_int_2}.png'
                
        plt.savefig(path_exit + output_name, format = 'png', dpi = 300) 
    
        plt.close(fig)        
        plt.gcf().clear()  
        
        #----------------------------------------------------------------------
        # Simple plot for visualisation of biomass evolution from COSMO-CLM simulation
        #----------------------------------------------------------------------       
        fig4 = plt.figure(figsize = (14,10))
        dx   = fig4.add_subplot(111) 
        # plot biomass field
        plot1par = vis.plots1(dx, df[clm_name[7]], 'BIOMASS')
        # settings for plot            
        plot1par = vis.lplots(dx, plot_label[3],
                                  axis_label[3],
                                  date_ind, nst, l_p,
                                  0.0, 0.21, 0.05   ,
                                  time_start[index] ,
                                  time_stop[index]  ) 
        # output name of plot 
        output_name = f'{clm_name[7]}_{time_int_1}_{time_int_2}.png'
                
        plt.savefig(path_exit + output_name, format = 'png', dpi = 300) 
    
        plt.close(fig)        
        plt.gcf().clear() 
    








