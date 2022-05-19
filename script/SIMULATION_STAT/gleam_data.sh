#!/bin/bash

#-------------------------------------------------------------------------------
#
# Current code owner:
#
#   Center Enviroment System Research (CESR)
#
# Authors:
#
#   CESR, 2021
#   Evgenii Churiulin, Vladimir Kopeykin, Merja Toelle
#   phone:  +49 170-261-51-04
#   email:  evgenychur@uni-kassel.de
#
#-------------------------------------------------------------------------------

#source ~/.profile


#========================= Settings ============================================
# Can be changed by user

#=============================== GLEAM =========================================

refer="GLEAM"

# The GLEAM version (v3.5a or v3.5b)
vrs_GLEAM=("v3.5a" "v3.5b")
vrs_count=1

# The GLEAM parameters for analysis
param_GLEAM=("Ep" "Et")

# The COSMO parameters for analysis
param_COSMO=("AEVAP_S" "ZVERBO")

# Count of parameters --> should be the same as param_GLEAM and param_COSMO
paramscount=1

# The research domain
#domain=("parc_domain" "linden_domain" "lindenberg_domain")
domain=("lindenberg_domain")
# Count of domain
domaincount=0
#============================== General settings ===============================

# The name of dataset: ds_1 - COSMO ORIGINAL
#                      ds_2 - COSMO v3.5 (b =  2000)
#                      ds_3 - COSMO v4.5 (b = 10000)
#                      ds_4 - COSMO v4.5 with changes in ztraleav
dataset="ds_4"

# Use dataset
if [ "${dataset}" == "ds_1" ]; then
    ds_folder="original"
    ds_name="COSMO_ORIG"

elif [ "${dataset}" == "ds_2" ]; then
    ds_folder="3.5"
    ds_name="COSMO_35"

elif [ "${dataset}" == "ds_3" ]; then
    ds_folder="4.5"
    ds_name="COSMO_45"

elif [ "${dataset}" == "ds_4" ]; then
    ds_folder="4.5e"
    ds_name="COSMO_45e"

else
    echo 'no datasets'
fi


#========================= Paths ===============================================
DIR=/work/bb1112/b381275
DIR_OUT=${DIR}/results/${refer}                                                 # Output path at WORK server of mistral
DIR_RESULT=/pf/b/b381275/COSMO_results/${refer}                                 # output path at PF   server of mistral

#========================= Generation operations ===============================
for (( z=0; z<=vrs_count ; z++ ));
do
    for (( i=0; i<=domaincount ; i++ ));
    do
        for (( j=0; j<=paramscount ; j++ ));
        do
            echo 'The GLEAM version: '${vrs_GLEAM[${z}]}
            echo 'The research territory: '${domain[${i}]}
            echo 'Parameters: '${param_COSMO[${j}]}

            #-----------------------------------------------------------------------
            # Define input and output GLEAM data
            #-----------------------------------------------------------------------

            # Input path for GLEAM data
            DIR_GLEAM=${DIR}/${refer}/${vrs_GLEAM[${z}]}/${domain}

            # Filename of observational dataset: GLEAM
            refer_fn="${param_GLEAM[${j}]}"_"${refer}"_"2010_2015".nc

            # Input and Output data:
            iGLEAM=${DIR_GLEAM}/"${refer_fn}"
            sGLEAM=${DIR_OUT}/"${refer}"_"${vrs_GLEAM[${z}]}"_"${param_COSMO[${j}]}"_"std"_"${domain[${i}]}"
            mGLEAM=${DIR_OUT}/"${refer}"_"${vrs_GLEAM[${z}]}"_"${param_COSMO[${j}]}"_"mean"_"${domain[${i}]}"
            dGLEAM=${DIR_OUT}/"${refer}"_"${vrs_GLEAM[${z}]}"_"${param_COSMO[${j}]}"_"mean_dav"_"${domain[${i}]}"

            #-----------------------------------------------------------------------
            # Define input and output COSMO data
            #-----------------------------------------------------------------------

            # Input paths for COSMO data
            DIR_COSMO=${DIR}/COSMO/"${domain[${i}]}"/"${ds_folder}"

            # Filename of model dataset: COSMO_ORIG, COSMO_v35, COSMO_v45, COSMO_v45e
            #model_ds="${param_COSMO[${j}]}"_"ts_1999_2015.nc"
            model_ds="${param_COSMO[${j}]}"_"ts_1999_2015.nc"
            # Paths for semiresults
            iCOSMO=${DIR_COSMO}/"${model_ds}"
            oCOSMO=${DIR_OUT}/"${ds_name}"_"${param_COSMO[${j}]}"_"${vrs_GLEAM[${z}]}"_"daily.nc"

            # Paths for COSMO mean data and std data
            sCOSMO=${DIR_OUT}/"${ds_name}"_"${param_COSMO[${j}]}"_"${vrs_GLEAM[${z}]}"_"std"             #  STD COSMO
            mCOSMO=${DIR_OUT}/"${ds_name}"_"${param_COSMO[${j}]}"_"${vrs_GLEAM[${z}]}"_"mean"            # Mean COSMO


            so_COSMO=${DIR_OUT}/"${ds_name}"_"${param_COSMO[${j}]}"_"${vrs_GLEAM[${z}]}"_"std"_"${domain[${i}]}"                                         # For COSMO STD
            mo_COSMO=${DIR_OUT}/"${ds_name}"_"${param_COSMO[${j}]}"_"${vrs_GLEAM[${z}]}"_"mean"_"${domain[${i}]}"                                        # For COSMO MEAN
            da_COSMO=${DIR_OUT}/"${ds_name}"_"${param_COSMO[${j}]}"_"${vrs_GLEAM[${z}]}"_"mean_dav"_"${domain[${i}]}"                                    # For COSMO MEAN DAV
            o_CORR=${DIR_OUT}/"Corr"_"${refer}"_"${vrs_GLEAM[${z}]}"_"${ds_name}"_"${param_COSMO[${j}]}"_"${domain[${i}]}"                               # For CORR


            #-------------------------------------------------------------------
            # Calculations for KGE and RMSD
            #-------------------------------------------------------------------

            # Convert MODEL data to GLEAM format
            if [ "${param_COSMO[${j}]}" == "AEVAP_S" ]; then
                cdo -daysum -mulc,-1.0 "${iCOSMO}" "${oCOSMO}"
            else
                cdo -daymean -mulc,-1.0 -mulc,100000 "${iCOSMO}" "${oCOSMO}"
            fi

            # Calcuate standart deviation (std) and mean for MODEL
            cdo timstd "${oCOSMO}" "${sCOSMO}"."nc"
            cdo timmean "${oCOSMO}" "${mCOSMO}"."nc"

            # Calcuate standart deviation (std) and mean for GLEAM
            cdo timstd "${iGLEAM}" "${sGLEAM}"."nc"
            cdo timmean "${iGLEAM}" "${mGLEAM}"."nc"

            #Calculation of correlation between reference and dataset
            cdo timcor "${iGLEAM}" "${oCOSMO}" "${o_CORR}"."nc"
            #cdo fldcor "${iGLEAM}" "${oCOSMO}" "${o_CORR}"."nc"
            #-----------------------------------------------------------------------

            # Output NetCDF > csv for MODEL dataset
            cdo -outputtab,date,lon,lat,value "${sCOSMO}"."nc" > "${so_COSMO}"."csv"
            cdo -outputtab,date,lon,lat,value "${mCOSMO}"."nc" > "${mo_COSMO}"."csv"

            # Output NetCDF > csv for GLEAM dataset
            cdo -outputtab,date,lon,lat,value "${sGLEAM}"."nc" > "${sGLEAM}"."csv"
            cdo -outputtab,date,lon,lat,value "${mGLEAM}"."nc" > "${mGLEAM}"."csv"

            # Output NetCDF > csv for correletion
            cdo -outputtab,date,lon,lat,value "${o_CORR}"."nc" > "${o_CORR}"."csv"

            #-----------------------------------------------------------------------
            # Calculations for DAV
            #-----------------------------------------------------------------------
            # Calculate std, mean for GLEAM dataset
            cdo fldmean "${iGLEAM}" "${dGLEAM}"."nc"
            # Calculate std, mean for MODEL dataset
            cdo fldmean "${oCOSMO}" "${da_COSMO}"."nc"
            #-----------------------------------------------------------------------
            # Output NetCDF > csv for GLEAM dataset
            cdo -outputts "${dGLEAM}"."nc" > "${dGLEAM}"."csv"
            # Output NetCDF > csv for MODEL dataset
            cdo -outputts "${da_COSMO}"."nc" > "${da_COSMO}"."csv"

        done
    done
done
cp -R ${DIR_OUT}/*."csv" ${DIR_RESULT}
rm -r ${DIR_OUT}/*






#"${refer}"_"${param_GLEAM[${j}]}"_"${domain[${i}]}".nc


        #cdo fldstd "${i_GLEAM_data}" "${s_GLEAM_data}"."nc"
        #cdo fldmean "${i_GLEAM_data}" "${m_GLEAM_data}"."nc"
        #cdo fldstd "${oCOSMO_s1}" "${s_cosmo_data}"."nc"
        #cdo fldmean "${oCOSMO_s1}" "${m_cosmo_data}"."nc"


        # Convert from NetCDF format > csv format for model dataset
        #cdo outputts "${s_cosmo_data}"."nc" > "${s_output_cosmo}"."csv"
        #cdo outputts "${m_cosmo_data}"."nc" > "${m_output_cosmo}"."csv"

                # NetCDF > csv for observations
        #cdo outputts "${s_GLEAM_data}"."nc" > "${s_GLEAM_data}"."csv"
        #cdo outputts "${m_GLEAM_data}"."nc" > "${m_GLEAM_data}"."csv"
