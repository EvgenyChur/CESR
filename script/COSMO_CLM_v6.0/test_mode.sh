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

# possible options: 1 - Parc; 2 - Linden; 3 - Lindenberg
mode="1"

# possible versions: clm6_balance
exp_version=("clm6_test")
ver_count=0

params=("AEVAP_S_ts"    "ALHFL_S_ts"    "ALHFL_BS_ts"    "ALHFL_PL_ts"     "ASHFL_S_ts" \
          "QV_2M_ts"       "QV_S_ts"        "T_2M_ts"     "TMAX_2M_ts"     "TMIN_2M_ts" \
            "T_S_ts"  "RELHUM_2M_ts"       "RSTOM_ts"          "PS_ts"        "W_SO_ts" \
       "TOT_PREC_ts" "SFLDIR_PAR_ts" "SFLDIFD_PAR_ts" "SFLDIFU_PAR_ts"    "ZTRALEAV_ts" \
        "ZTRANGS_ts"     "ZVERBO_ts"   "SUR_PANFM_ts"  "SUR_PANDAY_ts" "SUR_BIOMASS_ts" \
        "SUR_LAI_ts"       "COSZ_ts"         "GPP_ts"         "NPP_ts"     "RS_LEAF_ts" \
        "SLA_SUN_ts"    "SLA_SHA_ts"     "LAI_SUN_ts"     "LAI_SHA_ts"      "VC_SUN_ts" \
         "VC_SHA_ts"   "LEAF_PSN_ts"     "PAR_SUN_ts"     "PAR_SHA_ts"                  )

paramscount=38


fn_name="1999_2015"
fn_research="2010_2015"

#========================= Paths ===============================================
# Upload correct filenames
version="Parc"
domain="parc_domain"
f_out="PARC"

# Create paths:
for (( i=0; i<=ver_count ; i++ ));
do
    echo 'The experiment version: '${exp_version[${i}]}
    echo 'The research territory: '$f_out

    DIR=/work/bb1112/b381275
    DIR_COSMO_DATA=${DIR}/"${version}"_"${exp_version[${i}]}"/post
    DIR_OUT=${DIR}/results/COSMO_exp/"${f_out}"_"${exp_version[${i}]}"
    DIR_RESULT=/pf/b/b381275/COSMO_results/"${f_out}"/"${exp_version[${i}]}"
    #========================= Generation operations ===========================

    for (( j=0; j<=paramscount ; j++ ));
    do
        echo ' Actual parameter is: '${params[${j}]}
        # The name of COSMO-CLM parameter for analysis
        fn_COSMO=${params[${j}]}."nc"
        # The paths to COSMO-CLM data (control run)
        ts_COSMO="${DIR_COSMO_DATA}"/*/"${fn_COSMO}"

        # The full timeseries with COSMO data from 1999 to 2020
        fCOSMO=${DIR_OUT}/"${params[${j}]}"_"${fn_name}"
        # The research timeseries with COSMO data from 2010 to 2015
        rCOSMO=${DIR_OUT}/"${params[${j}]}"_"${fn_research}"

        # Output data from the first step
        COSMO_s1=${DIR_OUT}/"${params[${j}]}"_"vert"_"${fn_research}"
        # Second step
        COSMO_s2=${DIR_OUT}/"${params[${j}]}"_"vert_mean"_"${fn_name}"

        # Path to result
        res_COSMO=${DIR_OUT}/"${params[${j}]}"_"mean"_"${fn_name}"


        # Main calculations
        ncrcat -h ${ts_COSMO} "${fCOSMO}"."nc"
        cdo seldate,20100101,20153112 "${fCOSMO}"."nc" "${rCOSMO}"."nc"

        if [ "${params[${j}]}" == "ALHFL_PL_ts" ] || \
           [ "${params[${j}]}" == "ZTRANG_ts"   ]; then
            cdo -sellevel,1/4 "${rCOSMO}"."nc" "${COSMO_s1}"."nc"
            cdo vertsum  "${COSMO_s1}"."nc" "${COSMO_s2}"."nc"
            # Calculations
            cdo fldmean "${COSMO_s2}"."nc" "${res_COSMO}"."nc"
            cdo outputts "${res_COSMO}"."nc" > "${res_COSMO}"."csv"

        elif [ "${params[${j}]}" == "W_SO_ts" ]; then
            cdo -sellevidx,1/2 "${rCOSMO}"."nc" "${COSMO_s1}"."nc"
            cdo vertsum  "${COSMO_s1}"."nc" "${COSMO_s2}"."nc"
            # Calculations
            cdo fldmean "${COSMO_s2}"."nc" "${res_COSMO}"."nc"
            cdo outputts "${res_COSMO}"."nc" > "${res_COSMO}"."csv"
        else
            # Calculation
            cdo fldmean "${rCOSMO}"."nc" "${res_COSMO}"."nc"
            cdo outputts "${res_COSMO}"."nc" > "${res_COSMO}"."csv"
        fi

    done
    cp -R ${DIR_OUT}/*."csv" ${DIR_RESULT}
    rm -r ${DIR_OUT}/*
done

