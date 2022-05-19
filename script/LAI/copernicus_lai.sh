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

#============================ COPERNICUS =======================================

dataset="COPERNICUS"
fname="LAI_2010_2015_GERMANY.nc"
#============================== General settings ===============================

# The research domain
domain=("PARC" "LINDEN" "LINDENBERG")
grid=("grid_parc.txt" "grid_linden.txt" "grid_lindenberg.txt")
# Count of domain --> should be the same as domain and coordinates parameters
paramscount=2


#========================= Paths ===============================================
DIR=/work/bb1112/b381275
DIR_IN=${DIR}/LAI/${dataset}/DATA
DIR_OUT=${DIR}/LAI/${dataset}

DIR_RESULT=/pf/b/b381275/COSMO_results/${dataset}"_LAI"/
#========================= Generation operations ===============================

for (( i=0; i<=paramscount ; i++ ));
do
    DIR_COPERNICUS_OUT=${DIR_OUT}/${domain[${i}]}

    # Domain grid --> grid for interpolation
    mygrid=${DIR_IN}/"${grid[${i}]}"

    # Define input file
    IN_DATA=${DIR_IN}/${fname}

    # Define output file
    OUT_DATA=${DIR_COPERNICUS_OUT}/"LAI_"${domain[${i}]}"_2010_2015"

    # Interpolation EOBS data to COSMO grid
    cdo remapbil,${mygrid} ${IN_DATA} ${OUT_DATA}."nc"

    cdo -outputtab,name,date,lon,lat,value -select,name=LAI ${OUT_DATA}."nc" > ${OUT_DATA}."csv"

    cp -R ${OUT_DATA}."csv" ${DIR_RESULT}
    cp -R ${OUT_DATA}."nc"  ${DIR_RESULT}

    # Delete results from the work folder
    rm -r ${DIR_COPERNICUS_OUT}/*
done

echo 'All done MODIS'







