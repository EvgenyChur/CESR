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

#============================ MODIS ============================================
years=("2010" "2011" "2012" "2013" "2014" "2015")
yearcount=5

dataset="MODIS"
com_name="_LAI_FPAR_LPDAAC_"

# The name of dataset
#ifilename=${dataset}""${year}         #  EOBS dataset
#ofilename="tas_eobs_5_2010_2015.nc"

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

DIR_RESULT=/pf/b/b381275/COSMO_results/${dataset}"_LAI"
#========================= Generation operations ===============================

for (( i=0; i<=paramscount ; i++ ));
do
    DIR_MODIS_OUT=${DIR_OUT}/${domain[${i}]}

    for (( j=0; j<=yearcount ; j++ ));
    do
        # Domain grid --> grid for interpolation
        mygrid=${DIR_IN}/"${grid[${i}]}"

        # Define filename
        fname=${dataset}${com_name}${years[${j}]}.nc
 
        # Define input file
        inMODIS=${DIR_IN}/${fname}
        #echo ' MoDIS in: '${inMODIS}

        # Define output file
        outMODIS_lai=${DIR_MODIS_OUT}/"LAI_ONLY_"${fname}

        outMODIS=${DIR_MODIS_OUT}/${fname}
        #echo ' MoDIS out: '${outMODIS}

        cdo select,name=lai ${inMODIS} ${outMODIS_lai}

        # Interpolation EOBS data to COSMO grid
        cdo remapbil,${mygrid} ${outMODIS_lai} ${outMODIS}

        rm ${outMODIS_lai}
    done

    name_MODIS="LAI_"${domain[${i}]}"_2010_2015"

    res_MODIS_nc=${DIR_MODIS_OUT}/${name_MODIS}.nc
    res_MODIS_csv=${DIR_MODIS_OUT}/${name_MODIS}.csv

    cdo mergetime ${DIR_MODIS_OUT}/*.nc ${res_MODIS_nc}
    cdo -outputtab,name,date,lon,lat,value  ${res_MODIS_nc} > ${res_MODIS_csv}

    cp -R ${res_MODIS_csv} ${DIR_RESULT}
    cp -R ${res_MODIS_nc} ${DIR_RESULT}

    # Delete results from the work folder
    rm -r ${DIR_MODIS_OUT}/*
done



echo 'All done MODIS'







