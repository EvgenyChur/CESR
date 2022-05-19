!+ Data module for all parametric data in the vegetation modules
!------------------------------------------------------------------------------

MODULE sfc_phenology_data

!------------------------------------------------------------------------------
!
! Description:
!  This module declares and initializes all constant data for the new
!  stomatal resistance Ball-Berry approach, leaf photosynthesis and
!  evolution of leaf area index based on biomass evolution.
!
!  In module there are:
!  - general constants for vegetation and photosyntesis schemes:
!  - CLMv3.5 parameters for C3 and C4 grass (PFT parameters):
!
! Center for Environmental Systems Research, 2020 - 2022
! Evgeny Churiulin, Merja Toelle, Jürgen Helmert, Jean-Marie Bettems
! phone:  +49(170)261-51-04
! email:  evgenychur@uni-kassel.de, merja.toelle@uni-kassel.de,
!         Juergen.Helmert@dwd.de Jean-Marie.Bettems@meteoswiss.ch
!
! Acknowledge: Vladimir Kopeykin
!
! Current Code Owner: DWD, Juergen Helmert
!  phone:  +49  69  8062 2704
!  fax:    +49  69  8062 3721
!  email:  Juergen.Helmert@dwd.de
!
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! CLMv6.0    2022/03/02 Evgenii Churiulin
!  Initial release
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================
!
! Declarations:
!
! Modules used:
#ifdef __COSMO__
  USE kind_parameters, ONLY : vpp      ! KIND-type parameter for real variables (variable precision physics)
#endif

#ifdef __ICON__
  USE mo_kind        , ONLY : vpp      ! KIND-type parameter for real variables (variable precision physics)
#endif

!==============================================================================

IMPLICIT NONE

! All constants and variables in this module are public
PUBLIC

    !==========================================================================
    ! Declarations
    !
    ! Global (i.e. public) Declarations:
    ! 1. Constant variables for the new vegetation algorithms
    ! -------------------------------------------------------------------------
    REAL (KIND = vpp), PARAMETER ::               &
        ! Leaf boundary layer resistance algorithm
        TTC_CAN    = 0.01_vpp                   , & ! Turbulent transfer coefficient
                                                    ! between the canopy surface and canopy air       (m s**-0.5 )
        ! Available water in soil
        VEG_PAR    = 1.0_vpp                    , & ! Vegetation dependent PARAMETER
        ! Vcmax25 algorithm
        ALPHA_R_25 = 60.0_vpp                   , & ! Specific active of Rubisco                      (umol/gRubisco/s)
        ALPHA_VMAX = 2.4_vpp                    , & ! Temperature sensitive parameter
        F_NR       = 7.16_vpp                   , & ! Mass ration of total Rubisco molecular
                                                    ! mass to nitrogen in Rubisco                     (gRubisco/gN in Rubisco)
        ! Stomatal resistance algorithms
        RSMAX0     = 2.e4_vpp                   , & ! Maximum stomatal resistance                      (s/m)
        GLMIN      = 2000.0_vpp                 , & ! Minimum leaf conductance                         (umol/m2/s)
        ! Michaelis-Menten constants
        K_C25      = 30.0_vpp                   , & ! Concentration values of CO2 at 25 °C             (Pa)
        ALPHA_KC   = 2.1_vpp                    , & ! Temperature sensitive parameter for CO2
        K_O25      = 30000.0_vpp                , & ! Concentration values of O2 at 25 °C              (Pa)
        ALPHA_KO   = 1.2_vpp                    , & ! Temperature sensitive parameter for  O2
        ! Concentration of O2 and CO2 in air
        P_STD      = 101325.0_vpp               , & ! Standart pressure
        CON_O2     = 0.2_vpp * P_STD            , & ! O2 concentration                                 (Pa)
        CON_CO2    = 355_vpp * 1.e-6_vpp * P_STD, & ! CO2 concentration                                (Pa)
        ! Resperation algorithm
        B_RATE     = 2.525e-6_vpp               , & ! base rate for maintenance respiration
        Q_10       = 2.0_vpp                        ! temperature dependence

    ! 2. Parameters for C3 and C4 PFTs - (from CLM v3.5 model)
    !---------------------------------------------------------------------------
    ! Constants varibles depending on different PFT types:
    !---------------------------------------------------------------------------
    INTEGER, PARAMETER ::  &
        pft_type     = 2,  & ! current numbers of PFT (1 - C3 grass, 2 - C4 grass)
        param        = 48, & ! numbers of PFT parameters

        ! Parameters values for Table 1.
        R_ZOM        = 1,  & ! Ratio of momentum roughness length to canopy top height
        R_DISP       = 2,  & ! Ratio of displacement height to canopy top height
        D_LEAF       = 3,  & ! Characteristic dimension of the leaves in
                             ! the direction of wind flow
        C3_PSN       = 4,  & ! Photosynthetic pathway: 0. = c4, 1. = c3
        V_CMX25      = 5,  & ! Maximum rate of carboxylation at 25C                           (umol CO2/m**2/s)
        MP           = 6,  & ! Slope of conductance-to-photosynthesis relationship
        ALPHA        = 7,  & ! Quantum efficiency  at 25C                                     (μmol CO2 per μmol photon)
        ALPHA_vl     = 8,  & ! Weighted combination of the leaf reflectances                  (   VIS  )
        ALPHA_nl     = 9,  & ! Weighted combination of the leaf reflectances                  (   NIR  )
        ALPHA_vs     = 10, & ! Weighted combination of the stem reflectances                  (   VIS  )
        ALPHA_ns     = 11, & ! Weighted combination of the stem reflectances                  (   NIR  )
        TETA_vl      = 12, & ! Weighted combination of the leaf transmittances                (   VIS  )
        TETA_nl      = 13, & ! Weighted combination of the leaf transmittances                (   NIR  )
        TETA_vs      = 14, & ! Weighted combination of the stem transmittances                (   VIS  )
        TETA_ns      = 15, & ! Weighted combination of the stem transmittances                (   NIR  )
        X_l          = 16, & ! Departure of leaf angles from a random distribution
        ROOTA_PAR    = 17, & ! Rooting distribution parameter                                 (   1/m  )
        ROOTB_PAR    = 18, & ! Rooting distribution parameter                                 (   1/m  )
        SLA_o        = 19, & ! Value for SLA at the top of canopy                             (  m2/gC )
        SLA_m        = 20, & ! Linear slope coefficient
        CN_L         = 21, & ! Leaf carbon to nitrogen ratio (leaf C:N)                       (  gC/gN )
        F_LNR        = 22, & ! Fraction of leaf nitrogen in Rubisco                           ( no units)
        PSI_o        = 23, & ! Soil water potential when stomata fully open                   (    mm  )
        PSI_c        = 24, & ! Soil water potential when stomata fully close                  (    mm  )
        F_N          = 25, & ! Nitrogen availability factor                                   (  m2/gC )
        WOODY        = 26, & ! Binary flag for woody lifeform: 1. = woody, 0. = not woody
        LF_LIT_CN    = 27, & ! Leaf litter C:N                                                (  cG/gN )
        F_ROOT_CN    = 28, & ! Fine root C:N                                                  (  gC/gN )
        LIVEWD_CN    = 29, & ! Live wood (phloem and ray parenchyma) C:N                      (  gC/gN )
        DEADWD_CN    = 30, & ! Dead wood (xylem and heartwood) C:N                            (  gC/gN )
        F_ROOT_LEAF  = 31, & ! Allocation parameter: new fine root C per new leaf C           (  gC/gC )
        STEM_LEAF    = 32, & ! Allocation parameter: new stem c per new leaf C                (  gC/gC )
        C_ROOT_STEM  = 33, & ! Allocation parameter: new coarse root C per new stem C         (  gC/gC )
        F_LIVE_WD    = 34, & ! Allocation parameter: fraction of new wood that is live        ( no units)
        F_CUR        = 35, & ! Allocation parameter: fraction of allocation
                             ! that goes to currently displayed growth, remainder
                             ! to storage
        LF_FLAB      = 36, & ! Leaf litter labile fraction
        LF_FCEL      = 37, & ! Leaf litter cellulose fraction
        LF_FLIG      = 38, & ! Leaf litter lignin fraction
        FR_FLAB      = 39, & ! Fine root litter labile fraction
        FR_FCEL      = 40, & ! Fine root litter cellulose fraction
        FR_FLIG      = 41, & ! Fine root litter lignin fraction
        DW_FCEL      = 42, & ! Dead wood cellulose fraction
        DW_FLIG      = 43, & ! Dead wood lignin fraction
        LEAF_LONG    = 44, & ! Leaf longevity                                                  (  yrs   )
        EVERGREEN    = 45, & ! Binary flag for evergreen leaf habit                            ( 0 or 1 )
        STRESS_DECID = 46, & ! Binary flag for stress-deciduous leaf habit                     ( 0 or 1 )
        SEASON_DECID = 47, & ! Binary flag for seasonal-deciduous leaf habit                   ( 0 or 1 )
        RESIST       = 48    ! Fire resistance index                                           (unitless)

    ! Table 1: pft_aero_prop - Plant functional type aerodynamic parameters
    !                      C3_grass      C4_grass          name in CLMv3.5
    REAL (KIND=vpp), DIMENSION(pft_type, param), PARAMETER ::  &
        pft_CN_par = (/(/  0.120_vpp ,   0.120_vpp  /), & !  1 - z0mr
                       (/  0.68_vpp  ,   0.68_vpp   /), & !  2 - displar
                       (/  0.04_vpp  ,   0.04_vpp   /), & !  3 - dleaf
                       (/  1.0_vpp   ,   0.0_vpp    /), & !  4 - c3psn
                       (/ 52.0_vpp   ,  52.0_vpp    /), & !  5 - vcmx25
                       (/  9.0_vpp   ,   5.0_vpp    /), & !  6 - mp
                       (/  0.06_vpp  ,   0.04_vpp   /), & !  7 - qe25
                       (/  0.11_vpp  ,   0.11_vpp   /), & !  8 - rhol_vis
                       (/  0.35_vpp  ,   0.35_vpp   /), & !  9 - rhol_nir
                       (/  0.31_vpp  ,   0.31_vpp   /), & ! 10 - rhos_vis
                       (/  0.53_vpp  ,   0.53_vpp   /), & ! 11 - rhos_nir
                       (/  0.05_vpp  ,   0.05_vpp   /), & ! 12 - taul_vis
                       (/  0.34_vpp  ,   0.34_vpp   /), & ! 13 - taul_nir
                       (/  0.120_vpp ,   0.120_vpp  /), & ! 14 - taus_vis
                       (/  0.250_vpp ,   0.250_vpp  /), & ! 15 - taus_nir
                       (/ -0.30_vpp  ,  -0.30_vpp   /), & ! 16 - xl
                       (/ 11.0_vpp   ,  11.0_vpp    /), & ! 17 - roota_par
                       (/  2.0_vpp   ,   2.0_vpp    /), & ! 18 - rootb_par
                       (/  0.030_vpp ,   0.030_vpp  /), & ! 19 - slasun
                       (/  0.0_vpp   ,   0.0_vpp    /), & ! 20 - dsladlai
                       (/ 25.0_vpp   ,  25.0_vpp    /), & ! 21 - leafcn
                       (/  0.09_vpp  ,   0.09_vpp   /), & ! 22 - flnr
                       (/ -0.74e5_vpp,  -0.74e5_vpp /), & ! 23 - smpso
                       (/ -2.75e5_vpp,  -2.75e5_vpp /), & ! 24 - smpsc
                       (/  0.61_vpp  ,   0.64_vpp   /), & ! 25 - fnitr
                       (/  0.0_vpp   ,   0.0_vpp    /), & ! 26 - woody
                       (/ 50.0_vpp   ,  50.0_vpp    /), & ! 27 - lflitcn
                       (/ 42.0_vpp   ,  42.0_vpp    /), & ! 28 - frootcn
                       (/  0.0_vpp   ,   0.0_vpp    /), & ! 29 - livewdcn
                       (/  0.0_vpp   ,   0.0_vpp    /), & ! 30 - deadwdcn
                       (/  3.0_vpp   ,   3.0_vpp    /), & ! 31 - froot_leaf
                       (/  0.0_vpp   ,   0.0_vpp    /), & ! 32 - stem_leaf
                       (/  0.0_vpp   ,   0.0_vpp    /), & ! 33 - croot_stem
                       (/  0.0_vpp   ,   0.0_vpp    /), & ! 34 - flivewd
                       (/  0.5_vpp   ,   0.5_vpp    /), & ! 35 - fcur
                       (/  0.25_vpp  ,   0.25_vpp   /), & ! 36 - lf_flab
                       (/  0.50_vpp  ,   0.50_vpp   /), & ! 37 - lf_fcel
                       (/  0.25_vpp  ,   0.25_vpp   /), & ! 38 - lf_flig
                       (/  0.25_vpp  ,   0.25_vpp   /), & ! 39 - fr_flab
                       (/  0.50_vpp  ,   0.50_vpp   /), & ! 40 - fr_fcel
                       (/  0.25_vpp  ,   0.25_vpp   /), & ! 41 - fr_flig
                       (/  0.75_vpp  ,   0.75_vpp   /), & ! 42 - dw_fce
                       (/  0.25_vpp  ,   0.25_vpp   /), & ! 43 - dw_flig
                       (/  1.0_vpp   ,   1.0_vpp    /), & ! 44 - leaf_long
                       (/  0.0_vpp   ,   0.0_vpp    /), & ! 45 - evergreen
                       (/  1.0_vpp   ,   1.0_vpp    /), & ! 46 - stress_decid
                       (/  0.0_vpp   ,   0.0_vpp    /), & ! 47 - season_decid
                       (/  0.12_vpp  ,   0.12_vpp   /)/)  ! 48 - resist

    ! 3. Parameters for C3 and C4 grass from ISBA-Ag-s (SURFEX)
    ! -------------------------------------------------------------------------
    ! Constants varibles depending on different PFT types:
    !---------------------------------------------------------------------------

    INTEGER , PARAMETER ::  &
        pft_param = 2,      & ! current numbers of PFT (1 - C3 grass, 2 - C4 grass)
        par       = 7,      & ! numbers of PFT parameters

        ! Parameters values for Table 1.
        PSEFOLD   = 1,      & ! Maximum effective life expectancy                              (days )
        LAIMIN    = 2,      & ! Minimum leaf area index                                        (m2/m2)
        DMAX      = 3,      & ! Maximum specific humidity                                      (kg/kg)
        FZERO     = 4,      & ! Vegetation parameter                                           ( --- )
        GC        = 5,      & ! Cuticular conductance                                          ( m/s )
        GMES      = 6,      & ! Mesophyll conductance                                          ( m/s )
        PBSLAI    = 7         ! Ratio of biomass/lai

    ! Table 2: pft_aero_prop - Plant functional type aerodynamic parameters
    !                       C3_grass      C4_grass
    REAL (KIND=vpp), DIMENSION(pft_param,  par), PARAMETER :: &
        pft_SURFEX = (/(/ 150.0_vpp    , 150.0_vpp     /),    & ! 1 - teta_m
                       (/   0.3_vpp    ,   0.3_vpp     /),    & ! 2 - LAImin
                       (/   0.05_vpp   ,   0.033_vpp   /),    & ! 3 - Dmax;
                       (/   0.95_vpp   ,   0.6_vpp     /),    & ! 4 - f0
                       (/   0.00025_vpp,   0.00015_vpp /),    & ! 5 - Gc
                       (/   0.001_vpp  ,   0.009_vpp   /),    & ! 6 - Gm*
                       (/   0.06_vpp   ,   0.06_vpp    /)/)     ! 7 - alpha


    ! 4. Additional local parameters for new vegetation subroutines
    ! --------------------------------------------------------------------------

    ! Data for get_stomatal_grid and stomatal subroutines
    INTEGER                     :: &
        C3_grass = 1             , & ! Column index for C3 grass in tables with constant
        C4_grass = 2             , & ! Column index for C4 grass in tables with constant
        niter    = 3                 ! number of iterations for leaf photosynthesis

    ! Data for get_stomatal_data, stomatal, respiration subroutines
    REAL  (KIND=vpp), PARAMETER :: &
        MPE        = 1.e-6_vpp   , & ! prevents overflow error if division by zero  -> stomatal
        GPPFACT    = 0.5_vpp     , & ! time dicesion                                -> respiration
        H2OCAN     = 0.0_vpp     , & ! Temporal parameter for water content
        LEAFC      = 25.0_vpp    , & ! Leaf C                                       (kgC/m2)
        DEADSTEMC  = 5.5_vpp     , & ! Dead stem C                                  (kgC/m2)
        FORC_HGT_U = 30.0_vpp    , & ! Observational height of wind                 (   m  )
        TAPER      = 200.0_vpp   , & ! Ratio of height:radius breast height         (tree allometry)

        ! calculated parameters
        STOCKING   = 1000.0_vpp / 10000.0_vpp         , & ! Stocking density          (stems/m2)
        DWOOD      = 500.0_vpp * 1000.0_vpp * 0.5_vpp     ! Wood density              (  gC /m3)

    ! Data for biomass evolution
    REAL  (KIND=vpp), PARAMETER :: &
        XMC    = 12.0E-3_vpp     , & ! Molecular weights of carbon (Mc)
        XMCO2  = 44.0E-3_vpp     , & ! Molecular weights of CO2  (Mco2)
        XPCCO2 = 0.40_vpp        , & ! Proportion of carbon in the dry biomass
        M_CO2  = 44.0E-3_vpp     , & ! Molar mass of CO2                            (kg/mol)
        XDAY   = 86400_vpp       , & ! Day duration                                 (s)
        ! Recalculation coefficient for biomass processes
        zbmcoef = XMC / (XMCO2 * XPCCO2)

    REAL  (KIND=vpp), PARAMETER :: &
        BOLTZ  = 1.38065E-23_vpp  , & ! Boltzmann constant      (J/K)
        AVOGAD = 6.02214E26_vpp   , & ! Avogadro's number       (mol-1)
        RGAS   = avogad * boltz      ! Universal gas constant  (J/K/kmol)


    LOGICAL                    ::  &
        midnight                          ! logical parameter for midnight

#ifdef __COSMO__
    ! 4. vpp-equivalents for model parameters
    ! --------------------------------------------------------------------------

        ! These variables are from data_constants, where they are defined with KIND=wp
        ! For ease of use they have the same names, but are of type vpp, used in the
        ! surface schemes. They are set in sfc_interface, subroutine sfc_init

        REAL    (KIND = vpp) ::  &
                        pi     , & ! Circle constant
                        t0_melt    ! Freezing temperature of fresh water             ( K )
#endif

#ifdef ALLOC_WKARR
! Local arrays defined as allocatables here:
! ------------------------------------------------------------------------------

    REAL(KIND = vpp), ALLOCATABLE  ::  &
        ! variables for get_sun_data subroutine
        sun_decl        (:)          , & ! declanation of the sun
        dyl             (:)          , & ! actual daylength
        dyl_max         (:)          , & ! actual maximum daylength
        ! variables for get_stomatal_data subroutine
        elai            (:)          , & ! The one-sided leaf area index with burying by snow ( m2/m2)
        esai            (:)          , & ! One-sided stem area index with burying by snow     ( m2/m2)
        ! variables for two-big-leaf calclulations
        vai             (:)          , & ! Total leaf area index + stem area index            ( m2/m2)
        f_sun           (:)          , & ! Fraction of sunlit canopy                          (0 to 1)
        ! variables for biomass_evolution subroutine
        spanmax         (:)          , & ! maximum photosynthesis rate in optimal conditions
        spsn            (:)              ! values of photosynthesis (convert psn)             (kg CO2 m-2 s-1)

!==============================================================================

CONTAINS

SUBROUTINE phenology_wkarr_alloc (nproma, istat)

    INTEGER, INTENT(IN)  :: nproma
    INTEGER, INTENT(OUT) :: istat

    istat = 0

    ALLOCATE (                               &
              ! Get_sun_data subroutine - variables
              sun_decl        (nproma)     , & ! sun declanation
              dyl             (nproma)     , & ! actual daylength
              dyl_max         (nproma)     , & ! actual maximum daylength
              ! Get_stomatal_data subroutine - parameters
              elai            (nproma)     , & ! one-sided LAI with burying by snow
              esai            (nproma)     , & ! one-sided SAI (stem) with burying by snow
              vai             (nproma)     , & ! Total SAI + SAI
              f_sun           (nproma)     , & ! Fraction of sunlit canopy
              ! Biomass_evolution subroutine - parameters
              spanmax         (nproma)     , & ! maximum photosynthesis rate in optimal conditions
              spsn            (nproma)     , & ! values of photosynthesis (convert psn)
       STAT=istat)

    END SUBROUTINE phenology_wkarr_alloc


SUBROUTINE phenology_wkarr_dealloc (istat)

    INTEGER, INTENT(OUT) :: istat

    istat = 0

    DEALLOCATE (                                                           &
                ! Get_sun_data subroutine - parameters
                sun_decl, dyl     , dyl_max ,                              &
                ! Get_stomatal_data subroutine - parameters
                elai    , esai    , vai     , f_sun       ,                &
                ! Biomass_evolution subroutine - parameters
                spanmax , spsn    ,                                        &
            STAT=istat)

END SUBROUTINE terra_wkarr_dealloc

!===============================================================================
#endif
         !ALLOC_WKARR
!===============================================================================
END MODULE sfc_phenology_data

