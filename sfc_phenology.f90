!+ Source module  "src_phenology"
!------------------------------------------------------------------------------

MODULE sfc_phenology

!------------------------------------------------------------------------------
!! Description:
! This module contains the new vegetation algorithms for:
!      1. get_stomatal_grid - Creating a special grid with C3 and C4 grass
!
!      2. get_sun_data - Calculating additional solar parameters for stomata
!                        subroutine
!
!      3. get_stomatal_data - Calculating additional physical parameters for
!                             stomata subroutine. Also, in this subroutine there
!                             is algorithm for calculation two-big-leaf parameters
!
!      4. stomata - Calculating stomatal resistance and leaf photosynthesis
!                   parameters based on the Ball-Berry approach (stomatal
!                   resistance) and Farquhar and Collatz models (leaf
!                   photosynthesis)
!
!      5. respiration - Calculating GPP and NPP fluxes
!
!      6. biomass_evolution - Calculating changes in leaf area index (LAI)
!                             depending on biomass evolution

!
! Center for Environmental Systems Research, 2020 - 2022
! Evgeny Churiulin, Merja Toelle,
! phone:  +49(170)261-51-04
! email:  evgenychur@uni-kassel.de, merja.toelle@uni-kassel.de,
!         Juergen.Helmert@dwd.de Jean-Marie.Bettems@meteoswiss.ch
!
! Acknowledge: Vladimir Kopeykin , Juergen Helmert,
!              Jean-Marie Bettems, Markus Uebel
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
! CLM_6.0    2022/03/30 Evgenii Churiulin
!  Initial release
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================!
!
! Declarations:
!
! Modules used:
#ifdef __COSMO__
USE kind_parameters, ONLY :   &
    vpp           ! KIND-type parameter for real variables (variable precision physics)

! end of kind_parameters
!------------------------------------------------------------------------------!

USE data_constants,  ONLY :   &
    ! 1. physical constants and related variables
    ! -------------------------------------------
    pi_wp      => pi        , & ! circle constant
    t0_melt_wp => t0_melt       ! absolute zero for temperature

! end of data_constants
!------------------------------------------------------------------------------!

USE data_runcontrol, ONLY :   &
    ! 3. controlling the physics
    ! -------------------------------------------
    itype_root              , & ! type of root density distribution
    itype_tran              , & ! type of surface to tamospher transfer

    ! 5. additional control variables
    ! -------------------------------------------
    ntstep                  , & ! actual time step
    yakdat1                 , & ! actual date
    itype_calendar              ! for specifying the calendar used

! end of data_runcontrol
!------------------------------------------------------------------------------!

USE data_io,         ONLY :   &
    ydate_ini                   ! start of the forecast

! end of data_io
!------------------------------------------------------------------------------!

USE utilities,       ONLY:    &
    get_utc_date                 ! subroutine for calculating the actual day

! end of utilities
!------------------------------------------------------------------------------!

USE sfc_phenology_data  ! All variables from this data module are used by
                        ! this module.

!end of sfc_phenology_data
!------------------------------------------------------------------------------!
#endif

!==============================================================================

!------------------------------------------------------------------------------
! Declarations
!------------------------------------------------------------------------------

IMPLICIT NONE

PUBLIC             ! All constants and variables in this module are public

!==============================================================================

CONTAINS

    !===========================================================================
    ! The new subroutines for calculating vegetation parameters
    !---------------------------------------------------------------------------

    SUBROUTINE get_stomatal_grid(         &
                       nvec             , & ! array dimensions
                       ivstart          , & ! start index for computations in the parallel program
                       ivend            , & ! end index for computations in the parallel program
                       soiltyp_subs     , & ! type of the soil (keys 0-9)
                       zep_s            , & ! potential evaporation for t_s
                       pft                & ! Type of PFT
                                          )
        ! Creating a special grid with C3 and C4 grass
        !-----------------------------------------------------------------------

        ! Modules used:
        !-----------------------------------------------------------------------

        IMPLICIT NONE
        !----------------------------------------------------------------------!
        ! Declarations
        !----------------------------------------------------------------------!
        INTEGER, INTENT(IN)  ::      &
            nvec                   , & ! array dimensions
            ivstart                , & ! start index for computations in the parallel program
            ivend                      ! end index for computations in the parallel program

        INTEGER, DIMENSION(nvec), INTENT(IN)    :: &
            soiltyp_subs               ! type of the soil (keys 0-9)

        REAL (KIND = vpp), DIMENSION(nvec), INTENT(IN)    :: &
            zep_s                     ! potential evaporation for t_s

        INTEGER, DIMENSION(nvec), INTENT(INOUT) :: &
            pft                       ! type of PFT (1 - C3, 2 - C4)

        !======================================================================!
        ! Local scalars:
        ! ---------------------------------------------------------------------!
        INTEGER ::  &
        i           ! loop index in x-direction

        !- End of header
        !=======================================================================

        !-----------------------------------------------------------------------
        ! Begin subroutine: get_stomatal_grid
        !-----------------------------------------------------------------------

        ! Calculation of special grid for actual PFT grass type
        !-----------------------------------------------------------------------

        DO i = ivstart, ivend
#ifndef MESSY
            IF (soiltyp_subs(i) >= 3 .AND. zep_s(i) < 0.0_vpp) THEN  ! upwards directed potential evaporation
#else
            ! calculate this part always for MESSy
            IF (soiltyp_subs(i) >= 3) THEN ! < 3 is ice and rock
#endif
                pft(i) = C3_grass
            ENDIF
        ENDDO

        RETURN
    !===========================================================================
    END SUBROUTINE get_stomatal_grid





    SUBROUTINE get_sun_data(         &
                       nvec             , & ! array dimensions
                       ivstart          , & ! start index for computations in the parallel program
                       ivend            , & ! end index for computations in the parallel program
                       dt               , & ! time step
                       soiltyp_subs     , & ! type of the soil (keys 0-9)
                       zep_s            , & ! potential evaporation for t_s
                       rlat             , & ! geographical latitude                                  ( rad  )
                       sun_el           , & ! sun elevation angle                                    ( deg  )
                       day_length       , & ! daylength                                              (hours )
                       f_dyl              & ! daylength limitation function                          (0 to 1)
                                          )
        ! Calculating additional solar parameters for stomatal subroutine
        !-----------------------------------------------------------------------

        ! Modules used:
        !-----------------------------------------------------------------------

        IMPLICIT NONE
        !----------------------------------------------------------------------!
        ! Declarations
        !----------------------------------------------------------------------!
        INTEGER, INTENT(IN)  ::      &
            nvec                   , & ! array dimensions
            ivstart                , & ! start index for computations in the parallel program
            ivend                      ! end index for computations in the parallel program

        REAL (KIND = vpp), INTENT(IN)  ::  &
            dt                         ! time step

        INTEGER, DIMENSION(nvec), INTENT(IN)    :: &
            soiltyp_subs               ! type of the soil (keys 0-9)

        REAL (KIND = vpp), DIMENSION(nvec), INTENT(IN)    :: &
            zep_s                  , & ! potential evaporation for t_s
            rlat                   , & ! geographical latitude
            sun_el                     ! sun elevation angle

        REAL (KIND = vpp), DIMENSION(nvec), INTENT(INOUT) :: &
            day_length             , & ! daylength
            f_dyl                      ! daylength limitation function

        !======================================================================!
        ! Local scalars:
        ! ---------------------------------------------------------------------!
        INTEGER           ::         &
            i                      , & ! loop index in x-direction
            nactday                    ! day of the year

        REAL(KIND = vpp)  ::         &
            acthour                , & ! actual hour of the day
            jj                     , & ! actual year
            ztwo                   , & ! year regulation
            ztho                   , & ! angle of year
            ms_decl                , & ! maximum solar declination                         (  rad  )
            ! Temporal variables:
            temp_dyl               , & ! temporal daylength
            temp_dyl_max               ! temporal maximum daylength

        CHARACTER (LEN=14) yactdate1   ! actual date in the form yyyymmddhhmmss
        CHARACTER (LEN=28) yactdate2   ! actual date in the form &
                                       ! wd   dd.mm.yy  hh mm ss UTC

#ifndef ALLOC_WKARR
! Local (automatic) arrays:
! -------------------------
        REAL(KIND = vpp)  ::         &
            sun_decl        (nvec) , & ! declanation of the sun
            dyl             (nvec) , & ! actual daylength
            dyl_max         (nvec)     ! max temporal and actual daylength
#endif
        !- End of header
        !=======================================================================

        !-----------------------------------------------------------------------
        ! Begin subroutine: get_sun_data
        !-----------------------------------------------------------------------

#ifdef __COSMO__
! Double precision equivalents of module variables
        pi = REAL(pi_wp, vpp)
#endif

        ! Calculation of daylength and solar declination angle
        !-----------------------------------------------------------------------
        CALL get_utc_date(ntstep, ydate_ini, dt, itype_calendar, &
                          yactdate1, yactdate2, nactday, acthour )

        DO i = ivstart, ivend
#ifndef MESSY
            IF (soiltyp_subs(i) >= 3 .AND. zep_s(i) < 0.0_vpp) THEN  ! upwards directed potential evaporation
#else
            ! calculate this part always for MESSy
            IF (soiltyp_subs(i) >= 3) THEN ! < 3 is ice and rock
#endif
                READ (yactdate1(1:4),'(I4)') jj
                ztwo = 0.681_vpp + 0.2422_vpp * (jj - 1949) - (jj - 1949) / 4.0_vpp
                ztho = 2.0_vpp * pi * (REAL(nactday, vpp) - 1.0_vpp + ztwo) / 365.2422_vpp

                sun_decl(i) = 0.006918_vpp - 0.399912_vpp * COS(ztho) + 0.070257_vpp * SIN(ztho)           - &
                              0.006758_vpp * COS(2.0_vpp  * ztho)     + 0.000907_vpp * SIN(2.0_vpp * ztho) - &
                              0.002697_vpp * COS(3.0_vpp  * ztho)     + 0.001480_vpp * SIN(3.0_vpp * ztho)

                day_length(i) = 24 / pi * ACOS(-TAN(rlat(i)) * TAN(sun_decl(i)))

                !---------------------------------------------------------------
                ! Calculate daylength limitation function
                !---------------------------------------------------------------

                !IF (rlat >= 0.0) THEN
                ms_decl  = 23.4667_vpp * (pi / 180.0_vpp)
                !ELSE
                !    ms_decl  = (-1.0 * 23.4667) * D_to_R
                !ENDIF

                temp_dyl     = (-1.0_vpp * sin(rlat(i)) * sin(sun_decl(i))) / (cos(rlat(i)) * cos(sun_decl(i)))
                temp_dyl_max = (-1.0_vpp * sin(rlat(i)) * sin(ms_decl))     / (cos(rlat(i)) * cos(ms_decl))

                temp_dyl     = min(1.0_vpp, max(-1.0_vpp, temp_dyl))
                temp_dyl_max = min(1.0_vpp, max(-1.0_vpp, temp_dyl_max))

                dyl(i)       = 2.0_vpp * 13750.9871_vpp * acos(temp_dyl)
                dyl_max(i)   = 2.0_vpp * 13750.9871_vpp * acos(temp_dyl_max)

                ! calculate dayl_factor as the ratio of (current:max dayl)^2
                f_dyl(i) = min(1.0_vpp, max(0.01_vpp, (dyl(i) * dyl(i)) / (dyl_max(i) * dyl_max(i))))
            ENDIF
        ENDDO

        RETURN
    !===========================================================================
    END SUBROUTINE get_sun_data




    SUBROUTINE get_stomatal_data(         &
                       nvec             , & ! array dimensions
                       ivstart          , & ! start index for computations in the parallel program
                       ivend            , & ! end index for computations in the parallel program
                       soiltyp_subs     , & ! type of the soil                                        (keys 0-9)
                       pft              , & ! type of PFT                                             (1 - C3, 2 - C4)
                       zep_s            , & ! potential evaporation for t_s
                       zdz_snow         , & ! snow depth                                              ( m   )
                       ustar_fv         , & ! friction velocity (ustar)                               ( m/s )
                       ps               , & ! surface pressure                                        ( Pa  )
                       qv_s             , & ! specific humidity at the surface                        (kg/kg)
                       grad_zqvlow      , & ! specific humidity of lowest atmospheric layer           (kg/kg)
                       lai              , & ! leaf area index                                         (m2/m2)
                       cos_zen_ang      , & ! actual cosine of zenit angle                            (rad  )
                       sfldir_par       , & ! direct component of PAR at the ground                   (W/m2 )
                       sfldifd_par      , & ! diffuse downward component of PAR at the ground         (W/m2 )
                       sfldifu_par      , & ! diffuse upward component of PAR at the ground           (W/m2 )
                       fwet             , & ! Fraction of wet canopy                                  (0 to 1)
                       fdry             , & ! Fraction of dry canopy                                  (0 to 1)
                       wtl              , & ! Heat conductance for leaf                               (m/s  )
                       rs_leaf          , & ! canopy leaf resistance                                  (s/m  )
                       par_sun          , & ! photosynthetic active radiation for sunlit leaf         (W/m2 )
                       par_sha          , & ! photosynthetic active radiation for shaded leaf         (W/m2 )
                       ea               , & ! Vapor pressure of air in plant canopy                   (pa   )
                       lai_sun          , & ! sunlit leaf area index                                  (m2/m2)
                       lai_sha          , & ! shaded leaf area index                                  (m2/m2)
                       sla_sun          , & ! specific leaf area index for sunlit leaves              (m2/m2)
                       sla_sha            & ! specific leaf area index for shaded leaves              (m2/m2)
                                          )

        ! Calculating additional physical parameters for stomata subroutine
        ! also in this subroutine there is algorithm for calculation two-big-leaf
        ! approach
        !-----------------------------------------------------------------------

        ! Modules used:
        !-----------------------------------------------------------------------

        IMPLICIT NONE
        !----------------------------------------------------------------------!
        ! Declarations
        !----------------------------------------------------------------------!
        INTEGER, INTENT(IN)  ::      &
            nvec                   , & ! array dimensions
            ivstart                , & ! start index for computations in the parallel program
            ivend                      ! end index for computations in the parallel program

        INTEGER, DIMENSION(nvec), INTENT(IN)    :: &
            soiltyp_subs           , & ! type of the soil (keys 0-9)
            pft                        ! type of PFT

        REAL (KIND = vpp), DIMENSION(nvec), INTENT(IN)    :: &
            zep_s                  , & ! potential evaporation for t_s
            zdz_snow               , & ! snow depth
            ustar_fv               , & ! friction velocity
            ps                     , & ! surface pressure
            qv_s                   , & ! specific humidity at the surface
            grad_zqvlow            , & ! specific humidity of lowest atmospheric layer
            lai                    , & ! leaf area index
            cos_zen_ang            , & ! actual cosine of zenit angle
            sfldir_par             , & ! direct component of PAR at the ground
            sfldifd_par            , & ! diffuse downward component of PAR at the ground
            sfldifu_par                ! diffuse upward component of PAR at the ground

        REAL (KIND = vpp), DIMENSION(nvec), INTENT(INOUT) :: &
            fwet                   , & ! fraction of wet canopy
            fdry                   , & ! Fraction of dry canopy
            wtl                    , & ! heat conductance for leaf
            rs_leaf                , & ! leaf boundary layer resistance
            par_sun                , & ! PAR sunlit leaves
            par_sha                , & ! PAR shaded leaves
            ea                     , & ! vapor pressure of air in plant canopy
            lai_sun                , & ! LAI sunlit leaves
            lai_sha                , & ! LAI shaded leaves
            sla_sun                , & ! SLA sunlin leaves
            sla_sha                    ! SLA shaded leaves

        !======================================================================!
        ! Local scalars:
        ! ---------------------------------------------------------------------!
        INTEGER           ::         &
            i                          ! loop index in x-direction


        REAL(KIND = vpp)  ::         &
            ! variables for esai and elai calculations
            tlai                   , & ! one-sided leaf area index, no burying by snow      ( m2/m2)
            tsai                   , & ! one-sided stem area index, no burying by snow      ( m2/m2)
            htop                   , & ! canopy top                                         (  m   )
            hbot                   , & ! canopy bottom                                      (  m   )
            ol                     , & ! thickness of canopy layer covered by snow          (  m   )
            fb                     , & ! fraction of canopy layer covered by snow           (0 to 1)
            ! variables for two-big-leaf calculations
            cosz                   , & ! cosine zenith angle                                ( rad  )
            chil                   , & ! range -0.4 <= pft_CN_par(pft, X_l) <= 0.6
            G_mu                   , & ! leaf projection in solar direction                 (0 to 1)
            K_light                , & ! optical depth direct beam per unit LAI+SAI
            laifra                 , & ! fraction of canopy leaf area                       (0 to 1)
            sun_dir                , & ! total canopy absorbed indirect from direct         ( W/m2 )
            sun_up                 , & ! sun canopy absorbed indirect from indirect         ( W/m2 )
            sun_down               , & ! sun canopy absorbed indirect from direct           ( W/m2 )
            sun_alf                , & ! sun canopy total absorbed by leaves                ( W/m2 )
            sun_aperlai            , & ! sun canopy total absorbed per unit LAI             ( W/m2 )
            sha_up                 , & ! shade canopy absorbed indirect from indirect       ( W/m2 )
            sha_down               , & ! shade canopy absorbed indirect from direct         ( W/m2 )
            sun_atot               , & ! sun canopy total absorbed                          ( W/m2 )
            sha_atot               , & ! shade canopy total absorbed                        ( W/m2 )
            sha_alf                , & ! shade canopy total absored by leaves               ( W/m2 )
            sha_aperlai            , & ! shade canopy total absorbed per unit LAI           ( W/m2 )
            qaf                    , & ! humidity gradient                                  (kg/kg )
            ! temporal variables
            t1                     , & ! temporary variable for optical properties
            t2                     , & ! temporary variable for optical properties
            fau_1                  , & ! temporary variable for optical properties
            fau_2                      ! temporary variable for optical properties


#ifndef ALLOC_WKARR
! Local (automatic) arrays:
! -------------------------
        REAL(KIND = vpp)  ::         &
            elai        (nvec)     , & ! one-sided LAI with burying by snow                 ( m2/m2)
            esai        (nvec)     , & ! one-sided SAI with burying by snow                 ( m2/m2)
            vai         (nvec)     , & ! Total leaf area index + stem area index            ( m2/m2)
            f_sun       (nvec)         ! Fraction of sunlit canopy                          (0 to 1)

#endif

        LOGICAL :: lcontrol = .FALSE.
        !- End of header
        !=======================================================================

        !-----------------------------------------------------------------------
        ! Begin subroutine: get_stomatal_data
        !-----------------------------------------------------------------------

#ifdef __COSMO__
! Double precision equivalents of module variables
        pi = REAL(pi_wp, vpp)
#endif

        ! Calculation of esai and elai
        !-----------------------------------------------------------------------
        DO i = ivstart, ivend
#ifndef MESSY
            IF (soiltyp_subs(i) >= 3 .AND. zep_s(i) < 0.0_vpp) THEN  ! upwards directed potential evaporation
#else
            ! calculate this part always for MESSy
            IF (soiltyp_subs(i) >= 3) THEN ! < 3 is ice and rock
#endif
                IF (pft_CN_par(pft(i), SLA_m) > 0.0_vpp) THEN
                    tlai = (exp(LEAFC * pft_CN_par(pft(i), SLA_m)   + &
                                    log(pft_CN_par(pft(i), SLA_o))) - &
                                        pft_CN_par(pft(i), SLA_o))  / &
                                        pft_CN_par(pft(i), SLA_m)
                ELSE
                    tlai = pft_CN_par(pft(i), SLA_o) * LEAFC
                ENDIF

                tlai = max(0.0_vpp, tlai)


                IF (pft_CN_par(pft(i), WOODY) == 1.0_vpp) THEN
                    ! TREE
                    tsai = 0.25_vpp * tlai
                    htop = ((3.0_vpp * DEADSTEMC * TAPER * TAPER) / &
                            (pi * STOCKING * DWOOD))**(1.0_vpp / 3.0_vpp)

                    htop = min(htop, (FORC_HGT_U / (pft_CN_par(pft(i), R_DISP)  + &
                                                    pft_CN_par(pft(i), R_ZOM))) - 3.0_vpp)
                    htop = max(htop, 0.01_vpp)
                    hbot = max(0.0_vpp, min(3.0_vpp, htop - 1.0_vpp))
                ELSE
                    ! GRASS
                    tsai = 0.05_vpp * tlai
                    htop = max (0.25_vpp, tlai * 0.25_vpp)
                    htop = min(htop, (FORC_HGT_U / (pft_CN_par(pft(i), R_DISP) + &
                                                    pft_CN_par(pft(i), R_ZOM))) - 3.0_vpp)
                    htop = max(htop, 0.01_vpp)
                    hbot = max(0.0_vpp, min(0.05_vpp, htop - 0.20_vpp))
                ENDIF

                ol       = min(max(zdz_snow(i) - hbot, 0.0_vpp), htop - hbot)
                fb       = 1.0_vpp - ol / max(1.e-06_vpp  , htop - hbot)

                elai(i)  = max(tlai * fb, 0.0_vpp)
                esai(i)  = max(tsai * fb, 0.0_vpp)
            ENDIF
        ENDDO

        ! Calculation of: 1 - leaf boundary layer resistance [s/m]
        !                 2 - canopy air vapor pressure (Pa)
        !-----------------------------------------------------------------------
        DO i = ivstart, ivend
#ifndef MESSY
            IF (soiltyp_subs(i) >= 3 .AND. zep_s(i) < 0.0_vpp) THEN  ! upwards directed potential evaporation
#else
            ! calculate this part always for MESSy
            IF (soiltyp_subs(i) >= 3) THEN ! < 3 is ice and rock
#endif
                rs_leaf(i) = 1.0_vpp / TTC_CAN * &
                    ((ustar_fv(i) / pft_CN_par(pft(i), D_LEAF))**(-1.0_vpp / 2.0_vpp))

                qaf   = (qv_s(i) + grad_zqvlow(i)) / 2.0_vpp
                ea(i) = (  ps(i) * qaf) / 0.622_vpp
            ENDIF
        ENDDO


        ! Calculations of: 1. optical parameters based on work Sellers (1986);
        !                  2. sunlit and shaded parameters for v_cmax25;
        !                  3. SLA sunlit and shaded leaves
        !                     (Thornton and Zimmermann, 2007).
        !-----------------------------------------------------------------------
        DO i = ivstart, ivend
#ifndef MESSY
            IF (soiltyp_subs(i) >= 3 .AND. zep_s(i) < 0.0_vpp) THEN  ! upwards directed potential evaporation
#else
            ! calculate this part always for MESSy
            IF (soiltyp_subs(i) >= 3) THEN ! < 3 is ice and rock
#endif
                cosz = max(0.001_vpp, cos_zen_ang(i))
                chil = min(max(pft_CN_par(pft(i), X_l), -0.4_vpp), 0.6_vpp)
                IF (abs(chil) <= 0.01_vpp) THEN
                    chil = 0.01_vpp
                ENDIF

                fau_1 = 0.5_vpp - 0.633_vpp * chil - 0.330_vpp * chil**2.0_vpp
                fau_2 = 0.877_vpp * (1.0_vpp - 2.0_vpp * fau_1)
                G_mu  = fau_1 + fau_2 * cosz


                vai(i) = lai(i) + esai(i) ! stems = esai

                IF (cos_zen_ang(i) > 0.0_vpp .and. lai(i) > 0.0_vpp .and. G_mu > 0.0_vpp) THEN
                    K_light  = G_mu / cosz
                    t1       = min(K_light * lai(i), 40.0_vpp)
                    t2       = exp(-1.0_vpp * t1)
                    f_sun(i) = (1.0_vpp - t2) / t1
                    IF (lai(i) > 0.01_vpp) THEN
                        lai_sun(i) = lai(i) * f_sun(i)
                        lai_sha(i) = lai(i) * (1.0_vpp - f_sun(i))

                        ! calculate the average specific leaf area for sunlit and shaded
                        ! canopies, when effective LAI > 0
                        sla_sun(i) = (t2 * pft_CN_par(pft(i), SLA_m) * K_light * lai(i) + &
                                      t2 * pft_CN_par(pft(i), SLA_m) +                    &
                                      t2 * pft_CN_par(pft(i), SLA_o) * K_light -          &
                                           pft_CN_par(pft(i), SLA_m) -                    &
                                           pft_CN_par(pft(i), SLA_o) * K_light) /         &
                                     (K_light * (t2 - 1.0_vpp))

                        sla_sha(i) = ((pft_CN_par(pft(i), SLA_o) +                        &
                                     ((pft_CN_par(pft(i), SLA_m) * lai(i)) / 2.0_vpp)) *  &
                                         lai(i) - sla_sun(i) * lai_sun(i)) / lai_sha(i)
                    ELSE
                        f_sun(i)   = 1.0_vpp
                        lai_sun(i) = lai(i)
                        lai_sha(i) = 0.0_vpp
                        sla_sun(i) = pft_CN_par(pft(i), SLA_o)
                        sla_sha(i) = 0.0_vpp
                    ENDIF
                ELSE
                    f_sun(i)   = 0.0_vpp
                    lai_sun(i) = 0.0_vpp
                    lai_sha(i) = lai(i)
                    sla_sun(i) = 0.0_vpp
                    sla_sha(i) = pft_CN_par(pft(i), SLA_o) + &
                                (pft_CN_par(pft(i), SLA_m) * lai(i)) / 2.0_vpp
                ENDIF
            ENDIF
        ENDDO

        ! Calculation of PAR for sunlit and shaded leaves
        !-----------------------------------------------------------------------
        DO i = ivstart, ivend
#ifndef MESSY
            IF (soiltyp_subs(i) >= 3 .AND. zep_s(i) < 0.0_vpp) THEN  ! upwards directed potential evaporation
#else
            ! calculate this part always for MESSy
            IF (soiltyp_subs(i) >= 3) THEN ! < 3 is ice and rock
#endif
                IF (cos_zen_ang(i) > 0.0_vpp .and. lai(i) > 0.0_vpp) THEN

                    ! 1. calculation direct component of PAR
                    sun_dir = max(sfldir_par(i), 0.0_vpp)

                    ! 2. calculate the total indirect flux absorbed by the sunlit
                    !    and shaded canopy based on these fractions

                    sun_down = sfldifd_par(i) * f_sun(i)
                    sun_up   = sfldifu_par(i) * f_sun(i)
                    sha_down = sfldifd_par(i) * (1.0_vpp - f_sun(i))
                    sha_up   = sfldifu_par(i) * (1.0_vpp - f_sun(i))

                    ! 3. calculate the total flux absorbed in the sunlit and shaded
                    !    canopy as the sum of these terms

                    sun_atot = sun_dir + sun_down + sun_up
                    sha_atot = sha_down + sha_up

                    ! 4. calculate the total flux absorbed by leaves in the sunlit
                    !    and shaded canopies
                    laifra  = lai(i) / vai(i)
                    sun_alf = sun_atot * laifra
                    sha_alf = sha_atot * laifra

                    ! 5. calculate the fluxes per unit lai in the sunlit and shaded
                    !    canopies
                    IF (lai_sun(i) > 0.0_vpp) THEN
                        sun_aperlai = sun_alf / lai_sun(i)
                    ELSE
                        sun_aperlai = 0.0_vpp
                    ENDIF

                    IF (lai_sha(i) > 0.0_vpp) THEN
                        sha_aperlai = sha_alf / lai_sha(i)
                    ELSE
                        sha_aperlai = 0.0_vpp
                    ENDIF
                ELSE
                    sun_dir     = 0.0_vpp
                    sun_down    = 0.0_vpp
                    sun_up      = 0.0_vpp
                    sha_down    = 0.0_vpp
                    sha_up      = 0.0_vpp
                    sun_atot    = 0.0_vpp
                    sha_atot    = 0.0_vpp
                    sun_alf     = 0.0_vpp
                    sha_alf     = 0.0_vpp
                    sun_aperlai = 0.0_vpp
                    sha_aperlai = 0.0_vpp
                ENDIF

                ! PAR fluxes for sunlit and shaded leaves
                par_sun(i) = sun_aperlai
                par_sha(i) = sha_aperlai
            ENDIF
        ENDDO

        ! Calculations of additional parameters: 1. fraction of wet and dry leaves
        !                                        2. heat conductance for leaf
        !-----------------------------------------------------------------------
        DO i = ivstart, ivend
#ifndef MESSY
            IF (soiltyp_subs(i) >= 3 .AND. zep_s(i) < 0.0_vpp) THEN  ! upwards directed potential evaporation
#else
            ! calculate this part always for MESSy
            IF (soiltyp_subs(i) >= 3) THEN ! < 3 is ice and rock
#endif
                IF (H2OCAN > 0.2_vpp) THEN
                    ! look fracwet.
                    fwet(i) = 1.0_vpp
                    fwet(i) = min(fwet(i), 1.0_vpp)
                ELSE
                    fwet(i) = 0.0_vpp
                ENDIF

                fdry(i) = (1.0_vpp - fwet(i)) * lai(i) / vai(i)
                wtl(i)  = vai(i) / rs_leaf(i)
            ENDIF
        ENDDO

        ! Additional options for run control
        ! Write log file
        !-----------------------------------------------------------------------
        IF (lcontrol) THEN
            DO i = ivstart, ivend
                write(*,*) 'leaf area index         - lai     : ', lai(i)
                write(*,*) 'total LAI + SAI         - vai     : ', vai(i)
                write(*,*) 'LAI sunlit              - lai_sun : ', lai_sun(i)
                write(*,*) 'LAI shaded              - lai_sha : ', lai_sha(i)
                write(*,*) 'SLA sunlit              - sla_sun : ', sla_sun(i)
                write(*,*) 'SLA shaded              - sla_sha : ', sla_sha(i)
                write(*,*) 'PAR sunlit              - par_sun : ', par_sun(i)
                write(*,*) 'PAR shaded              - par_sha : ', par_sha(i)
                write(*,*) 'Vapor pressure          - ea      : ', ea(i)
                write(*,*) 'Leaf resistance         - rs_leaf : ', ea(i)
                write(*,*) 'Fraction of wet canopy  - fwet    : ', fwet(i)
                write(*,*) 'Fraction of dry canopy  - fdry    : ', fdry(i)
                write(*,*) 'Heat conductance (leaf) - wtl     : ', wtl(i)
            ENDDO
        ENDIF

        RETURN
    !===========================================================================
    END SUBROUTINE get_stomatal_data



    !  Calculating stomatal resistance and leaf photosynthesis parameters
    ! based on the Ball-Berry approach (stomatal resistance) and Farquhar
    ! and Collatz models (leaf photosynthesis)
    !--------------------------------------------------------------------------!
    SUBROUTINE stomata(                   &
                       nvec             , & ! array dimensions
                       ivstart          , & ! start index for computations in the parallel program
                       ivend            , & ! end index for computations in the parallel program
                       soiltyp_subs     , & ! type of the soil                                        (keys 0-9)
                       pft              , & ! type of PFT                                             (1 - C3, 2 - C4)
                       zep_s            , & ! potential evaporation for t_s
                       ps               , & ! surface pressure                                        ( Pa  )
                       zth_low          , & ! potential temperature of lowest layer                   ( K   )
                       ztsk             , & ! skin temperature                                        ( K   )
                       zepsat           , & ! saturation vapour pressure at near surface temperature
                       ea               , & ! vapor pressure of air in the plant canopy               ( pa  )
                       sla              , & ! specific leaf area index      (sunlit or shaded)        (m2/m2)
                       apar             , & ! PAR                           (sunlit or shaded)        (W/m2 )
                       zf_wat           , & ! soil water function for stomatal resistance
                       f_dyl            , & ! daylength limitation function                           (0 to 1)
                       rs_leaf          , & ! canopy leaf resistance                                  (s/m)
                       rs               , & ! stomatal resistance           (sunlit or shaded)        (s/m)
                       leaf_psn         , & ! leaf photosynthesis           (sunlit or shaded)        (umol CO2 m^-2 s^-1)
                       v_cmax             & ! maximum rate of carboxylation (sunlit or shaded)        (umol CO2 m^-2 s^-1)
                                          )
        ! Modules used:
        !-----------------------------------------------------------------------

        IMPLICIT NONE
        !----------------------------------------------------------------------!
        ! Declarations
        !----------------------------------------------------------------------!
        INTEGER, INTENT(IN)  ::      &
            nvec                   , & ! array dimensions
            ivstart                , & ! start index for computations in the parallel program
            ivend                      ! end index for computations in the parallel program

        INTEGER, DIMENSION(nvec), INTENT(IN)    :: &
            soiltyp_subs           , & ! type of the soil (keys 0-9)
            pft                        ! type of PFT

        REAL (KIND = vpp), DIMENSION(nvec), INTENT(IN)    :: &
            zep_s                  , & ! potential evaporation for t_s
            ps                     , & ! surface pressure
            zth_low                , & ! potential temperature of lowest layer
            ztsk                   , & ! skin temperature
            zepsat                 , & ! saturation vapour pressure at near surface temperature
            ea                     , & ! vapor pressure of air in the plant canopy
            sla                    , & ! specific leaf area index                      (sunlit or shaded)
            apar                   , & ! PAR                                           (sunlit or shaded)
            zf_wat                 , & ! soil water function for stomatal resistance
            f_dyl                  , & ! daylength limitation function
            rs_leaf                    ! canopy leaf resistance

        REAL (KIND = vpp), DIMENSION(nvec), INTENT(INOUT) :: &
            rs                     , & ! stomatal resistance                           (sunlit or shaded)
            leaf_psn               , & ! leaf photosynthesis                           (sunlit or shaded)
            v_cmax                     ! maximum rate of carboxylation                 (sunlit or shaded)


        !======================================================================!
        ! Local scalars:
        ! ---------------------------------------------------------------------!
        INTEGER           ::         &
            i                      , & ! loop index in x-direction
            iter                       ! iteration index

        REAL(KIND = vpp)  ::         &
            zf_water               , & ! soil water content  (limitation function)          (0 to 1)
            f_t_veg                , & ! mimics thermal function                            (0 to 1)
            slamax                 , & ! spatial leaf area                                  (m2/m2)
            ppf                    , & ! absorb photosynthetic photon flux                  (umol photons/m2/s)
            lnc                    , & ! leaf N concentration per unit projected LAI        (gN leaf/m2)
            act                    , & ! rubisco activity                                   (umol/mgRubisco/min)
            k_c                    , & ! michaelis-Menten constant for CO2                  (Pa)
            k_o                    , & ! michaelis-Menten constant for  O2                  (Pa)
            wc                     , & ! rubisco limited photosynthesis                     (umol co2/m2/s)
            wj                     , & ! light limited photosynthesis                       (umol co2/m2/s)
            we                     , & ! export limited photosynthesis                      (umol co2/m2/s)
            ci                     , & ! intracellular leaf CO2                             (Pa)
            cea                    , & ! constrain ea or else model blows up                (Pa)
            cp                     , & ! CO2 compensation point                             (Pa)
            cs                     , & ! CO2 concentration at leaf surface                  (Pa)
            ! Additional parameters for calculations
            recal_rs_coef          , & ! stomatal resistance recalculation
                                       ! coefficient from (s m2/umol) to (s/m)
            tc                     , & ! skin temperature                                   (degree C)
            can_leaf               , & ! canopy leaf resistance                             (s m**2 / umol)
            awc                    , & ! intermediate calcuation for wc
            a_tmp                  , & ! intermediate calculations for rs                   ( s/m )
            b_tmp                  , & ! intermediate calculations for rs                   ( s/m )
            c_tmp                  , & ! intermediate calculations for rs                   ( s/m )
            q                      , & ! intermediate calculations for rs                   ( s/m )
            r_s1                   , & ! roots for rs                                       (s/m)
            r_s2                       ! roots for rs                                       (s/m)

        LOGICAL :: lcontrol = .FALSE.
        !- End of header
        !=======================================================================

        !-----------------------------------------------------------------------
        ! Begin subroutine: stomata
        !-----------------------------------------------------------------------

#ifdef __COSMO__
! Double precision equivalents of module variables
        pi      = REAL(pi_wp     , vpp)
        t0_melt = REAL(t0_melt_wp, vpp)
#endif

        DO i = ivstart, ivend
#ifndef MESSY
            IF (soiltyp_subs(i) >= 3 .AND. zep_s(i) < 0.0_vpp) THEN  ! upwards directed potential evaporation
#else
            ! calculate this part always for MESSy
            IF (soiltyp_subs(i) >= 3) THEN ! < 3 is ice and rock
#endif
                ! Convert variables to appropriate units
                !-------------------------------------------------------------------
                recal_rs_coef = ps(i) / (RGAS * 0.001_vpp * zth_low(i)) * 1.e06_vpp                ! (s/m) to (s m**2 / umol)
                tc            = ztsk(i) - t0_melt                                                  ! (K) to (C)
                can_leaf      = rs_leaf(i) / recal_rs_coef                                         ! (s/m) to (s m**2 / umol)

                ! Calculation of a function that mimics thermal breakdown
                !           of metabolic process
                !-------------------------------------------------------------------
                f_t_veg = 1.0_vpp + exp((-2.2e5_vpp + 710.0_vpp * (tc + t0_melt)) / &
                                        (RGAS       * 0.001_vpp * (tc + t0_melt))   )

                ! Calculation of area based leaf nitrogen concentration for sunlit
                ! and shaded leaves. (Thornton and Zimmermann, 2007)
                !--------------------------------------------------------------------
                slamax    = max(1.e-5_vpp, sla(i))
                lnc       = 1.0_vpp / (pft_CN_par(pft(i), CN_L) * slamax)
                act       = ALPHA_R_25 * ALPHA_VMAX**((tc - 25.0_vpp) / 10.0_vpp)
                v_cmax(i) = lnc * pft_CN_par(pft(i), F_LNR) * F_NR * act / f_t_veg * &
                      zf_wat(i) * pft_CN_par(pft(i), F_N  ) * f_dyl(i)

                IF (zf_wat(i) == 0.0_vpp) THEN
                    zf_water = 0.2_vpp
                ELSE
                    zf_water = zf_wat(i)
                ENDIF

                IF (apar(i) <= 0.0_vpp) THEN ! night time
                    rs(i)  = min(RSMAX0, 1.0_vpp / GLMIN * recal_rs_coef)
                    leaf_psn(i) = 0.0_vpp

                ELSE
                    ! Calculation of the absorb photosynthetic photon flux
                    ppf = 4.6_vpp * apar(i)
                    ! Calculation of the Michaelis-Menten constants
                    !                    and CO2 compensation point
                    k_c = K_C25 * ALPHA_KC**(( tc - 25.0_vpp) / 10.0_vpp)
                    k_o = K_O25 * ALPHA_KO**(( tc - 25.0_vpp) / 10.0_vpp)
                    awc = k_c * (1.0_vpp + CON_O2 / k_o)
                    cp  = 0.5_vpp * k_c / k_o * 0.21_vpp * CON_O2

                    ! Calculation of the intracellular leaf CO2 and the constrain
                    ! ea or else model blows up
                    IF (REAL(pft(i), vpp) == 1.0_vpp) THEN
                        ci  = 0.7_vpp * CON_CO2
                        cea = max(min(ea(i), zepsat(i)), 0.25_vpp * zepsat(i))
                    ELSE
                        ci  = 0.4_vpp * CON_CO2
                        cea = max(min(ea(i), zepsat(i)), 0.4_vpp * zepsat(i))
                    ENDIF

                    ! Calculation of stomatal resistance and leaf photosynthesis
                    DO iter = 1, niter
                        IF (REAL(pft(i), vpp) == 1.0_vpp) THEN
                            ! Calculation of photosynthesis in C3 plants is based on
                            ! the model of Farquhar et al., (1980).
                            !-------------------------------------------------------
                            wc = max(ci - cp, 0.0_vpp) * v_cmax(i) / (ci + awc)
                            wj = max(ci - cp, 0.0_vpp) * (pft_CN_par(pft(i), ALPHA) * ppf) / (ci + 2.0_vpp * cp)
                            we = 0.5_vpp * v_cmax(i)
                        ELSE
                            ! Calculation of photosynthesis in C4 plants is based on
                            ! the model of Collatz et al., (1992).
                            !-------------------------------------------------------
                            wc = v_cmax(i)
                            wj = pft_CN_par(pft(i), ALPHA) * ppf
                            we = 4000.0_vpp * v_cmax(i) * ci / ps(i)
                        ENDIF

                        leaf_psn(i) = min(wc, wj, we)
                        cs     = max(CON_CO2 - 1.37_vpp * can_leaf * ps(i) * leaf_psn(i), MPE)

                        !-----------------------------------------------------------
                        ! Calculation of stomatal resistance is based on
                        ! the Ball-Berry approach (Ball and Berry, 1988)
                        !-----------------------------------------------------------
                        a_tmp =  pft_CN_par(pft(i), MP) * leaf_psn(i) * ps(i) * cea  / (cs * zepsat(i)) + GLMIN
                        b_tmp = (pft_CN_par(pft(i), MP) * leaf_psn(i) * ps(i) / cs + GLMIN) * can_leaf - 1.0_vpp
                        c_tmp = - can_leaf

                        IF (b_tmp >= 0.0_vpp) THEN
                            q = -0.5_vpp * (b_tmp + sqrt(b_tmp * b_tmp - 4.0_vpp * a_tmp * c_tmp))
                        ELSE
                            q = -0.5_vpp * (b_tmp - sqrt(b_tmp * b_tmp - 4.0_vpp * a_tmp * c_tmp))
                        ENDIF
                        r_s1 = q / a_tmp
                        r_s2 = c_tmp / q
                        rs(i) = max(r_s1, r_s2)
                        ci = max(cs - leaf_psn(i) * ps(i) * 1.65_vpp * rs(i), 0.0_vpp)
                    ENDDO

                    rs(i) = min(RSMAX0, rs(i) * recal_rs_coef)
                ENDIF
            ENDIF

        ENDDO

        ! Additional options for run control
        ! Write log file
        !-----------------------------------------------------------------------
        IF (lcontrol) THEN
            DO i = ivstart, ivend
                write(*,*) 'wc                                       : ', wc
                write(*,*) 'wj                                       : ', wj
                write(*,*) 'we                                       : ', we
                write(*,*) 'maximum rate of carboxylation - v_cmax   : ', v_cmax(i)
                write(*,*) 'stomatal resistance           - rs       : ', rs(i)
                write(*,*) 'leaf photosynthesis           - leaf_psn : ', leaf_psn(i)
            ENDDO
        ENDIF


        RETURN
    !===========================================================================
    END SUBROUTINE stomata




     ! Calculating GPP and NPP fluxes
     !-------------------------------------------------------------------------!
    SUBROUTINE respiration(                   &
                           nvec             , & ! array dimensions
                           ivstart          , & ! start index for computations in the parallel program
                           ivend            , & ! end index for computations in the parallel program
                           soiltyp_subs     , & ! type of the soil (keys 0-9)
                           zep_s            , & ! potential evaporation for t_s
                           ps               , & ! surface pressure                                        ( Pa  )
                           zth_low          , & ! potential temperature of lowest layer                   ( K   )
                           rs_leaf          , & ! canopy leaf resistance                                  ( s/m )
                           par_sun          , & ! photosynthetic active radiation for sunlit leaf         ( W/m2)
                           par_sha          , & ! photosynthetic active radiation for shaded leaf         ( W/m2)
                           vcmax_sun        , & ! maximum rate of carboxylation for sunlit leaves         (umol CO2 m^-2 s^-1)
                           vcmax_sha        , & ! maximum rate of carboxylation for shaded leaves         (umol CO2 m^-2 s^-1)
                           lai_sun          , & ! LAI                 (sunlit leaves)                     (m2/m2)
                           lai_sha          , & ! LAI                 (shaded leaves)                     (m2/m2)
                           lai_can          , & ! LAI                 (    canopy   )                     (m2/m2)
                           psn_sun          , & ! leaf photosynthesis (sunlit leaves)                     (umol CO2 m^-2 s^-1)
                           psn_sha          , & ! leaf photosynthesis (shaded leaves)                     (umol CO2 m^-2 s^-1)
                           psn_can          , & ! leaf photosynthesis (    canopy   )                     (umol CO2 m^-2 s^-1)
                           fplres           , & ! total plant respiration                                 (umol CO2 m^-2 s^-1)
                           cleaf            , & ! CO2 partial pressure at leaf surface                    ( Pa  )
                           gpp              , & ! gross primary production                                (gC/m2/s)
                           fnpp               & ! net primary production                                  (gC/m2/s)
                                              )
        ! Modules used:
        !-----------------------------------------------------------------------

        IMPLICIT NONE
        !----------------------------------------------------------------------!
        ! Declarations
        !----------------------------------------------------------------------!
        INTEGER, INTENT(IN)  ::      &
            nvec                   , & ! array dimensions
            ivstart                , & ! start index for computations in the parallel program
            ivend                      ! end index for computations in the parallel program

        INTEGER, DIMENSION(nvec), INTENT(IN)    :: &
            soiltyp_subs               ! type of the soil (keys 0-9)

        REAL (KIND = vpp), DIMENSION(nvec), INTENT(IN)    :: &
            zep_s                  , & ! potential evaporation for t_s
            ps                     , & ! surface pressure
            zth_low                , & ! potential temperature of lowest layer
            rs_leaf                , & ! canopy leaf resistance
            par_sun                , & ! photosynthetic active radiation for sunlit leaf
            par_sha                , & ! photosynthetic active radiation for shaded leaf
            vcmax_sun              , & ! maximum rate of carboxylation for sunlit leaves
            vcmax_sha              , & ! maximum rate of carboxylation for shaded leaves
            lai_sun                , & ! LAI                 (sunlit leaves)
            lai_sha                , & ! LAI                 (shaded leaves)
            lai_can                , & ! LAI                 (    canopy   )
            psn_sun                , & ! Leaf photosynthesis (sunlit leaves)
            psn_sha                , & ! Leaf photosynthesis (shaded leaves)
            psn_can                    ! Leaf photosynthesis (    canopy   )

        REAL (KIND = vpp), DIMENSION(nvec), INTENT(INOUT) :: &
            fplres                 , & ! total plant respiration
            cleaf                  , & ! CO2 partial pressure at leaf surface
            gpp                    , & ! gross primary production
            fnpp                       ! net primary production


        !======================================================================!
        ! Local scalars:
        ! ---------------------------------------------------------------------!
        INTEGER           ::         &
            i                          ! loop index in x-direction

        REAL(KIND = vpp)  ::         &
            rb                     , & ! leaf boundary resistance                           (s m2/umol)
            coef                   , & ! Stomatal resistance recalculation
                                       ! coefficient from                                   (s m2/umol) to (s/m)
            psnsun_to_cpool        , & ! gross primary production sunlit leaf               (gC/m2/s)
            psnsha_to_cpool            ! gross primary production shaded leaf               (gC/m2/s)
        !- End of header
        !=======================================================================

        !-----------------------------------------------------------------------
        ! Begin subroutine: respiration
        !-----------------------------------------------------------------------

        ! Section 1. Calculation total plant respiration and CO2 partial pressure
        !-----------------------------------------------------------------------
        DO i = ivstart, ivend
#ifndef MESSY
            IF (soiltyp_subs(i) >= 3 .AND. zep_s(i) < 0.0_vpp) THEN  ! upwards directed potential evaporation
#else
            ! calculate this part always for MESSy
            IF (soiltyp_subs(i) >= 3) THEN ! < 3 is ice and rock
#endif
                IF ((par_sha(i) > 0.0_vpp) .or. (par_sun(i) > 0.0_vpp )) THEN
                    fplres(i) = 0.015_vpp * (vcmax_sun(i) * lai_sun(i) + &
                                             vcmax_sha(i) * lai_sha(i)   )
                ELSE
                    fplres(i) = 0.015_vpp * vcmax_sha(i) * lai_can(i)
                ENDIF

                coef     = ps(i) / (RGAS * 0.001_vpp * zth_low(i)) * 1.e06_vpp
                rb       = rs_leaf(i) / coef
                cleaf(i) = CON_CO2 - 1.37_vpp * rb * ps(i) * psn_can(i)

                ! Convert psn from umol/m2/s -> gC/m2/s
                psnsun_to_cpool = psn_sun(i) * lai_sun(i) * 12.011e-6_vpp
                psnsha_to_cpool = psn_sha(i) * lai_sha(i) * 12.011e-6_vpp

                gpp(i)  = psnsun_to_cpool + psnsha_to_cpool
                fnpp(i) = GPPFACT * psn_can(i) * 12.011e-6_vpp                                        ! umolC/m2/sec to gC/m2/sec

                !npp = gpp - (fsres_a * 12.011e-6)
            ENDIF
        ENDDO

        RETURN
    !===========================================================================
    END SUBROUTINE respiration




    SUBROUTINE biomass_evolution(                   &
                                 nvec             , & ! array dimensions
                                 ivstart          , & ! start index for computations in the parallel program
                                 ivend            , & ! end index for computations in the parallel program
                                 dt               , & ! time step
                                 pft              , & ! type of PFT                                             (1 - C3, 2 - C4)
                                 soiltyp_subs     , & ! type of the soil (keys 0-9)
                                 zep_s            , & ! potential evaporation for t_s
                                 plcov            , & ! fraction of plant cover                                  --
                                 psn_can          , & ! leaf photosynthesis (canopy)                            (umol CO2 m^-2 s^-1)
                                 sur_lai          , & ! leaf area index depending on biomass evolution          (m2/m2)
                                 sur_panday       , & ! daily values of photosyntesis                           ()
                                 sur_biomass      , & ! leaf biomass                                            ()
                                 sur_panfm          & ! maximum leaf assimilation                               ()
                                                    )
        ! Calculating changes in leaf area index (LAI) depending on biomass
        ! evolution
        !-----------------------------------------------------------------------

        ! Modules used:
        !-----------------------------------------------------------------------

        IMPLICIT NONE
        !----------------------------------------------------------------------!
        ! Declarations
        !----------------------------------------------------------------------!
        INTEGER, INTENT(IN)  ::      &
            nvec                   , & ! array dimensions
            ivstart                , & ! start index for computations in the parallel program
            ivend                      ! end index for computations in the parallel program

        REAL (KIND = vpp), INTENT(IN)  ::  &
            dt                         ! time step

        INTEGER, DIMENSION(nvec), INTENT(IN)    :: &
            pft                    , & ! type of PFT
            soiltyp_subs               ! type of the soil (keys 0-9)

        REAL (KIND = vpp), DIMENSION(nvec), INTENT(IN)    :: &
            zep_s                  , & ! potential evaporation for t_s
            plcov                  , & ! fraction of plant cover                         --
            psn_can                    ! leaf photosynthesis (canopy)

        REAL (KIND = vpp), DIMENSION(nvec), INTENT(INOUT) :: &
            sur_lai                , & ! leaf area index depending on biomass evolution          (m2/m2)
            sur_panday             , & ! daily values of photosyntesis                           (kg CO2 m-2 s-1)
            sur_biomass            , & ! leaf biomass                                            ()
            sur_panfm                  ! maximum leaf assimilation                               (kg CO2 m-2 s-1)

        !======================================================================!
        ! Local scalars:
        ! ---------------------------------------------------------------------!
        INTEGER           ::         &
            i                      , & ! loop index in x-direction
            hour, minut, sec       , & ! actual hour, min, sec
            nactday                    ! day of the year

        REAL(KIND = vpp)  ::         &
            acthour                , & ! actual hour of the day
            zxsefold               , & ! leaf life expectancy
            zxm                        ! senesence of active biomass

        CHARACTER (LEN=14) yactdate1   ! actual date in the form yyyymmddhhmmss
        CHARACTER (LEN=28) yactdate2   ! actual date in the form &
                                       ! wd   dd.mm.yy  hh mm ss UTC
#ifndef ALLOC_WKARR
! Local (automatic) arrays:
! -------------------------
        REAL(KIND = vpp)  ::         &
            spanmax      (nvec)     , & ! maximum photosynthesis rate in optimal conditions
            spsn         (nvec)         ! values of photosynthesis                                   (kg CO2 m-2 s-1)
#endif
        !- End of header
        !=======================================================================

        !-----------------------------------------------------------------------
        ! Begin subroutine: biomass_evolution
        !-----------------------------------------------------------------------

        ! Read current date
        !-----------------------------------------------------------------------
        CALL get_utc_date(ntstep, ydate_ini, dt, itype_calendar, &
                          yactdate1, yactdate2, nactday, acthour)

        ! Get actual information about month and day
        !-----------------------------------------------------------------------
        READ (yactdate1(9:10) ,'(I2)') hour
        READ (yactdate1(11:12),'(I2)') minut
        READ (yactdate1(13:14),'(I2)') sec

        ! Define: midnight True or False
        !midnight = (hour + minut + sec == 0)
        DO i = ivstart, ivend
#ifndef MESSY
            IF (soiltyp_subs(i) >= 3 .AND. zep_s(i) < 0.0_vpp) THEN  ! upwards directed potential evaporation
#else
            ! calculate this part always for MESSy
            IF (soiltyp_subs(i) >= 3) THEN ! < 3 is ice and rock
#endif
                ! Initialization:
                spsn(i)    = (psn_can(i)   * M_CO2) / 1000.0_vpp                             ! Convert [umolCO2 m-2 s-1] to kg CO2 m-2 s-1
                spanmax(i) = (2.854717_vpp * M_CO2) / 1000.0_vpp                             ! Set maximum values of An
            ENDIF
        ENDDO

        ! Calculation of the biomass evolution
        !-----------------------------------------------------------------------
        DO i = ivstart, ivend
#ifndef MESSY
            IF (soiltyp_subs(i) >= 3 .AND. zep_s(i) < 0.0_vpp) THEN  ! upwards directed potential evaporation
#else
            ! calculate this part always for MESSy
            IF (soiltyp_subs(i) >= 3) THEN ! < 3 is ice and rock
#endif
                IF ((REAL(hour, vpp) + REAL(minut, vpp) + REAL(sec, vpp)) == 0.0_vpp) THEN  ! Biomass evolution
                    IF (plcov(i) > 0.0_vpp) THEN

                        sur_panday(i) = (sur_panday(i) / 3312.0_vpp) * 28.0_vpp             ! temporary

                        ! Calculates the time change in LAI due to biomass
                        ! evolution Calvet at al (1997, 1998)
                        !-----------------------------------------------------------
                        ! Leaf life expectancy
                        zxsefold = (pft_SURFEX(pft(i), PSEFOLD) * XDAY) * &
                                    min(1.0_vpp, sur_panfm(i) / spanmax(i)) / XDAY

                        ! Avoid possible but unlikely division by zero
                        zxsefold = max(1.0e-8_vpp, zxsefold)
                        ! Limitation of leaf life expectancy
                        zxsefold = max(5.0_vpp   , zxsefold)
                        ! Senesence of active biomass
                        zxm = sur_biomass(i) * (1.0_vpp - exp(-1.0_vpp / zxsefold))
                        ! Decrease biomass:
                        sur_biomass(i) = sur_biomass(i) - zxm
                        ! Growth biomass:
                        sur_biomass(i) = sur_biomass(i) + sur_panday(i) * zbmcoef

                        ! Biomass evolution
                        sur_biomass(i) = max(pft_SURFEX(pft(i), LAIMIN) * &
                                             pft_SURFEX(pft(i), PBSLAI), sur_biomass(i))

                        ! Re-calculation of LAI due to biomass evolution
                        sur_lai(i) = sur_biomass(i) / pft_SURFEX(pft(i), PBSLAI)

                        ! Reset: max and daily net assimilation values for next day (kgCO2 kgAir-1 m s-1):
                        sur_panfm(i)  = 0.0_vpp
                        sur_panday(i) = 0.0_vpp

                    ENDIF
                ELSE
                    ! Define the max values of An for one day
                    sur_panfm(i)  = max(spsn(i), sur_panfm(i))
                    ! Accumulation of panday values
                    sur_panday(i) = sur_panday(i) + spsn(i)
                ENDIF
            ENDIF
        ENDDO

        RETURN
    !===========================================================================
    END SUBROUTINE biomass_evolution



END MODULE sfc_phenology














