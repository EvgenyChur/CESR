1. Added to data_fields:
#===============================================================================
    ! V6_0 CLM     2021-02-28 Churiulin Evgenii
    !  Added 30 fields for the new vegetation scheme

    ! Lines  1088 - 1041
    ! 14. fields for the vegetation parameterization                            ! Chur 2022 - new parameters for vegetation
    ! -----------------------------------------
    #ifndef MESSY
      REAL  (KIND=wp), TARGET, ALLOCATABLE :: &
    #else
      REAL  (KIND=wp), POINTER             :: &
    #endif
        ! All parameters can be activated if lphenology is TRUE. Parameters has tiles.
        ! Radiation
        ! Additional fields for the radiation correction scheme
        ! these are actual values
        sfldir_par   (:,:),       & ! direct component of PAR at the ground            ( W/m2 )
        sfldifd_par  (:,:),       & ! diffuse downward component PAR at the ground     ( W/m2 )
        sfldifu_par  (:,:),       & ! diffuse upward component of PAR at the ground    ( W/m2 )
        sfltrdir_par (:,:),       & ! direct comp. of PAR at surface                   ( W/m2 )
        sfltrdifd_par(:,:),       & ! diffuse downward comp. of PAR at surface         ( W/m2 )
        sfltrdifu_par(:,:),       & ! diffuse upward comp. of PAR at surface           ( W/m2 )
        cos_zen_ang  (:,:),       & ! actual cosine of zenit angle                     ( rad  )
        ! these are accumulated values for aditional fields
        asfldir_par  (:,:),       & ! direct component of PAR at the ground            ( W/m2)
        asfldifd_par (:,:),       & ! diffuse downward component of PAR at the ground  ( W/m2)
        asfldifu_par (:,:),       & ! diffuse upward component of PAR at the ground    ( W/m2)

        ! Fields that are computed in the dynamics and / or physics - TERRA-ML
        ! Field for new vegetation
        sur_lai      (:,:,:),     & ! leaf area index depending on biomass evolution   ( m2/m2)
        sur_panday   (:,:,:),     & ! parameter for daily values of photosyntesis      ( kg CO2 m-2 s-1)
        sur_biomass  (:,:,:),     & ! parameter for surfex biomass                     ( kg of dry matter m^2)
        sur_panfm    (:,:,:),     & ! parameter for surfex maximum leaf assimilation   ( kg CO2 m-2 s-1)
        ztraleav     (:,:,:),     & ! transpiration rate of dry leaves surface         ( mm/hour)
        zverbo       (:,:,:),     & ! total evapotranspiration                         ( mm/hour)
        lai_sun      (:,:,:),     & ! sunlit leaf area index                           ( m2/m2  )
        lai_sha      (:,:,:),     & ! shaded leaf area index                           ( m2/m2  )
        sla_sun      (:,:,:),     & ! specific leaf area index for sunlit leaves       ( m2/m2  )
        sla_sha      (:,:,:),     & ! specific leaf area index for shaded leaves       ( m2/m2  )
        vcmax_sun    (:,:,:),     & ! maximum rate of carboxylation for sunlit leaves  (umol CO2 m^-2 s^-1)
        vcmax_sha    (:,:,:),     & ! maximum rate of carboxylation for shaded leaves  (umol CO2 m^-2 s^-1)
        par_sun      (:,:,:),     & ! photosynthetic active radiation for sunlit leaf  ( W/m2)
        par_sha      (:,:,:),     & ! photosynthetic active radiation for shaded leaf  ( W/m2)
        rs_leaf      (:,:,:),     & ! canopy leaf resistance                           ( s/m )
        psn          (:,:,:),     & ! leaf photosynthesis                              (umol CO2 m^-2 s^-1)
        gpp_flux     (:,:,:),     & ! gross primary production                         (gC/m2/s )
        npp_flux     (:,:,:),     & ! net primary production                           (gC/m2/s )
        !  Field that are computed for vegetation (average)
        aztraleav    (:,:),       & ! average transpiration rate of dry leaves         ( mm/day)
        azverbo      (:,:)          ! average total evapotranspiration                 ( mm/day)
    !===========================================================================

! Total: 30 new parameters
!         5 - accumulated and average fields (3 accumulated and 2 average)
!        25 - usual
!==============================================================================!





2. Added to data_runcontrol.f90
#===============================================================================
! V6_0 CLM     2021-05-15 Evgenii Churiulin
!  Added logical switch lphenology to compute new vegetation scheme

! Lines 552
LOGICAL                          ::           &
    ...,                  &
    lphenology,           & ! calculations of new vegetation scheme                                                                 ! Chur 2022
!==============================================================================!





3. Added data to organize_physics
#===============================================================================
    ! V6_0 CLM     2022-04-05 Evgenii Churiulin
    !  Added additional variables for radiation: sfldir_par, sfldifd_par, sfldifu_par
    !  Added logical switch lphenology to compute new vegetation scheme

    ! Lines 395
    USE data_runcontrol,    ONLY:
        lphenology

    ! Lines 444 - 445
    ! need additionally for updating the host data
    USE data_fields,        ONLY:                                               &
            sfldir_par, sfldifd_par,   sfldifu_par                              &

    ! Lines 3490
    lphenology_d,     & ! forecast with new vegetation scheme                                                                       ! Chur 2022

    ! Lines 3579
    NAMELIST /phyctl/ ..., lphenology, ...

    ! Lines 3645
    lphenology_d        = .FALSE.                                                                                                   ! Chur 2022

    ! Lines 3804
    lphenology        = lphenology_d                                                                                                ! Chur 2022

    ! Lines 4518 - 4521
    IF (lphenology) THEN                                                                                                            ! Chur 2022
        PRINT  *, ' ERROR  *** GPU: lphenology=.TRUE. not supported on GPU'
        ierrstat = 1002
    END IF

    ! Lines 4652
    logbuf (58) = lphenology                                                                                                        ! Chur 2022

    ! Lines 4678
    CALL distribute_values (logbuf, 58, 0, imp_logical,  icomm_world, ierr)                                                         ! Chur 2022

    ! Lines 4798
    lphenology               = logbuf (58)

    ! Lines 4946 - 4947
    WRITE (nuspecif, '(T8,A,T33,L12  ,T52,L12  ,T71,A3)')                      &                                                    ! Chur 2022
                                  'lphenology',lphenology,lphenology_d,' L '

!==============================================================================!





4. Added data to src_allocation
#===============================================================================
    ! V6_0 CLM     2021-02-28 Churiulin Evgenii
    !  Added 30 parameters for the new vegetation scheme
    !  Added logical parameter lphenology for the new vegetation parameterisation scheme.

    ! Lines 334
    ! 3. controlling the physics
    USE data_runcontrol , ONLY :   &
        lphenology,   & ! forecast with new vegetation scheme
    ! --------------------------

    ! Lines 1010 - 1056
    USE data_fields     , ONLY :   &
    ! 14. fields for the vegetation parameterization                                                                                    ! Chur 2022 - new parameters for vegetation
    ! -----------------------------------------
        !   Additional fields for the radiation correction scheme
        sfldir_par    , & ! direct component of PAR at the ground            ( W/m2)
        sfldifd_par   , & ! diffuse downward component of PAR at the ground  ( W/m2)
        sfldifu_par   , & ! diffuse upward component of PAR at the ground    ( W/m2)
        sfltrdir_par  , & ! direct comp. of PAR at surface                   ( W/m2)
        sfltrdifd_par , & ! diffuse downward comp. of PAR at surface         ( W/m2)
        sfltrdifu_par , & ! diffuse upward comp. of PAR at surface           ( W/m2)
        cos_zen_ang   , & ! actual cosine of zenit angle                     ( rad )
        ! these are accumulated values
        asfldir_par   , & ! direct component of PAR at the ground            ( W/m2)
        asfldifd_par  , & ! diffuse downward component of PAR at the ground  ( W/m2)
        asfldifu_par  , & ! diffuse upward component of PAR at the ground    ( W/m2)
        ! Additional fiels for vegetation scheme are computed in the dynamics and / or physics
        ztraleav      , & ! transpiration rate of dry leaves surface         ( mm/hour)
        zverbo        , & ! total evapotranspiration                         ( mm/hour)
        lai_sun       , & ! sunlit leaf area index                           ( m2/m2)
        lai_sha       , & ! shaded leaf area index                           ( m2/m2)
        sla_sun       , & ! specific leaf area index for sunlit leaves       ( m2/m2)
        sla_sha       , & ! specific leaf area index for shaded leaves       ( m2/m2)
        vcmax_sun     , & ! maximum rate of carboxylation for sunlit leaves  ( umol CO2 m^-2 s^-1)
        vcmax_sha     , & ! maximum rate of carboxylation for shaded leaves  ( umol CO2 m^-2 s^-1)
        par_sun       , & ! photosynthetic active radiation for sunlit leaf  ( W/m2)
        par_sha       , & ! photosynthetic active radiation for shaded leaf  ( W/m2)
        rs_leaf       , & ! canopy leaf resistance                           ( s/m)
        psn           , & ! leaf photosynthesis                              ( umol CO2 m^-2 s^-1)
        sur_lai       , & ! leaf area index depending on biomass evolution   ( m2/m2)
        sur_panday    , & ! parameter for daily values of photosyntesis      ( kg CO2 m-2 s-1)
        sur_biomass   , & ! parameter for surfex biomass                     ( kg of dry matter m^2)
        sur_panfm     , & ! parameter for surfex maximum leaf assimilation   ( kg CO2 m-2 s-1)
        gpp_flux      , & ! gross primary production                         ( gC/m2/s)
        npp_flux      , & ! net primary production                           ( gC/m2/s)
        !   field that are computed for vegetation (average)
        aztraleav     , & ! average transpiration rate of dry leaves         ( mm/day)
        azverbo           ! average total evapotranspiration                 ( mm/day)
    ! end of data_fields
    !------------------------------------------------------------------------------


    ! Lines 1726 - 1760
    ! from the new vegetation parameterisation scheme - Churiulin, CESR
    ALLOCATE ( ztraleav     (ie,je,0:ntiles)     , STAT=izl ) ; ztraleav     = 0.0_wp  ; ist = ist + izl
    ALLOCATE ( zverbo       (ie,je,0:ntiles)     , STAT=izl ) ; zverbo       = 0.0_wp  ; ist = ist + izl
    ALLOCATE ( aztraleav    (ie,je)              , STAT=izl ) ; aztraleav    = 0.0_wp  ; ist = ist + izl
    ALLOCATE ( azverbo      (ie,je)              , STAT=izl ) ; azverbo      = 0.0_wp  ; ist = ist + izl
    IF (lphenology) THEN
        ! Additional parameters for radiation scheme
        ALLOCATE ( sfldir_par   (ie,je)            , STAT=izl ) ; sfldir_par   = 0.0_wp  ; ist = ist + izl
        ALLOCATE ( sfldifd_par  (ie,je)            , STAT=izl ) ; sfldifd_par  = 0.0_wp  ; ist = ist + izl
        ALLOCATE ( sfldifu_par  (ie,je)            , STAT=izl ) ; sfldifu_par  = 0.0_wp  ; ist = ist + izl
        ALLOCATE ( sfltrdir_par (ie,je)            , STAT=izl ) ; sfltrdir_par = 0.0_wp  ; ist = ist + izl
        ALLOCATE ( sfltrdifd_par(ie,je)            , STAT=izl ) ; sfltrdifd_par= 0.0_wp  ; ist = ist + izl
        ALLOCATE ( sfltrdifu_par(ie,je)            , STAT=izl ) ; sfltrdifu_par= 0.0_wp  ; ist = ist + izl
        ALLOCATE ( asfldir_par  (ie,je)            , STAT=izl ) ; asfldir_par  = 0.0_wp  ; ist = ist + izl
        ALLOCATE ( asfldifd_par (ie,je)            , STAT=izl ) ; asfldifd_par = 0.0_wp  ; ist = ist + izl
        ALLOCATE ( asfldifu_par (ie,je)            , STAT=izl ) ; asfldifu_par = 0.0_wp  ; ist = ist + izl
        ALLOCATE ( cos_zen_ang  (ie,je)            , STAT=izl ) ; cos_zen_ang  = 0.0_wp  ; ist = ist + izl
        ! Parameters for the new vegetation scheme
        ALLOCATE ( sur_lai      (ie,je,0:ntiles)   , STAT=izl ) ; sur_lai      = 0.0_wp  ; ist = ist + izl
        ALLOCATE ( sur_panday   (ie,je,0:ntiles)   , STAT=izl ) ; sur_panday   = 0.0_wp  ; ist = ist + izl
        ALLOCATE ( sur_biomass  (ie,je,0:ntiles)   , STAT=izl ) ; sur_biomass  = 0.0_wp  ; ist = ist + izl
        ALLOCATE ( sur_panfm    (ie,je,0:ntiles)   , STAT=izl ) ; sur_panfm    = 0.0_wp  ; ist = ist + izl
        ALLOCATE ( lai_sun      (ie,je,0:ntiles)   , STAT=izl ) ; lai_sun      = 0.0_wp  ; ist = ist + izl
        ALLOCATE ( lai_sha      (ie,je,0:ntiles)   , STAT=izl ) ; lai_sha      = 0.0_wp  ; ist = ist + izl
        ALLOCATE ( sla_sun      (ie,je,0:ntiles)   , STAT=izl ) ; sla_sun      = 0.0_wp  ; ist = ist + izl
        ALLOCATE ( sla_sha      (ie,je,0:ntiles)   , STAT=izl ) ; sla_sha      = 0.0_wp  ; ist = ist + izl
        ALLOCATE ( vcmax_sun    (ie,je,0:ntiles)   , STAT=izl ) ; vcmax_sun    = 0.0_wp  ; ist = ist + izl
        ALLOCATE ( vcmax_sha    (ie,je,0:ntiles)   , STAT=izl ) ; vcmax_sha    = 0.0_wp  ; ist = ist + izl
        ALLOCATE ( par_sun      (ie,je,0:ntiles)   , STAT=izl ) ; par_sun      = 0.0_wp  ; ist = ist + izl
        ALLOCATE ( par_sha      (ie,je,0:ntiles)   , STAT=izl ) ; par_sha      = 0.0_wp  ; ist = ist + izl
        ALLOCATE ( rs_leaf      (ie,je,0:ntiles)   , STAT=izl ) ; rs_leaf      = 0.0_wp  ; ist = ist + izl
        ALLOCATE ( psn          (ie,je,0:ntiles)   , STAT=izl ) ; psn          = 0.0_wp  ; ist = ist + izl
        ALLOCATE ( gpp_flux     (ie,je,0:ntiles)   , STAT=izl ) ; gpp_flux     = 0.0_wp  ; ist = ist + izl
        ALLOCATE ( npp_flux     (ie,je,0:ntiles)   , STAT=izl ) ; npp_flux     = 0.0_wp  ; ist = ist + izl
    ENDIF

    ! Lines 2573 - 2607
    ! from the new vegetation parameterisation scheme - Churiulin, CESR
    DEALLOCATE ( ztraleav    , STAT=istat )
    DEALLOCATE ( zverbo      , STAT=istat )
    DEALLOCATE ( aztraleav   , STAT=istat )
    DEALLOCATE ( azverbo     , STAT=istat )
    IF (lphenology) THEN
        ! Additional parameters for radiation scheme
        DEALLOCATE ( sfldir_par   , STAT=istat )
        DEALLOCATE ( sfldifd_par  , STAT=istat )
        DEALLOCATE ( sfldifu_par  , STAT=istat )
        DEALLOCATE ( sfltrdir_par , STAT=istat )
        DEALLOCATE ( sfltrdifd_par, STAT=istat )
        DEALLOCATE ( sfltrdifu_par, STAT=istat )
        DEALLOCATE ( asfldir_par  , STAT=istat )
        DEALLOCATE ( asfldifd_par , STAT=istat )
        DEALLOCATE ( asfldifu_par , STAT=istat )
        DEALLOCATE ( cos_zen_ang  , STAT=istat )
        ! Parameters for the new vegetation scheme
        DEALLOCATE ( sur_lai     , STAT=istat )
        DEALLOCATE ( sur_panday  , STAT=istat )
        DEALLOCATE ( sur_biomass , STAT=istat )
        DEALLOCATE ( sur_panfm   , STAT=istat )
        DEALLOCATE ( lai_sun     , STAT=istat )
        DEALLOCATE ( lai_sha     , STAT=istat )
        DEALLOCATE ( sla_sun     , STAT=istat )
        DEALLOCATE ( sla_sha     , STAT=istat )
        DEALLOCATE ( vcmax_sun   , STAT=istat )
        DEALLOCATE ( vcmax_sha   , STAT=istat )
        DEALLOCATE ( par_sun     , STAT=istat )
        DEALLOCATE ( par_sha     , STAT=istat )
        DEALLOCATE ( rs_leaf     , STAT=istat )
        DEALLOCATE ( psn         , STAT=istat )
        DEALLOCATE ( gpp_flux    , STAT=istat )
        DEALLOCATE ( npp_flux    , STAT=istat )
    ENDIF


!==============================================================================!





5. Added data to src_setup_vartab:
#===============================================================================
    ! V6_0 CLM     2021-02-28 Churiulin Evgenii
    !  Added 27 parameters for the new vegetation scheme
    !  Added logical parameter lphenology for the new vegetation parameterisation scheme.

    ! Lines 239
    USE data_runcontrol , ONLY :   &
        lphenology,   & ! forecast with new vegetation scheme


    ! Lines 667 - 699
    USE data_fields,        ONLY:                              &
    ! 14. fields for the new vegetation parameterization scheme                                                                     ! Chur 2022
    ! -----------------------------------------
        !   Additional fields for the radiation correction scheme
        sfldir_par    , & ! direct comp. of PAR at the ground                ( W/m2)
        sfldifd_par   , & ! diffuse downward comp. of PAR at the ground      ( W/m2)
        sfldifu_par   , & ! diffuse upward comp. of PAR at the ground        ( W/m2)
        cos_zen_ang   , & ! actual cosine of zenit angle                     ( rad )
        ! these are accumulated values
        asfldir_par   , & ! direct comp. of PAR at the ground                ( W/m2)
        asfldifd_par  , & ! diffuse downward comp. of PAR at the ground      ( W/m2)
        asfldifu_par  , & ! diffuse upward comp. of PAR at the ground        ( W/m2)
        ! Parameters for the new vegetation scheme
        sur_lai       , & ! leaf area index depending on biomass evolution   ( m2/m2)
        sur_panday    , & ! parameter for daily values of photosyntesis      ( kg CO2 m-2 s-1)
        sur_biomass   , & ! parameter for surfex biomass                     ( kg of dry matter m^2)
        sur_panfm     , & ! parameter for surfex maximum leaf assimilation   ( kg CO2 m-2 s-1)
        lai_sun       , & ! sunlit leaf area index                           ( m2/m2)
        lai_sha       , & ! shaded leaf area index                           ( m2/m2)
        sla_sun       , & ! specific leaf area index for sunlit leaves       ( m2/m2)
        sla_sha       , & ! specific leaf area index for shaded leaves       ( m2/m2)
        vcmax_sun     , & ! maximum rate of carboxylation for sunlit leaves  ( umol CO2 m^-2 s^-1)
        vcmax_sha     , & ! maximum rate of carboxylation for shaded leaves  ( umol CO2 m^-2 s^-1)
        par_sun       , & ! photosynthetic active radiation for sunlit leaf  ( W/m2)
        par_sha       , & ! photosynthetic active radiation for shaded leaf  ( W/m2)
        rs_leaf       , & ! canopy leaf resistance                           ( s/m)
        psn           , & ! leaf photosynthesis                              ( umol CO2 m^-2 s^-1)
        gpp_flux      , & ! gross primary production                         ( gC/m2/s)
        npp_flux      , & ! net primary production                           ( gC/m2/s)
        ztraleav      , & ! transpiration rate of dry leaves surface         ( mm/hour)
        zverbo        , & ! total evapotranspiration                         ( mm/hour)
        !   field that are computed for vegetation (average)
        aztraleav     , & ! average transpiration rate of dry leaves         ( mm/day)
        azverbo           ! average total evapotranspiration                 ( mm/day)

    ! Lines 1384 - 1451
    ! from the new vegetation parameterisation scheme - Churiulin, CESR
    var(1,202,3)= ar_des('ZTRALEAV      ',      1, 0,  0,   1.0_wp    , 0.0_wp,   0,   3,  -1,        dum5, dum4,  dum4,  ztraleav   , dum3, dum2   ,1,&
          'mm/hour        ','-',                                          'transpiration from dry leaf surface'       , ' ',  .TRUE.,.FALSE.)
    var(1,236,3)= ar_des('ZVERBO        ',      1, 0,  0,   1.0_wp    , 0.0_wp,   0,   3,  -1,        dum5, dum4,  dum4,  zverbo     , dum3, dum2   ,1,&
          'mm/hour        ','-',                                                     'total evapotranspiration'       , ' ',  .TRUE.,.FALSE.)
    !  Accumulated fields for vegetation parameters
    var(1,237,3)= ar_des('AZTRALEAV     ',      1, 0,  0,   1.0_wp    , 0.0_wp,   3,   2,  -1,        dum5, dum4,  dum4,  dum3,  dum3, aztraleav    ,1,&
          'mm/day         ','-',                              'accumulated transpiration rate of dry leaves'          , ' ', .FALSE.,.FALSE.)
    var(1,238,3)= ar_des('AZVERBO       ',      1, 0,  0,   1.0_wp    , 0.0_wp,   3,   2,  -1,        dum5, dum4,  dum4,  dum3,  dum3, azverbo      ,1,&
          'mm/day         ','-',                                      'accumulated total evapotranspiration'          , ' ', .FALSE.,.FALSE.)
    IF     (lphenology) THEN                                                                                                            ! Chur 2022
    ! Parameters from radiation scheme
    var(1,239,3)= ar_des('SFLDIR_PAR    ',      1, 0,  0,   1.0_wp    , 0.0_wp,   0,   2,  -1,        dum5, dum4,  dum4,  dum3,   dum3, sfldir_par  ,1,&
          'W m-2          ','-',                                'direct component of PAR flux at the ground'          , ' ', .FALSE.,.FALSE.)
    var(2,240,1)= ar_des('SFLDIFD_PAR   ',      1, 0,  0,   1.0_wp    , 0.0_wp,   0,   2,  -1,        dum5, dum4,  dum4,  dum3,   dum3, sfldifd_par ,1,&
          'W m-2          ','-',                         'diffuse downward component of PAR flux at the ground'       , ' ', .FALSE.,.FALSE.)
    var(2,240,2)= ar_des('SFLDIFU_PAR   ',      1, 0,  0,   1.0_wp    , 0.0_wp,   0,   2,  -1,        dum5, dum4,  dum4,  dum3,   dum3, sfldifu_par ,1,&
          'W m-2          ','-',                                'diffuse upward component of PAR at the ground'       , ' ', .FALSE.,.FALSE.)
    var(2,241,1)= ar_des('COSZ          ',      1, 0,  0,   1.0_wp    , 0.0_wp,   0,   2,  -1,        dum5, dum4,  dum4,  dum3,   dum3, cos_zen_ang ,1,&
          'rad            ','-',                                                 'actual cosine of zenit angle'       , ' ', .FALSE.,.FALSE.)
    ! Average fields
    var(2,242,3)= ar_des('ASFLDIR_PAR   ',      1, 0,  0,   1.0_wp    , 0.0_wp,   3,   2,  -1,        dum5, dum4,  dum4,  dum3,   dum3, asfldir_par ,1,&
          'W m-2          ','-',                      'average direct component of PAR flux at the ground'            , ' ', .FALSE.,.FALSE.)
    var(2,243,1)= ar_des('ASFLDIFD_PAR  ',      1, 0,  0,   1.0_wp    , 0.0_wp,   3,   2,  -1,        dum5, dum4,  dum4,  dum3,   dum3, asfldir_par ,1,&
          'W m-2          ','-',            'average diffuse downward component of PAR flux at the ground'            , ' ', .FALSE.,.FALSE.)
    var(2,243,2)= ar_des('ASFLDIFU_PAR  ',      1, 0,  0,   1.0_wp    , 0.0_wp,   3,   2,  -1,        dum5, dum4,  dum4,  dum3,   dum3, asfldifu_par,1,&
          'W m-2          ','-',                   'average diffuse upward component of PAR at the ground'            , ' ', .FALSE.,.FALSE.)
    ! Parameters from the new vegetation scheme - tiles parameters
    var(2,243,3)= ar_des('SUR_LAI       ',      1, 0,  0,   1.0_wp    , 0.0_wp,   0,   3,  -1,        dum5, dum4,  dum4,  sur_lai    , dum3, dum2   ,1,&
           'm2/m2         ','-',                                                     'leaf area index isba'           , ' ',  .TRUE.,.FALSE.)
    var(1,244,1)= ar_des('SUR_PANFM     ',      1, 0,  0,   1.0_wp    , 0.0_wp,   0,   3,  -1,        dum5, dum4,  dum4,  sur_panfm  , dum3, dum2   ,1,&
          'kg CO2 m-2 s-1 ','-',                                       'maximum rate of leaf assimilation'            , ' ',  .TRUE.,.FALSE.)
    var(2,244,1)= ar_des('SUR_PANDAY    ',      1, 0,  0,   1.0_wp    , 0.0_wp,   0,   3,  -1,        dum5, dum4,  dum4,  sur_panday , dum3, dum2   ,1,&
          'kg CO2 m-2 s-1 ','-',                                           'daily leaf photosynthesis isba'           , ' ',  .TRUE.,.FALSE.)
    var(1,244,2)= ar_des('SUR_BIOMASS   ',      1, 0,  0,   1.0_wp    , 0.0_wp,   0,   3,  -1,        dum5, dum4,  dum4,  sur_biomass, dum3, dum2   ,1,&
          'kg m2          ','-',                                                            'biomass isba'            , ' ',  .TRUE.,.FALSE.)
    var(1,244,3)= ar_des('SLA_SUN       ',      1, 0,  0,   1.0_wp    , 0.0_wp,   0,   3,  -1,        dum5, dum4,  dum4,  sla_sun    , dum3, dum2   ,1,&
           'm2/gC         ','-',                                  'specific leaf area index for sunlit leaves'        , ' ',  .TRUE.,.FALSE.)
    var(1,245,1)= ar_des('SLA_SHA       ',      1, 0,  0,   1.0_wp    , 0.0_wp,   0,   3,  -1,        dum5, dum4,  dum4,  sla_sha    , dum3, dum2   ,1,&
           'm2/gC         ','-',                                  'specific leaf area index for shaded leaves'        , ' ',  .TRUE.,.FALSE.)
    var(1,245,2)= ar_des('LAI_SUN       ',      1, 0,  0,   1.0_wp    , 0.0_wp,   0,   3,  -1,        dum5, dum4,  dum4,  lai_sun    , dum3, dum2   ,1,&
           'm2/m2         ','-',                                                       'sunlit leaf area index'       , ' ',  .TRUE.,.FALSE.)
    var(1,245,3)= ar_des('LAI_SHA       ',      1, 0,  0,   1.0_wp    , 0.0_wp,   0,   3,  -1,        dum5, dum4,  dum4,  lai_sha    , dum3, dum2   ,1,&
           'm2/m2         ','-',                                                      'shaded leaf area index'        , ' ',  .TRUE.,.FALSE.)
    var(1,246,1)= ar_des('VC_SUN        ',      1, 0,  0,   1.0_wp    , 0.0_wp,   0,   3,  -1,        dum5, dum4,  dum4,  vcmax_sun  , dum3, dum2   ,1,&
           'umol CO2/m2/s ','-',                           'maximum rate of carboxylation for sunlit leaves'          , ' ',  .TRUE.,.FALSE.)
    var(1,246,2)= ar_des('VC_SHA        ',      1, 0,  0,   1.0_wp    , 0.0_wp,   0,   3,  -1,        dum5, dum4,  dum4,  vcmax_sha  , dum3, dum2   ,1,&
           'umol CO2/m2/s ','-',                           'maximum rate of carboxylation for shaded leaves'          , ' ',  .TRUE.,.FALSE.)
    var(1,246,3)= ar_des('PAR_SUN       ',      1, 0,  0,   1.0_wp    , 0.0_wp,   0,   3,  -1,        dum5, dum4,  dum4,  par_sun    , dum3, dum2   ,1,&
           'W m-2         ','surface_downwelling_photosynthetic_radiative_flux_in_air',                                          &
                                                                  'surface photosynthetic active radiation'           , ' ',  .TRUE.,.FALSE.)
    var(1,247,1)= ar_des('PAR_SHA       ',      1, 0,  0,   1.0_wp    , 0.0_wp,   0,   3,  -1,        dum5, dum4,  dum4,  par_sha    , dum3, dum2   ,1,&
           'W m-2         ','surface_downwelling_photosynthetic_radiative_flux_in_air',                                          &
                                                                  'surface photosynthetic active radiation'           , ' ',  .TRUE.,.FALSE.)
    var(1,247,2)= ar_des('RS_LEAF       ',      1, 0,  0,   1.0_wp    , 0.0_wp,   0,   3,  -1,        dum5, dum4,  dum4,  rs_leaf    , dum3, dum2   ,1,&
           's/m           ','-',                                                   'canopy layer resistance'          , ' ',  .TRUE.,.FALSE.)
    var(1,247,3)= ar_des('LEAF_PSN      ',      1, 0,  0,   1.0_wp    , 0.0_wp,   0,   3,  -1,        dum5, dum4,  dum4,  psn        , dum3, dum2   ,1,&
           'umol CO2 /m2/s','-',                                                      'leaf photosynthesis'           , ' ',  .TRUE.,.FALSE.)
    var(1,248,1)= ar_des('GPP           ',      1, 0,  0,   1.0_wp    , 0.0_wp,   0,   3,  -1,        dum5, dum4,  dum4,  gpp_flux   , dum3, dum2   ,1,&
           'gC/m2/s       ','-',                                                'gross primary production'            , ' ',  .TRUE.,.FALSE.)
    var(1,248,2)= ar_des('NPP           ',      1, 0,  0,   1.0_wp    , 0.0_wp,   0,   3,  -1,        dum5, dum4,  dum4,  npp_flux   , dum3, dum2   ,1,&
           'gC/m2         ','-',                                                  'net primary production'            , ' ',  .TRUE.,.FALSE.)
    ENDIF





6. Add to data_block_fields:
#===============================================================================
    ! V6.0 CLM     2022-03-07 Evgenii Churiulin
    !  Added tile-dimension to variables for the new vegetation scheme:
    !      (lai_b      , sfldir_par_b, sfldifd_par_b, sfldifu_par_b, cos_zen_ang_b,
    !       sun_el     , sur_lai_b   , sur_panday_b , sur_biomass_b, sur_panfm_b  ,
    !       ztraleav_b , zverbo_b    , lai_sun_b    , lai_sha_b    , sla_sun_b    ,
    !       sla_sha_b  , vcmax_sun_b , vcmax_sha_b  , par_sun_b    , par_sha_b    ,
    !       rs_leaf_b  , psn_b       , gpp_flux_b   , npp_flux_b                  )

    ! Lines 594 - 622
    ! 11. fields for the vegetation parameterization  (lphenology = .TRUE.)
    ! -----------------------------------------
      REAL (KIND=wp), ALLOCATABLE, TARGET     ::   &
        ! Additional external parameter fields
        lai_b        (:,:),     & ! leaf area index                                          --
        !   Additional fields for the radiation correction scheme
        sfldir_par_b (:),         & ! direct component of PAR at the ground            ( W/m2 )
        sfldifd_par_b(:),         & ! diffuse downward component of PAR at the ground  ( W/m2 )
        sfldifu_par_b(:),         & ! diffuse upward component of PAR at the ground    ( W/m2 )
        cos_zen_ang_b(:),         & ! actual cosine of zenit angle                     ( rad  )
        sun_el_b     (:),         & ! sun declination                                         (deg)
        ! Additional field for TERRA
        sur_lai_b    (:,:),       & ! leaf area index depending on biomass evolution   ( m2/m2  )
        sur_panday_b (:,:),       & ! parameter for daily values of photosyntesis      ( kg CO2 m-2 s-1)
        sur_biomass_b(:,:),       & ! parameter for surfex biomass                     ( kg of dry matter m^2)
        sur_panfm_b  (:,:),       & ! parameter for surfex maximum leaf assimilation   ( kg CO2 m-2 s-1)
        ztraleav_b   (:,:),       & ! transpiration rate of dry leaves surface         ( mm/hour)
        zverbo_b     (:,:),       & ! total evapotranspiration
        lai_sun_b    (:,:),       & ! sunlit leaf area index                           ( m2/m2  )
        lai_sha_b    (:,:),       & ! shaded leaf area index                           ( m2/m2  )
        sla_sun_b    (:,:),       & ! specific leaf area index for sunlit leaves       ( m2/m2  )
        sla_sha_b    (:,:),       & ! specific leaf area index for shaded leaves       ( m2/m2  )
        vcmax_sun_b  (:,:),       & ! maximum rate of carboxylation for sunlit leaves  ( umol CO2 m^-2 s^-1)
        vcmax_sha_b  (:,:),       & ! maximum rate of carboxylation for shaded leaves  ( umol CO2 m^-2 s^-1)
        par_sun_b    (:,:),       & ! photosynthetic active radiation for sunlit leaf  ( W/m2)
        par_sha_b    (:,:),       & ! photosynthetic active radiation for shaded leaf  ( W/m2)
        rs_leaf_b    (:,:),       & ! canopy leaf resistance                           ( s/m)
        psn_b        (:,:),       & ! leaf photosynthesis                              ( umol CO2 m^-2 s^-1)
        gpp_flux_b   (:,:),       & ! gross primary production                         ( gC/m2/s )
        npp_flux_b   (:,:)          ! net primary production                           ( gC/m2/s )

!==============================================================================!





7. Added data to src_block_fields_org
#===============================================================================
    ! V6.0 CLM     2022-03-07 Evgenii Churiulin
    !  Added tile-dimension to variables for the new vegetation scheme:
    !      (lai_b      , sfldir_par_b, sfldifd_par_b, sfldifu_par_b, cos_zen_ang_b,
    !       sun_el     , sur_lai_b   , sur_panday_b , sur_biomass_b, sur_panfm_b  ,
    !       ztraleav_b , zverbo_b    , lai_sun_b    , lai_sha_b    , sla_sun_b    ,
    !       sla_sha_b  , vcmax_sun_b , vcmax_sha_b  , par_sun_b    , par_sha_b    ,
    !       rs_leaf_b  , psn_b       , gpp_flux_b   , npp_flux_b                  )
    !  Added logical parameter lphenology for the new vegetation parameterisation scheme.

    ! Lines 124
    USE data_runcontrol,  ONLY : ..., lphenology, ...

    ! Lines 209 - 214
    ! from new vegetation - Chur 2022
    USE data_fields,      ONLY : &
        lai       , sfldir_par, sfldifd_par, sfldifu_par, cos_zen_ang,     &
        sun_el    , sur_lai   , sur_panday , sur_biomass, sur_panfm  ,     &
        ztraleav  , zverbo    , lai_sun    , lai_sha    , sla_sun    ,     &
        sla_sha   , vcmax_sun , vcmax_sha  , par_sun    , par_sha    ,     &
        rs_leaf   , psn       , gpp_flux   , npp_flux   ,                  &

    ! Lines 1114 - 1143
    ! from new vegetation parameterisation scheme - Churiulin, CESR
    ALLOCATE(ztraleav_b     (nproma,0:ntiles) , STAT=izl); ztraleav_b       = r_init_val; ist=ist+izl
    ALLOCATE(zverbo_b       (nproma,0:ntiles) , STAT=izl); zverbo_b         = r_init_val; ist=ist+izl
    IF (lphenology) THEN
        ! external parameters, needed in several parameterizations
        ALLOCATE(lai_b        (nproma,0:ntiles) , STAT=izl); lai_b          = r_init_val; ist=ist+izl
        ! from radiation scheme
        ALLOCATE(sfldir_par_b  (nproma)         , STAT=izl); sfldir_par_b   = r_init_val; ist=ist+izl
        ALLOCATE(sfldifd_par_b (nproma)         , STAT=izl); sfldifd_par_b  = r_init_val; ist=ist+izl
        ALLOCATE(sfldifu_par_b (nproma)         , STAT=izl); sfldifu_par_b  = r_init_val; ist=ist+izl
        ALLOCATE(cos_zen_ang_b (nproma)         , STAT=izl); cos_zen_ang_b  = r_init_val; ist=ist+izl
        ALLOCATE(sun_el_b      (nproma)         , STAT=izl); sun_el_b       = r_init_val; ist=ist+izl
        ! from vegetation scheme
        ALLOCATE(sur_lai_b     (nproma,0:ntiles), STAT=izl); sur_lai_b      = r_init_val; ist=ist+izl
        ALLOCATE(sur_panday_b  (nproma,0:ntiles), STAT=izl); sur_panday_b   = r_init_val; ist=ist+izl
        ALLOCATE(sur_biomass_b (nproma,0:ntiles), STAT=izl); sur_biomass_b  = r_init_val; ist=ist+izl
        ALLOCATE(sur_panfm_b   (nproma,0:ntiles), STAT=izl); sur_panfm_b    = r_init_val; ist=ist+izl
        ALLOCATE(lai_sun_b     (nproma,0:ntiles), STAT=izl); lai_sun_b      = r_init_val; ist=ist+izl
        ALLOCATE(lai_sha_b     (nproma,0:ntiles), STAT=izl); lai_sha_b      = r_init_val; ist=ist+izl
        ALLOCATE(sla_sun_b     (nproma,0:ntiles), STAT=izl); sla_sun_b      = r_init_val; ist=ist+izl
        ALLOCATE(sla_sha_b     (nproma,0:ntiles), STAT=izl); sla_sha_b      = r_init_val; ist=ist+izl
        ALLOCATE(vcmax_sun_b   (nproma,0:ntiles), STAT=izl); vcmax_sun_b    = r_init_val; ist=ist+izl
        ALLOCATE(vcmax_sha_b   (nproma,0:ntiles), STAT=izl); vcmax_sha_b    = r_init_val; ist=ist+izl
        ALLOCATE(par_sun_b     (nproma,0:ntiles), STAT=izl); par_sun_b      = r_init_val; ist=ist+izl
        ALLOCATE(par_sha_b     (nproma,0:ntiles), STAT=izl); par_sha_b      = r_init_val; ist=ist+izl
        ALLOCATE(rs_leaf_b     (nproma,0:ntiles), STAT=izl); rs_leaf_b      = r_init_val; ist=ist+izl
        ALLOCATE(psn_b         (nproma,0:ntiles), STAT=izl); psn_b          = r_init_val; ist=ist+izl
        ALLOCATE(gpp_flux_b    (nproma,0:ntiles), STAT=izl); gpp_flux_b     = r_init_val; ist=ist+izl
        ALLOCATE(npp_flux_b    (nproma,0:ntiles), STAT=izl); npp_flux_b     = r_init_val; ist=ist+izl
    ENDIF

    ! Lines 1783 - 1810
    ! from new vegetation parameterisation scheme - Churiulin, CESR
    DEALLOCATE(ztraleav_b     , STAT=izl);     ist=ist+izl
    DEALLOCATE(zverbo_b       , STAT=izl);     ist=ist+izl
    IF (lphenology) THEN
        ! external parameters, needed in several parameterizations
        DEALLOCATE(lai_b          , STAT=izl);     ist=ist+izl
        ! from radiation scheme
        DEALLOCATE(sfldir_par_b   , STAT=izl);     ist=ist+izl
        DEALLOCATE(sfldifd_par_b  , STAT=izl);     ist=ist+izl
        DEALLOCATE(sfldifu_par_b  , STAT=izl);     ist=ist+izl
        DEALLOCATE(cos_zen_ang_b  , STAT=izl);     ist=ist+izl
        DEALLOCATE(sun_el_b       , STAT=izl);     ist=ist+izl
        ! from vegetation scheme
        DEALLOCATE(sur_lai_b      , STAT=izl);     ist=ist+izl
        DEALLOCATE(sur_panday_b   , STAT=izl);     ist=ist+izl
        DEALLOCATE(sur_biomass_b  , STAT=izl);     ist=ist+izl
        DEALLOCATE(sur_panfm_b    , STAT=izl);     ist=ist+izl
        DEALLOCATE(lai_sun_b      , STAT=izl);     ist=ist+izl
        DEALLOCATE(lai_sha_b      , STAT=izl);     ist=ist+izl
        DEALLOCATE(sla_sun_b      , STAT=izl);     ist=ist+izl
        DEALLOCATE(sla_sha_b      , STAT=izl);     ist=ist+izl
        DEALLOCATE(vcmax_sun_b    , STAT=izl);     ist=ist+izl
        DEALLOCATE(vcmax_sha_b    , STAT=izl);     ist=ist+izl
        DEALLOCATE(par_sun_b      , STAT=izl);     ist=ist+izl
        DEALLOCATE(par_sha_b      , STAT=izl);     ist=ist+izl
        DEALLOCATE(rs_leaf_b      , STAT=izl);     ist=ist+izl
        DEALLOCATE(psn_b          , STAT=izl);     ist=ist+izl
        DEALLOCATE(gpp_flux_b     , STAT=izl);     ist=ist+izl
        DEALLOCATE(npp_flux_b     , STAT=izl);     ist=ist+izl
    ENDIF

    ! Lines 2238 - 2269
    ! from new vegetation parameterisation scheme - Churiulin, CESR
    CALL register_block_field ("ztraleav"      , ztraleav       , ztraleav_b       )
    CALL register_block_field ("zverbo"        , zverbo         , zverbo_b         )
    IF (lphenology) THEN
        ! external parameters, needed in several parameterizations
        CALL register_block_field ("lai"         , lai            , lai_b          )
        ! from radiation scheme
        CALL register_block_field ("sfldir_par"  , sfldir_par     , sfldir_par_b   )
        CALL register_block_field ("sfldifd_par" , sfldifd_par    , sfldifd_par_b  )
        CALL register_block_field ("sfldifu_par" , sfldifu_par    , sfldifu_par_b  )
        CALL register_block_field ("cos_zen_ang" , cos_zen_ang    , cos_zen_ang_b  )
        CALL register_block_field ("sun_el"      , sun_el         , sun_el_b       )
        ! from vegetation scheme
        CALL register_block_field ("sur_lai"     , sur_lai        , sur_lai_b      )
        CALL register_block_field ("sur_panday"  , sur_panday     , sur_panday_b   )
        CALL register_block_field ("sur_biomass" , sur_biomass    , sur_biomass_b  )
        CALL register_block_field ("sur_panfm"   , sur_panfm      , sur_panfm_b    )
        CALL register_block_field ("lai_sun"     , lai_sun        , lai_sun_b      )
        CALL register_block_field ("lai_sha"     , lai_sha        , lai_sha_b      )
        CALL register_block_field ("sla_sun"     , sla_sun        , sla_sun_b      )
        CALL register_block_field ("sla_sha"     , sla_sha        , sla_sha_b      )
        CALL register_block_field ("vcmax_sun"   , vcmax_sun      , vcmax_sun_b    )
        CALL register_block_field ("vcmax_sha"   , vcmax_sha      , vcmax_sha_b    )
        CALL register_block_field ("par_sun"     , par_sun        , par_sun_b      )
        CALL register_block_field ("par_sha"     , par_sha        , par_sha_b      )
        CALL register_block_field ("rs_leaf"     , rs_leaf        , rs_leaf_b      )
        CALL register_block_field ("psn"         , psn            , psn_b          )
        CALL register_block_field ("gpp_flux"    , gpp_flux       , gpp_flux_b     )
        CALL register_block_field ("npp_flux"    , npp_flux       , npp_flux_b     )
    ENDIF
!==============================================================================!






8. Added to sfc_tile_approach
#===============================================================================
    ! V6.0 CLM     2022-03-07 Evgenii Churiulin
    !  Added tile-dimension to variables for the new vegetation scheme:
    !      (sur_lai_t   , sur_panday_t , sur_biomass_t, sur_panfm_t  , ztraleav_t ,
    !       zverbo_t    , lai_sun_t    , lai_sha_t    , sla_sun_t    , sla_sha_t  ,
    !       vcmax_sun_t , vcmax_sha_t  , par_sun_t    , par_sha_t    , rs_leaf_t  ,
    !       psn_t       , gpp_flux_t   , npp_flux_t                               )
    !  Added logical parameter lphenology for the new vegetation parameterisation scheme.

    ! Lines 155
    USE data_runcontrol,  ONLY:   &
        lphenology,   & ! forecast with new vegetation scheme

    ! Lines 417 - 421
    SUBROUTINE tile_average_ground                                            &
        (...
         ztraleav_t         , zverbo_t           , sur_lai_t    , sur_panday_t,&
         sur_biomass_t      , sur_panfm_t        , lai_sun_t    , lai_sha_t  , &
         sla_sun_t          , sla_sha_t          , vcmax_sun_t  , vcmax_sha_t, &
         par_sun_t          , par_sha_t          , rs_leaf_t    , psn_t      , &
         gpp_flux_t         , npp_flux_t                                       )

    ! Lines 508 - 509
    REAL(KIND=wp), INTENT(INOUT) ::                                            &
        ztraleav_t          (nvec,0:ntiles)  ,                                 &  ! Chur 2022
        zverbo_t            (nvec,0:ntiles)                                       ! Chur 2022

    ! Lines 511 - 528
    ! from new vegetation parameterisation scheme - Churiulin, CESR
    REAL(KIND=wp), OPTIONAL,  INTENT(INOUT) ::                                 &
        sur_lai_t           (nvec,0:ntiles)  ,                                 &
        sur_panday_t        (nvec,0:ntiles)  ,                                 &
        sur_biomass_t       (nvec,0:ntiles)  ,                                 &
        sur_panfm_t         (nvec,0:ntiles)  ,                                 &
        lai_sun_t           (nvec,0:ntiles)  ,                                 &
        lai_sha_t           (nvec,0:ntiles)  ,                                 &
        sla_sun_t           (nvec,0:ntiles)  ,                                 &
        sla_sha_t           (nvec,0:ntiles)  ,                                 &
        vcmax_sun_t         (nvec,0:ntiles)  ,                                 &
        vcmax_sha_t         (nvec,0:ntiles)  ,                                 &
        par_sun_t           (nvec,0:ntiles)  ,                                 &
        par_sha_t           (nvec,0:ntiles)  ,                                 &
        rs_leaf_t           (nvec,0:ntiles)  ,                                 &
        psn_t               (nvec,0:ntiles)  ,                                 &
        gpp_flux_t          (nvec,0:ntiles)  ,                                 &
        npp_flux_t          (nvec,0:ntiles)

    ! Lines 601 - 621
    ! from new vegetation parameterisation scheme - Churiulin, CESR
    ztraleav_t          (iv,0)  = 0.0_wp
    zverbo_t            (iv,0)  = 0.0_wp
    IF (lphenology) THEN
        sur_lai_t           (iv,0)  = 0.0_wp
        sur_panday_t        (iv,0)  = 0.0_wp
        sur_biomass_t       (iv,0)  = 0.0_wp
        sur_panfm_t         (iv,0)  = 0.0_wp
        lai_sun_t           (iv,0)  = 0.0_wp
        lai_sha_t           (iv,0)  = 0.0_wp
        sla_sun_t           (iv,0)  = 0.0_wp
        sla_sha_t           (iv,0)  = 0.0_wp
        vcmax_sun_t         (iv,0)  = 0.0_wp
        vcmax_sha_t         (iv,0)  = 0.0_wp
        par_sun_t           (iv,0)  = 0.0_wp
        par_sha_t           (iv,0)  = 0.0_wp
        rs_leaf_t           (iv,0)  = 0.0_wp
        psn_t               (iv,0)  = 0.0_wp
        gpp_flux_t          (iv,0)  = 0.0_wp
        npp_flux_t          (iv,0)  = 0.0_wp
    ENDIF

    ! Lines 722 - 744
    ! from new vegetation parameterisation scheme - Churiulin, CESR
    ztraleav_t        (iv,0) = ztraleav_t   (iv,0) + frc_t(iv,ntls) * ztraleav_t   (iv,ntls)
    zverbo_t          (iv,0) = zverbo_t     (iv,0) + frc_t(iv,ntls) * zverbo_t     (iv,ntls)
    IF (lphenology) THEN
        sur_lai_t     (iv,0) = sur_lai_t    (iv,0) + frc_t(iv,ntls) * sur_lai_t    (iv,ntls)
        sur_panday_t  (iv,0) = sur_panday_t (iv,0) + frc_t(iv,ntls) * sur_panday_t (iv,ntls)
        sur_biomass_t (iv,0) = sur_biomass_t(iv,0) + frc_t(iv,ntls) * sur_biomass_t(iv,ntls)
        sur_panfm_t   (iv,0) = sur_panfm_t  (iv,0) + frc_t(iv,ntls) * sur_panfm_t  (iv,ntls)
        lai_sun_t     (iv,0) = lai_sun_t    (iv,0) + frc_t(iv,ntls) * lai_sun_t    (iv,ntls)
        lai_sha_t     (iv,0) = lai_sha_t    (iv,0) + frc_t(iv,ntls) * lai_sha_t    (iv,ntls)
        sla_sun_t     (iv,0) = sla_sun_t    (iv,0) + frc_t(iv,ntls) * sla_sun_t    (iv,ntls)
        sla_sha_t     (iv,0) = sla_sha_t    (iv,0) + frc_t(iv,ntls) * sla_sha_t    (iv,ntls)
        vcmax_sun_t   (iv,0) = vcmax_sun_t  (iv,0) + frc_t(iv,ntls) * vcmax_sun_t  (iv,ntls)
        vcmax_sha_t   (iv,0) = vcmax_sha_t  (iv,0) + frc_t(iv,ntls) * vcmax_sha_t  (iv,ntls)
        par_sun_t     (iv,0) = par_sun_t    (iv,0) + frc_t(iv,ntls) * par_sun_t    (iv,ntls)
        par_sha_t     (iv,0) = par_sha_t    (iv,0) + frc_t(iv,ntls) * par_sha_t    (iv,ntls)
        rs_leaf_t     (iv,0) = rs_leaf_t    (iv,0) + frc_t(iv,ntls) * rs_leaf_t    (iv,ntls)
        psn_t         (iv,0) = psn_t        (iv,0) + frc_t(iv,ntls) * psn_t        (iv,ntls)
        gpp_flux_t    (iv,0) = gpp_flux_t   (iv,0) + frc_t(iv,ntls) * gpp_flux_t   (iv,ntls)
        npp_flux_t    (iv,0) = npp_flux_t   (iv,0) + frc_t(iv,ntls) * npp_flux_t   (iv,ntls)
    ENDIF

    ! LLines 910 - 930
    ! from new vegetation parameterisation scheme - Churiulin, CESR
    WRITE(*,'(A,3F28.16)') ' ztraleav           :  ', ztraleav_t     (iv,0), ztraleav_t     (iv,1), ztraleav_t     (iv,2)
    WRITE(*,'(A,3F28.16)') ' zverbo             :  ', zverbo_t       (iv,0), zverbo_t       (iv,1), zverbo_t       (iv,2)
    IF (lphenology) THEN
    WRITE(*,'(A,3F28.16)') ' sur_lai            :  ', sur_lai_t      (iv,0), sur_lai_t      (iv,1), sur_lai_t      (iv,2)
    WRITE(*,'(A,3F28.16)') ' sur_panday         :  ', sur_panday_t   (iv,0), sur_panday_t   (iv,1), sur_panday_t   (iv,2)
    WRITE(*,'(A,3F28.16)') ' sur_biomass        :  ', sur_biomass_t  (iv,0), sur_biomass_t  (iv,1), sur_biomass_t  (iv,2)
    WRITE(*,'(A,3F28.16)') ' sur_panfm          :  ', sur_panfm_t    (iv,0), sur_panfm_t    (iv,1), sur_panfm_t    (iv,2)
    WRITE(*,'(A,3F28.16)') ' lai_sun            :  ', lai_sun_t      (iv,0), lai_sun_t      (iv,1), lai_sun_t      (iv,2)
    WRITE(*,'(A,3F28.16)') ' lai_sha            :  ', lai_sha_t      (iv,0), lai_sha_t      (iv,1), lai_sha_t      (iv,2)
    WRITE(*,'(A,3F28.16)') ' sla_sun            :  ', sla_sun_t      (iv,0), sla_sun_t      (iv,1), sla_sun_t      (iv,2)
    WRITE(*,'(A,3F28.16)') ' sla_sha            :  ', sla_sha_t      (iv,0), sla_sha_t      (iv,1), sla_sha_t      (iv,2)
    WRITE(*,'(A,3F28.16)') ' vcmax_sun          :  ', vcmax_sun_t    (iv,0), vcmax_sun_t    (iv,1), vcmax_sun_t    (iv,2)
    WRITE(*,'(A,3F28.16)') ' vcmax_sha          :  ', vcmax_sha_t    (iv,0), vcmax_sha_t    (iv,1), vcmax_sha_t    (iv,2)
    WRITE(*,'(A,3F28.16)') ' par_sun            :  ', par_sun_t      (iv,0), par_sun_t      (iv,1), par_sun_t      (iv,2)
    WRITE(*,'(A,3F28.16)') ' par_sha            :  ', par_sha_t      (iv,0), par_sha_t      (iv,1), par_sha_t      (iv,2)
    WRITE(*,'(A,3F28.16)') ' rs_leaf            :  ', rs_leaf_t      (iv,0), rs_leaf_t      (iv,1), rs_leaf_t      (iv,2)
    WRITE(*,'(A,3F28.16)') ' psn                :  ', psn_t          (iv,0), psn_t          (iv,1), psn_t          (iv,2)
    WRITE(*,'(A,3F28.16)') ' gpp_flux           :  ', gpp_flux_t     (iv,0), gpp_flux_t     (iv,1), gpp_flux_t     (iv,2)
    WRITE(*,'(A,3F28.16)') ' npp_flux           :  ', npp_flux_t     (iv,0), npp_flux_t     (iv,1), npp_flux_t     (iv,2)
    ENDIF
!==============================================================================!





9. Add to sfc_interface
#===============================================================================
    ! V6.0 CLM     2022-03-14 Evgenii Churiulin
    !  Added new parameters for vegetation
    !  Added logical parameter lphenology for the new vegetation parameterisation scheme.

    ! Lines 177 - 180 and 198 - 202
    USE data_block_fields,ONLY :
        ! for vegetation - Churiulin, CESR ---> in and inout variables
        lai_b           , sfldir_par_b    , sfldifd_par_b   , sfldifu_par_b   , &
        cos_zen_ang_b   , rlat_b          , sun_el_b        , sur_lai_b       , &
        sur_panday_b    , sur_biomass_b   , sur_panfm_b                       , &
        ! for vegetation - Churiulin, CESR ---> ! output variables
        ztraleav_b      , zverbo_b        , lai_sun_b       , lai_sha_b       , &
        sla_sun_b       , sla_sha_b       , vcmax_sun_b     , vcmax_sha_b     , &
        par_sun_b       , par_sha_b       , rs_leaf_b       , psn_b           , &
        gpp_flux_b      , npp_flux_b      ,                                     &

    ! Lines 272
    USE data_runcontrol,  ONLY : ..., lphenology

    ! Lines 289
    USE sfc_phenology_data, ONLY:  phenology_wkarr_alloc, phenology_wkarr_dealloc

    ! Lines 418 - 426
    ! from new vegetation parameterisation scheme - Churiulin, CESR
    REAL (KIND=vpp), DIMENSION(:), ALLOCATABLE                 ::
        lai_t              , sfldir_par_t       , sfldifd_par_t      , sfldifu_par_t      ,  &
        cos_zen_ang_t      , rlat_t             , sun_el_t           ,  sur_lai_t         ,  &
        sur_panday_t       , sur_biomass_t      , sur_panfm_t        , ztraleav_t         ,  &
        zverbo_t           , lai_sun_t          , lai_sha_t          , sla_sun_t          ,  &
        sla_sha_t          , vcmax_sun_t        , vcmax_sha_t        , par_sun_t          ,  &
        par_sha_t          , rs_leaf_t          , psn_t              , gpp_flux_t         ,  &
        npp_flux_t

    ! Lines 653 - 656
    ! from new vegetation parameterisation scheme - Churiulin, CESR
    IF (lphenology) THEN
        CALL phenology_wkarr_alloc(nproma, istat)
    ENDIF

    ! Lines 827 - 857
    ! from new vegetation parameterisation scheme - Churiulin, CESR
    CALL register_copy (ztraleav_b        , sfcCopyList, copyToBlockF)! print *, 'ztraleav        ok'
    CALL register_copy (zverbo_b          , sfcCopyList, copyToBlockF)! print *, 'zverbo          ok'
    IF (lphenology) THEN
        ! Variables with in intent IN
        CALL register_copy (lai_b         , sfcCopyList, copyToBlockF)! print *, 'lai             ok'
        CALL register_copy (rlat_b        , sfcCopyList, copyToBlockF)! print *, 'rlat            ok'
        CALL register_copy (sun_el_b      , sfcCopyList, copyToBlockF)! print *, 'sun_el          ok'
        CALL register_copy (sfldir_par_b  , sfcCopyList, copyToBlockF)! print *, 'sfldir_par      ok'
        CALL register_copy (sfldifd_par_b , sfcCopyList, copyToBlockF)! print *, 'sfldifd_par     ok'
        CALL register_copy (sfldifu_par_b , sfcCopyList, copyToBlockF)! print *, 'sfldifu_par     ok'
        CALL register_copy (cos_zen_ang_b , sfcCopyList, copyToBlockF)! print *, 'cos_zen_ang     ok'
        ! Variables with in intent INOUT
        CALL register_copy (sur_lai_b     , sfcCopyList, copyToBlockF)! print *, 'sur_lai         ok'
        CALL register_copy (sur_panday_b  , sfcCopyList, copyToBlockF)! print *, 'sur_panday      ok'
        CALL register_copy (sur_biomass_b , sfcCopyList, copyToBlockF)! print *, 'sur_biomass     ok'
        CALL register_copy (sur_panfm_b   , sfcCopyList, copyToBlockF)! print *, 'sur_panfm       ok'
        ! Variables with in intent OUT
        CALL register_copy (lai_sun_b     , sfcCopyList, copyToBlockF)! print *, 'lai_sun         ok'
        CALL register_copy (lai_sha_b     , sfcCopyList, copyToBlockF)! print *, 'lai_sha         ok'
        CALL register_copy (sla_sun_b     , sfcCopyList, copyToBlockF)! print *, 'sla_sun         ok'
        CALL register_copy (sla_sha_b     , sfcCopyList, copyToBlockF)! print *, 'sla_sha         ok'
        CALL register_copy (vcmax_sun_b   , sfcCopyList, copyToBlockF)! print *, 'vcmax_sun       ok'
        CALL register_copy (vcmax_sha_b   , sfcCopyList, copyToBlockF)! print *, 'vcmax_sha       ok'
        CALL register_copy (par_sun_b     , sfcCopyList, copyToBlockF)! print *, 'par_sun         ok'
        CALL register_copy (par_sha_b     , sfcCopyList, copyToBlockF)! print *, 'par_sha         ok'
        CALL register_copy (rs_leaf_b     , sfcCopyList, copyToBlockF)! print *, 'rs_leaf         ok'
        CALL register_copy (psn_b         , sfcCopyList, copyToBlockF)! print *, 'psn             ok'
        CALL register_copy (gpp_flux_b    , sfcCopyList, copyToBlockF)! print *, 'gpp_flux        ok'
        CALL register_copy (npp_flux_b    , sfcCopyList, copyToBlockF)! print *, 'npp_flux        ok'
    ENDIF

    ! Lines 1025 - 1049
    ! from new vegetation parameterisation scheme - Churiulin, CESR
    CALL register_copy (ztraleav_b          , sfcCopyList, copyFromBlockF)! print *, 'ztraleav      ok'
    CALL register_copy (zverbo_b            , sfcCopyList, copyFromBlockF)! print *, 'zverbo        ok'
    IF (lphenology) THEN
        ! Variables with in intent INOUT
        CALL register_copy (sur_lai_b         , sfcCopyList, copyFromBlockF)! print *, 'sur_lai       ok'
        CALL register_copy (sur_panday_b      , sfcCopyList, copyFromBlockF)! print *, 'sur_panday    ok'
        CALL register_copy (sur_biomass_b     , sfcCopyList, copyFromBlockF)! print *, 'sur_biomass   ok'
        CALL register_copy (sur_panfm_b       , sfcCopyList, copyFromBlockF)! print *, 'sur_panfm     ok'
        ! Variables with in intent OUT
        CALL register_copy (lai_sun_b         , sfcCopyList, copyFromBlockF)! print *, 'lai_sun       ok'
        CALL register_copy (lai_sha_b         , sfcCopyList, copyFromBlockF)! print *, 'lai_sha       ok'
        CALL register_copy (sla_sun_b         , sfcCopyList, copyFromBlockF)! print *, 'sla_sun       ok'
        CALL register_copy (sla_sha_b         , sfcCopyList, copyFromBlockF)! print *, 'sla_sha       ok'
        CALL register_copy (vcmax_sun_b       , sfcCopyList, copyFromBlockF)! print *, 'vcmax_sun     ok'
        CALL register_copy (vcmax_sha_b       , sfcCopyList, copyFromBlockF)! print *, 'vcmax_sha     ok'
        CALL register_copy (par_sun_b         , sfcCopyList, copyFromBlockF)! print *, 'par_sun       ok'
        CALL register_copy (par_sha_b         , sfcCopyList, copyFromBlockF)! print *, 'par_sha       ok'
        CALL register_copy (rs_leaf_b         , sfcCopyList, copyFromBlockF)! print *, 'rs_leaf       ok'
        CALL register_copy (psn_b             , sfcCopyList, copyFromBlockF)! print *, 'psn           ok'
        CALL register_copy (gpp_flux_b        , sfcCopyList, copyFromBlockF)! print *, 'gpp_flux      ok'
        CALL register_copy (npp_flux_b        , sfcCopyList, copyFromBlockF)! print *, 'npp_flux      ok'
    ENDIF

    ! Lines 1206 - 1212
    ! from new vegetation parameterisation scheme - Churiulin, CESR
    REAL (KIND=vpp), DIMENSION(nproma)                         ::                              &
          lai_t          , sfldir_par_t   , sfldifd_par_t  , sfldifu_par_t  , cos_zen_ang_t  , &
          rlat_t         , sun_el_t       , sur_lai_t      , sur_panday_t   , sur_biomass_t  , &
          sur_panfm_t    , ztraleav_t     , zverbo_t       , lai_sun_t      , lai_sha_t      , &
          sla_sun_t      , sla_sha_t      , vcmax_sun_t    , vcmax_sha_t    , par_sun_t      , &
          par_sha_t      , rs_leaf_t      , psn_t          , gpp_flux_t     , npp_flux_t

    ! Lines 1571 - 1584
    ! from new vegetation parameterisation scheme - Churiulin, CESR
    IF (lphenology) THEN  
        lai_t        (il)    = REAL( lai_b           (ifull,ntls), vpp)
        rlat_t       (il)    = REAL( rlat_b          (ifull), vpp)
        sun_el_t     (il)    = REAL( sun_el_b        (ifull), vpp)
        sfldir_par_t (il)    = REAL( sfldir_par_b    (ifull), vpp)
        sfldifd_par_t(il)    = REAL( sfldifd_par_b   (ifull), vpp)
        sfldifu_par_t(il)    = REAL( sfldifu_par_b   (ifull), vpp)
        cos_zen_ang_t(il)    = REAL( cos_zen_ang_b   (ifull), vpp)
        sur_lai_t    (il)    = REAL( sur_lai_b       (ifull,ntls), vpp)
        sur_panday_t (il)    = REAL( sur_panday_b    (ifull,ntls), vpp)
        sur_biomass_t(il)    = REAL( sur_biomass_b   (ifull,ntls), vpp)
        sur_panfm_t  (il)    = REAL( sur_panfm_b     (ifull,ntls), vpp)
    ENDIF

    ! Lines 1919 - 1945
    ztraleav          = ztraleav_t    (:)       , & !OUT transpiration rate of dry leaves surface (mm/hour)                    ! Chur 2022
    zverbo            = zverbo_t      (:)       , & !OUT total evapotranspiration                 (mm/hour)                    ! Chur 2022
    ! from new vegetation parameterisation scheme - Churiulin, CESR
    ! optional parameters
    lai               = lai_t         (:)       , & ! IN leaf area index                                     --
    sfldir_par        = sfldir_par_t  (:)       , & ! IN direct component of PAR at the ground             (W/m2 )
    sfldifd_par       = sfldifd_par_t (:)       , & ! IN diffuse downward component of PAR at the ground   (W/m2 )
    sfldifu_par       = sfldifu_par_t (:)       , & ! IN diffuse upward component of PAR at the ground     (W/m2 )
    cos_zen_ang       = cos_zen_ang_t (:)       , & ! IN actual cosine of zenit angle                      (rad  )
    rlat              = rlat_t        (:)       , & ! IN latitide                                          (rad  )
    sun_el            = sun_el_t      (:)       , & ! IN sun declination                                   (deg  )
    sur_lai           = sur_lai_t     (:)       , & ! INOUT leaf area index depending on biomass evolution (m2/m2)
    sur_panday        = sur_panday_t  (:)       , & ! INOUT daily values of photosyntesis                  (kg CO2 m-2 s-1)
    sur_biomass       = sur_biomass_t (:)       , & ! INOUT leaf biomass                                   (kg of dry matter m^2)
    sur_panfm         = sur_panfm_t   (:)       , & ! INOUT maximum leaf assimilation                      (kg CO2 m-2 s-1)
    lai_sun           = lai_sun_t     (:)       , & ! OUT sunlit leaf area index                           (m2/m2)
    lai_sha           = lai_sha_t     (:)       , & ! OUT shaded leaf area index                           (m2/m2)
    sla_sun           = sla_sun_t     (:)       , & ! OUT specific leaf area index for sunlit leaves       (m2/m2)
    sla_sha           = sla_sha_t     (:)       , & ! OUT specific leaf area index for shaded leaves       (m2/m2)
    vcmax_sun         = vcmax_sun_t   (:)       , & ! OUT maximum rate of carboxylation for sunlit leaves  (umol co2/m**2/s)
    vcmax_sha         = vcmax_sha_t   (:)       , & ! OUT maximum rate of carboxylation for shaded leaves  (umol co2/m**2/s)
    par_sun           = par_sun_t     (:)       , & ! OUT photosynthetic active radiation for sunlit leaf  (W/m2)
    par_sha           = par_sha_t     (:)       , & ! OUT photosynthetic active radiation for shaded leaf  (W/m2)
    rs_leaf           = rs_leaf_t     (:)       , & ! OUT canopy leaf resistance                           (s/m)
    psn               = psn_t         (:)       , & ! OUT leaf photosynthesis
    gpp_flux          = gpp_flux_t    (:)       , & ! OUT gross primary production                         (gC/m2/s)
    npp_flux          = npp_flux_t    (:)       , & ! OUT net primary production                           (gC/m2/s)

    ! Lines 2048 - 2071
    ! from new vegetation parameterisation scheme - Churiulin, CESR
    ztraleav_b      (ifull,ntls)  = REAL( ztraleav_t     (il), wp)
    zverbo_b        (ifull,ntls)  = REAL( zverbo_t       (il), wp)
    IF (lphenology) THEN
        ! intent INOUT - 4 variables
        sur_lai_b     (ifull,ntls)  = REAL( sur_lai_t      (il), wp)
        sur_panday_b  (ifull,ntls)  = REAL( sur_panday_t   (il), wp)
        sur_biomass_b (ifull,ntls)  = REAL( sur_biomass_t  (il), wp)
        sur_panfm_b   (ifull,ntls)  = REAL( sur_panfm_t    (il), wp)
        ! intent OUT - 12 variables
        lai_sun_b     (ifull,ntls)  = REAL( lai_sun_t      (il), wp)
        lai_sha_b     (ifull,ntls)  = REAL( lai_sha_t      (il), wp)
        sla_sun_b     (ifull,ntls)  = REAL( sla_sun_t      (il), wp)
        sla_sha_b     (ifull,ntls)  = REAL( sla_sha_t      (il), wp)
        vcmax_sun_b   (ifull,ntls)  = REAL( vcmax_sun_t    (il), wp)
        vcmax_sha_b   (ifull,ntls)  = REAL( vcmax_sha_t    (il), wp)
        par_sun_b     (ifull,ntls)  = REAL( par_sun_t      (il), wp)
        par_sha_b     (ifull,ntls)  = REAL( par_sha_t      (il), wp)
        rs_leaf_b     (ifull,ntls)  = REAL( rs_leaf_t      (il), wp)
        psn_b         (ifull,ntls)  = REAL( psn_t          (il), wp)
        gpp_flux_b    (ifull,ntls)  = REAL( gpp_flux_t     (il), wp)
        npp_flux_b    (ifull,ntls)  = REAL( npp_flux_t     (il), wp)
    ENDIF

    ! Lines 2236 - 2241
    CALL tile_average_ground                                                  &
        ztraleav_b         , zverbo_b           ,                             &
        ! Optional parameters for new vegetation - Churiulin, CESR
        sur_lai_b          , sur_panday_b       , sur_biomass_b, sur_panfm_b, &
        lai_sun_b          , lai_sha_b          , sla_sun_b    , sla_sha_b  , &
        vcmax_sun_b        , vcmax_sha_b        , par_sun_b    , par_sha_b  , &
        rs_leaf_b          , psn_b              , gpp_flux_b   , npp_flux_b   )

    ! Lines 3011 - 3014
    ! from new vegetation parameterisation scheme - Churiulin, CESR
    IF (lphenology) THEN
        CALL phenology_wkarr_dealloc (izerror)
    ENDIF

    ! Lines 3204 - 3234
    ! from new vegetation parameterisation scheme - Churiulin, CESR
    ALLOCATE ( ztraleav_t           (nproma          ) ,STAT=izl ); ztraleav_t    = 0.0_vpp; ist=ist+izl
    ALLOCATE ( zverbo_t             (nproma          ) ,STAT=izl ); zverbo_t      = 0.0_vpp; ist=ist+izl
    IF (lphenology) THEN
        ! external
        ALLOCATE ( lai_t            (nproma          ), STAT=izl ); lai_t         = 0.0_vpp; ist=ist+izl
        ! radiation
        ALLOCATE ( rlat_t           (nproma          ), STAT=izl ); rlat_t        = 0.0_vpp; ist=ist+izl
        ALLOCATE ( sun_el_t         (nproma          ), STAT=izl ); sun_el_t      = 0.0_vpp; ist=ist+izl
        ALLOCATE ( sfldir_par_t     (nproma          ), STAT=izl ); sfldir_par_t  = 0.0_vpp; ist=ist+izl
        ALLOCATE ( sfldifd_par_t    (nproma          ), STAT=izl ); sfldifd_par_t = 0.0_vpp; ist=ist+izl
        ALLOCATE ( sfldifu_par_t    (nproma          ), STAT=izl ); sfldifu_par_t = 0.0_vpp; ist=ist+izl
        ALLOCATE ( cos_zen_ang_t    (nproma          ), STAT=izl ); cos_zen_ang_t = 0.0_vpp; ist=ist+izl
        ! vegetation
        ALLOCATE ( sur_lai_t        (nproma          ), STAT=izl ); sur_lai_t     = 0.0_vpp; ist=ist+izl
        ALLOCATE ( sur_panday_t     (nproma          ), STAT=izl ); sur_panday_t  = 0.0_vpp; ist=ist+izl
        ALLOCATE ( sur_biomass_t    (nproma          ), STAT=izl ); sur_biomass_t = 0.0_vpp; ist=ist+izl
        ALLOCATE ( sur_panfm_t      (nproma          ), STAT=izl ); sur_panfm_t   = 0.0_vpp; ist=ist+izl
        ALLOCATE ( lai_sun_t        (nproma          ), STAT=izl ); lai_sun_t     = 0.0_vpp; ist=ist+izl
        ALLOCATE ( lai_sha_t        (nproma          ), STAT=izl ); lai_sha_t     = 0.0_vpp; ist=ist+izl
        ALLOCATE ( sla_sun_t        (nproma          ), STAT=izl ); sla_sun_t     = 0.0_vpp; ist=ist+izl
        ALLOCATE ( sla_sha_t        (nproma          ), STAT=izl ); sla_sha_t     = 0.0_vpp; ist=ist+izl
        ALLOCATE ( vcmax_sun_t      (nproma          ), STAT=izl ); vcmax_sun_t   = 0.0_vpp; ist=ist+izl
        ALLOCATE ( vcmax_sha_t      (nproma          ), STAT=izl ); vcmax_sha_t   = 0.0_vpp; ist=ist+izl
        ALLOCATE ( par_sun_t        (nproma          ), STAT=izl ); par_sun_t     = 0.0_vpp; ist=ist+izl
        ALLOCATE ( par_sha_t        (nproma          ), STAT=izl ); par_sha_t     = 0.0_vpp; ist=ist+izl
        ALLOCATE ( rs_leaf_t        (nproma          ), STAT=izl ); rs_leaf_t     = 0.0_vpp; ist=ist+izl
        ALLOCATE ( psn_t            (nproma          ), STAT=izl ); psn_t         = 0.0_vpp; ist=ist+izl
        ALLOCATE ( gpp_flux_t       (nproma          ), STAT=izl ); gpp_flux_t    = 0.0_vpp; ist=ist+izl
        ALLOCATE ( npp_flux_t       (nproma          ), STAT=izl ); npp_flux_t    = 0.0_vpp; ist=ist+izl
    ENDIF

    ! Lines 3457 - 3486
    ! from new vegetation parameterisation scheme - Churiulin, CESR
    DEALLOCATE ( ztraleav_t         , STAT=izl ); ist=ist+izl
    DEALLOCATE ( zverbo_t           , STAT=izl ); ist=ist+izl
    IF (lphenology) THEN
        ! external
        DEALLOCATE ( lai_t            , STAT=izl ); ist=ist+izl
        ! radiation
        DEALLOCATE ( sfldir_par_t     , STAT=izl ); ist=ist+izl
        DEALLOCATE ( sfldifd_par_t    , STAT=izl ); ist=ist+izl
        DEALLOCATE ( sfldifu_par_t    , STAT=izl ); ist=ist+izl
        DEALLOCATE ( cos_zen_ang_t    , STAT=izl ); ist=ist+izl
        DEALLOCATE ( rlat_t           , STAT=izl ); ist=ist+izl
        ! vegetation
        DEALLOCATE ( sun_el_t         , STAT=izl ); ist=ist+izl
        DEALLOCATE ( sur_lai_t        , STAT=izl ); ist=ist+izl
        DEALLOCATE ( sur_panday_t     , STAT=izl ); ist=ist+izl
        DEALLOCATE ( sur_biomass_t    , STAT=izl ); ist=ist+izl
        DEALLOCATE ( sur_panfm_t      , STAT=izl ); ist=ist+izl
        DEALLOCATE ( lai_sun_t        , STAT=izl ); ist=ist+izl
        DEALLOCATE ( lai_sha_t        , STAT=izl ); ist=ist+izl
        DEALLOCATE ( sla_sun_t        , STAT=izl ); ist=ist+izl
        DEALLOCATE ( sla_sha_t        , STAT=izl ); ist=ist+izl
        DEALLOCATE ( vcmax_sun_t      , STAT=izl ); ist=ist+izl
        DEALLOCATE ( vcmax_sha_t      , STAT=izl ); ist=ist+izl
        DEALLOCATE ( par_sun_t        , STAT=izl ); ist=ist+izl
        DEALLOCATE ( par_sha_t        , STAT=izl ); ist=ist+izl
        DEALLOCATE ( rs_leaf_t        , STAT=izl ); ist=ist+izl
        DEALLOCATE ( psn_t            , STAT=izl ); ist=ist+izl
        DEALLOCATE ( gpp_flux_t       , STAT=izl ); ist=ist+izl
        DEALLOCATE ( npp_flux_t       , STAT=izl ); ist=ist+izl
    ENDIF
!==============================================================================!





10. Added data to sfc_terra_data
#===============================================================================
    ! V6.0 CLM     2022-03-15 Evgenii Churiulin
    !  Added 20 new variables for vegetation scheme
    !  Commented 2 fields (ztraleav, zverbo) -> change to global

    ! Lines 570, 691
    !zverbo      (:)              , & ! total evapotranspiration                                                                    ! Chur 2022
    !ztraleav    (:)              , & ! transpiration rate of dry leaves                                                            ! Chur 2022

    ! Lines 774 - 803
    ! from new vegetation parameterisation scheme - Churiulin, CESR
    INTEGER,          ALLOCATABLE   ::  &
        ! Arrays for get_stomatal_grid subroutine
        pft           (:)                   ! type of PFT                                          (1 - C3, 2 - C4)

    REAL(KIND = vpp), ALLOCATABLE   ::  &
        ! Additional variables for TERRA
        pzdqs         (:)               , & ! derivative of zqs with respect to T_s
        pzqvlow       (:)               , & ! specific humidity of lowest atmospheric layer        (kg/kg)
        pustar_fv     (:)               , & ! friction velocity (ustar)                            ( m/s )
        pzcatm        (:)               , & ! transfer function CA
        pzrla         (:)               , & ! atmospheric resistance                               (s/m)
        pzf_wat       (:)               , & ! soil water function for stomatal resistance          (0 - 1)
        pzepsat       (:)               , & ! saturation vapour pressure at near surface temperatures
        ! Variables for get_sun_data subroutine
        day_length    (:)               , & ! daylength                                            (hours )
        f_dyl         (:)               , & ! daylength effect on Vcmax                            (0 to 1)
        ! Variables for get_stomatal_data subroutine
        fwet          (:)               , & ! fraction of canopy that is wet                       (0 to 1)
        fdry          (:)               , & ! fraction of canopy that is dry                       (0 to 1)
        wtl           (:)               , & ! heat conductance for leaf                            ( m/s  )
        eah           (:)               , & ! vapor pressure of air in the plant canopy            ( Pa   )
        ! Variables for stomata subroutine
        rs_sun        (:)               , & ! stomatal resistance (sunlit leaves)                  (  s/m )
        rs_sha        (:)               , & ! stomatal resistance (shaded leaves)                  (  s/m )
        psn_sun       (:)               , & ! leaf photosynthesis (sunlit leaves)                  (umol co2 m^-2 s^-1)
        psn_sha       (:)               , & ! leaf photosynthesis (shaded leaves)                  (umol co2 m^-2 s^-1)
        ! Variables for respiration subroutine
        fplres        (:)               , & ! total plant respiration                              (umol CO2 m^-2 s^-1)
        cleaf         (:)                   ! CO2 partial pressure at leaf surface                 ( Pa  )

    ! Lines 880, 996
    !zverbo      (nproma)         , & ! total evapotranspiration                                                                    ! Chur 2022
    !ztraleav    (nproma)         , & ! transpiration rate of dry leaves                                                            ! Chur 2022

    ! Lines 1070 - 1097
    ! from new vegetation parameterisation scheme - Churiulin, CESR
    ALLOCATE (                                 &
            ! Variables for get_stomatal_grid subroutine
            pft           (nproma)         , & ! type of PFT
            ! Additional variables for TERRA
            pzdqs         (nproma)         , & ! derivative of zqs with respect to T_s
            pzqvlow       (nproma)         , & ! specific humidity of lowest atmospheric layer        (kg/kg)
            pustar_fv     (nproma)         , & ! friction velocity (ustar)                            ( m/s )
            pzcatm        (nproma)         , & ! transfer function CA
            pzrla         (nproma)         , & ! atmospheric resistance                               (s/m)
            pzf_wat       (nproma)         , & ! soil water function for stomatal resistance          (0 - 1)
            pzepsat       (nproma)         , & ! saturation vapour pressure at near surface temperatures
            ! Variables for get_sun_data subroutine
            day_length    (nproma)         , & ! daylength
            f_dyl         (nproma)         , & ! daylength effect on Vcmax
            ! Variables for get_stomatal_data subroutine
            fwet          (nproma)         , & ! fraction of canopy that is wet
            fdry          (nproma)         , & ! fraction of canopy that is dry
            wtl           (nproma)         , & ! heat conductance for leaf
            eah           (nproma)         , & ! vapor pressure of air in the plant canopy
            ! Variables for stomata subroutine
            rs_sun        (nproma)         , & ! stomatal resistance (sunlit leaves)
            rs_sha        (nproma)         , & ! stomatal resistance (shaded leaves)
            psn_sun       (nproma)         , & ! leaf photosynthesis (sunlit leaves)
            psn_sha       (nproma)         , & ! leaf photosynthesis (shaded leaves)
            ! Variables for respiration subroutine
            fplres        (nproma)         , & ! total plant respiration
            cleaf         (nproma)             ! CO2 partial pressure at leaf surface
        STAT=istat)

    ! Lines 1285 - 1292
    ! from new vegetation parameterisation scheme - Churiulin, CESR
    DEALLOCATE (                                                                  &
              pft            , pzdqs          , pzqvlow        , pustar_fv      , &
              pzcatm         , pzrla          , pzf_wat        , pzepsat        , &
              day_length     , f_dyl          , fwet           , fdry           , &
              wtl            , eah            , rs_sun         , rs_sha         , &
              psn_sun        , psn_sha        , fplres         , cleaf          , &
           STAT=istat)

!==============================================================================!





11. Added data to sfc_terra
#===============================================================================
    ! V6.0 CLM     2022-03-14 Evgenii Churiulin
    !  Added 25 new global parameters for vegetation:
    !       IN   : lai     , sfldir_par, sfldifd_par, sfldifu_par, cos_zen_ang, rlat, sun_el (7)
    !       INOUT: sur_lai , sur_panday, sur_biomass, sur_panfm                              (4)
    !       OUT  : ztraleav, zverbo    , lai_sun    , lai_sha
    !              sla_sun , sla_sha   , vcmax_sun  , vcmax_sha  , par_sun    , par_sha
    !              rs_leaf , psn       , gpp_flux   , npp_flux                               (14)
    !  Added 20 new local parameters for vegetation
    ! Added 6 new subroutines for vegetation
    ! Added logical parameter lphenology for the new vegetation parameterisation scheme.

    ! Lines 190
    USE data_runcontrol, ONLY :   &
        lphenology,   & ! forecast with new vegetation scheme        

    ! Lines 220 - 222
    USE sfc_phenology,         ONLY: get_stomatal_grid, get_sun_data, &             ! Chur 2022
                                     get_stomatal_data, stomata,      &             ! Chur 2022
                                     respiration, biomass_evolution                 ! Chur 2022
    !---------------------------------------------------------------------------

    ! Lines 450 - 478
    ! from new vegetation parameterisation scheme - Churiulin, CESR
        ztraleav         , & ! transpiration rate of dry leaves surface      (mm/hour)
        zverbo           , & ! total evapotranspiration                      (mm/hour)
    ! optional parameters
        lai              , & ! leaf area index                                 --                  ! Chur 2022
        sfldir_par       , & ! direct component of PAR at the ground           (W/m2 )
        sfldifd_par      , & ! diffuse downward component of PAR at the ground (W/m2 )
        sfldifu_par      , & ! diffuse upward component of PAR at the ground   (W/m2 )
        cos_zen_ang      , & ! actual cosine of zenit angle                    (rad  )
        rlat             , & ! latitide                                        (rad  )
        sun_el           , & ! sun elevation                                   (deg  )
        !
        sur_lai          , & ! leaf area index depending on biomass evolution  (m2/m2)
        sur_panday       , & ! daily values of photosyntesis                   (kg CO2 m-2 s-1)
        sur_biomass      , & ! leaf biomass                                    (kg of dry matter m^2)
        sur_panfm        , & ! maximum leaf assimilation                       (kg CO2 m-2 s-1)
        !
        lai_sun          , & ! sunlit leaf area index                          (m2/m2)
        lai_sha          , & ! shaded leaf area index                          (m2/m2)
        sla_sun          , & ! specific leaf area index for sunlit leaves      (m2/m2)
        sla_sha          , & ! specific leaf area index for shaded leaves      (m2/m2)
        vcmax_sun        , & ! maximum rate of carboxylation for sunlit leaves (umol CO2/m2/s)
        vcmax_sha        , & ! maximum rate of carboxylation for shaded leaves (umol CO2/m2/s)
        par_sun          , & ! photosynthetic active radiation for sunlit leaf (W/m2)
        par_sha          , & ! photosynthetic active radiation for shaded leaf (W/m2)
        rs_leaf          , & ! canopy leaf resistance                          (s/m)
        psn              , & ! leaf photosynthesis                             (umol CO2/m2/s)
        gpp_flux         , & ! gross primary production                        (gC/m2/s)
        npp_flux         , & ! net primary production                          (gC/m2/s)


    ! Lines 613 - 614
    REAL    (KIND = vpp), DIMENSION(nvec), INTENT(OUT) :: &
        ! from new vegetation parameterisation scheme - Churiulin, CESR
        ztraleav         , & ! transpiration rate of dry leaves surface      (mm/hour)
        zverbo               ! total evapotranspiration                      (mm/hour)

    ! Lines 650 - 678
    ! from new vegetation parameterisation scheme - Churiulin, CESR
    REAL    (KIND = vpp), DIMENSION(nvec), OPTIONAL, INTENT(IN)  :: &
        lai              , & ! leaf area index                                 --
        sfldir_par       , & ! direct component of PAR at the ground           (W/m2)
        sfldifd_par      , & ! diffuse downward component of PAR at the ground (W/m2)
        sfldifu_par      , & ! diffuse upward component of PAR at the ground   (W/m2)
        cos_zen_ang      , & ! actual cosine of zenit angle                    (rad )
        rlat             , & ! latitide                                        (rad )
        sun_el               ! sun elevation                                   (deg )

    REAL    (KIND = vpp), DIMENSION(nvec), OPTIONAL, INTENT(INOUT) :: &
        sur_lai          , & ! leaf area index depending on biomass evolution  (m2/m2)
        sur_panday       , & ! daily values of photosyntesis                   (kg CO2 m-2 s-1)
        sur_biomass      , & ! leaf biomass                                    (kg of dry matter m^2)
        sur_panfm            ! maximum leaf assimilation                       (kg CO2 m-2 s-1)

    REAL    (KIND = vpp), DIMENSION(nvec), OPTIONAL, INTENT(OUT) :: &
        lai_sun          , & ! sunlit leaf area index                          (m2/m2)
        lai_sha          , & ! shaded leaf area index                          (m2/m2)
        sla_sun          , & ! specific leaf area index for sunlit leaves      (m2/m2)
        sla_sha          , & ! specific leaf area index for shaded leaves      (m2/m2)
        vcmax_sun        , & ! maximum rate of carboxylation for sunlit leaves (umol co2/m**2/s)
        vcmax_sha        , & ! maximum rate of carboxylation for shaded leaves (umol co2/m**2/s)
        par_sun          , & ! photosynthetic active radiation for sunlit leaf (W/m2)
        par_sha          , & ! photosynthetic active radiation for shaded leaf (W/m2)
        rs_leaf          , & ! canopy leaf resistance                          (s/m)
        psn              , & ! leaf photosynthesis                             (umol CO2/m2/s)
        gpp_flux         , & ! gross primary production                        (gC/m2/s)
        npp_flux             ! net primary production                          (gC/m2/s)

    ! Lines 1031, 1153
    !zverbo      (nvec)             , & ! total evapotranspiration                                                                  ! Chur 2022
    !ztraleav    (nvec)             , & ! transpiration rate of dry leaves                                                          ! Chur 2022

    ! Lines 1220 - 1229
    ! from new vegetation parameterisation scheme - Churiulin, CESR
    INTEGER                        ::  &
        ! Arrays for get_stomatal_grid subroutine
        pft         (nvec)                 ! type of PFT                                          (1 - C3, 2 - C4)

    REAL    (KIND=vpp) ::              &
        ! Additional variables for TERRA
        pzdqs       (nvec)             , & ! derivative of zqs with respect to T_s
        pzqvlow     (nvec)             , & ! specific humidity of lowest atmospheric layer        (kg/kg)
        pustar_fv   (nvec)             , & ! friction velocity (ustar)                            ( m/s )
        pzcatm      (nvec)             , & ! transfer function CA
        pzrla       (nvec)             , & ! atmospheric resistance                               (s/m)
        pzf_wat     (nvec)             , & ! soil water function for stomatal resistance          (0 - 1)
        pzepsat     (nvec)             , & ! saturation vapour pressure at near surface temperatures
        ! Variables for get_sun_data subroutine
        day_length  (nvec)             , & ! daylength                                            (hours )
        f_dyl       (nvec)             , & ! daylength effect on Vcmax                            (0 to 1)
        ! Variables for get_stomatal_data subroutine
        fwet        (nvec)             , & ! fraction of canopy that is wet                       (0 to 1)
        fdry        (nvec)             , & ! fraction of canopy that is dry                       (0 to 1)
        wtl         (nvec)             , & ! heat conductance for leaf                            ( m/s  )
        eah         (nvec)             , & ! vapor pressure of air in the plant canopy            ( Pa   )
        ! Variables for stomata subroutine
        rs_sun      (nvec)             , & ! stomatal resistance (sunlit leaves)                  (  s/m )
        rs_sha      (nvec)             , & ! stomatal resistance (shaded leaves)                  (  s/m )
        psn_sun     (nvec)             , & ! leaf photosynthesis (sunlit leaves)                  (umol co2 m^-2 s^-1)
        psn_sha     (nvec)             , & ! leaf photosynthesis (shaded leaves)                  (umol co2 m^-2 s^-1)
        ! Variables for respiration subroutine
        fplres      (nvec)             , & ! total plant respiration                              (umol CO2 m^-2 s^-1)
        cleaf       (nvec)                 ! CO2 partial pressure at leaf surface                 ( Pa  )

    ! Lines 1434 - 1447
    ! from new vegetation parameterisation scheme - Churiulin, CESR
    IF (lphenology) THEN
        WRITE(*,'(A,F28.16)') '   sfldir_par       :  ', sfldir_par  (i)
        WRITE(*,'(A,F28.16)') '   sfldifd_par      :  ', sfldifd_par (i)
        WRITE(*,'(A,F28.16)') '   sfldifu_par      :  ', sfldifu_par (i)
        WRITE(*,'(A,F28.16)') '   cos_zen_ang      :  ', cos_zen_ang (i)
        WRITE(*,'(A,F28.16)') '   rlat             :  ', rlat        (i)
        WRITE(*,'(A,F28.16)') '   sun_el           :  ', sun_el      (i)
        WRITE(*,'(A,F28.16)') '   sur_lai (in)     :  ', sur_lai     (i)
        WRITE(*,'(A,F28.16)') '   sur_panday (in)  :  ', sur_panday  (i)
        WRITE(*,'(A,F28.16)') '   sur_biomass (in) :  ', sur_biomass (i)
        WRITE(*,'(A,F28.16)') '   sur_panfm (in)   :  ', sur_panfm   (i)
    ENDIF

    ! Lines 1813 - 1829
    ! from new vegetation parameterisation scheme - Churiulin, CESR
    ztraleav     (i)        = 0.0_vpp
    zverbo       (i)        = 0.0_vpp
    IF (lphenology) THEN
        lai_sun  (i)        = 0.0_vpp
        lai_sha  (i)        = 0.0_vpp
        sla_sun  (i)        = 0.0_vpp
        sla_sha  (i)        = 0.0_vpp
        vcmax_sun(i)        = 0.0_vpp
        vcmax_sha(i)        = 0.0_vpp
        par_sun  (i)        = 0.0_vpp
        par_sha  (i)        = 0.0_vpp
        rs_leaf  (i)        = 0.0_vpp
        psn      (i)        = 0.0_vpp
        gpp_flux (i)        = 0.0_vpp
        npp_flux (i)        = 0.0_vpp
    ENDIF

    ! Lines 2519 - 2523
    ! from new vegetation parameterisation scheme - Churiulin, CESR
    IF (lphenology) THEN
      pzqvlow(i) = zqvlow
      pzdqs  (i) = zdqs
    ENDIF

    ! Lines 3016 - 3028
    ! from new vegetation parameterisation scheme - Churiulin, CESR
    IF (lphenology) THEN
        ! Saturation deficit function
        z2iw       = ztsk_pm(i)*b2w + (1.0_vpp - ztsk_pm(i))*b2i
        z4iw       = ztsk_pm(i)*b4w + (1.0_vpp - ztsk_pm(i))*b4i

        ! Chur 2022 - add arrays for variables: zcatm, zustar, zrla, zf_wat, zepsat
        pustar_fv(i) = zuv * SQRT(tcm(i))
        pzcatm(i)    = zcatm
        pzrla(i)     = zrla                                                     ! Doesn't work with  icant = 2
        pzf_wat(i)   = zf_wat
        pzepsat(i)   = zsf_psat_iw(t(i), z2iw, z4iw)
    ENDIF


    ! Lines 3030 - 3062
    ! Define the current algorithm for calculations of stomatal resistance - Churiulin, CESR
    ! In case of standart run the parameter lphenology is .FALSE.
        IF (.NOT. lphenology) THEN
          IF (lstomata) THEN
            IF (itype_trvg == 3) THEN
              ! Modification of rsmin depending on accumulated plant evaporation; the z0 dependency
              ! is used to get a stronger effect for trees than for low vegetation
#ifdef __COSMO__
              IF (z0(i) <= 0.4_vpp) zxx = MIN(1.25_vpp, zxx)
              zzz = MAX(0.5_vpp, EXP(SQRT(z0(i))*LOG(zxx)) )
              ! limit reduction of rsmin-factor below 1 at low temperatures
              zzz = MAX(zzz, MIN(1.0_vpp,(t0_melt+15.0_vpp-t(i))/15.0_vpp))
#endif
#ifdef __ICON__
              zzz = MAX(0.5_vpp+MIN(0.5_vpp,1.0_vpp-laifac(i)), EXP(SQRT(z0(i))*LOG(zxx)) )
#endif
            ELSE
              zzz = 1.0_vpp
            ENDIF
            zedrstom   = 1.0_vpp/crsmax + (1.0_vpp/MAX(40.0_vpp,zzz*rsmin2d(i)) - &
                         1.0_vpp/crsmax) * zf_rad(i) * zf_wat * zf_tem * zf_sat
          ELSE
            zedrstom   = 1.0_vpp/crsmax + (1.0_vpp/crsmin -                     &
                         1.0_vpp/crsmax) * zf_rad(i) * zf_wat * zf_tem * zf_sat
          ENDIF
          zrstom   = 1.0_vpp/zedrstom              ! stomatal resistance
          rstom(i) = zrstom
          zrveg    = zrla + zrstom

          ! Transpiration rate of dry leaves:
          ztraleav(i)=zep_s(i)*tai(i)/(sai(i)+zrveg*zcatm)
        END IF  ! upwards directed potential evaporation only


    ! Lines 3066 - 3149
    ! New algorithm for calculations of stomatal resistance - Churiulin, CESR
    ! In case of standart run the parameter lphenology is .TRUE.
    IF (lphenology) THEN

      ! Run get_stomatal grid subrotine
      !-------------------------------------
      ! Need to define the type of vegetation at PFT grid
      CALL get_stomatal_grid(nvec, ivstart, ivend, soiltyp_subs, zep_s, pft)

      ! Run get_sun_data subroutine
      !-------------------------------------
      ! Need to define additional solar parameters and daylenth limitation function
      CALL get_sun_data(nvec, ivstart, ivend, dt, soiltyp_subs, zep_s,     &
                        rlat, sun_el , day_length, f_dyl                   )

      ! Run get_stomatal_data subroutine
      !-------------------------------------
      ! Need to defeni parameteres from two-big leaf approach
      CALL get_stomatal_data(nvec       , ivstart   , ivend      , soiltyp_subs, &
                             pft        , zep_s     , zdz_snow   , pustar_fv   , &
                             ps         , qv_s      , pzqvlow    , lai         , &
                             cos_zen_ang, sfldir_par, sfldifd_par, sfldifu_par , &
                             fwet       , fdry      , wtl        , rs_leaf     , &
                             par_sun    , par_sha   , eah        , lai_sun     , &
                             lai_sha    , sla_sun   , sla_sha                    )

      ! Run stomata subroutine - sunlit and shaded leaves
      !-------------------------------------
      ! Calculations of stomatal resistance for sunlit leaves
      CALL stomata(nvec     , ivstart, ivend, soiltyp_subs, pft    , zep_s  , &
                   ps       , zth_low, ztsk , pzepsat     , eah    , sla_sun, &
                   par_sun  , pzf_wat, f_dyl, rs_leaf     , rs_sun , psn_sun, &
                   vcmax_sun                                                  )
      ! Calculations of stomatal resistance for shaded leaves
      CALL stomata(nvec     , ivstart, ivend, soiltyp_subs, pft    , zep_s  , &
                   ps       , zth_low, ztsk , pzepsat     , eah    , sla_sha, &
                   par_sha  , pzf_wat, f_dyl, rs_leaf     , rs_sha , psn_sha, &
                   vcmax_sha                                                  )

      ! Start final calculations
      DO i = ivstart, ivend
#ifndef MESSY
        IF (soiltyp_subs(i) >= 3 .AND. zep_s(i) < 0.0_vpp) THEN  ! upwards directed potential evaporation
#else
        ! calculate this part always for MESSy
        IF (soiltyp_subs(i) >= 3) THEN ! < 3 is ice and rock
#endif

          ! Final calculations of stomatal resistance
          IF ((lai(i) <= 1.0_vpp).AND.(cos_zen_ang(i) <= 0.01_vpp)) THEN
            zedrstom = 1.0_vpp / rs_sha(i)
          ELSE
            zedrstom = (1.0_vpp / rs_sun(i)) * lai_sun(i) + &
                       (1.0_vpp / rs_sha(i)) * lai_sha(i)
          ENDIF

          ! stomatal resistance
          zrstom   = 1.0_vpp / zedrstom
          rstom(i) = zrstom
          zrveg    = pzrla(i) + zrstom

          ! Transpiration rate of dry leaves:
          ztraleav(i) = zep_s  (i) * tai(i) / (sai(i) + zrveg * pzcatm(i))
          psn     (i) = psn_sun(i) * lai_sun(i) + &
                        psn_sha(i) * lai_sha(i)
        ENDIF  ! upwards directed potential evaporation only
      ENDDO

      ! Run respiration subroutine
      !-------------------------------------
      ! Calculations of plant respiration
      CALL respiration(nvec     , ivstart  , ivend  , soiltyp_subs, zep_s  , &
                       ps       , zth_low  , rs_leaf, par_sun     , par_sha, &
                       vcmax_sun, vcmax_sha, lai_sun, lai_sha     , lai    , &
                       psn_sun  , psn_sha  , psn    , fplres      , cleaf  , &
                       gpp_flux , npp_flux                                   )

      ! Run biomass_evolution subroutine - Chur 2022
      !-------------------------------------
      ! Calculations of LAI based on biomass evolution
      CALL biomass_evolution(nvec   , ivstart, ivend, dt, pft, soiltyp_subs, &
                             zep_s  , plcov  , psn  , sur_lai, sur_panday  , &
                             sur_biomass     , sur_panfm                     )
    ENDIF ! end of new vegetation calculations

    ! Lines 7262 - 7284
    ! from new vegetation parameterisation scheme - Churiulin, CESR
    WRITE(*,'(A, F28.16)') '   ztraleav             :  ', ztraleav(i)
    WRITE(*,'(A, F28.16)') '   zverbo               :  ', zverbo(i)
    IF (lphenology) THEN
        ! intent INOUT
        WRITE(*,'(A, F28.16)') '   sur_lai (out)      :  ', sur_lai(i)
        WRITE(*,'(A, F28.16)') '   sur_panday (out)   :  ', sur_panday(i)
        WRITE(*,'(A, F28.16)') '   sur_biomass (out)  :  ', sur_biomass(i)
        WRITE(*,'(A, F28.16)') '   sur_panfm (out)    :  ', sur_panfm(i)
        ! intent OUT
        WRITE(*,'(A, F28.16)') '   lai_sun            :  ', lai_sun(i)
        WRITE(*,'(A, F28.16)') '   lai_sha            :  ', lai_sha(i)
        WRITE(*,'(A, F28.16)') '   sla_sun            :  ', sla_sun(i)
        WRITE(*,'(A, F28.16)') '   sla_sha            :  ', sla_sha(i)
        WRITE(*,'(A, F28.16)') '   vcmax_sun          :  ', vcmax_sun(i)
        WRITE(*,'(A, F28.16)') '   vcmax_sha          :  ', vcmax_sha(i)
        WRITE(*,'(A, F28.16)') '   par_sun            :  ', par_sun(i)
        WRITE(*,'(A, F28.16)') '   par_sha            :  ', par_sha(i)
        WRITE(*,'(A, F28.16)') '   rs_leaf            :  ', rs_leaf(i)
        WRITE(*,'(A, F28.16)') '   psn                :  ', psn(i)
        WRITE(*,'(A, F28.16)') '   gpp_flux           :  ', gpp_flux(i)
        WRITE(*,'(A, F28.16)') '   npp_flux           :  ', npp_flux(i)
    ENDIF

!==============================================================================!




12. Added new module sfc_phenology
#===============================================================================
!-------------------------------------------------------------------------------




13. Added new module sfc_phenology_data
#===============================================================================
!-------------------------------------------------------------------------------




14. Added data to dfi_initialization
#===============================================================================
    ! V6_0 CLM     2021-03-16 Evgenii Churiulin
    ! Added vegetation fields (aztraleav, azverbo)
    ! Added logical parameter lphenology for the new vegetation parameterisation scheme.

    ! Lines 296 - 300
    USE data_fields   , ONLY :   &
    ! 14. fields for the vegetation parameterization - Churiulin, CESR
    ! -----------------------------------------
        aztraleav     , & ! average transpiration rate of dry leaves      (mm/day)
        azverbo           ! average total evapotranspiration              (mm/day)

    ! Lines 409
    USE data_runcontrol , ONLY :   &
        lphenology,   & ! forecast with new vegetation scheme 

    ! Lines 2334 - 2335
    aztraleav= 0.0_wp                                                                 ! Chur 2022
    azverbo  = 0.0_wp                                                                 ! Chur 2022
!==============================================================================!





15. Added data to near_surface
#===============================================================================
    ! V6_0 CLM     2022-03-15 Evgenii Churiulin
    !  Added new accumulated variables for solar radiation:
    !                             (asfldir_par, asfldifd_par, asfldifu_par)
    !  Added tile 0 to ztraleav, zverbo

    ! Lines 149 - 150
    USE data_fields,        ONLY:                                        &
       ztraleav, zverbo, aztraleav, azverbo, sfldir_par, sfldifd_par,    &          ! Chur 2022
       sfldifu_par, asfldir_par, asfldifd_par, asfldifu_par,             &          ! Chur 2022

    ! Lines 169
    USE data_runcontrol,    ONLY:
        ..., lphenology

    ! Lines 1441 - 1443
    ! from new vegetation parameterisation scheme - Churiulin, CESR
    aztraleav(i,j)= aztraleav(i,j)+ ztraleav(i,j,0)
    azverbo(i,j)  = azverbo(i,j)  +   zverbo(i,j,0)

    ! from new vegetation parameterisation scheme - Churiulin, CESR
    IF (lphenology) THEN
        asfldir_par(i,j)  = asfldir_par(i,j)  + sfldir_par(i,j)
        asfldifd_par(i,j) = asfldifd_par(i,j) + sfldifd_par(i,j)
        asfldifu_par(i,j) = asfldifu_par(i,j) + sfldifu_par(i,j)
    ENDIF
!==============================================================================!




16. Added data to organize_data
#===============================================================================
    ! V6_0 CLM     2021-05-12 Churiulin Evgenii
    !  Added new vegetation parameters for restart files:
    !          1. vegetation fields               : sur_lai, sur_panday, sur_biomass, sur_panfm
    !          2. vegetation (average fields)     : aztraleav, azverbo,
    !          3. radiation  (accumulation fields): asfldir_par, asfldifd_par, asfldifu_par)
    !  Added logical parameter lphenology for the new vegetation parameterisation scheme.

    ! Lines 709
    USE data_runcontrol,    ONLY:                               &
        lphenology,        & ! forecast with new vegetation scheme                                                                      ! Chur 2022

    ! Lines 1870 - 1877
    ! from new vegetation parameterisation scheme - Churiulin, CESR
    IF (lphenology) THEN
        yvarini(nyvar_i + 1) = 'SUR_LAI    '
        yvarini(nyvar_i + 2) = 'SUR_BIOMASS'
        yvarini(nyvar_i + 3) = 'SUR_PANMAX '
        yvarini(nyvar_i + 4) = 'SUR_PANDAY '
        nyvar_i = nyvar_i + 4
    ENDIF

    ! Lines 2078
    'SUR_LAI   ','SUR_BIOMASS','SUR_PANFM ','SUR_PANDAY'   )                                                        ! Chur 2022

    ! Lines 2510 - 2523
    ! from new vegetation parameterisation scheme - Churiulin, CESR
    IF (lphenology) THEN
        ! from TERRA
        pp_restart%yvarml(pp_restart%nyvar_m + 1)  = 'SUR_LAI    '
        pp_restart%yvarml(pp_restart%nyvar_m + 2)  = 'SUR_BIOMASS'
        pp_restart%yvarml(pp_restart%nyvar_m + 3)  = 'SUR_PANMAX '
        pp_restart%yvarml(pp_restart%nyvar_m + 4)  = 'SUR_PANDAY '
        ! from radiation
        pp_restart%yvarml(pp_restart%nyvar_m + 5)  = 'ASFLDIR_PAR  '
        pp_restart%yvarml(pp_restart%nyvar_m + 6)  = 'ASFLDIFD_PAR '
        pp_restart%yvarml(pp_restart%nyvar_m + 7)  = 'ASFLDIFU_PAR '

        pp_restart%nyvar_m = pp_restart%nyvar_m + 7
    ENDIF

    ! Lines 3023
    ('SUR_LAI   ','SUR_BIOMASS','SUR_PANFM  ','SUR_PANDAY'   )

    ! Lines 3210
    ('SUR_LAI   ','SUR_BIOMASS','SUR_PANFM  ','SUR_PANDAY'   )

    ! Lines 8241 - 8242
    CASE(...
         'ASFLDIR_PAR'  , 'ASFLDIFD_PAR' , 'ASFLDIFU_PAR ', 'AZTRALEAV'    , &
         'AZVERBO'                                                         , &
         ...
        )
!==============================================================================!





17. Added data to src_input
#===============================================================================22222222
    ! V6_0 CLM     2021-05-12 Churiulin Evgenii
    !  Added new vegetation parameters for restart files:
    !          1. vegetation fields: sur_lai, sur_panday, sur_biomass, sur_panfm
    !  Added logical parameter lphenology for the new vegetation parameterisation scheme.

    ! Lines 560
    USE data_runcontrol, ONLY : &
        lphenology,   & ! forecast with new vegetation scheme                                                                           ! Chur 2022


    ! Lines 5039
    ('SUR_LAI   ','SUR_BIOMASS','SUR_PANFM  ','SUR_PANDAY' )
!==============================================================================!





18. Added data to src_gridpoints
#===============================================================================
    ! V6_0         2022-04-05 Evgenii Churiulin
    !  Added additional parameters for radiation: sfldir_par, sfldifd_par, sfldifu_par

    ! Lines 351 - 353
    USE data_fields     , ONLY :   &
        sfldir_par ,    & ! direct comp. of PAR at the ground                (W/m2)              ! Chur 2022
        sfldifd_par,    & ! diffuse downward comp. of PAR at the ground      (W/m2)              ! Chur 2022
        sfldifu_par,    & ! diffuse upward comp. of PAR at the ground        (W/m2)              ! Chur 2022

    ! Lines 397
    USE data_runcontrol , ONLY :   &
        lphenology,   & ! forecast with new vegetation scheme                                    ! Chur 2022

    ! Lines 558 - 560
    ! from new vegetation parameterisation scheme 3 parameters - Churiulin, CESR
    gpsfldir_par(:,:), & ! direct comp. of PAR at the ground            ( W/m2)
    gpsfldifd_par(:,:),& ! diffuse downward comp. of PAR at the ground  ( W/m2)
    gpsfldifu_par(:,:),& ! diffuse upward comp. of PAR at the ground    ( W/m2)

    ! Lines 1482 - 1487
    ! from new vegetation parameterisation scheme - Churiulin, CESR
    IF (lphenology) THEN
        gpsfldir_par (mzg,nzpa) = sfldir_par (izgp,jzgp)
        gpsfldifd_par(mzg,nzpa) = sfldifd_par(izgp,jzgp)
        gpsfldifu_par(mzg,nzpa) = sfldifu_par(izgp,jzgp)
    ENDIF

    ! Lines 1856 - 1867
    ! from new vegetation parameterisation scheme - Churiulin, CESR
    IF (lphenology) THEN
    WRITE (yline(k),'(A,F10.3)')                                             &
        '              (w/m**2)   SFLDIR_PAR : ', gpsfldir_par(mzg,nzpa)
    k = k+1
    WRITE (yline(k),'(A,F10.3)')                                             &
        '              (w/m**2)   SFLDIFD_PAR: ', gpsfldifd_par(mzg,nzpa)
    k = k+1
    WRITE (yline(k),'(A,F10.3)')                                             &
        '              (w/m**2)   SFLDIFU_PAR: ', gpsfldifu_par(mzg,nzpa)
    k = k+1
    ENDIF

    ! Lines 2496 - 2501
    ! from new vegetation parameterisation scheme - Churiulin, CESR
    IF (lphenology) THEN
      ALLOCATE ( gpsfldir_par  (ngp, nstepsgp), STAT=istat );  gpsfldir_par  = 0.0_wp
      ALLOCATE ( gpsfldifd_par (ngp, nstepsgp), STAT=istat );  gpsfldifd_par = 0.0_wp
      ALLOCATE ( gpsfldifu_par (ngp, nstepsgp), STAT=istat );  gpsfldifu_par = 0.0_wp
    ENDIF

    ! Lines 2677 - 2682
    ! from new vegetation parameterisation scheme - Churiulin, CESR
    IF (lphenology) THEN
      DEALLOCATE ( gpsfldir_par , STAT=istat)
      DEALLOCATE ( gpsfldifd_par, STAT=istat)
      DEALLOCATE ( gpsfldifu_par, STAT=istat)
    ENDIF
!==============================================================================!





19. Added data to radiation_interface
#===============================================================================
    ! V6_0         2022-04-05 Evgenii Churiulin
    !  - Added new variables for photosynthetic active radiation for sfc_phenology
    !    module: sfldir_par, sfldifd_par, sfldifu_par, cos_zen_ang
    !  - Added tile variables:  sfldir_par_b, sfldifd_par_b,
    !                          sfldifu_par_b, cos_zen_ang_b
    !  - Added logical parameter lphenology for the new vegetation parameterisation scheme.

    ! Lines 249 - 256
    ! new variables for PAR and cosine zenit angle for src_phenology - Churiulin, CESR
    sfldir_par   ,  & ! direct comp. of PAR flux at the ground           ( W/m2)
    sfldifd_par  ,  & ! diffuse downward comp. of PAR flux at the ground ( W/m2)
    sfldifu_par  ,  & ! diffuse upward comp. of PAR flux at the ground   ( W/m2)
    sfltrdir_par ,  & ! direct comp. of PAR at surface                   ( W/m2)
    sfltrdifd_par,  & ! diffuse downward comp. of PAR at surface         ( W/m2)
    sfltrdifu_par,  & ! diffuse upward comp. of PAR at surface           ( W/m2)
    cos_zen_ang  ,  & ! cosine of solar zenith angle for src_phenology   ( rad )

    ! Lines 296
    USE data_runcontrol, ONLY : &
        lphenology,   & ! forecast with new vegetation scheme                                                                       ! Chur 2022

    ! Lines 523
    USE data_block_fields,  ONLY :
        sfldir_par_b  , sfldifd_par_b  , sfldifu_par_b, cos_zen_ang_b                                                               ! Chur 2022

    ! Lines 1470 - 1476
    ! from new vegetation parameterisation scheme - Churiulin, CESR
    IF (lphenology) THEN
        CALL register_copy (sfldir_par_b , radCopyList, copyFromBlockF)! print *, 'sfldir_par      ok'
        CALL register_copy (sfldifd_par_b, radCopyList, copyFromBlockF)! print *, 'sfldifd_par     ok'
        CALL register_copy (sfldifu_par_b, radCopyList, copyFromBlockF)! print *, 'sfldifu_par_b   ok'
        CALL register_copy (cos_zen_ang_b, radCopyList, copyFromBlockF)! print *, 'cos_zen_ang_b   ok'
    ENDIF

    ! Lines 1494 - 1500
    ! from new vegetation parameterisation scheme - Churiulin, CESR
    IF (lphenology) THEN
        CALL register_copy (sfldir_par_b , radsoCopyList, copyFromBlockF)! print *, 'sfldir_par   ok'
        CALL register_copy (sfldifd_par_b, radsoCopyList, copyFromBlockF)! print *, 'sfldifd_par  ok'
        CALL register_copy (sfldifu_par_b, radsoCopyList, copyFromBlockF)! print *, 'sfldifu_par  ok'
        CALL register_copy (cos_zen_ang_b, radsoCopyList, copyFromBlockF)! print *, 'cos_zen_ang  ok'
    ENDIF

    ! Lines 2116 - 2121
    ! from new vegetation parameterisation scheme - Churiulin, CESR
    IF (lphenology) THEN
        sfltrdir_par (i,j) = 0.0_wp
        sfltrdifd_par(i,j) = 0.0_wp
        sfltrdifu_par(i,j) = 0.0_wp
    ENDIF

    ! Lines 2133 - 2138
    ! from new vegetation parameterisation scheme - Churiulin, CESR
    IF (lphenology) THEN
      sfltrdir_par(i,j)   = zflpar_s_dir (ip) / sod_t(i,j)
      sfltrdifd_par(i,j)  = zflpar_s_difd(ip) / sod_t(i,j)
      sfltrdifu_par(i,j)  = zflpar_s_difu(ip) / sod_t(i,j)
    ENDIF

    ! Lines 2250 - 2260
    ! from new vegetation parameterisation scheme - Churiulin, CESR
    IF (lphenology) THEN
      sfldir_par   (i,j) = 0.0_wp
      sfldifd_par  (i,j) = 0.0_wp
      sfldifu_par  (i,j) = 0.0_wp
      cos_zen_ang  (i,j) = 0.0_wp
      sfldir_par_b (ip)  = 0.0_wp
      sfldifd_par_b(ip)  = 0.0_wp
      sfldifu_par_b(ip)  = 0.0_wp
      cos_zen_ang_b(ip)  = 0.0_wp
    ENDIF

    ! Lines 2289 - 2299
    ! from new vegetation parameterisation scheme - Churiulin, CESR
    IF (lphenology) THEN
        sfldir_par  (i,j) = sfltrdir_par (i,j) * sod_t(i,j)
        sfldifd_par (i,j) = sfltrdifd_par(i,j) * sod_t(i,j)
        sfldifu_par (i,j) = sfltrdifu_par(i,j) * sod_t(i,j)
        cos_zen_ang (i,j) = zsmu0_fg(ip)
        sfldir_par_b (ip) = sfltrdir_par (i,j) * sod_t(i,j)
        sfldifd_par_b(ip) = sfltrdifd_par(i,j) * sod_t(i,j)
        sfldifu_par_b(ip) = sfltrdifu_par(i,j) * sod_t(i,j)
        cos_zen_ang_b(ip) = zsmu0_fg(ip)
    ENDIF
!==============================================================================!