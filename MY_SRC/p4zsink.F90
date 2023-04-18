MODULE p4zsink
   !!======================================================================
   !!                         ***  MODULE p4zsink  ***
   !! TOP :  PISCES  vertical flux of particulate matter due to 
   !!        gravitational sinking
   !!        This module is the same for both PISCES and PISCES-QUOTA
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Change aggregation formula
   !!             3.5  !  2012-07  (O. Aumont) Introduce potential time-splitting
   !!             4.0  !  2019     (O. Aumont) an external subroutine is called
   !!                                          to compute the impact of sinking
   !!----------------------------------------------------------------------
   !!   p4z_sink       :  Compute vertical flux of particulate matter due to gravitational sinking
   !!   p4z_sink_init  :  Unitialisation of sinking speed parameters
   !!   p4z_sink_alloc :  Allocate sinking speed variables
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE trcsink         !  General routine to compute sedimentation
   USE prtctl_trc      !  print control for debugging
   USE iom             !  I/O manager
   USE lib_mpp

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_sink         ! called in p4zbio.F90
   PUBLIC   p4z_sink_init    ! called in trcini_pisces.F90
   PUBLIC   p4z_sink_alloc   ! called in trcini_pisces.F90

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sinking, sinking2  !: POC sinking fluxes 
   !                                                          !  (different meanings depending on the parameterization)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sinkingn, sinking2n  !: PON sinking fluxes 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sinkingp, sinking2p  !: POP sinking fluxes 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sinkcal, sinksil   !: CaCO3 and BSi sinking fluxes
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sinkfer            !: Small BFe sinking fluxes
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sinkfer2           !: Big iron sinking fluxes
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sinklfe
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sinklfa
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sinkafs
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   sinkafb
   INTEGER  :: ik100

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zsink.F90 13233 2020-07-02 18:34:16Z aumont $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_sink ( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_sink  ***
      !!
      !! ** Purpose :   Compute vertical flux of particulate matter due to 
      !!                gravitational sinking. 
      !!
      !! ** Method  : - An external advection subroutine is called to compute
      !!                the impact of sinking on the particles. The tracers
      !!                concentrations are updated in this subroutine which
      !!                is mandatory to deal with negative concentrations
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) :: kt, knt
      INTEGER  ::   ji, jj, jk
      CHARACTER (len=25) :: charout
      REAL(wp) :: zmax, zfact
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zw3d
      REAL(wp), ALLOCATABLE, DIMENSION(:,:  ) :: zw2d
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('p4z_sink')

      ! Initialization of some global variables
      ! ---------------------------------------
      prodpoc(:,:,:) = 0.
      conspoc(:,:,:) = 0.
      prodgoc(:,:,:) = 0.
      consgoc(:,:,:) = 0.

      ! Sinking speeds of big detritus is increased with depth as shown
      ! by data and from the coagulation theory. This is controled by
      ! wsbio2max and wsbio2scale. If wsbio2max is set to wsbio2, then
      ! sinking speed is constant with depth.
      ! CaCO3 and bSi are supposed to sink at the big particles speed 
      ! due to their high density
      ! ---------------------------------------------------------------
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1,jpi
               zmax  = MAX( heup_01(ji,jj), hmld(ji,jj) )
               zfact = MAX( 0., gdepw_n(ji,jj,jk+1) - zmax ) / wsbio2scale
               wsbio4(ji,jj,jk) = wsbio2 + MAX(0., ( wsbio2max - wsbio2 )) * zfact
            END DO
         END DO
      END DO

      ! Sinking speed of the small particles is always constant
      wsbio3(:,:,:) = wsbio
       
      ! Initialize to zero all the sinking arrays 
      ! -----------------------------------------
      sinking (:,:,:) = 0.e0
      sinking2(:,:,:) = 0.e0
      sinkcal (:,:,:) = 0.e0
      sinkfer (:,:,:) = 0.e0
      sinksil (:,:,:) = 0.e0
      sinkfer2(:,:,:) = 0.e0

      ! Compute the sedimentation term using trc_sink for all the sinking particles
      ! ---------------------------------------------------------------------------
      CALL trc_sink( kt, wsbio3, sinking , jppoc, rfact2 )
      CALL trc_sink( kt, wsbio3, sinkfer , jpsfe, rfact2 )
      CALL trc_sink( kt, wsbio4, sinking2, jpgoc, rfact2 )
      CALL trc_sink( kt, wsbio4, sinkfer2, jpbfe, rfact2 )
      CALL trc_sink( kt, wsbio4, sinksil , jpgsi, rfact2 )
      CALL trc_sink( kt, wsbio4, sinkcal , jpcal, rfact2 )

      ! sinking of lithogenic Fe
      IF( ln_bait ) THEN

      wsbio5(:,:,:) = wslfe

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1,jpi
               zmax  = MAX( heup_01(ji,jj), hmld(ji,jj) )
               zfact = MAX( 0., gdepw_n(ji,jj,jk+1) - zmax ) / wsbio2scale
               wsbio6(ji,jj,jk) = wslfea + MAX(0., ( wslfeamax - wslfea )) * zfact
            END DO
         END DO
      END DO
      sinklfe (:,:,:) = 0.e0
      sinklfa (:,:,:) = 0.e0
      IF ( ln_feauth ) THEN
      sinkafs (:,:,:) = 0.e0
      sinkafb (:,:,:) = 0.e0
      ENDIF ! IF ( ln_feauth ) THEN

      ! Compute the sedimentation term using trc_sink for all the sinking
      ! particles
      ! ---------------------------------------------------------------------------
      CALL trc_sink( kt, wsbio5, sinklfe , jplfe, rfact2 )
      CALL trc_sink( kt, wsbio6, sinklfa , jplfa, rfact2 )
      IF ( ln_feauth ) THEN
      wsbio7(:,:,:) = wsafes

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1,jpi
               zmax  = MAX( heup_01(ji,jj), hmld(ji,jj) )
               zfact = MAX( 0., gdepw_n(ji,jj,jk+1) - zmax ) / wsbio2scale
               wsbio8(ji,jj,jk) = wsafeb + MAX(0., ( wsafebmax - wsafeb )) * zfact
            END DO
         END DO
      END DO

      CALL trc_sink( kt, wsbio7, sinkafs , jpafs, rfact2 )
      CALL trc_sink( kt, wsbio8, sinkafb , jpafb, rfact2 )
      ENDIF ! IF ( ln_feauth ) THEN
      ENDIF ! IF( ln_bait ) THEN 

      ! PISCES-QUOTA part
      IF( ln_p5z ) THEN
         sinkingn (:,:,:) = 0.e0
         sinking2n(:,:,:) = 0.e0
         sinkingp (:,:,:) = 0.e0
         sinking2p(:,:,:) = 0.e0

         ! Compute the sedimentation term using trc_sink for all the sinking particles
         ! ---------------------------------------------------------------------------
         CALL trc_sink( kt, wsbio3, sinkingn , jppon, rfact2 )
         CALL trc_sink( kt, wsbio3, sinkingp , jppop, rfact2 )
         CALL trc_sink( kt, wsbio4, sinking2n, jpgon, rfact2 )
         CALL trc_sink( kt, wsbio4, sinking2p, jpgop, rfact2 )
      ENDIF

     ! Total carbon export per year
     IF( iom_use( "tcexp" ) .OR. ( ln_check_mass .AND. kt == nitend .AND. knt == nrdttrc )  )  &
        &   t_oce_co2_exp = glob_sum( 'p4zsink', ( sinking(:,:,ik100) + sinking2(:,:,ik100) ) * e1e2t(:,:) * tmask(:,:,1) )
     !
     IF( lk_iomput ) THEN
       IF( knt == nrdttrc ) THEN
          ALLOCATE( zw2d(jpi,jpj), zw3d(jpi,jpj,jpk) )
          zfact = 1.e+3 * rfact2r  !  conversion from mol/l/kt to  mol/m3/s
          !
          IF( iom_use( "EPC100" ) )  THEN
              zw2d(:,:) = ( sinking(:,:,ik100) + sinking2(:,:,ik100) ) * zfact * tmask(:,:,1) ! Export of carbon at 100m
              CALL iom_put( "EPC100"  , zw2d )
          ENDIF
          IF( iom_use( "EPFE100" ) )  THEN
              zw2d(:,:) = ( sinkfer(:,:,ik100) + sinkfer2(:,:,ik100) ) * zfact * tmask(:,:,1) ! Export of iron at 100m
              CALL iom_put( "EPFE100"  , zw2d )
          ENDIF
          IF( iom_use( "EPCAL100" ) )  THEN
              zw2d(:,:) = sinkcal(:,:,ik100) * zfact * tmask(:,:,1) ! Export of calcite at 100m
              CALL iom_put( "EPCAL100"  , zw2d )
          ENDIF
          IF( iom_use( "EPSI100" ) )  THEN
              zw2d(:,:) =  sinksil(:,:,ik100) * zfact * tmask(:,:,1) ! Export of bigenic silica at 100m
              CALL iom_put( "EPSI100"  , zw2d )
          ENDIF
          IF( iom_use( "EXPC" ) )  THEN
              zw3d(:,:,:) = ( sinking(:,:,:) + sinking2(:,:,:) ) * zfact * tmask(:,:,:) ! Export of carbon in the water column
              CALL iom_put( "EXPC"  , zw3d )
          ENDIF
          IF( iom_use( "EXPFE" ) )  THEN
              zw3d(:,:,:) = ( sinkfer(:,:,:) + sinkfer2(:,:,:) ) * zfact * tmask(:,:,:) ! Export of iron 
              CALL iom_put( "EXPFE"  , zw3d )
          ENDIF
          IF( iom_use( "EXPCAL" ) )  THEN
              zw3d(:,:,:) = sinkcal(:,:,:) * zfact * tmask(:,:,:) ! Export of calcite 
              CALL iom_put( "EXPCAL"  , zw3d )
          ENDIF
          IF( iom_use( "EXPSI" ) )  THEN
              zw3d(:,:,:) = sinksil(:,:,:) * zfact * tmask(:,:,:) ! Export of bigenic silica
              CALL iom_put( "EXPSI"  , zw3d )
          ENDIF
          IF( iom_use( "EXPLFE" ) )  THEN
              zw3d(:,:,:) = sinklfe(:,:,:) * zfact * tmask(:,:,:) ! Export of bigenic silica
              CALL iom_put( "EXPLFE"  , zw3d )
          ENDIF
          IF( iom_use( "EXPLFEA" ) )  THEN
              zw3d(:,:,:) = sinklfa(:,:,:) * zfact * tmask(:,:,:) ! Export of
              CALL iom_put( "EXPLFEA"  , zw3d )
          ENDIF
          IF( iom_use( "EXPAFES" ) )  THEN
              zw3d(:,:,:) = sinkafs(:,:,:) * zfact * tmask(:,:,:) ! Export of
              CALL iom_put( "EXPAFES"  , zw3d )
          ENDIF
          IF( iom_use( "EXPAFEB" ) )  THEN
              zw3d(:,:,:) = sinkafb(:,:,:) * zfact * tmask(:,:,:) ! Export of
              CALL iom_put( "EXPAFEB"  , zw3d )
          ENDIF
          IF( iom_use( "tcexp" ) )  CALL iom_put( "tcexp" , t_oce_co2_exp * zfact )   ! molC/s
          ! 
          DEALLOCATE( zw2d, zw3d )
        ENDIF
      ENDIF
      !
      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('sink')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( ln_timing )   CALL timing_stop('p4z_sink')
      !
   END SUBROUTINE p4z_sink


   SUBROUTINE p4z_sink_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_sink_init  ***
      !!
      !! ** Purpose :   Initialization of sinking parameters
      !!
      !! ** Method  :   
      !!
      !! ** input   :   
      !!----------------------------------------------------------------------
      INTEGER :: jk
      !!----------------------------------------------------------------------
      !
      ik100 = 10        !  last level where depth less than 100 m
      DO jk = jpkm1, 1, -1
         IF( gdept_1d(jk) > 100. )  ik100 = jk - 1
      END DO
      IF (lwp) WRITE(numout,*)
      IF (lwp) WRITE(numout,*) ' Level corresponding to 100m depth ',  ik100 + 1
      IF (lwp) WRITE(numout,*)
      !
      t_oce_co2_exp = 0._wp
      !
   END SUBROUTINE p4z_sink_init

   INTEGER FUNCTION p4z_sink_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_sink_alloc  ***
      !!----------------------------------------------------------------------
      INTEGER :: ierr(2)
      !!----------------------------------------------------------------------
      !
      ierr(:) = 0
      !
      ALLOCATE( sinking(jpi,jpj,jpk) , sinking2(jpi,jpj,jpk)                    ,     &                
         &      sinkcal(jpi,jpj,jpk) , sinksil (jpi,jpj,jpk)                    ,     &                
         &      sinkfer2(jpi,jpj,jpk)                                           ,     &                
         &      sinkfer(jpi,jpj,jpk)                                            , STAT=ierr(1) )                
         !
      IF( ln_p5z    ) ALLOCATE( sinkingn(jpi,jpj,jpk), sinking2n(jpi,jpj,jpk)   ,     &
         &                      sinkingp(jpi,jpj,jpk), sinking2p(jpi,jpj,jpk)   , STAT=ierr(2) )
      !
      IF( ln_bait  ) ALLOCATE( sinklfe(jpi,jpj,jpk) , sinklfa(jpi,jpj,jpk)       ,   &
         &                     sinkafs(jpi,jpj,jpk) , sinkafb(jpi,jpj,jpk)      , STAT=ierr(3) )
      p4z_sink_alloc = MAXVAL( ierr )
      IF( p4z_sink_alloc /= 0 ) CALL ctl_stop( 'STOP', 'p4z_sink_alloc : failed to allocate arrays.' )
      !
   END FUNCTION p4z_sink_alloc
   
   !!======================================================================
END MODULE p4zsink
