MODULE p4zfechem
   !!======================================================================
   !!                         ***  MODULE p4zfechem  ***
   !! TOP :   PISCES Compute iron chemistry and scavenging
   !!======================================================================
   !! History :   3.5  !  2012-07 (O. Aumont, A. Tagliabue, C. Ethe) Original code
   !!             3.6  !  2015-05  (O. Aumont) PISCES quota
   !!----------------------------------------------------------------------
   !!   p4z_fechem       : Compute remineralization/scavenging of iron
   !!   p4z_fechem_init  : Initialisation of parameters for remineralisation
   !!   p4z_fechem_alloc : Allocate remineralisation variables
   !!----------------------------------------------------------------------
   USE oce_trc         ! shared variables between ocean and passive tracers
   USE trc             ! passive tracers common variables 
   USE sms_pisces      ! PISCES Source Minus Sink variables
   USE p4zche          ! chemical model
   USE p4zsbc          ! Boundary conditions from sediments
   USE prtctl_trc      ! print control for debugging
   USE iom             ! I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_fechem        ! called in p4zbio.F90
   PUBLIC   p4z_fechem_init   ! called in trcsms_pisces.F90

   LOGICAL          ::   ln_ligvar    !: boolean for variable ligand concentration following Tagliabue and voelker
   REAL(wp), PUBLIC ::   xlam1        !: scavenging rate of Iron 
   REAL(wp), PUBLIC ::   xlamdust     !: scavenging rate of Iron by dust 
   REAL(wp), PUBLIC ::   ligand       !: ligand concentration in the ocean 
   REAL(wp), PUBLIC ::   kfep         !: rate constant for nanoparticle formation
   REAL(wp), PUBLIC ::   scaveff      !: Fraction of scavenged iron that is considered as being subject to solubilization
   REAL(wp), PUBLIC ::   xlamdust1     !: scavenging rate of Iron by lfe
   REAL(wp), PUBLIC ::   xlamdust2     !: scavenging rate of Iron by lfe aggregates
   LOGICAL          ::   ln_agglfe    !: boolean for colloidal aggregation via lithogenic Fe
   LOGICAL          ::   ln_dyncol    !: 
   REAL(wp), PUBLIC ::   mincolfe     !: 
   LOGICAL          ::   ln_l1l2      !:
   REAL(wp), PUBLIC ::   l2_min       !:
   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zfechem.F90 14276 2021-01-07 22:09:56Z aumont $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_fechem( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_fechem  ***
      !!
      !! ** Purpose :   Compute remineralization/scavenging of iron
      !!
      !! ** Method  :   A simple chemistry model of iron from Aumont and Bopp (2006)
      !!                based on one ligand and one inorganic form
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, knt   ! ocean time step
      !
      INTEGER  ::   ji, jj, jk, jic, jn
      REAL(wp) ::   zlam1a, zlam1b
      REAL(wp) ::   zkeq, zfesatur, zfecoll, fe3sol, zligco
      REAL(wp) ::   zscave, zaggdfea, zaggdfeb, ztrc, zdust, zklight
      REAL(wp) ::   ztfe, zhplus, zxlam, zaggliga, zaggligb
      REAL(wp) ::   zrfact2
      REAL(wp) :: zaggdfec, zaggdfed, zlam1c, zlam1d
      CHARACTER (len=25) :: charout
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zTL1, zFe3, ztotlig, precip, precipno3, zFeL1
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zcoll3d, zscav3d, zlcoll3d, zprecip3d
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zcolll3d, zscavl3d, zlcolll3d, zcfe
      REAL(wp) :: zscave1, zscave2, zdust2, zl1, zl2, zk1, zk2, zdoct, zpht, fzl1
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('p4z_fechem')
      !
      zFe3 (:,:,:) = 0.   ;   zFeL1(:,:,:) = 0.
      zTL1 (:,:,:) = 0.   ;   zcfe(:,:,:) = 0.

      ! Total ligand concentration : Ligands can be chosen to be constant or variable
      ! Parameterization from Pham and Ito (2018)
      ! -------------------------------------------------
      IF( ln_ligvar ) THEN
         ztotlig(:,:,:) =  0.09 * 0.667 * trb(:,:,:,jpdoc) * 1E6 + ligand * 1E9 + MAX(0., chemo2(:,:,:) - trb(:,:,:,jpoxy) ) / 400.E-6
         ztotlig(:,:,:) =  MIN( ztotlig(:,:,:), 10. )
      ELSE
        IF( ln_ligand ) THEN  ;   ztotlig(:,:,:) = trb(:,:,:,jplgw) * 1E9
        ELSE                  ;   ztotlig(:,:,:) = ligand * 1E9 
        ENDIF
      ENDIF

      IF( ln_bait ) THEN
      zcolll3d(:,:,:) = 0. ; zscavl3d(:,:,:) = 0. ; zlcolll3d(:,:,:) = 0.
      ENDIF
      ! simple param to account for L1 and L2 ligands as well as effects of 
      ! pH and DOC based on Ye et al., 2020 GBC doi: 10.1029/2019GB006425 
      IF( ln_bait .AND. ln_l1l2 ) THEN
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi

         zl2 = min ( trb(ji,jj,jk, jplgw) , l2_min )
         zl1 = max( 0., trb(ji,jj,jk, jplgw) - zl2 )
         zdoct = max( trb(ji,jj,jk, jpdoc) + 40.E-6, rtrn )
         zpht = -1. * LOG10( MAX( hi(ji,jj,jk), rtrn ) )
         zk2 = ((-1)*2.E-4*zdoct + 0.034 )  * zdoct + ( zpht * (-1) * 1.67 ) + 24.36
         zk1 = zk2 +  2.67
         fzl1 = (zl1/(zl1+zl2 + rtrn ) )

         fekeq(ji,jj,jk) = 10**( (zk1 * fzl1 ) + ( zk2 * (1. - fzl1) ) )

       ENDDO
       ENDDO
       ENDDO
       ENDIF
      ! ------------------------------------------------------------
      !  from Aumont and Bopp (2006)
      ! This model is based on one ligand, Fe2+ and Fe3+ 
      ! Chemistry is supposed to be fast enough to be at equilibrium
      ! ------------------------------------------------------------
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               zTL1(ji,jj,jk)  = ztotlig(ji,jj,jk)
               zkeq            = fekeq(ji,jj,jk)
               zklight         = 4.77E-7 * etot(ji,jj,jk) * 0.5 / 10**-6.3
               zfesatur        = zTL1(ji,jj,jk) * 1E-9
               ztfe            = (1.0 + zklight) * trb(ji,jj,jk,jpfer) 
               ! Fe' is the root of a 2nd order polynom
               zFe3 (ji,jj,jk) = ( -( 1. + zfesatur * zkeq + zklight + consfe3(ji,jj,jk)/10**-6.3 - zkeq * trb(ji,jj,jk,jpfer) )               &
                  &              + SQRT( ( 1. + zfesatur * zkeq + zklight + consfe3(ji,jj,jk)/10**-6.3 - zkeq * trb(ji,jj,jk,jpfer) )**2       &
                  &              + 4. * ztfe * zkeq) ) / ( 2. * zkeq )
               zFeL1(ji,jj,jk) = MAX( 0., trb(ji,jj,jk,jpfer) - zFe3(ji,jj,jk) )
           END DO
         END DO
      END DO
      !
      plig(:,:,:) =  MAX( 0., ( zFeL1(:,:,:) / ( trb(:,:,:,jpfer) + rtrn ) ) )
      !
      zdust = 0.         ! if no dust available
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               ! Scavenging rate of iron. This scavenging rate depends on the load of particles of sea water. 
               ! This parameterization assumes a simple second order kinetics (k[Particles][Fe]).
               ! Scavenging onto dust is also included as evidenced from the DUNE experiments.
               ! --------------------------------------------------------------------------------------
               zhplus  = max( rtrn, hi(ji,jj,jk) )
               fe3sol  = fesol(ji,jj,jk,1) * ( zhplus**3 + fesol(ji,jj,jk,2) * zhplus**2  &
               &         + fesol(ji,jj,jk,3) * zhplus + fesol(ji,jj,jk,4)     &
               &         + fesol(ji,jj,jk,5) / zhplus )
               !
               ! option to calculate the colloidal fraction dynamically to vary
               ! between 50-90% depending of Fe solublity, operates only on FeL
               IF ( ln_bait .AND. ln_dyncol ) THEN
               zfecoll = MIN( MAX( ( zFeL1(ji,jj,jk) - fe3sol ) , mincolfe * zFeL1(ji,jj,jk) ) , ( 0.9 * zFeL1(ji,jj,jk) ) )
               ELSE
               zfecoll = 0.5 * zFeL1(ji,jj,jk)
               ENDIF
               zcfe(ji,jj,jk) = zfecoll
               ! precipitation of Fe3+, creation of nanoparticles, operates on
               ! uncomplexed Fe3+
               precip(ji,jj,jk) = MAX( 0., ( zFe3(ji,jj,jk) - fe3sol ) ) * kfep * xstep * ( 1.0 - nitrfac(ji,jj,jk) ) 
               ! Precipitation of Fe2+ due to oxidation by NO3 (Croot et al., 2019)
               ! This occurs in anoxic waters only
               precipno3(ji,jj,jk) = 2.0 * 130.0 * trb(ji,jj,jk,jpno3) * nitrfac(ji,jj,jk) * xstep * zFe3(ji,jj,jk)
               !
               ztrc   = ( trb(ji,jj,jk,jppoc) + trb(ji,jj,jk,jpgoc) + trb(ji,jj,jk,jpcal) + trb(ji,jj,jk,jpgsi) ) * 1.e6 
               ztrc = MAX( rtrn, ztrc )

               IF( ln_bait .AND. ln_dust ) THEN 
                   zdust  = trb(ji,jj,jk,jplfe) / 0.035 * 55.85 ! converted to g/L 
                   zdust2 = trb(ji,jj,jk,jplfa) / 0.035 * 55.85 ! converted to g/L
               ELSE
               IF( ln_dust )  zdust  = dust(ji,jj) / ( wdust / rday ) * tmask(ji,jj,jk) ! g/L
               ENDIF

               IF( ln_bait) THEN
               ! Scavenging of Fe3+ by organic and lithogenic particles
               zxlam  = MAX( 1.E-3, (1. - EXP(-2 * trb(ji,jj,jk,jpoxy) / 100.E-6 ) ))
               zlam1b = 3.e-5 + ( xlam1 * ztrc ) * zxlam ! organic particles
               zscave = zFe3(ji,jj,jk) * zlam1b * xstep
               ! lithogenic particles
               zscave1 = zFe3(ji,jj,jk) * ( xlamdust1 * zdust ) * zxlam * xstep
               zscave2 = zFe3(ji,jj,jk) * ( xlamdust2 * zdust2 ) * zxlam * xstep
               ELSE
               ! Scavenging of Fe3+ by organic particles only
               zxlam  = MAX( 1.E-3, (1. - EXP(-2 * trb(ji,jj,jk,jpoxy) / 100.E-6 ) ))
               zlam1b = 3.e-5 + ( xlamdust * zdust + xlam1 * ztrc ) * zxlam
               zscave = zFe3(ji,jj,jk) * zlam1b * xstep
               ENDIF
               ! Compute the coagulation of colloidal iron. This parameterization 
               ! could be thought as an equivalent of colloidal pumping.
               ! It requires certainly some more work as it is very poorly constrained.
               ! ----------------------------------------------------------------
               zlam1a   = ( 12.0  * 0.3 * trb(ji,jj,jk,jpdoc) + 9.05  * trb(ji,jj,jk,jppoc) ) * xdiss(ji,jj,jk)    &
                   &    + ( 2.49  * trb(ji,jj,jk,jppoc) )     &
                   &    + ( 127.8 * 0.3 * trb(ji,jj,jk,jpdoc) + 725.7 * trb(ji,jj,jk,jppoc) )
               zaggdfea = zlam1a * xstep * zfecoll
               !
               zlam1b   = ( 1.94 * xdiss(ji,jj,jk) + 1.37 ) * trb(ji,jj,jk,jpgoc)
               zaggdfeb = zlam1b * xstep * zfecoll
               !
               IF( ln_bait .AND. ln_agglfe ) THEN
               ! aggegation of colloidal Fe with LFe and LFea
               ! First term represents aggregation with Lfe and second with Lfea
               ! converts both lithogenic pools to carbon equivalents and, for
               ! now, assumes same aggregation kernals as for POC and GOC
               zdust  = trb(ji,jj,jk,jplfe) / 0.035 * 55.85 ! converted to g/L
               zdust2 = trb(ji,jj,jk,jplfa) / 0.035 * 55.85 ! converted to g/L
!               zlam1c = ( 9.05  * xdiss(ji,jj,jk) * ( zdust / 24. ) ) + ( 2.49 * ( zdust / 24. ) ) &
!               &        + ( 725.7 * ( zdust / 24. ) )
!               zlam1d = ( ( 1.94 * xdiss(ji,jj,jk) + 1.37 ) * ( zdust2 / 24. ) )
               zlam1c = ( 18.1  * xdiss(ji,jj,jk) * ( zdust / 24. ) ) + ( 4.98 * ( zdust / 24. ) ) &
               &        + ( 1451.4 * ( zdust / 24. ) )
               zlam1d = ( ( 3.88 * xdiss(ji,jj,jk) + 2.74 ) * ( zdust2 / 24. ) )
               zaggdfec = zlam1c * xstep * zfecoll
               zaggdfed = zlam1d * xstep * zfecoll

               tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) - zscave - zaggdfea - zaggdfeb &
               &                     - precip(ji,jj,jk) - precipno3(ji,jj,jk)           &
               &                     - zscave1 - zscave2 - zaggdfec - zaggdfed
               tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) + zscave * scaveff * trb(ji,jj,jk,jppoc) / ztrc
               tra(ji,jj,jk,jpbfe) = tra(ji,jj,jk,jpbfe) + zscave * scaveff * trb(ji,jj,jk,jpgoc) / ztrc
               tra(ji,jj,jk,jplfe) = tra(ji,jj,jk,jplfe) + zscave1 + zaggdfec 
               tra(ji,jj,jk,jplfa) = tra(ji,jj,jk,jplfa) + zscave2 + zaggdfed
               zscavl3d(ji,jj,jk)   = zscave1 + zscave2
               zcolll3d(ji,jj,jk)   = zaggdfec + zaggdfed
               ELSE
               tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) - zscave - zaggdfea - zaggdfeb &
               &                     - precip(ji,jj,jk) - precipno3(ji,jj,jk)
               tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) + zscave * scaveff * trb(ji,jj,jk,jppoc) / ztrc
               tra(ji,jj,jk,jpbfe) = tra(ji,jj,jk,jpbfe) + zscave * scaveff * trb(ji,jj,jk,jpgoc) / ztrc
               ENDIF
               ! Precipitated iron is supposed to be permanently lost.
               ! Scavenged iron is supposed to be released back to seawater
               ! when POM is solubilized. This is highly uncertain as probably
               ! a significant part of it may be rescavenged back onto 
               ! the particles. An efficiency factor is applied that is read
               ! in the namelist. 
               ! See for instance Tagliabue et al. (2019).
               ! Aggregated FeL is considered as biogenic Fe as it 
               ! probably remains  complexed when the particle is solubilized.
               ! -------------------------------------------------------------
               tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) + zaggdfea
               tra(ji,jj,jk,jpbfe) = tra(ji,jj,jk,jpbfe) + zaggdfeb
               zscav3d(ji,jj,jk)   = zscave 
               zcoll3d(ji,jj,jk)   = zaggdfea + zaggdfeb
               zprecip3d(ji,jj,jk) = precip(ji,jj,jk) + precipno3(ji,jj,jk)
               !
            END DO
         END DO
      END DO
      !
      !  Define the bioavailable fraction of iron
      !  ----------------------------------------
      biron(:,:,:) = trb(:,:,:,jpfer) 
      !
      IF( ln_ligand ) THEN
         !
               IF( ln_bait .AND. ln_agglfe ) THEN
         DO jk = 1, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
               zhplus  = max( rtrn, hi(ji,jj,jk) )
               fe3sol  = fesol(ji,jj,jk,1) * ( zhplus**3 + fesol(ji,jj,jk,2) * zhplus**2  &
               &         + fesol(ji,jj,jk,3) * zhplus + fesol(ji,jj,jk,4)     &
               &         + fesol(ji,jj,jk,5) / zhplus )
               ! aggegation of colloidal FeL with LFe and LFea
               ! First term represents aggregation with Lfe and second with Lfea
               ! converts both lithogenic pools to carbobn equivalents and, for
               ! now, assumes same aggregation kernals as for POC and GOC
               zdust  = trb(ji,jj,jk,jplfe) / 0.035 * 55.85 ! converted to g/L
               zdust2 = trb(ji,jj,jk,jplfa) / 0.035 * 55.85 ! converted to g/L
!              zlam1c = ( 9.05  * xdiss(ji,jj,jk) * ( zdust / 24. ) ) + ( 2.49 * ( zdust / 24. ) )
!               zlam1d = ( ( 1.94 * xdiss(ji,jj,jk) + 1.37 ) * ( zdust2 / 24. ) )
               zlam1c = ( 18.1  * xdiss(ji,jj,jk) * ( zdust / 24. ) ) + ( 4.98 * ( zdust / 24. ) ) &
               &        + ( 1451.4 * ( zdust / 24. ) )
               zlam1d = ( ( 3.88 * xdiss(ji,jj,jk) + 2.74 ) * ( zdust2 / 24. ) )
                  IF ( ln_dyncol ) THEN
                  zligco   = MIN( MAX( ( zFeL1(ji,jj,jk) - fe3sol ) , mincolfe * zFeL1(ji,jj,jk) ) , ( 0.9 * zFeL1(ji,jj,jk) ) ) 
                  ELSE
                  zligco   = 0.5 * trb(ji,jj,jk,jplgw)
                  ENDIF
               zaggdfec = zlam1c * xstep * zfecoll
               zaggdfed = zlam1d * xstep * zfecoll
                  tra(ji,jj,jk,jplgw) = tra(ji,jj,jk,jplgw) - zaggdfec - zaggdfed
                  zlcolll3d(ji,jj,jk) = zaggdfec + zaggdfed
               END DO
            END DO
         END DO
         ENDIF


         DO jk = 1, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  ! Coagulation of ligands due to various processes (Brownian, shear, diff. sedimentation
                  ! Coefficients are taken from p4zagg
                  ! -------------------------------------------------------------------------------------
                  zlam1a   = ( 12.0  * 0.3 * trb(ji,jj,jk,jpdoc) + 9.05  * trb(ji,jj,jk,jppoc) ) * xdiss(ji,jj,jk)    &
                      &    + ( 2.49  * trb(ji,jj,jk,jppoc) )     &
                      &    + ( 127.8 * 0.3 * trb(ji,jj,jk,jpdoc) + 725.7 * trb(ji,jj,jk,jppoc) )
                  !
                  zlam1b   = ( 1.94 * xdiss(ji,jj,jk) + 1.37 ) * trb(ji,jj,jk,jpgoc)
                  ! 50% of the ligands are supposed to be in the colloidal size fraction
                  ! as for FeL, or can be set dynamically between 50-90%
                  ! depending on solubility of Fe
                  IF ( ln_dyncol ) THEN
                  zligco   = MIN( MAX( ( zFeL1(ji,jj,jk) - fe3sol ) , mincolfe * zFeL1(ji,jj,jk) ) , ( 0.9 * zFeL1(ji,jj,jk) ) )
                  ELSE
                  zligco   = 0.5 * trb(ji,jj,jk,jplgw)
                  ENDIF
                  zaggliga = zlam1a * xstep * zligco 
                  zaggligb = zlam1b * xstep * zligco
                  tra(ji,jj,jk,jplgw) = tra(ji,jj,jk,jplgw) - zaggliga - zaggligb
                  zlcoll3d(ji,jj,jk)  = zaggliga + zaggligb
               END DO
            END DO
         END DO
         !
      ENDIF
      !  Output of some diagnostics variables
      !     ---------------------------------
      IF( lk_iomput .AND. knt == nrdttrc ) THEN
         zrfact2 = 1.e3 * rfact2r  ! conversion from mol/L/timestep into mol/m3/s
         IF( iom_use("Fe3")    )  CALL iom_put("Fe3"    , zFe3   (:,:,:)       * tmask(:,:,:) )   ! Fe3+
         IF( iom_use("FeL1")   )  CALL iom_put("FeL1"   , zFeL1  (:,:,:)       * tmask(:,:,:) )   ! FeL1
         IF( iom_use("TL1")    )  CALL iom_put("TL1"    , zTL1   (:,:,:)       * tmask(:,:,:) )   ! TL1
         IF( iom_use("Totlig") )  CALL iom_put("Totlig" , ztotlig(:,:,:)       * tmask(:,:,:) )   ! TL
         IF( iom_use("Biron")  )  CALL iom_put("Biron"  , biron  (:,:,:)  * 1e9 * tmask(:,:,:) )   ! biron
         IF( iom_use("FESCAV") )  CALL iom_put("FESCAV" , zscav3d(:,:,:)  * 1e9 * tmask(:,:,:) * zrfact2 )
         IF( iom_use("FECOLL") )  CALL iom_put("FECOLL" , zcoll3d(:,:,:)  * 1e9 * tmask(:,:,:) * zrfact2 )
         IF( iom_use("FEPREC") )  CALL iom_put("FEPREC" , zprecip3d(:,:,:)  * 1e9 * tmask(:,:,:) * zrfact2 )
         IF( iom_use("LGWCOLL"))  CALL iom_put("LGWCOLL", zlcoll3d(:,:,:) * 1e9 * tmask(:,:,:) * zrfact2 )
         IF( iom_use("FESCAV2") )  CALL iom_put("FESCAV2" , zscavl3d(:,:,:)  * 1e9 * tmask(:,:,:) * zrfact2 )
         IF( iom_use("FECOLL2") )  CALL iom_put("FECOLL2" , zcolll3d(:,:,:)  * 1e9 * tmask(:,:,:) * zrfact2 )
         IF( iom_use("LGWCOLL2"))  CALL iom_put("LGWCOLL2", zlcolll3d(:,:,:) * 1e9 * tmask(:,:,:) * zrfact2 )
         IF( iom_use("cFe")    )  CALL iom_put("cFe"    , zcfe   (:,:,:)       * tmask(:,:,:) )   ! Fe3+
         IF( iom_use("LogK")   )  CALL iom_put("LogK"   , LOG10(fekeq(:,:,:))       * tmask(:,:,:) )   !
      ENDIF

      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('fechem')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( ln_timing )   CALL timing_stop('p4z_fechem')
      !
   END SUBROUTINE p4z_fechem


   SUBROUTINE p4z_fechem_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_fechem_init  ***
      !!
      !! ** Purpose :   Initialization of iron chemistry parameters
      !!
      !! ** Method  :   Read the nampisfer namelist and check the parameters
      !!      called at the first timestep
      !!
      !! ** input   :   Namelist nampisfer
      !!
      !!----------------------------------------------------------------------
      INTEGER ::   ios   ! Local integer 
      !!
      NAMELIST/nampisfer/ ln_ligvar, xlam1, xlamdust, ligand, kfep, scaveff, xlamdust1, xlamdust2, & 
      &                   ln_agglfe, ln_dyncol, mincolfe, ln_l1l2, l2_min 
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'p4z_rem_init : Initialization of iron chemistry parameters'
         WRITE(numout,*) '~~~~~~~~~~~~'
      ENDIF
      !
      REWIND( numnatp_ref )            ! Namelist nampisfer in reference namelist : Pisces iron chemistry
      READ  ( numnatp_ref, nampisfer, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'nampisfer in reference namelist' )
      REWIND( numnatp_cfg )            ! Namelist nampisfer in configuration namelist : Pisces iron chemistry
      READ  ( numnatp_cfg, nampisfer, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'nampisfer in configuration namelist' )
      IF(lwm) WRITE( numonp, nampisfer )

      IF(lwp) THEN                     ! control print
         WRITE(numout,*) '   Namelist : nampisfer'
         WRITE(numout,*) '      variable concentration of ligand          ln_ligvar    =', ln_ligvar
         WRITE(numout,*) '      scavenging rate of Iron                   xlam1        =', xlam1
         WRITE(numout,*) '      scavenging rate of Iron by dust           xlamdust     =', xlamdust
         WRITE(numout,*) '      ligand concentration in the ocean         ligand       =', ligand
         WRITE(numout,*) '      rate constant for nanoparticle formation  kfep         =', kfep
         WRITE(numout,*) '      Scavenged iron that is added to POFe      scaveff      =', scaveff
      IF( ln_bait) THEN
         WRITE(numout,*) '      scavenging rate of Iron by lfe            xlamdust1     =', xlamdust1
         WRITE(numout,*) '      scavenging rate of Iron by lfe aggregates xlamdust2     =', xlamdust2
         WRITE(numout,*) '      dynamic computation of colloidal Fe        ln_dyncol =', ln_dyncol
         WRITE(numout,*) ' minimum colloidal Fe fraction mincolfe=',mincolfe
         WRITE(numout,*) ' fancy l1 and l2 ligands?          ln_l1l2',ln_l1l2
         WRITE(numout,*) ' set value for l2                  l2_min',l2_min
      ENDIF
      ENDIF
      ! 
   END SUBROUTINE p4z_fechem_init
   
   !!======================================================================
END MODULE p4zfechem
