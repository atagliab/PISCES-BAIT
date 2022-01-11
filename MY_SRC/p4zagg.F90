MODULE p4zagg
   !!======================================================================
   !!                         ***  MODULE p4zagg  ***
   !! TOP :  PISCES  aggregation of particles (DOC, POC, GOC)
   !!        This module is the same for both PISCES and PISCES-QUOTA
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Change aggregation formula
   !!             3.5  !  2012-07  (O. Aumont) Introduce potential time-splitting
   !!             3.6  !  2015-05  (O. Aumont) PISCES quota
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   p4z_agg       :  Compute aggregation of particles
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE prtctl_trc      !  print control for debugging

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_agg         ! called in p4zbio.F90

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zagg.F90 14276 2021-01-07 22:09:56Z aumont $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_agg ( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_agg  ***
      !!
      !! ** Purpose :   Compute aggregation of particle. Aggregation by 
      !!                brownian motion, differential settling and shear
      !!                are considered.
      !!
      !! ** Method  : - Aggregation rates are computed assuming a fixed and 
      !!                constant size spectrum in the different particulate 
      !!                pools. The coagulation rates have been computed 
      !!                externally using dedicated programs (O. Aumont). They 
      !!                are hard-coded because they can't be changed 
      !!                independently of each other. 
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, knt   !
      !
      INTEGER  ::   ji, jj, jk
      REAL(wp) ::   zagg, zagg1, zagg2, zagg3, zagg4
      REAL(wp) ::   zaggpoc1, zaggpoc2, zaggpoc3, zaggpoc4
      REAL(wp) ::   zaggpoc , zaggfe, zaggdoc, zaggdoc2, zaggdoc3
      REAL(wp) ::   zaggpon , zaggdon, zaggdon2, zaggdon3
      REAL(wp) ::   zaggpop, zaggdop, zaggdop2, zaggdop3
      REAL(wp) ::   zaggtmp, zfact, zmax
      REAL(wp) ::   r1_rday    ! 1 / rday
      REAL(wp) ::   zlfe, zagglfe, zagglfe2, zdislfea
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('p4z_agg')

      !
      !  Exchange between organic matter compartments due to 
      !  coagulation/disaggregation
      !  ---------------------------------------------------
      r1_rday   = 1._wp / rday
      ! PISCES part
      IF( ln_p4z ) THEN
         !
         DO jk = 1, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  !
                  zfact = xstep * xdiss(ji,jj,jk)
                  ! Part I : Coagulation dependent on turbulence
                  ! The stickiness has been assumed to be 0.1
                  zagg1 = 12.5  * zfact * trb(ji,jj,jk,jppoc) * trb(ji,jj,jk,jppoc)
                  zagg2 = 169.7 * zfact * trb(ji,jj,jk,jppoc) * trb(ji,jj,jk,jpgoc)

                  ! Part II : Differential settling
                  ! Aggregation of small into large particles
                  ! The stickiness has been assumed to be 0.1
                  zagg3 =  8.63  * xstep * trb(ji,jj,jk,jppoc) * trb(ji,jj,jk,jppoc)
                  zagg4 =  132.8 * xstep * trb(ji,jj,jk,jppoc) * trb(ji,jj,jk,jpgoc)

                  zagg   = zagg1 + zagg2 + zagg3 + zagg4
                  zaggfe = zagg * trb(ji,jj,jk,jpsfe) / ( trb(ji,jj,jk,jppoc) + rtrn )

                  ! Aggregation of DOC to POC : 
                  ! 1st term is shear aggregation of DOC-DOC
                  ! 2nd term is shear aggregation of DOC-POC
                  ! 3rd term is differential settling of DOC-POC
                  ! 1/3 of DOC is supposed to experience aggregation (HMW)
                  zaggdoc  = ( ( 12.0 * 0.3 * trb(ji,jj,jk,jpdoc) + 9.05 * trb(ji,jj,jk,jppoc) ) * zfact       &
                  &            + 2.49 * xstep * trb(ji,jj,jk,jppoc) ) * 0.3 * trb(ji,jj,jk,jpdoc)
                  ! transfer of DOC to GOC : 
                  ! 1st term is shear aggregation
                  ! 1/3 of DOC is supposed to experience aggregation (HMW)
                  zaggdoc2 = ( 1.94 * zfact + 1.37 * xstep ) * trb(ji,jj,jk,jpgoc) * 0.3 * trb(ji,jj,jk,jpdoc)
                  ! tranfer of DOC to POC due to brownian motion
                  ! The temperature dependency has been omitted.
                  zaggdoc3 =  ( 127.8 * 0.3 * trb(ji,jj,jk,jpdoc) + 725.7 * trb(ji,jj,jk,jppoc) ) * xstep * 0.3 * trb(ji,jj,jk,jpdoc)

                  !  Update the trends
                  tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) - zagg + zaggdoc + zaggdoc3
                  tra(ji,jj,jk,jpgoc) = tra(ji,jj,jk,jpgoc) + zagg + zaggdoc2
                  tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) - zaggfe
                  tra(ji,jj,jk,jpbfe) = tra(ji,jj,jk,jpbfe) + zaggfe
                  tra(ji,jj,jk,jpdoc) = tra(ji,jj,jk,jpdoc) - zaggdoc - zaggdoc2 - zaggdoc3
                  !
                  IF (ln_bait .AND. ln_feauth ) THEN
      ! aggregation of small authigenic Fe into big authigenic Fe
                  tra(ji,jj,jk,jpafs) = tra(ji,jj,jk,jpafs) - zagg * ( trb(ji,jj,jk,jpafs) / ( trb(ji,jj,jk,jppoc) + rtrn ) )
                  tra(ji,jj,jk,jpafb) = tra(ji,jj,jk,jpafb) + zagg * ( trb(ji,jj,jk,jpafs) / ( trb(ji,jj,jk,jppoc) + rtrn ) )
      ! (2) autocatalytic
!                  tra(ji,jj,jk,jpafs) = tra(ji,jj,jk,jpafs) - cfeagg * trb(ji,jj,jk,jpafs) * zfact
!                  tra(ji,jj,jk,jpafb) = tra(ji,jj,jk,jpafb) + cfeagg * trb(ji,jj,jk,jpafs) * zfact
                  tra(ji,jj,jk,jpafs) = tra(ji,jj,jk,jpafs) - ( cfeagg * trb(ji,jj,jk,jpafs) * zfact ) &
                  &                   - ( cfeagg2 * trb(ji,jj,jk,jpafb) * xstep )
                  tra(ji,jj,jk,jpafb) = tra(ji,jj,jk,jpafb) + ( cfeagg * trb(ji,jj,jk,jpafs) * zfact ) &
                  &                   + ( cfeagg2 * trb(ji,jj,jk,jpafb) * xstep ) 
                 ENDIF
                  conspoc(ji,jj,jk) = conspoc(ji,jj,jk) - zagg + zaggdoc + zaggdoc3
                  prodgoc(ji,jj,jk) = prodgoc(ji,jj,jk) + zagg + zaggdoc2
                  !
               END DO
            END DO
         END DO
      ELSE    ! ln_p5z
        !
        ! PISCES-QUOTA part
         DO jk = 1, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  !
                  zfact = xstep * xdiss(ji,jj,jk)
                  ! Part I : Coagulation dependent on turbulence
                  ! The stickiness has been assumed to be 0.1
                  zaggtmp = 25.9  * zfact * trb(ji,jj,jk,jppoc)
                  zaggpoc1 = zaggtmp * trb(ji,jj,jk,jppoc)
                  zaggtmp = 4452. * zfact * trb(ji,jj,jk,jpgoc)
                  zaggpoc2 = zaggtmp * trb(ji,jj,jk,jppoc)

                  ! Part II : Differential settling
                  ! The stickiness has been assumed to be 0.1
   
                  !  Aggregation of small into large particles
                  zaggtmp =  47.1 * xstep * trb(ji,jj,jk,jpgoc)
                  zaggpoc3 = zaggtmp * trb(ji,jj,jk,jppoc)
                  zaggtmp =  3.3  * xstep * trb(ji,jj,jk,jppoc)
                  zaggpoc4 = zaggtmp * trb(ji,jj,jk,jppoc)

                  zaggpoc   = zaggpoc1 + zaggpoc2 + zaggpoc3 + zaggpoc4
                  zaggpon = zaggpoc * trb(ji,jj,jk,jppon) / ( trb(ji,jj,jk,jppoc) + rtrn)
                  zaggpop = zaggpoc * trb(ji,jj,jk,jppop) / ( trb(ji,jj,jk,jppoc) + rtrn)
                  zaggfe = zaggpoc * trb(ji,jj,jk,jpsfe) / ( trb(ji,jj,jk,jppoc)  + rtrn )

                  ! Aggregation of DOC to POC : 
                  ! 1st term is shear aggregation of DOC-DOC
                  ! 2nd term is shear aggregation of DOC-POC
                  ! 3rd term is differential settling of DOC-POC
                  ! 1/3 of DOC is supposed to experience aggregation (HMW)
                  zaggtmp = ( ( 0.37 * 0.3 * trb(ji,jj,jk,jpdoc) + 20.5 * trb(ji,jj,jk,jppoc) ) * zfact       &
                  &            + 0.15 * xstep * trb(ji,jj,jk,jppoc) )
                  zaggdoc  = zaggtmp * 0.3 * trb(ji,jj,jk,jpdoc)
                  zaggdon  = zaggtmp * 0.3 * trb(ji,jj,jk,jpdon)
                  zaggdop  = zaggtmp * 0.3 * trb(ji,jj,jk,jpdop)

                  ! transfer of DOC to GOC : 
                  ! 1st term is shear aggregation
                  ! 2nd term is differential settling 
                  ! 1/3 of DOC is supposed to experience aggregation (HMW)
                  zaggtmp = 655.4 * zfact * trb(ji,jj,jk,jpgoc)
                  zaggdoc2 = zaggtmp * 0.3 * trb(ji,jj,jk,jpdoc)
                  zaggdon2 = zaggtmp * 0.3 * trb(ji,jj,jk,jpdon)
                  zaggdop2 = zaggtmp * 0.3 * trb(ji,jj,jk,jpdop)

                  ! tranfer of DOC to POC due to brownian motion
                  ! 1/3 of DOC is supposed to experience aggregation (HMW)
                  zaggtmp = ( 260.2 * 0.3 * trb(ji,jj,jk,jpdoc) +  418.5 * trb(ji,jj,jk,jppoc) ) * xstep
                  zaggdoc3 =  zaggtmp * 0.3 * trb(ji,jj,jk,jpdoc)
                  zaggdon3 =  zaggtmp * 0.3 * trb(ji,jj,jk,jpdon)
                  zaggdop3 =  zaggtmp * 0.3 * trb(ji,jj,jk,jpdop)

                  !  Update the trends
                  tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) - zaggpoc + zaggdoc + zaggdoc3
                  tra(ji,jj,jk,jppon) = tra(ji,jj,jk,jppon) - zaggpon + zaggdon + zaggdon3
                  tra(ji,jj,jk,jppop) = tra(ji,jj,jk,jppop) - zaggpop + zaggdop + zaggdop3
                  tra(ji,jj,jk,jpgoc) = tra(ji,jj,jk,jpgoc) + zaggpoc + zaggdoc2
                  tra(ji,jj,jk,jpgon) = tra(ji,jj,jk,jpgon) + zaggpon + zaggdon2
                  tra(ji,jj,jk,jpgop) = tra(ji,jj,jk,jpgop) + zaggpop + zaggdop2
                  tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) - zaggfe
                  tra(ji,jj,jk,jpbfe) = tra(ji,jj,jk,jpbfe) + zaggfe
                  tra(ji,jj,jk,jpdoc) = tra(ji,jj,jk,jpdoc) - zaggdoc - zaggdoc2 - zaggdoc3
                  tra(ji,jj,jk,jpdon) = tra(ji,jj,jk,jpdon) - zaggdon - zaggdon2 - zaggdon3
                  tra(ji,jj,jk,jpdop) = tra(ji,jj,jk,jpdop) - zaggdop - zaggdop2 - zaggdop3
                  IF (ln_bait .AND. ln_feauth ) THEN
      ! aggregation of small authigenic Fe into big authigenic Fe
      ! (1) catalysed by organic matter
                  tra(ji,jj,jk,jpafs) = tra(ji,jj,jk,jpafs) - zaggpoc * ( trb(ji,jj,jk,jpafs) / ( trb(ji,jj,jk,jppoc) + rtrn ) )
                  tra(ji,jj,jk,jpafb) = tra(ji,jj,jk,jpafb) + zaggpoc * ( trb(ji,jj,jk,jpafs) / ( trb(ji,jj,jk,jppoc) + rtrn ) )
      ! (2) autocatalytic
!                  tra(ji,jj,jk,jpafs) = tra(ji,jj,jk,jpafs) - cfeagg * trb(ji,jj,jk,jpafs) * zfact
!                  tra(ji,jj,jk,jpafb) = tra(ji,jj,jk,jpafb) + cfeagg * trb(ji,jj,jk,jpafs) * zfact
                  tra(ji,jj,jk,jpafs) = tra(ji,jj,jk,jpafs) - ( cfeagg * trb(ji,jj,jk,jpafs) * zfact ) &
                  &                   - ( cfeagg2  * trb(ji,jj,jk,jpafb) * xstep )
                  tra(ji,jj,jk,jpafb) = tra(ji,jj,jk,jpafb) + ( cfeagg * trb(ji,jj,jk,jpafs) * zfact ) &
                  &                   + ( cfeagg2  * trb(ji,jj,jk,jpafb) * xstep )
                  ENDIF
                  !
                  conspoc(ji,jj,jk) = conspoc(ji,jj,jk) - zaggpoc + zaggdoc + zaggdoc3
                  prodgoc(ji,jj,jk) = prodgoc(ji,jj,jk) + zaggpoc + zaggdoc2
                  !
               END DO
            END DO
         END DO
         !
      ENDIF

      IF(ln_bait) THEN ! Bait iron submodule
         DO jk = 1, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
      ! aggregation of lithogenic Fe (converted from molFe/L ->  mg(dust)/m3
      ! assuming 3.5% fixed Fe content of dust and atomic mass of Fe)
      ! organic matter converted from mol/L -> mg/m3 using 24 g per mol C
      ! constants are in (mg/m3)-1 d-1 and come from Ye and Voelker (2019, GBC)
      ! As constants are per day, use zfact with shear and xstep without
      ! option
      ! to include some brownian aggregation and extra term is added to zagglfe
                 zlfe = ( trb(ji,jj,jk,jplfe) / 0.035 * 55.85 ) * 1e6 ! mg/m3 dust
                 ! turbulent aggregation accounts for shear
                 zagglfe = zlfe * ( ragglfe * zfact ) ! produes a per kt rate
                 zagglfe = trb(ji,jj,jk,jplfe) * zagglfe ! molLFe/L/kt
                 !! aggregation with lfe and poc/goc not due to shear
                 zagglfe2 = ( zlfe + ( ( trb(ji,jj,jk,jppoc) + trb(ji,jj,jk,jpgoc) ) * 24.E6 ) ) &
                 &        *  ( ragglfe * aggrat ) * xstep ! produces a per kt rate
                 zagglfe2 = trb(ji,jj,jk,jplfe) * zagglfe2 ! molLFe/L/kt
      ! turbulent dissaggration of LFe aggregates
                 zdislfea = trb(ji,jj,jk,jplfa) * rdislfea * zfact ! molLFea/L/kt
      ! update SMS
                 tra(ji,jj,jk,jplfe) = tra(ji,jj,jk,jplfe) - zagglfe - zagglfe2 + zdislfea 
                 tra(ji,jj,jk,jplfa) = tra(ji,jj,jk,jplfa) + zagglfe + zagglfe2 - zdislfea
               END DO
            END DO
         END DO
      ENDIF
      !
      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('agg')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( ln_timing )   CALL timing_stop('p4z_agg')
      !
   END SUBROUTINE p4z_agg

   !!======================================================================
END MODULE p4zagg
