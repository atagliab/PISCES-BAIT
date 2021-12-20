MODULE p4zsed
   !!======================================================================
   !!                         ***  MODULE p4sed  ***
   !! TOP :   PISCES Compute loss of biogenic matter in the sediments
   !!                Compute gain/loss of tracers from dust, rivers and 
   !!                sediments 
   !!                This module is used both by PISCES and PISCES-QUOTA
   !!======================================================================
   !! History :   1.0  !  2004-03 (O. Aumont) Original code
   !!             2.0  !  2007-12 (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06 (C. Ethe) USE of fldread
   !!             3.5  !  2012-07 (O. Aumont) improvment of river input of nutrients 
   !!-----------------------------------------------------------------------
   !!   p4z_sed        :  Compute loss of organic matter in the sediments
   !!                  :  Compute gain/loss of tracers from dust, rivers and 
   !!                     sediments 
   !!-----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE p4zlim          !  Co-limitations of differents nutrients
   USE p4zsbc          !  External source of nutrients 
   USE p4zint          !  interpolation and computation of various fields
   USE sed             !  Sediment module
   USE iom             !  I/O manager
   USE prtctl_trc      !  print control for debugging

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_sed  
   PUBLIC   p4z_sed_alloc
 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: nitrpot    !: Nitrogen fixation 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:  ) :: sdenit     !: Nitrate reduction in the sediments
   REAL(wp) :: r1_rday                  !: inverse of rday
   LOGICAL, SAVE :: lk_sed              !: Explicit sediment module

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zsed.F90 14276 2021-01-07 22:09:56Z aumont $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_sed( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_sed  ***
      !!
      !! ** Purpose : Compute the loss of biogenic matter in the sediments. This
      !!              is by no way a real sediment model. The loss is simply 
      !!              computed from metamodels.
      !!              Loss/gain of tracers are also computed here for 
      !!              dust, rivers, sediments and hydrothermal vents (Fe) 
      !!              N2 fixation is modeled using an implicit approach
      !!
      !! ** Method  : - Fluxes with the sediments are mainly modeled using
      !!                statiscal metamodels.
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt, knt ! ocean time step
      INTEGER  ::  ji, jj, jk, ikt
      REAL(wp) ::  zrivalk, zrivsil, zrivno3
      REAL(wp) ::  zwflux, zlim, zfact, zfactcal
      REAL(wp) ::  zo2, zno3, zflx, zpdenit, z1pdenit, zolimit
      REAL(wp) ::  zsiloss, zcaloss, zws3, zws4, zwsc, zdep
      REAL(wp) ::  zwstpoc, zwstpon, zwstpop
      REAL(wp) ::  ztrpo4s, ztrdp, zwdust, zmudia, ztemp
      REAL(wp) ::  xdiano3, xdianh4
      !
      CHARACTER (len=25) :: charout
      REAL(wp), DIMENSION(jpi,jpj    ) :: zdenit2d, zbureff, zwork
      REAL(wp), DIMENSION(jpi,jpj    ) :: zwsbio3, zwsbio4
      REAL(wp), DIMENSION(jpi,jpj    ) :: zwsbio5, zwsbio6,  zwsbio7, zwsbio8
      REAL(wp), DIMENSION(jpi,jpj    ) :: zsedcal, zsedsi, zsedc
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zsoufer, zlight, ztrfer
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: ztrpo4, ztrdop, zirondep, zpdep, zsidep
      REAL(wp), ALLOCATABLE, DIMENSION(:,:  ) :: zironice
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )  CALL timing_start('p4z_sed')
      !
      IF( kt == nittrc000 .AND. knt == 1 )   THEN
          r1_rday  = 1. / rday
          ! Configuration with an active two-way sediment module 
          IF (ln_sediment .AND. ln_sed_2way) THEN
             lk_sed = .TRUE.
          ELSE
             lk_sed = .FALSE.
          ENDIF
      ENDIF
      !
      IF( kt == nittrc000 .AND. knt == 1 )   r1_rday  = 1. / rday
      !
      ! Allocate temporary workspace
      ALLOCATE( ztrpo4(jpi,jpj,jpk) )
      IF( ln_p5z )    ALLOCATE( ztrdop(jpi,jpj,jpk) )

      zdenit2d(:,:) = 0.e0
      zbureff (:,:) = 0.e0
      zwork   (:,:) = 0.e0
      zsedsi  (:,:) = 0.e0
      zsedcal (:,:) = 0.e0
      zsedc   (:,:) = 0.e0

      ! Iron input/uptake due to sea ice : Crude parameterization based on 
      ! Lancelot et al. Iron concentration in sea-ice is constant and set 
      ! in the namelist_pisces (icefeinput). ln_ironice is forced to false
      ! when nn_ice_tr = 1
      ! ----------------------------------------------------
      IF( ln_ironice ) THEN  
         !                                              
         ALLOCATE( zironice(jpi,jpj) )

         ! Compute the iron flux between sea ice and sea water
         ! Simple parameterization assuming a fixed constant concentration in
         ! sea-ice (icefeinput)
         ! ------------------------------------------------------------------
         DO jj = 1, jpj
            DO ji = 1, jpi
               zdep    = rfact2 / e3t_n(ji,jj,1)
               zwflux  = fmmflx(ji,jj) / 1000._wp
               zironice(ji,jj) =  MAX( -0.99 * trb(ji,jj,1,jpfer), -zwflux * icefeinput * zdep )
            END DO
         END DO
         ! Update of the TRA array
         tra(:,:,1,jpfer) = tra(:,:,1,jpfer) + zironice(:,:) 
         ! 
         IF( lk_iomput .AND. knt == nrdttrc .AND. iom_use( "Ironice" ) )   &
            &   CALL iom_put( "Ironice", zironice(:,:) * 1.e+3 * rfact2r * e3t_n(:,:,1) * tmask(:,:,1) ) ! iron flux from ice
         !
         DEALLOCATE( zironice )
         !                                              
      ENDIF

      ! Add the external input of nutrients from dust deposition
      ! ----------------------------------------------------------
      IF( ln_dust ) THEN
         !                                              
         ALLOCATE( zsidep(jpi,jpj,jpk), zpdep(jpi,jpj,jpk), zirondep(jpi,jpj,jpk) )
         IF( ln_bait ) zirondep(:,:,:) = 0.
         ! Iron, P and Si deposition at the surface
         ! Iron flux at the surface due to dust deposition. Solubility can be 
         ! be variable if ln_solub is set to true. In that case, solubility 
         ! has to be provided in the specific input file (read in p4zsbc)
         ! mfrac is the mean iron relative weight content of dust
         ! ------------------------------------------------------------------
         IF( ln_bait ) THEN
         ! The bait model uses a combination of soluble Fe and dust inputs fromn
         ! CESM2CAM6 model in kg/m2/s
!          zirondep(:,:,1) = solfe(:,:) * rfact2 / e3t_n(:,:,1) / 55.85 + 3.e-10 * r1_ryyss 
          zirondep(:,:,1) = solfe(:,:) * dustsolub * rfact2 / e3t_n(:,:,1) / 55.85 + 3.e-10 * r1_ryyss
         ELSE
         IF( ln_solub ) THEN
            zirondep(:,:,1) = solub(:,:) * dust(:,:) * mfrac * rfact2 / e3t_n(:,:,1) / 55.85 + 3.e-10 * r1_ryyss 
         ELSE
            zirondep(:,:,1) = dustsolub  * dust(:,:) * mfrac * rfact2 / e3t_n(:,:,1) / 55.85 + 3.e-10 * r1_ryyss 
         ENDIF
         ENDIF ! ln_bait

         ! Si and P flux at the surface due to dust deposition. The content 
         ! and the solubility are hard coded
         ! ----------------------------------------------------------------
         ! Si crustal abundance is about 26.9% by mass and 7.5% is soluble
         ! (see for instance, Tegen and Kohfeld, 2006).
         zsidep(:,:,1) = 0.269 * 0.075 * dust(:,:) * rfact2 / e3t_n(:,:,1) / 28.1 
         ! P Crustal abundance is about 0.1% by mass and about 10% of it soluble
         ! (Paytan and McLaughlin, 2007).
         zpdep (:,:,1) = 0.1 * 1.e-3 * dust(:,:) * rfact2 / e3t_n(:,:,1) / 31. / po4r 

         IF( ln_bait ) THEN
         ! account for supply of lithogenic particulate Fe from dust input
         ! a specific fraction goes directly into fine lfe and the other already
         ! into the aggregated lfe
         ! We model directly the Fe content of the dust
         ! We can convert to particulate mass in g/L (where needed) by dividing by
         ! mfrac and multiplying by 55.85
         tra(:,:,1,jplfe) = tra(:,:,1,jplfe) + flfe * ( dust(:,:) * &
         &                  mfrac * rfact2 / e3t_n(:,:,1) / 55.85 + 3.e-10 * r1_ryyss )
         tra(:,:,1,jplfa) = tra(:,:,1,jplfa) + (1. - flfe) * ( dust(:,:) * &
         &                  mfrac * rfact2 / e3t_n(:,:,1) / 55.85 + 3.e-10 * r1_ryyss )
         ! account for solubilisation of lithogenic dust associated Fe
         DO jk = 1, jpkm1
            zirondep(:,:,jk) = zirondep(:,:,jk) + trb(:,:,jk,jplfe) * ( sollfe / tsollfe ) * r1_rday * rfact2
            zpdep   (:,:,jk) = zirondep(:,:,jk) * 1.e-3 / mfrac * 55.85 / 31.
            zsidep  (:,:,jk) = zirondep(:,:,jk) * 0.269 / mfrac * 55.85 / 28.1
            tra(:,:,jk,jplfe) = tra(:,:,jk,jplfe) - zirondep(:,:,jk)
         END DO
         ELSE ! no ln_bait
         ! Iron solubilization of particles in the water column
         ! dust in kg/m2/s ---> 1/55.85 to put in mol/Fe ;  wdust in m/d
         ! Dust are supposed to sink at wdust sinking speed. 3% of the iron 
         ! in dust is hypothesized to be soluble at a dissolution rate set to 
         ! 1/(250 days). The vertical distribution of iron in dust is computed 
         ! from a steady state assumption. Parameters are very uncertain and 
         ! are estimated from the literature quoted in Raiswell et al. (2011) 
         ! ------------------------------------------------------------------- 
         zwdust = 0.03 * rday / ( wdust * 55.85 ) / ( 250. * rday )
         DO jk = 2, jpkm1
            zirondep(:,:,jk) = dust(:,:) * mfrac * zwdust * rfact2 * EXP( -gdept_n(:,:,jk) / (250. * wdust) )
            zpdep   (:,:,jk) = zirondep(:,:,jk) * 1.e-3 / mfrac * 55.85 / 31.
            zsidep  (:,:,jk) = zirondep(:,:,jk) * 0.269 / mfrac * 55.85 / 28.1
         END DO
         ENDIF ! ln_bait

         ! Solubilization of particles in the water column (Si, P, Fe)
         DO jk = 1, jpkm1
            tra(:,:,jk,jppo4) = tra(:,:,jk,jppo4) + zpdep   (:,:,jk)
            tra(:,:,jk,jpfer) = tra(:,:,jk,jpfer) + zirondep(:,:,jk) 
            tra(:,:,jk,jpsil) = tra(:,:,jk,jpsil) + zsidep  (:,:,jk)
         ENDDO

         DO jk = 1, jpkm1
            tra(:,:,jk,jplgw) = tra(:,:,jk,jplgw) + zirondep(:,:,jk ) * lgw_ratd
         ENDDO

         ! 
         IF( lk_iomput ) THEN
            IF( knt == nrdttrc ) THEN
                IF( iom_use( "Irondep" ) )   &
                &  CALL iom_put( "Irondep", zirondep(:,:,1) * 1.e+3 * rfact2r * e3t_n(:,:,1) * tmask(:,:,1) ) ! surface downward dust depo of iron
                IF( iom_use( "pdust" ) )   &
                &  CALL iom_put( "pdust"  , dust(:,:) / ( wdust / rday )  * tmask(:,:,1) ) ! dust concentration at surface
            ENDIF
         ENDIF
         DEALLOCATE( zsidep, zpdep, zirondep )
         !                                              
      ENDIF
     
      ! Add the external input of nutrients from river
      ! ----------------------------------------------
      IF( ln_river ) THEN
         DO jj = 1, jpj
            DO ji = 1, jpi
               DO jk = 1, nk_rnf(ji,jj)
                  tra(ji,jj,jk,jppo4) = tra(ji,jj,jk,jppo4) +  rivdip(ji,jj) * rfact2
                  tra(ji,jj,jk,jpno3) = tra(ji,jj,jk,jpno3) +  rivdin(ji,jj) * rfact2
                  tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) +  rivdic(ji,jj) * 5.e-5 * rfact2
                  tra(ji,jj,jk,jpsil) = tra(ji,jj,jk,jpsil) +  rivdsi(ji,jj) * rfact2
                  tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) +  rivdic(ji,jj) * rfact2
                  tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) +  ( rivalk(ji,jj) - rno3 * rivdin(ji,jj) ) * rfact2
                  tra(ji,jj,jk,jpdoc) = tra(ji,jj,jk,jpdoc) +  rivdoc(ji,jj) * rfact2
               ENDDO
            ENDDO
         ENDDO

         ! When prognostic ligands are activated, ligands are supplied 
         ! to the ocean by rivers. We assume that the amount of ligands
         ! is equal to that of iron (iron is completely complexed)
         ! ------------------------------------------------------------
         IF (ln_ligand) THEN
            DO jj = 1, jpj
               DO ji = 1, jpi
                  DO jk = 1, nk_rnf(ji,jj)
                     tra(ji,jj,jk,jplgw) = tra(ji,jj,jk,jplgw) +  rivdic(ji,jj) * 5.e-5 * rfact2
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
         ! PISCES-QUOTA part
         IF( ln_p5z ) THEN
            DO jj = 1, jpj
               DO ji = 1, jpi
                  DO jk = 1, nk_rnf(ji,jj)
                     tra(ji,jj,jk,jpdop) = tra(ji,jj,jk,jpdop) + rivdop(ji,jj) * rfact2
                     tra(ji,jj,jk,jpdon) = tra(ji,jj,jk,jpdon) + rivdon(ji,jj) * rfact2
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
      ENDIF
      
      ! Add the external input of nutrients from nitrogen deposition
      ! ----------------------------------------------------------
      IF( ln_ndepo ) THEN
         tra(:,:,1,jpno3) = tra(:,:,1,jpno3) + nitdep(:,:) * rfact2
         tra(:,:,1,jptal) = tra(:,:,1,jptal) - rno3 * nitdep(:,:) * rfact2
      ENDIF

      ! Add the external input of iron from hydrothermal vents
      ! Please refer to Tagliabue et al. (2010) for more information
      ! ------------------------------------------------------------
      IF( ln_hydrofe ) THEN
            tra(:,:,:,jpfer) = tra(:,:,:,jpfer) + hydrofe(:,:,:) * rfact2
         IF( ln_ligand ) THEN
            tra(:,:,:,jplgw) = tra(:,:,:,jplgw) + ( hydrofe(:,:,:) * lgw_rath ) * rfact2
         ENDIF
         !
         IF( lk_iomput .AND. knt == nrdttrc .AND. iom_use( "HYDR" ) )   &
            &   CALL iom_put( "HYDR", hydrofe(:,:,:) * 1.e+3 * tmask(:,:,:) ) ! hydrothermal iron input
      ENDIF

      ! OA: Warning, the following part is necessary to avoid CFL problems 
      ! above the sediments. Vertical sinking speed is limited using the 
      ! typical CFL criterion
      ! --------------------------------------------------------------------
      DO jj = 1, jpj
         DO ji = 1, jpi
            ikt  = mbkt(ji,jj)
            zdep = e3t_n(ji,jj,ikt) / xstep
            zwsbio4(ji,jj) = MIN( 0.99 * zdep, wsbio4(ji,jj,ikt) )
            zwsbio3(ji,jj) = MIN( 0.99 * zdep, wsbio3(ji,jj,ikt) )
         END DO
      END DO

         IF( ln_bait ) THEN
      DO jj = 1, jpj
         DO ji = 1, jpi
            ikt  = mbkt(ji,jj)
            zdep = e3t_n(ji,jj,ikt) / xstep
            zwsbio6(ji,jj) = MIN( 0.99 * zdep, wsbio6(ji,jj,ikt) )
            zwsbio5(ji,jj) = MIN( 0.99 * zdep, wsbio5(ji,jj,ikt) )
         END DO
      END DO
          ENDIF

         IF( ln_feauth ) THEN
      DO jj = 1, jpj
         DO ji = 1, jpi
            ikt  = mbkt(ji,jj)
            zdep = e3t_n(ji,jj,ikt) / xstep
            zwsbio8(ji,jj) = MIN( 0.99 * zdep, wsbio8(ji,jj,ikt) )
            zwsbio7(ji,jj) = MIN( 0.99 * zdep, wsbio7(ji,jj,ikt) )
         END DO
      END DO
          ENDIF
      !
      ! No sediment module activated
      IF( .NOT.lk_sed ) THEN
!
         ! Add the external input of iron from sediment mobilization
         ! ------------------------------------------------------
         IF( ln_ironsed ) THEN
            tra(:,:,:,jpfer) = tra(:,:,:,jpfer) + ironsed(:,:,:) * rfact2
            !
            IF( lk_iomput .AND. knt == nrdttrc .AND. iom_use( "Ironsed" ) )   &
               &   CALL iom_put( "Ironsed", ironsed(:,:,:) * 1.e+3 * tmask(:,:,:) ) ! iron inputs from sediments
         ENDIF

         ! Computation of the sediment denitrification proportion: The metamodel 
         ! from Middleburg (2006) is used
         ! Computation of the fraction of organic matter that is permanently 
         ! buried from Dunne's model (2007)
         ! -------------------------------------------------------
         DO jj = 1, jpj
            DO ji = 1, jpi
              IF( tmask(ji,jj,1) == 1 ) THEN
                 ikt = mbkt(ji,jj)
                 zflx = (  trb(ji,jj,ikt,jpgoc) * zwsbio4(ji,jj)   &
                   &     + trb(ji,jj,ikt,jppoc) * zwsbio3(ji,jj) )  * 1E3 * 1E6 / 1E4
                 zflx  = LOG10( MAX( 1E-3, zflx ) )
                 zo2   = LOG10( MAX( 10. , trb(ji,jj,ikt,jpoxy) * 1E6 ) )
                 zno3  = LOG10( MAX( 1.  , trb(ji,jj,ikt,jpno3) * 1E6 * rno3 ) )
                 zdep  = LOG10( gdepw_n(ji,jj,ikt+1) )
                 zdenit2d(ji,jj) = -2.2567 - 1.185 * zflx - 0.221 * zflx**2 - 0.3995 * zno3 * zo2 + 1.25 * zno3    &
                   &                + 0.4721 * zo2 - 0.0996 * zdep + 0.4256 * zflx * zo2
                 zdenit2d(ji,jj) = 10.0**( zdenit2d(ji,jj) )
                   !
                 zflx = (  trb(ji,jj,ikt,jpgoc) * zwsbio4(ji,jj)   &
                   &     + trb(ji,jj,ikt,jppoc) * zwsbio3(ji,jj) ) * 1E6
                 zbureff(ji,jj) = 0.013 + 0.13 * zflx**2 / ( 7.0 + zflx )**2
              ENDIF
            END DO
         END DO 
         !
      ENDIF

      ! Fraction of dSi that is dissolved in the sediments. This fraction is  
      ! set to a constant value in p4zsbc
      ! --------------------------------------------------------------------
      IF( .NOT.lk_sed )  zrivsil = 1._wp - sedsilfrac

      ! Loss of bSi and CaCO3 to the sediments
      DO jj = 1, jpj
         DO ji = 1, jpi
            ikt  = mbkt(ji,jj)
            zdep = xstep / e3t_n(ji,jj,ikt) 
            zwsc = zwsbio4(ji,jj) * zdep
            zsiloss = trb(ji,jj,ikt,jpgsi) * zwsc
            zcaloss = trb(ji,jj,ikt,jpcal) * zwsc
            !
            tra(ji,jj,ikt,jpgsi) = tra(ji,jj,ikt,jpgsi) - zsiloss
            tra(ji,jj,ikt,jpcal) = tra(ji,jj,ikt,jpcal) - zcaloss
         END DO
      END DO
      !
      IF( .NOT.lk_sed ) THEN
         ! Dissolution of CaCO3 and bSi in the sediments. This is 
         ! instantaneous since here sediments are not explicitly 
         ! modeled. The amount of CaCO3 that dissolves in the sediments
         ! is computed using a metamodel constructed from Archer (1996)
         ! A minimum set to sedcalfrac is preserved. This value is defined
         ! in p4zsbc
         ! ---------------------------------------------------------------
         DO jj = 1, jpj
            DO ji = 1, jpi
               ikt  = mbkt(ji,jj)
               zdep = xstep / e3t_n(ji,jj,ikt) 
               zwsc = zwsbio4(ji,jj) * zdep
               zsiloss = trb(ji,jj,ikt,jpgsi) * zwsc
               zcaloss = trb(ji,jj,ikt,jpcal) * zwsc
               tra(ji,jj,ikt,jpsil) = tra(ji,jj,ikt,jpsil) + zsiloss * zrivsil 
               !
               zfactcal = MAX(-0.1, MIN( excess(ji,jj,ikt), 0.2 ) )
               zfactcal = 0.3 + 0.7 * MIN( 1., (0.1 + zfactcal) / ( 0.5 - zfactcal ) )
               zrivalk  = sedcalfrac * zfactcal
               tra(ji,jj,ikt,jptal) =  tra(ji,jj,ikt,jptal) + zcaloss * zrivalk * 2.0
               tra(ji,jj,ikt,jpdic) =  tra(ji,jj,ikt,jpdic) + zcaloss * zrivalk
               zsedcal(ji,jj) = (1.0 - zrivalk) * zcaloss * e3t_n(ji,jj,ikt) 
               zsedsi (ji,jj) = (1.0 - zrivsil) * zsiloss * e3t_n(ji,jj,ikt) 
            END DO
         END DO
      ENDIF
      !
      ! Loss of particulate organic carbon and Fe to the sediments

       DO jj = 1, jpj
         DO ji = 1, jpi
            ikt  = mbkt(ji,jj)
            zdep = xstep / e3t_n(ji,jj,ikt) 
            zws4 = zwsbio4(ji,jj) * zdep
            zws3 = zwsbio3(ji,jj) * zdep
            tra(ji,jj,ikt,jpgoc) = tra(ji,jj,ikt,jpgoc) - trb(ji,jj,ikt,jpgoc) * zws4 
            tra(ji,jj,ikt,jppoc) = tra(ji,jj,ikt,jppoc) - trb(ji,jj,ikt,jppoc) * zws3
            tra(ji,jj,ikt,jpbfe) = tra(ji,jj,ikt,jpbfe) - trb(ji,jj,ikt,jpbfe) * zws4
            tra(ji,jj,ikt,jpsfe) = tra(ji,jj,ikt,jpsfe) - trb(ji,jj,ikt,jpsfe) * zws3
         END DO
      END DO
      !
      ! Loss of particulate organic N and P to the sediments (p5z)
      IF( ln_p5z ) THEN
         DO jj = 1, jpj
            DO ji = 1, jpi
               ikt  = mbkt(ji,jj)
               zdep = xstep / e3t_n(ji,jj,ikt) 
               zws4 = zwsbio4(ji,jj) * zdep
               zws3 = zwsbio3(ji,jj) * zdep
               tra(ji,jj,ikt,jpgon) = tra(ji,jj,ikt,jpgon) - trb(ji,jj,ikt,jpgon) * zws4
               tra(ji,jj,ikt,jppon) = tra(ji,jj,ikt,jppon) - trb(ji,jj,ikt,jppon) * zws3
               tra(ji,jj,ikt,jpgop) = tra(ji,jj,ikt,jpgop) - trb(ji,jj,ikt,jpgop) * zws4
               tra(ji,jj,ikt,jppop) = tra(ji,jj,ikt,jppop) - trb(ji,jj,ikt,jppop) * zws3
            END DO
         END DO
      ENDIF

      IF ( ln_bait ) THEN
      DO jj = 1, jpj
         DO ji = 1, jpi
            ikt  = mbkt(ji,jj)
            zdep = xstep / e3t_n(ji,jj,ikt)
            zws4 = zwsbio6(ji,jj) * zdep
            zws3 = zwsbio5(ji,jj) * zdep
            tra(ji,jj,ikt,jplfe) = tra(ji,jj,ikt,jplfe) - trb(ji,jj,ikt,jplfe) * zws3
            tra(ji,jj,ikt,jplfa) = tra(ji,jj,ikt,jplfa) - trb(ji,jj,ikt,jplfa) * zws4
      IF ( ln_feauth ) THEN
            zws4 = zwsbio8(ji,jj) * zdep
            zws3 = zwsbio7(ji,jj) * zdep
            tra(ji,jj,ikt,jpafs) = tra(ji,jj,ikt,jpafs) - trb(ji,jj,ikt,jpafs) *zws3
            tra(ji,jj,ikt,jpafb) = tra(ji,jj,ikt,jpafb) - trb(ji,jj,ikt,jpafb) * zws4
      ENDIF
         END DO
      END DO
      ENDIF

      IF( .NOT.lk_sed ) THEN
         ! Degradation of organic matter in the sediments. The metamodel of 
         ! Middleburg (2006) is used here to mimic the diagenetic reactions. 
         ! The 0.5 factor in zpdenit is to avoid negative NO3 concentration after
         ! denitrification in the sediments. Not very clever, but simpliest option.
         ! ------------------------------------------------------------------------
         DO jj = 1, jpj
            DO ji = 1, jpi
               ikt  = mbkt(ji,jj)
               zdep = xstep / e3t_n(ji,jj,ikt) 
               zws4 = zwsbio4(ji,jj) * zdep
               zws3 = zwsbio3(ji,jj) * zdep
               ! Fraction that is permanently buried in the sediments
               zrivno3 = 1. - zbureff(ji,jj)
               zwstpoc = trb(ji,jj,ikt,jpgoc) * zws4 + trb(ji,jj,ikt,jppoc) * zws3
               ! Denitrification in the sediments
               zpdenit  = MIN( 0.5 * ( trb(ji,jj,ikt,jpno3) - rtrn ) / rdenit, zdenit2d(ji,jj) * zwstpoc * zrivno3 )
               ! Fraction that is not denitrified
               z1pdenit = zwstpoc * zrivno3 - zpdenit
               ! Oxic remineralization of organic matter in the sediments
               zolimit = MIN( ( trb(ji,jj,ikt,jpoxy) - rtrn ) / o2ut, z1pdenit * ( 1.- nitrfac(ji,jj,ikt) ) )
               ! The fraction that cannot be denitrified nor oxidized by O2
               ! is released back to the water column as DOC
               tra(ji,jj,ikt,jpdoc) = tra(ji,jj,ikt,jpdoc) + z1pdenit - zolimit
               ! Update of the tracers concentrations
               tra(ji,jj,ikt,jppo4) = tra(ji,jj,ikt,jppo4) + zpdenit + zolimit
               tra(ji,jj,ikt,jpnh4) = tra(ji,jj,ikt,jpnh4) + zpdenit + zolimit
               tra(ji,jj,ikt,jpno3) = tra(ji,jj,ikt,jpno3) - rdenit * zpdenit
               tra(ji,jj,ikt,jpoxy) = tra(ji,jj,ikt,jpoxy) - zolimit * o2ut
               tra(ji,jj,ikt,jptal) = tra(ji,jj,ikt,jptal) + rno3 * (zolimit + (1.+rdenit) * zpdenit )
               tra(ji,jj,ikt,jpdic) = tra(ji,jj,ikt,jpdic) + zpdenit + zolimit 
               sdenit(ji,jj) = rdenit * zpdenit * e3t_n(ji,jj,ikt)
               zsedc(ji,jj)   = (1. - zrivno3) * zwstpoc * e3t_n(ji,jj,ikt)
               ! PISCES-QUOTA (p5z)
               IF( ln_p5z ) THEN
                  zwstpop              = trb(ji,jj,ikt,jpgop) * zws4 + trb(ji,jj,ikt,jppop) * zws3
                  zwstpon              = trb(ji,jj,ikt,jpgon) * zws4 + trb(ji,jj,ikt,jppon) * zws3
                  tra(ji,jj,ikt,jpdon) = tra(ji,jj,ikt,jpdon) + ( z1pdenit - zolimit ) * zwstpon / (zwstpoc + rtrn)
                  tra(ji,jj,ikt,jpdop) = tra(ji,jj,ikt,jpdop) + ( z1pdenit - zolimit ) * zwstpop / (zwstpoc + rtrn)
               ENDIF
            END DO
         END DO
       ENDIF


      ! Nitrogen fixation process : light limitation of diazotrophy
      ! Small source of iron from particulate inorganic iron (photochemistry)
      ! This is a purely adhoc param.
      !----------------------------------------------------------------------
      DO jk = 1, jpkm1
         zlight (:,:,jk) =  ( 1.- EXP( -etot_ndcy(:,:,jk) / diazolight ) ) * ( 1. - fr_i(:,:) ) 
         zsoufer(:,:,jk) = zlight(:,:,jk) * 1E-10 / ( 1E-10 + biron(:,:,jk) )
      ENDDO

      ! Diazotrophy (nitrogen fixation) is modeled according to an empirical
      ! formulation. This is described in Aumont et al. (2015). Limitation 
      ! by P and Fe is computed. Inhibition by high N concentrations is imposed.
      ! Diazotrophy sensitivity to temperature is parameterized as in 
      ! Ye et al. (2012)  
      ! ------------------------------------------------------------------------
      IF( ln_p4z ) THEN
         ! PISCES part
         DO jk = 1, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  ! Potential nitrogen fixation dependant on temperature
                  ztemp = tsn(ji,jj,jk,jp_tem)
                  zmudia = MAX( 0.,-0.001096*ztemp**2 + 0.057*ztemp -0.637 ) / rno3
                  ! Nitrogen fixation is inhibited when enough NO3 and/or NH4
                  zlim = ( 1.- xnanonh4(ji,jj,jk) - xnanono3(ji,jj,jk) )
                  zfact = zlim * rfact2
                  ! Nitrogen fixation limitation by PO4 and Fe
                  ztrfer(ji,jj,jk) = biron(ji,jj,jk) / ( concfediaz + biron(ji,jj,jk) )
                  ztrpo4(ji,jj,jk) = trb(ji,jj,jk,jppo4) / ( concnnh4 + trb(ji,jj,jk,jppo4) )
                  ztrdp = ztrpo4(ji,jj,jk)
                  nitrpot(ji,jj,jk) =  zmudia * r1_rday * zfact * MIN( ztrfer(ji,jj,jk), ztrdp ) * zlight(ji,jj,jk)
               END DO
            END DO
         END DO
      ELSE       ! p5z
         ! PISCES-QUOTA part
         DO jk = 1, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  ! Potential nitrogen fixation dependant on temperature
                  ztemp = tsn(ji,jj,jk,jp_tem)
                  zmudia = MAX( 0.,-0.001096*ztemp**2 + 0.057*ztemp -0.637 ) * 7.625
                  ! Nitrogen fixation is inhibited when enough NO3 and/or NH4
                  xdianh4 = trb(ji,jj,jk,jpnh4) / ( concnnh4 + trb(ji,jj,jk,jpnh4) )
                  xdiano3 = trb(ji,jj,jk,jpno3) / ( concnno3 + trb(ji,jj,jk,jpno3) ) * (1. - xdianh4)
                  zlim = ( 1.- xdiano3 - xdianh4 )
                  zfact = zlim * rfact2
                  ! Nitrogen fixation limitation by PO4/DOP and Fe
                  ztrfer(ji,jj,jk) = biron(ji,jj,jk) / ( concfediaz + biron(ji,jj,jk) )
                  ztrpo4(ji,jj,jk) = trb(ji,jj,jk,jppo4) / ( 1E-6 + trb(ji,jj,jk,jppo4) )
                  ztrdop(ji,jj,jk) = trb(ji,jj,jk,jpdop) / ( 1E-6 + trb(ji,jj,jk,jpdop) ) * (1. - ztrpo4(ji,jj,jk))
                  ztrdp = ztrpo4(ji,jj,jk) + ztrdop(ji,jj,jk)
                  nitrpot(ji,jj,jk) =  zmudia * r1_rday * zfact * MIN( ztrfer(ji,jj,jk), ztrdp ) * zlight(ji,jj,jk)
               END DO
            END DO
         END DO
      ENDIF

      ! Update of the TRA arrays due to nitrogen fixation
      ! -------------------------------------------------
      IF( ln_p4z ) THEN
         ! PISCES part
         DO jk = 1, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zfact = nitrpot(ji,jj,jk) * nitrfix
                  ! 1/3 of the diazotrophs growth is supposed to be excreted
                  ! as NH4. 1/3 as DOC and the rest is routed to POC/GOC as 
                  ! a result of mortality by predation. Completely adhoc param 
                  tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) + zfact / 3.0
                  tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) + rno3 * zfact / 3.0
                  tra(ji,jj,jk,jppo4) = tra(ji,jj,jk,jppo4) - zfact * 2.0 / 3.0
                  tra(ji,jj,jk,jpdoc) = tra(ji,jj,jk,jpdoc) + zfact * 1.0 / 3.0
                  tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) + zfact * 1.0 / 3.0 * 2.0 / 3.0
                  tra(ji,jj,jk,jpgoc) = tra(ji,jj,jk,jpgoc) + zfact * 1.0 / 3.0 * 1.0 / 3.0
                  tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy) + ( o2ut + o2nit ) * zfact * 2.0 / 3.0 + o2nit * zfact / 3.0
                  ! Fe/c of diazotrophs is assumed to be 30umol Fe/mol C at max
                  tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) - 30E-6 * ztrfer(ji,jj,jk) * zfact * 1.0 / 3.0
                  tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) + 30E-6 * ztrfer(ji,jj,jk) * zfact * 1.0 / 3.0 * 2.0 / 3.0
                  tra(ji,jj,jk,jpbfe) = tra(ji,jj,jk,jpbfe) + 30E-6 * ztrfer(ji,jj,jk) * zfact * 1.0 / 3.0 * 1.0 / 3.0
                  tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) + 0.005 * 4E-10 * zsoufer(ji,jj,jk) * rfact2 / rday
                  tra(ji,jj,jk,jppo4) = tra(ji,jj,jk,jppo4) + concdnh4 / ( concdnh4 + trb(ji,jj,jk,jppo4) ) &
                  &                     * 0.001 * trb(ji,jj,jk,jpdoc) * xstep
              END DO
            END DO 
         END DO
      ELSE    ! p5z
         ! PISCES-QUOTA part
         DO jk = 1, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zfact = nitrpot(ji,jj,jk) * nitrfix
                  ! 1/3 of the diazotrophs growth is supposed to be excreted
                  ! as NH4. 1/3 as DOC and the rest is routed POC and GOC as 
                  ! a result of mortality by predation. Completely adhoc param 
                  tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) + zfact / 3.0
                  tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) + rno3 * zfact / 3.0
                  ! N/P ratio of diazotrophs is supposed to be 46
                  tra(ji,jj,jk,jppo4) = tra(ji,jj,jk,jppo4) - 16.0 / 46.0 * zfact * ( 1.0 - 1.0 / 3.0 ) &
                  &                     * ztrpo4(ji,jj,jk) / (ztrpo4(ji,jj,jk) + ztrdop(ji,jj,jk) + rtrn)
                  tra(ji,jj,jk,jpdon) = tra(ji,jj,jk,jpdon) + zfact * 1.0 / 3.0
                  tra(ji,jj,jk,jpdoc) = tra(ji,jj,jk,jpdoc) + zfact * 1.0 / 3.0
                  tra(ji,jj,jk,jpdop) = tra(ji,jj,jk,jpdop) + 16.0 / 46.0 * zfact / 3.0  &
                  &                     - 16.0 / 46.0 * zfact * ztrdop(ji,jj,jk)   &
                  &                     / (ztrpo4(ji,jj,jk) + ztrdop(ji,jj,jk) + rtrn)
                  tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) + zfact * 1.0 / 3.0 * 2.0 / 3.0
                  tra(ji,jj,jk,jppon) = tra(ji,jj,jk,jppon) + zfact * 1.0 / 3.0 * 2.0 /3.0
                  tra(ji,jj,jk,jppop) = tra(ji,jj,jk,jppop) + 16.0 / 46.0 * zfact * 1.0 / 3.0 * 2.0 /3.0
                  tra(ji,jj,jk,jpgoc) = tra(ji,jj,jk,jpgoc) + zfact * 1.0 / 3.0 * 1.0 / 3.0
                  tra(ji,jj,jk,jpgon) = tra(ji,jj,jk,jpgon) + zfact * 1.0 / 3.0 * 1.0 /3.0
                  tra(ji,jj,jk,jpgop) = tra(ji,jj,jk,jpgop) + 16.0 / 46.0 * zfact * 1.0 / 3.0 * 1.0 /3.0
                  tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy) + ( o2ut + o2nit ) * zfact * 2.0 / 3.0 + o2nit * zfact / 3.0
                  ! Fe/c of diazotrophs is assumed to be 30umol Fe/mol C at max
                  tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) - 30E-6 * ztrfer(ji,jj,jk) * zfact * 1.0 / 3.0 
                  tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) + 30E-6 * ztrfer(ji,jj,jk) * zfact * 1.0 / 3.0 * 2.0 / 3.0
                  tra(ji,jj,jk,jpbfe) = tra(ji,jj,jk,jpbfe) + 30E-6 * ztrfer(ji,jj,jk) * zfact * 1.0 / 3.0 * 1.0 / 3.0
                  tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) + 0.002 * 4E-10 * zsoufer(ji,jj,jk) * rfact2 / rday
              END DO
            END DO 
         END DO
         !
      ENDIF

      IF( lk_iomput ) THEN
         IF( knt == nrdttrc ) THEN
            zfact = 1.e+3 * rfact2r !  conversion from molC/l/kt  to molN/m3/s
            IF( iom_use("Nfix"   ) ) CALL iom_put( "Nfix", nitrpot(:,:,:) * nitrfix * rno3 * zfact * tmask(:,:,:) )  ! nitrogen fixation 
            IF( iom_use("INTNFIX") ) THEN   ! nitrogen fixation rate in ocean ( vertically integrated )
               zwork(:,:) = 0.
               DO jk = 1, jpkm1
                 zwork(:,:) = zwork(:,:) + nitrpot(:,:,jk) * nitrfix * rno3 * zfact * e3t_n(:,:,jk) * tmask(:,:,jk)
               ENDDO
               CALL iom_put( "INTNFIX" , zwork ) 
            ENDIF
            IF( iom_use("SedCal" ) ) CALL iom_put( "SedCal", zsedcal(:,:) * zfact ) ! Permanent burial of CaCO3 in sediments
            IF( iom_use("SedSi" ) )  CALL iom_put( "SedSi",  zsedsi (:,:) * zfact ) ! Permanent burial of bSi in sediments
            IF( iom_use("SedC" ) )   CALL iom_put( "SedC",   zsedc  (:,:) * zfact ) ! Permanent burial of OC in sediments
!            IF( iom_use("Sdenit" ) ) CALL iom_put( "Sdenit", sdenit (:,:) * zfact * rno3 ) ! Denitrification in the sediments
            IF( iom_use("Sdenit" ) ) CALL iom_put( "Sdenit", zdenit2d(:,:) ) ! Denitrification in the sediments
         ENDIF
      ENDIF
      !
      IF(ln_ctl) THEN  ! print mean trends (USEd for debugging)
         WRITE(charout, fmt="('sed ')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( ln_p5z )    DEALLOCATE( ztrpo4, ztrdop )
      !
      IF( ln_timing )  CALL timing_stop('p4z_sed')
      !
   END SUBROUTINE p4z_sed


   INTEGER FUNCTION p4z_sed_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_sed_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( nitrpot(jpi,jpj,jpk), sdenit(jpi,jpj), STAT=p4z_sed_alloc )
      !
      IF( p4z_sed_alloc /= 0 )   CALL ctl_stop( 'STOP', 'p4z_sed_alloc: failed to allocate arrays' )
      !
   END FUNCTION p4z_sed_alloc

   !!======================================================================
END MODULE p4zsed
