MODULE p4zprod
   !!======================================================================
   !!                         ***  MODULE p4zprod  ***
   !! TOP :  Growth Rate of the two phytoplankton groups of PISCES 
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-05  (O. Aumont, C. Ethe) New parameterization of light limitation
   !!----------------------------------------------------------------------
   !!   p4z_prod       : Compute the growth Rate of the two phytoplanktons groups
   !!   p4z_prod_init  : Initialization of the parameters for growth
   !!   p4z_prod_alloc : Allocate variables for growth
   !!----------------------------------------------------------------------
   USE oce_trc         ! shared variables between ocean and passive tracers
   USE trc             ! passive tracers common variables 
   USE sms_pisces      ! PISCES Source Minus Sink variables
   USE p4zlim          ! Co-limitations of differents nutrients
   USE prtctl_trc      ! print control for debugging
   USE iom             ! I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_prod         ! called in p4zbio.F90
   PUBLIC   p4z_prod_init    ! called in trcsms_pisces.F90
   PUBLIC   p4z_prod_alloc   ! called in trcini_pisces.F90

   REAL(wp), PUBLIC ::   pislopen     !:  P-I slope of nanophytoplankton
   REAL(wp), PUBLIC ::   pisloped     !:  P-I slope of diatoms
   REAL(wp), PUBLIC ::   xadap        !:  Adaptation factor to low light 
   REAL(wp), PUBLIC ::   excretn      !:  Excretion ratio of nanophyto
   REAL(wp), PUBLIC ::   excretd      !:  Excretion ratio of diatoms
   REAL(wp), PUBLIC ::   bresp        !:  Basal respiration rate
   REAL(wp), PUBLIC ::   chlcnm       !:  Maximum Chl/C ratio of nano
   REAL(wp), PUBLIC ::   chlcdm       !:  Maximum Chl/C ratio of diatoms
   REAL(wp), PUBLIC ::   chlcmin      !:  Minimum Chl/C ratio of phytoplankton
   REAL(wp), PUBLIC ::   fecnm        !:  Maximum Fe/C ratio of nano
   REAL(wp), PUBLIC ::   fecdm        !:  Maximum Fe/C ratio of diatoms
   REAL(wp), PUBLIC ::   grosip       !:  Mean Si/C ratio of diatoms

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   quotan   !: proxy of N quota in Nanophyto
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   quotad   !: proxy of N quota in diatoms
   
   REAL(wp) ::   r1_rday    ! 1 / rday
   REAL(wp) ::   texcretn   ! 1 - excretn 
   REAL(wp) ::   texcretd   ! 1 - excretd        

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zprod.F90 14276 2021-01-07 22:09:56Z aumont $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_prod( kt , knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_prod  ***
      !!
      !! ** Purpose :   Computes phytoplankton production depending on
      !!                light, temperature and nutrient availability
      !!                Computes also the uptake of Iron and Si as well 
      !!                as the chlorophyll content of the cells
      !!                PISCES relies on a mixed Monod-Quota formalism 
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, knt   !
      !
      INTEGER  ::   ji, jj, jk
      REAL(wp) ::   zsilfac, znanotot, zdiattot, zconctemp, zconctemp2
      REAL(wp) ::   zratio, zmax, zsilim, ztn, zadap, zlim, zsiborn
      REAL(wp) ::   zprod, zproreg, zproreg2, zprochln, zprochld
      REAL(wp) ::   zdocprod, zpislopen, zpisloped, zfact
      REAL(wp) ::   zratiosi, zmaxsi, zlimfac, zsizetmp, zfecnm, zfecdm
      REAL(wp) ::   zrum, zcodel, zargu, zval, zfeup
      CHARACTER (len=25) :: charout
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: zw2d
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zw3d
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zprmaxn,zprmaxd
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zpislopeadn, zpislopeadd, zysopt  
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zprdia, zprbio, zprchld, zprchln   
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zprorcan, zprorcad, zprofed, zprofen
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zpronewn, zpronewd
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zmxl_fac, zmxl_chl
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zpligprod1, zpligprod2
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('p4z_prod')
      !
      !  Allocate temporary workspace
      !
      zprorcan(:,:,:) = 0._wp ; zprorcad(:,:,:) = 0._wp ; zprofed (:,:,:) = 0._wp
      zprofen (:,:,:) = 0._wp ; zysopt  (:,:,:) = 0._wp
      zpronewn(:,:,:) = 0._wp ; zpronewd(:,:,:) = 0._wp ; zprdia  (:,:,:) = 0._wp
      zprbio  (:,:,:) = 0._wp ; zprchld (:,:,:) = 0._wp ; zprchln (:,:,:) = 0._wp 
      zmxl_fac(:,:,:) = 0._wp ; zmxl_chl(:,:,:) = 0._wp 
      consfe3 (:,:,:) = 0._wp

      ! Computation of the maximimum production. Based on a Q10 description
      ! of the thermal dependency
      ! Parameters are taken from Bissinger et al. (2008)
      zprmaxn(:,:,:) = 0.8_wp * r1_rday * tgfunc(:,:,:)
      zprmaxd(:,:,:) = zprmaxn(:,:,:)

      ! Impact of the day duration and light intermittency on phytoplankton growth
      ! Intermittency is supposed to have a similar effect on production as 
      ! day length (Shatwell et al., 2012). The correcting factor is zmxl_fac. 
      ! zmxl_chl is the fractional day length and is used to compute the mean
      ! PAR during daytime. The effect of mixing is computed using the 
      ! absolute light level definition of the euphotic zone
      ! ------------------------------------------------------------------------- 
      DO jk = 1, jpkm1
         DO jj = 1 ,jpj
            DO ji = 1, jpi
               IF( etot_ndcy(ji,jj,jk) > 1.E-3 ) THEN
                  zval = MAX( 1., strn(ji,jj) )
                  IF( gdepw_n(ji,jj,jk+1) <= hmld(ji,jj) ) THEN
                     zval = zval * MIN(1., heup_01(ji,jj) / ( hmld(ji,jj) + rtrn ))
                  ENDIF
                  zmxl_chl(ji,jj,jk) = zval / 24.
                  zmxl_fac(ji,jj,jk) = 1.0 - exp( -0.26 * zval )
               ENDIF
            END DO
         END DO
      END DO

      zprbio(:,:,:) = zprmaxn(:,:,:) * zmxl_fac(:,:,:)
      zprdia(:,:,:) = zprmaxd(:,:,:) * zmxl_fac(:,:,:)

      ! Computation of the P-I slope for nanos and diatoms
      ! The formulation proposed by Geider et al. (1997) has been modified 
      ! to exclude the effect of nutrient limitation and temperature in the PI
      ! curve following Vichi et al. (2007)
      ! -----------------------------------------------------------------------
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( etot_ndcy(ji,jj,jk) > 1.E-3 ) THEN
                  ztn         = MAX( 0., tsn(ji,jj,jk,jp_tem) - 15. )
                  zadap       = xadap * ztn / ( 2.+ ztn )
                  zconctemp   = MAX( 0.e0 , trb(ji,jj,jk,jpdia) - xsizedia )
                  zconctemp2  = trb(ji,jj,jk,jpdia) - zconctemp

                  ! The initial slope of the PI curve can be increased for nano
                  ! to account for photadaptation, for instance in the DCM
                  ! This parameterization is adhoc and should be either 
                  ! improved or removed in future versions of the model

                  ! Nanophytoplankton
                  zpislopeadn(ji,jj,jk) = pislopen * ( 1.+ zadap  * EXP( -0.25 * enano(ji,jj,jk) ) )  &
                  &                   * trb(ji,jj,jk,jpnch) /( trb(ji,jj,jk,jpphy) * 12. + rtrn)

                  ! Diatoms
                  zpislopeadd(ji,jj,jk) = (pislopen * zconctemp2 + pisloped * zconctemp) / ( trb(ji,jj,jk,jpdia) + rtrn )   &
                  &                   * trb(ji,jj,jk,jpdch) /( trb(ji,jj,jk,jpdia) * 12. + rtrn)
               ENDIF
            END DO
         END DO
      END DO

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( etot_ndcy(ji,jj,jk) > 1.E-3 ) THEN
                   ! Computation of production function for Carbon
                   ! Actual light levels are used here 
                   ! ----------------------------------------------
                   zpislopen = zpislopeadn(ji,jj,jk) / ( ( r1_rday + bresp * r1_rday ) &
                   &            * zmxl_fac(ji,jj,jk) * rday + rtrn)
                   zpisloped = zpislopeadd(ji,jj,jk) / ( ( r1_rday + bresp * r1_rday ) &
                   &            * zmxl_fac(ji,jj,jk) * rday + rtrn)
                   zprbio(ji,jj,jk) = zprbio(ji,jj,jk) * ( 1.- EXP( -zpislopen * enano(ji,jj,jk) )  )
                   zprdia(ji,jj,jk) = zprdia(ji,jj,jk) * ( 1.- EXP( -zpisloped * ediat(ji,jj,jk) )  )

                   !  Computation of production function for Chlorophyll
                   !  Mean light level in the mixed layer (when appropriate)
                   !  is used here (acclimation is in general slower than 
                   !  the characteristic time scales of vertical mixing)
                   !  ------------------------------------------------------
                   zpislopen = zpislopeadn(ji,jj,jk) / ( zprmaxn(ji,jj,jk) * zmxl_chl(ji,jj,jk) * rday + rtrn )
                   zpisloped = zpislopeadd(ji,jj,jk) / ( zprmaxd(ji,jj,jk) * zmxl_chl(ji,jj,jk) * rday + rtrn )
                   zprchln(ji,jj,jk) = zprmaxn(ji,jj,jk) * ( 1.- EXP( -zpislopen * enanom(ji,jj,jk) ) )
                   zprchld(ji,jj,jk) = zprmaxd(ji,jj,jk) * ( 1.- EXP( -zpisloped * ediatm(ji,jj,jk) ) )
               ENDIF
            END DO
         END DO
      END DO

      !  Computation of a proxy of the N/C quota from nutrient limitation 
      !  and light limitation. Steady state is assumed to allow the computation
      !  ----------------------------------------------------------------------
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
                zval = MIN( xnanopo4(ji,jj,jk), ( xnanonh4(ji,jj,jk) + xnanono3(ji,jj,jk) ) )   &
                &      * zprmaxn(ji,jj,jk) / ( zprbio(ji,jj,jk) + rtrn )
                quotan(ji,jj,jk) = MIN( 1., 0.3 + 0.7 * zval )
                zval = MIN( xdiatpo4(ji,jj,jk), ( xdiatnh4(ji,jj,jk) + xdiatno3(ji,jj,jk) ) )   &
                &      * zprmaxd(ji,jj,jk) / ( zprdia(ji,jj,jk) + rtrn )
                quotad(ji,jj,jk) = MIN( 1., 0.3 + 0.7 * zval )
            END DO
         END DO
      END DO


      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi

               IF( etot_ndcy(ji,jj,jk) > 1.E-3 ) THEN

                   ! Si/C of diatoms
                   ! ------------------------
                   ! Si/C increases with iron stress and silicate availability
                   ! Si/C is arbitrariliy increased for very high Si concentrations
                   ! to mimic the very high ratios observed in the Southern Ocean (zsilfac)
                   ! A parameterization derived from Flynn (2003) is used for the control
                   ! when Si is not limiting which is similar to the parameterisation
                   ! proposed by Gurney and Davidson (1999).
                   ! -----------------------------------------------------------------------
                 zlim  = trb(ji,jj,jk,jpsil) / ( trb(ji,jj,jk,jpsil) + xksi1 )
                 zsilim = xlimdia(ji,jj,jk) * zprdia(ji,jj,jk) / ( zprmaxd(ji,jj,jk) + rtrn )
                 zsiborn = trb(ji,jj,1,jpsil)**3
                 IF (gphit(ji,jj) < -30.0 ) THEN
                   zsilfac = 1. + 2. * zsiborn / ( zsiborn + xksi2**3 )
                 ELSE
                   zsilfac = 1. +      zsiborn / ( zsiborn + xksi2**3 )
                 ENDIF
                 zratiosi = 1.0 - trb(ji,jj,jk,jpdsi) / ( trb(ji,jj,jk,jpdia) + rtrn ) / ( zsilfac * grosip * 3.0 + rtrn )
                 zratiosi = MAX(0., MIN(1.0, zratiosi) )
                 zmaxsi  = (1.0 + 0.1**4) * zratiosi**4 / ( zratiosi**4 + 0.1**4 )
                 IF ( xlimsi(ji,jj,jk) /= xlimdia(ji,jj,jk) ) THEN
                    zysopt(ji,jj,jk) = zlim * zsilfac * grosip * 1.0 * zmaxsi
                 ELSE
                    zysopt(ji,jj,jk) = zlim * zsilfac * grosip * 1.0 * zsilim**0.7 * zmaxsi
                 ENDIF
              ENDIF
            END DO
         END DO
      END DO

      ! Sea-ice effect on production
      ! No production is assumed below sea ice
      ! -------------------------------------- 
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               zprbio(ji,jj,jk) = zprbio(ji,jj,jk) * ( 1. - fr_i(ji,jj) )
               zprdia(ji,jj,jk) = zprdia(ji,jj,jk) * ( 1. - fr_i(ji,jj) )
            END DO
         END DO
      END DO

      ! Computation of the various production  and nutrient uptake terms
      ! ---------------------------------------------------------------
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( etot_ndcy(ji,jj,jk) > 1.E-3 ) THEN
                  !  production term of nanophyto. (C)
                  zprorcan(ji,jj,jk) = zprbio(ji,jj,jk)  * xlimphy(ji,jj,jk) * trb(ji,jj,jk,jpphy) * rfact2

                  !  New production (uptake of NO3)
                  zpronewn(ji,jj,jk)  = zprorcan(ji,jj,jk)* xnanono3(ji,jj,jk) / ( xnanono3(ji,jj,jk) + xnanonh4(ji,jj,jk) + rtrn )

                  ! Size computation
                  ! Size is made a function of the limitation of of phytoplankton growth
                  ! Strongly limited cells are supposed to be smaller. sizena is the 
                  ! size at time step t+1 and is thus updated at the end of the 
                  ! current time step
                  ! --------------------------------------------------------------------
                  zlimfac = xlimphy(ji,jj,jk) * zprchln(ji,jj,jk) / ( zprmaxn(ji,jj,jk) + rtrn )
                  zsizetmp = 1.0 + 1.3 * ( xsizern - 1.0 ) * zlimfac**3/(0.3 + zlimfac**3)
                  sizena(ji,jj,jk) = min(xsizern, max( sizena(ji,jj,jk), zsizetmp ) )

                  ! Iron uptake rates of nanophytoplankton. Upregulation is  
                  ! not parameterized at low iron concentrations as observations
                  ! do not suggest it for accimated cells. Uptake is
                  ! downregulated when the quota is close to the maximum quota
                  zfecnm = xqfuncfecn(ji,jj,jk) + ( fecnm - xqfuncfecn(ji,jj,jk) ) * ( xnanono3(ji,jj,jk) + xnanonh4(ji,jj,jk) )
                  zratio = 1.0 - MIN(1.0,trb(ji,jj,jk,jpnfe) / ( trb(ji,jj,jk,jpphy) * zfecnm + rtrn ) )
                  zmax   = MAX( 0., MIN( 1.0, zratio**2/ (0.05**2+zratio**2) ) ) 
                  zprofen(ji,jj,jk) = zfecnm * zprmaxn(ji,jj,jk) * ( 1.0 - fr_i(ji,jj) )  &
                  &          * (1. + 0.8 * xnanono3(ji,jj,jk) / ( rtrn + xnanono3(ji,jj,jk)  &
                  &          + xnanonh4(ji,jj,jk) ) * (1. - xnanofer(ji,jj,jk) ) )   &
                  &          * xnanofer(ji,jj,jk) * zmax * trb(ji,jj,jk,jpphy) * rfact2
                  ! production terms of diatoms (C)
                  zprorcad(ji,jj,jk) = zprdia(ji,jj,jk) * xlimdia(ji,jj,jk) * trb(ji,jj,jk,jpdia) * rfact2

                  ! New production (uptake of NO3)
                  zpronewd(ji,jj,jk) = zprorcad(ji,jj,jk) * xdiatno3(ji,jj,jk) / ( xdiatno3(ji,jj,jk) + xdiatnh4(ji,jj,jk) + rtrn )

                  ! Size computation
                  ! Size is made a function of the limitation of of phytoplankton growth
                  ! Strongly limited cells are supposed to be smaller. sizeda is
                  ! size at time step t+1 and is thus updated at the end of the 
                  ! current time step. 
                  ! --------------------------------------------------------------------
                  zlimfac = zprchld(ji,jj,jk) * xlimdia(ji,jj,jk) / ( zprmaxd(ji,jj,jk) + rtrn )
                  zsizetmp = 1.0 + 1.3 * ( xsizerd - 1.0 ) * zlimfac**3/(0.3 + zlimfac**3)
                  sizeda(ji,jj,jk) = min(xsizerd, max( sizeda(ji,jj,jk), zsizetmp ) )

                  ! Iron uptake rates of diatoms. Upregulation is  
                  ! not parameterized at low iron concentrations as observations
                  ! do not suggest it for accimated cells. Uptake is
                  ! downregulated when the quota is close to the maximum quota
                  zfecdm = xqfuncfecd(ji,jj,jk) + ( fecdm - xqfuncfecd(ji,jj,jk) ) * ( xdiatno3(ji,jj,jk) + xdiatnh4(ji,jj,jk) )
                  zratio = 1.0 - MIN(1.0, trb(ji,jj,jk,jpdfe) / ( trb(ji,jj,jk,jpdia) * zfecdm + rtrn ) )
                  zmax   = MAX( 0., MIN( 1.0, zratio**2/ (0.05**2+zratio**2) ) ) 
                  zprofed(ji,jj,jk) = zfecdm * zprmaxd(ji,jj,jk) * (1.0 - fr_i(ji,jj) )  &
                  &          * (1. + 0.8 * xdiatno3(ji,jj,jk) / ( rtrn + xdiatno3(ji,jj,jk)  &
                  &          + xdiatnh4(ji,jj,jk) ) * (1. - xdiatfer(ji,jj,jk) ) )   &
                  &          * xdiatfer(ji,jj,jk) * zmax * trb(ji,jj,jk,jpdia) * rfact2
               ENDIF
            END DO
         END DO
      END DO

      ! Computation of the chlorophyll production terms
      ! The parameterization is taken from Geider et al. (1997)
      ! -------------------------------------------------------
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( etot_ndcy(ji,jj,jk) > 1.E-3 ) THEN
                  !  production term for nanophyto. ( chlorophyll )
                  znanotot = enanom(ji,jj,jk) / ( zmxl_chl(ji,jj,jk) + rtrn )
                  zprod    = rday * zprorcan(ji,jj,jk) * zprchln(ji,jj,jk) * xlimphy(ji,jj,jk)
                  zprochln = chlcmin * 12. * zprorcan (ji,jj,jk)
                  zprochln = zprochln + (chlcnm - chlcmin) * 12. * zprod / &
                                        & (  zpislopeadn(ji,jj,jk) * znanotot +rtrn)

                  ! production terms of diatoms ( chlorophyll )
                  zdiattot = ediatm(ji,jj,jk) / ( zmxl_chl(ji,jj,jk) + rtrn )
                  zprod    = rday * zprorcad(ji,jj,jk) * zprchld(ji,jj,jk) * xlimdia(ji,jj,jk)
                  zprochld = chlcmin * 12. * zprorcad(ji,jj,jk)
                  zprochld = zprochld + (chlcdm - chlcmin) * 12. * zprod / &
                                        & ( zpislopeadd(ji,jj,jk) * zdiattot +rtrn )

                  ! Update the arrays TRA which contain the Chla sources and sinks
                  tra(ji,jj,jk,jpnch) = tra(ji,jj,jk,jpnch) + zprochln * texcretn
                  tra(ji,jj,jk,jpdch) = tra(ji,jj,jk,jpdch) + zprochld * texcretd
               ENDIF
            END DO
         END DO
      END DO

      !   Update the arrays TRA which contain the biological sources and sinks
      DO jk = 1, jpkm1
         DO jj = 1, jpj
           DO ji =1 ,jpi
              IF( etot_ndcy(ji,jj,jk) > 1.E-3 ) THEN
                 zproreg  = zprorcan(ji,jj,jk) - zpronewn(ji,jj,jk)
                 zproreg2 = zprorcad(ji,jj,jk) - zpronewd(ji,jj,jk)
                 zdocprod = excretd * zprorcad(ji,jj,jk) + excretn * zprorcan(ji,jj,jk)
                 tra(ji,jj,jk,jppo4) = tra(ji,jj,jk,jppo4) - zprorcan(ji,jj,jk) - zprorcad(ji,jj,jk)
                 tra(ji,jj,jk,jpno3) = tra(ji,jj,jk,jpno3) - zpronewn(ji,jj,jk) - zpronewd(ji,jj,jk)
                 tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) - zproreg - zproreg2
                 tra(ji,jj,jk,jpphy) = tra(ji,jj,jk,jpphy) + zprorcan(ji,jj,jk) * texcretn
                 tra(ji,jj,jk,jpnfe) = tra(ji,jj,jk,jpnfe) + zprofen(ji,jj,jk) * texcretn
                 tra(ji,jj,jk,jpdia) = tra(ji,jj,jk,jpdia) + zprorcad(ji,jj,jk) * texcretd
                 tra(ji,jj,jk,jpdfe) = tra(ji,jj,jk,jpdfe) + zprofed(ji,jj,jk) * texcretd
                 tra(ji,jj,jk,jpdsi) = tra(ji,jj,jk,jpdsi) + zprmaxd(ji,jj,jk) * zysopt(ji,jj,jk)   &
                 &                   * rfact2 * trb(ji,jj,jk,jpdia)
                 tra(ji,jj,jk,jpdoc) = tra(ji,jj,jk,jpdoc) + zdocprod
                 tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy) + o2ut * ( zproreg + zproreg2) &
                 &                   + ( o2ut + o2nit ) * ( zpronewn(ji,jj,jk) + zpronewd(ji,jj,jk) )
                 !
                 tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) - ( texcretn * zprofen(ji,jj,jk)    &
                 &                   + texcretd * zprofed(ji,jj,jk) )
                 consfe3(ji,jj,jk) = ( texcretn * zprofen(ji,jj,jk) + texcretd * zprofed(ji,jj,jk) )     &
                 &                    * 75.0 / ( rtrn + ( plig(ji,jj,jk) + 75.0 * (1.0 - plig(ji,jj,jk) ) ) * trb(ji,jj,jk,jpfer) ) / rfact2
                 tra(ji,jj,jk,jpsil) = tra(ji,jj,jk,jpsil) - zprmaxd(ji,jj,jk) * zysopt(ji,jj,jk)   &
                 &                   * rfact2 * trb(ji,jj,jk,jpdia)
                 tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) - zprorcan(ji,jj,jk) - zprorcad(ji,jj,jk)
                 tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) + rno3 * ( zpronewn(ji,jj,jk) + zpronewd(ji,jj,jk) ) &
                 &                                         - rno3 * ( zproreg + zproreg2 )
              ENDIF
           END DO
        END DO
     END DO

     ! Production and uptake of ligands by phytoplankton. This part is activated 
     ! when ln_ligand is set to .true. in the namelist. Ligand uptake is small 
     ! and based on the FeL model by Morel et al. (2008) and on the study of
     ! Shaked et al. (2020)
     ! -------------------------------------------------------------------------
     IF( ln_ligand ) THEN
         zpligprod1(:,:,:) = 0._wp    ;    zpligprod2(:,:,:) = 0._wp
         DO jk = 1, jpkm1
            DO jj = 1, jpj
              DO ji =1 ,jpi
                 IF( etot_ndcy(ji,jj,jk) > 1.E-3 ) THEN
                    zdocprod = excretd * zprorcad(ji,jj,jk) + excretn * zprorcan(ji,jj,jk)
                    zfeup    = texcretn * zprofen(ji,jj,jk) + texcretd * zprofed(ji,jj,jk)
                    tra(ji,jj,jk,jplgw) = tra(ji,jj,jk,jplgw) + zdocprod * ldocp    &
                    &       - zfeup * plig(ji,jj,jk) / ( rtrn + plig(ji,jj,jk) + 75.0 * (1.0 - plig(ji,jj,jk) ) ) * lthet
                    zpligprod1(ji,jj,jk) = zdocprod * ldocp
                    zpligprod2(ji,jj,jk) = zfeup * plig(ji,jj,jk) / ( rtrn + plig(ji,jj,jk) &
                    &                      + 75.0 * (1.0 - plig(ji,jj,jk) ) ) * lthet
                 ENDIF
              END DO
           END DO
        END DO
     ENDIF


    ! Output of the diagnostics
    ! Total primary production per year
    IF( iom_use( "tintpp" ) .OR. ( ln_check_mass .AND. kt == nitend .AND. knt == nrdttrc )  )  &
         & tpp = glob_sum( 'p4zprod', ( zprorcan(:,:,:) + zprorcad(:,:,:) ) * cvol(:,:,:) )

    IF( lk_iomput ) THEN
       IF( knt == nrdttrc ) THEN
          ALLOCATE( zw2d(jpi,jpj), zw3d(jpi,jpj,jpk) )
          zfact = 1.e+3 * rfact2r  !  conversion from mol/l/kt to  mol/m3/s
          !
          IF( iom_use( "PPPHYN" ) .OR. iom_use( "PPPHYD" ) )  THEN
              zw3d(:,:,:) = zprorcan(:,:,:) * zfact * tmask(:,:,:)  ! primary production by nanophyto
              CALL iom_put( "PPPHYN"  , zw3d )
              !
              zw3d(:,:,:) = zprorcad(:,:,:) * zfact * tmask(:,:,:)  ! primary production by diatomes
              CALL iom_put( "PPPHYD"  , zw3d )
          ENDIF
          IF( iom_use( "PPNEWN" ) .OR. iom_use( "PPNEWD" ) )  THEN
              zw3d(:,:,:) = zpronewn(:,:,:) * zfact * tmask(:,:,:)  ! new primary production by nanophyto
              CALL iom_put( "PPNEWN"  , zw3d )
              !
              zw3d(:,:,:) = zpronewd(:,:,:) * zfact * tmask(:,:,:)  ! new primary production by diatoms
              CALL iom_put( "PPNEWD"  , zw3d )
          ENDIF
          IF( iom_use( "PBSi" ) )  THEN
              zw3d(:,:,:) = zprmaxd(:,:,:) * 1.E3 * tmask(:,:,:) * zysopt(:,:,:) * trb(:,:,:,jpdia) ! biogenic silica production
              CALL iom_put( "PBSi"  , zw3d )
          ENDIF
          IF( iom_use( "PFeN" ) .OR. iom_use( "PFeD" ) )  THEN
              zw3d(:,:,:) = zprofen(:,:,:) * zfact * tmask(:,:,:)  ! biogenic iron production by nanophyto
              CALL iom_put( "PFeN"  , zw3d )
              !
              zw3d(:,:,:) = zprofed(:,:,:) * zfact * tmask(:,:,:)  ! biogenic iron production by  diatoms
              CALL iom_put( "PFeD"  , zw3d )
          ENDIF
          IF( iom_use( "LPRODP" ) )  THEN
              zw3d(:,:,:) = zpligprod1(:,:,:) * 1e9 * zfact * tmask(:,:,:)  ! Ligand production by phytoplankton
              CALL iom_put( "LPRODP"  , zw3d )
          ENDIF
          IF( iom_use( "QFEMAXN" ) )  THEN
              zw3d(:,:,:) = xqfuncfecn(:,:,:) + ( fecnm - xqfuncfecn(:,:,:) ) * ( xnanono3(:,:,:) + xnanonh4(:,:,:) ) * tmask(:,:,:)  ! Ligand production by phytoplankton
              CALL iom_put( "QFEMAXD"  , zw3d )
          ENDIF
          IF( iom_use( "QFEMAXN" ) )  THEN
              zw3d(:,:,:) = xqfuncfecd(:,:,:) + ( fecdm - xqfuncfecn(:,:,:) ) * ( xdiatno3(:,:,:) + xdiatnh4(:,:,:) ) * tmask(:,:,:)  ! Ligand production by phytoplankton
              CALL iom_put( "QFEMAXD"  , zw3d )
          ENDIF

          IF( iom_use( "LDETP" ) )  THEN
              zw3d(:,:,:) = zpligprod2(:,:,:) * 1e9 * zfact * tmask(:,:,:)  ! Uptake of ligands by phytoplankton
              CALL iom_put( "LDETP"  , zw3d )
          ENDIF
          IF( iom_use( "Mumax" ) )  THEN
              zw3d(:,:,:) = zprmaxn(:,:,:) * tmask(:,:,:)   ! Maximum growth rate
              CALL iom_put( "Mumax"  , zw3d )
          ENDIF
          IF( iom_use( "MuN" ) .OR. iom_use( "MuD" ) )  THEN
              zw3d(:,:,:) = zprbio(:,:,:) * xlimphy(:,:,:) * tmask(:,:,:)  ! Realized growth rate for nanophyto
              CALL iom_put( "MuN"  , zw3d )
              !
              zw3d(:,:,:) =  zprdia(:,:,:) * xlimdia(:,:,:) * tmask(:,:,:)  ! Realized growth rate for diatoms
              CALL iom_put( "MuD"  , zw3d )
          ENDIF
          IF( iom_use( "LNlight" ) .OR. iom_use( "LDlight" ) )  THEN
              zw3d(:,:,:) = zprbio (:,:,:) / (zprmaxn(:,:,:) + rtrn) * tmask(:,:,:) ! light limitation term of nanophytoplankton
              CALL iom_put( "LNlight"  , zw3d )
              !
              zw3d(:,:,:) = zprdia (:,:,:) / (zprmaxd(:,:,:) + rtrn) * tmask(:,:,:)  ! light limitation term of diatoms
              CALL iom_put( "LDlight"  , zw3d )
          ENDIF
          IF( iom_use( "TPP" ) )  THEN
              zw3d(:,:,:) = ( zprorcan(:,:,:) + zprorcad(:,:,:) ) * zfact * tmask(:,:,:)  ! total primary production
              CALL iom_put( "TPP"  , zw3d )
          ENDIF
          IF( iom_use( "TPNEW" ) )  THEN
              zw3d(:,:,:) = ( zpronewn(:,:,:) + zpronewd(:,:,:) ) * zfact * tmask(:,:,:)  ! total new production
              CALL iom_put( "TPNEW"  , zw3d )
          ENDIF
          IF( iom_use( "TPBFE" ) )  THEN
              zw3d(:,:,:) = ( zprofen(:,:,:) + zprofed(:,:,:) ) * zfact * tmask(:,:,:)  ! total biogenic iron production
              CALL iom_put( "TPBFE"  , zw3d )
          ENDIF
          IF( iom_use( "INTPPPHYN" ) .OR. iom_use( "INTPPPHYD" ) ) THEN  
             zw2d(:,:) = 0.
             DO jk = 1, jpkm1
               zw2d(:,:) = zw2d(:,:) + zprorcan(:,:,jk) * e3t_n(:,:,jk) * zfact * tmask(:,:,jk)  ! vert. integrated  primary produc. by nano
             ENDDO
             CALL iom_put( "INTPPPHYN" , zw2d )
             !
             zw2d(:,:) = 0.
             DO jk = 1, jpkm1
                zw2d(:,:) = zw2d(:,:) + zprorcad(:,:,jk) * e3t_n(:,:,jk) * zfact * tmask(:,:,jk) ! vert. integrated  primary produc. by diatom
             ENDDO
             CALL iom_put( "INTPPPHYD" , zw2d )
          ENDIF
          IF( iom_use( "INTPP" ) ) THEN   
             zw2d(:,:) = 0.
             DO jk = 1, jpkm1
                zw2d(:,:) = zw2d(:,:) + ( zprorcan(:,:,jk) + zprorcad(:,:,jk) ) * e3t_n(:,:,jk) * zfact * tmask(:,:,jk) ! vert. integrated pp
             ENDDO
             CALL iom_put( "INTPP" , zw2d )
          ENDIF
          IF( iom_use( "INTPNEW" ) ) THEN    
             zw2d(:,:) = 0.
             DO jk = 1, jpkm1
                zw2d(:,:) = zw2d(:,:) + ( zpronewn(:,:,jk) + zpronewd(:,:,jk) ) * e3t_n(:,:,jk) * zfact * tmask(:,:,jk)  ! vert. integrated new prod
             ENDDO
             CALL iom_put( "INTPNEW" , zw2d )
          ENDIF
          IF( iom_use( "INTPBFE" ) ) THEN           !   total biogenic iron production  ( vertically integrated )
             zw2d(:,:) = 0.
             DO jk = 1, jpkm1
                zw2d(:,:) = zw2d(:,:) + ( zprofen(:,:,jk) + zprofed(:,:,jk) ) * e3t_n(:,:,jk) * zfact * tmask(:,:,jk) ! vert integr. bfe prod
             ENDDO
            CALL iom_put( "INTPBFE" , zw2d )
          ENDIF
          IF( iom_use( "INTPBSI" ) ) THEN           !   total biogenic silica production  ( vertically integrated )
             zw2d(:,:) = 0.
             DO jk = 1, jpkm1
                zw2d(:,:) = zw2d(:,:) + zprorcad(:,:,jk) * zysopt(:,:,jk) * e3t_n(:,:,jk) * zfact * tmask(:,:,jk)  ! vert integr. bsi prod
             ENDDO
             CALL iom_put( "INTPBSI" , zw2d )
          ENDIF
          IF( iom_use( "tintpp" ) )  CALL iom_put( "tintpp" , tpp * zfact )  !  global total integrated primary production molC/s
          !
          DEALLOCATE( zw2d, zw3d )
       ENDIF
     ENDIF

     IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('prod')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
     ENDIF
      !
      IF( ln_timing )  CALL timing_stop('p4z_prod')
      !
   END SUBROUTINE p4z_prod


   SUBROUTINE p4z_prod_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_prod_init  ***
      !!
      !! ** Purpose :   Initialization of phytoplankton production parameters
      !!
      !! ** Method  :   Read the namp4zprod namelist and check the parameters
      !!      called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist namp4zprod
      !!----------------------------------------------------------------------
      INTEGER ::   ios   ! Local integer
      !
      ! Namelist block
      NAMELIST/namp4zprod/ pislopen, pisloped, xadap, bresp, excretn, excretd,  &
         &                 chlcnm, chlcdm, chlcmin, fecnm, fecdm, grosip
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN                         ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'p4z_prod_init : phytoplankton growth'
         WRITE(numout,*) '~~~~~~~~~~~~~'
      ENDIF
      !
      REWIND( numnatp_ref )              ! Namelist namp4zprod in reference namelist : Pisces phytoplankton production
      READ  ( numnatp_ref, namp4zprod, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namp4zprod in reference namelist' )
      REWIND( numnatp_cfg )              ! Namelist namp4zprod in configuration namelist : Pisces phytoplankton production
      READ  ( numnatp_cfg, namp4zprod, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namp4zprod in configuration namelist' )
      IF(lwm) WRITE( numonp, namp4zprod )

      IF(lwp) THEN                         ! control print
         WRITE(numout,*) '   Namelist : namp4zprod'
         WRITE(numout,*) '      mean Si/C ratio                           grosip       =', grosip
         WRITE(numout,*) '      P-I slope                                 pislopen     =', pislopen
         WRITE(numout,*) '      Acclimation factor to low light           xadap        =', xadap
         WRITE(numout,*) '      excretion ratio of nanophytoplankton      excretn      =', excretn
         WRITE(numout,*) '      excretion ratio of diatoms                excretd      =', excretd
         WRITE(numout,*) '      basal respiration in phytoplankton        bresp        =', bresp
         WRITE(numout,*) '      Maximum Chl/C in phytoplankton            chlcmin      =', chlcmin
         WRITE(numout,*) '      P-I slope  for diatoms                    pisloped     =', pisloped
         WRITE(numout,*) '      Minimum Chl/C in nanophytoplankton        chlcnm       =', chlcnm
         WRITE(numout,*) '      Minimum Chl/C in diatoms                  chlcdm       =', chlcdm
         WRITE(numout,*) '      Maximum Fe/C in nanophytoplankton         fecnm        =', fecnm
         WRITE(numout,*) '      Minimum Fe/C in diatoms                   fecdm        =', fecdm
      ENDIF
      !
      r1_rday   = 1._wp / rday 
      texcretn  = 1._wp - excretn
      texcretd  = 1._wp - excretd
      tpp       = 0._wp
      !
   END SUBROUTINE p4z_prod_init


   INTEGER FUNCTION p4z_prod_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_prod_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( quotan(jpi,jpj,jpk), quotad(jpi,jpj,jpk), STAT = p4z_prod_alloc )
      !
      IF( p4z_prod_alloc /= 0 ) CALL ctl_stop( 'STOP', 'p4z_prod_alloc : failed to allocate arrays.' )
      !
   END FUNCTION p4z_prod_alloc

   !!======================================================================
END MODULE p4zprod
