MODULE usrdef_istate
   !!======================================================================
   !!                     ***  MODULE usrdef_istate   ***
   !!
   !!                      ===  BENCH configuration  ===
   !!
   !! User defined : set the initial state of a user configuration
   !!======================================================================
   !! History :  NEMO !
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!  usr_def_istate : initial state in Temperature and salinity
   !!----------------------------------------------------------------------
   USE par_oce        ! ocean space and time domain
   USE dom_oce        
   USE phycst         ! physical constants
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE lbclnk          ! lateral boundary conditions - mpp exchanges
   !   
   USE usrdef_nam
   
   IMPLICIT NONE
   PRIVATE

   PUBLIC   usr_def_istate   ! called by istate.F90

   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2016)
   !! $Id$ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS
  
   SUBROUTINE usr_def_istate( pdept, ptmask, pts, pu, pv, pssh )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE usr_def_istate  ***
      !! 
      !! ** Purpose :   Initialization of the dynamics and tracers
      !!                Here BENCH configuration 
      !!
      !! ** Method  :   Set a gaussian anomaly of pressure and associated
      !!                geostrophic velocities
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(in   ) ::   pdept   ! depth of t-point               [m]
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(in   ) ::   ptmask  ! t-point ocean mask             [m]
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts), INTENT(  out) ::   pts     ! T & S fields      [Celsius ; g/kg]
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(  out) ::   pu      ! i-component of the velocity  [m/s] 
      REAL(wp), DIMENSION(jpi,jpj,jpk)     , INTENT(  out) ::   pv      ! j-component of the velocity  [m/s] 
      REAL(wp), DIMENSION(jpi,jpj)         , INTENT(  out) ::   pssh    ! sea-surface height
      !
      REAL(wp), DIMENSION(jpi,jpj) ::   z2d   ! 2D workspace
      REAL(wp) ::   zfact
      INTEGER  ::   ji, jj, jk
      INTEGER  ::   igloi, igloj   ! to be removed in the future, see comment bellow
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'usr_def_istate : BENCH configuration, analytical definition of initial state'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~   '
      !
      ! define unique value on each point of the inner global domain. z2d ranging from 0.05 to -0.05
      !
      ! WARNING: to keep compatibility with the trunk that was including periodocity into the input data, 
      ! we must define z2d as bellow.
      ! Once we decide to forget trunk compatibility, we must simply define z2d as:
!!$      DO_2D( 0, 0, 0, 0 )
!!$         z2d(ji,jj) = 0.1 * ( 0.5 - REAL( mig0_oldcmp(ji) + (mjg0_oldcmp(jj)-1) * Ni0glo, wp ) / REAL( Ni0glo * Nj0glo, wp ) )
!!$      END_2D
      igloi = Ni0glo + 2 * COUNT( (/ jperio == 1 .OR. jperio == 4 .OR. jperio == 6 .OR. jperio == 7 /) )
      igloj = Nj0glo + 2 * COUNT( (/ jperio == 2 .OR. jperio == 7 /) ) + 1 * COUNT( (/ jperio >= 4 .AND. jperio <= 6 /) )
      DO_2D( 0, 0, 0, 0 )
         z2d(ji,jj) = 0.1 * ( 0.5 - REAL( mig0_oldcmp(ji) + (mjg0_oldcmp(jj)-1) * igloi, wp ) / REAL( igloi * igloj, wp ) )
      END_2D
      !
      ! sea level:
      pssh(:,:) = z2d(:,:)                                                ! +/- 0.05 m
      !
      DO_3D( 0, 0, 0, 0, 1, jpkm1 )
         zfact = REAL(jk-1,wp) / REAL(jpk-1,wp)   ! 0 to 1 to add a basic stratification
         ! temperature choosen to lead to ~50% ice at the beginning if rn_thres_sst = 0.5
         pts(:,:,jk,jp_tem) = 20._wp*z2d(:,:) - 1._wp - 0.5_wp * zfact    ! -1 to -1.5 +/-1.0 degG
         ! salinity:  
         pts(:,:,jk,jp_sal) = 30._wp + 1._wp * zfact + z2d(:,:)           ! 30 to 31 +/- 0.05 psu
         ! velocities:
         pu(:,:,jk) = z2d(:,:) *  0.1_wp * umask(:,:,jk)                  ! +/- 0.005  m/s
         pv(:,:,jk) = z2d(:,:) * 0.01_wp * vmask(:,:,jk)                  ! +/- 0.0005 m/s
      END_3D
      pts(:,:,jpk,:) = 0._wp
      pu( :,:,jpk  ) = 0._wp
      pv( :,:,jpk  ) = 0._wp
      !
      CALL lbc_lnk('usrdef_istate', pssh, 'T',  1. )            ! apply boundary conditions
      CALL lbc_lnk('usrdef_istate',  pts, 'T',  1. )            ! apply boundary conditions
      CALL lbc_lnk('usrdef_istate',   pu, 'U', -1. )            ! apply boundary conditions
      CALL lbc_lnk('usrdef_istate',   pv, 'V', -1. )            ! apply boundary conditions
      
   END SUBROUTINE usr_def_istate

   !!======================================================================
END MODULE usrdef_istate
