MODULE gpu_manager
    USE lib_mpp        ! distributed memory computing library
    USE cudafor
    USE par_oce
    !!======================================================================
    !!                       ***  MODULE  gpu_manager  ***
    !! GPU CUDA utilities : Defines run parameters for GPGPU routines
    !!=====================================================================
    !! History :   1.0  !  2019-10  (M. Faria)   original code
    !!----------------------------------------------------------------------

    !!----------------------------------------------------------------------
    IMPLICIT NONE
    PUBLIC
    !!----------------------------------------------------------------------
    !!                   gpu configuration parameters
    !!----------------------------------------------------------------------
    LOGICAL          ::   gpu      = .TRUE.     ! to control Cuda use
    TYPE (cudaEvent) :: startEvent, stopEvent
    REAL             :: cudatime
    INTEGER          :: cudaistat, cudadevice

   !pinned dom_oce.f90
    !   INTEGER , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)           ::   mbkt, mbku, mbkv                        !: bottom last wet T-, U- and V-level
!   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, TARGET, DIMENSION(:,:)   ::   e1v   , e2u  , r1_e1v, r1_e2v , r1_e1u  !: horizontal scale factors at v-point [m]
!   REAL(wp), PUBLIC, ALLOCATABLE, SAVE        , DIMENSION(:,:)   ::   e1e2t , r1_e1e2t                        !: associated metrics at t-point
!   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)         ::   e3t_n ,   e3u_n ,   e3v_n ,  e3w_n      !: t- vert. scale factor [m]
   !pinned oce.f90
!   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)           ::           sshn           !: sea surface height at t-point [m]
!   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:)       ::           tsn            !: 4D T-S fields (page-locked)    [Celsius,psu]
!   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)         ::           un             !: i-horizontal velocity        [m/s]
!   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)         ::           wn             !: vertical velocity            [m/s]
!   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)         ::           vn             !: j-horizontal velocity        [m/s]
!   !pinned zdf_oce.f90
!   REAL(wp), PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:,:)         ::   avm, avt, avs  !: vertical mixing coefficients (w-point) [m2/s]
!   !pinned zdfdrg
!   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:), PUBLIC           ::  rCdU_bot   !: top/bottom drag coeff. at t-point (<0)  [m/s]

!!----------------------------------------------------------------------
!! NEMO/OCE 4.0 , NEMO Consortium (2018)
!! $Id:
!! Software governed by the CeCILL license (see ./LICENSE)
!!=====================================================================
CONTAINS
    SUBROUTINE setdevice()
!         cudaistat = cudaSetDevice( MOD( mpprank /  ( 20 / 2 )  , 4 )  ) ! Set which GPU use. Ex for 4 GPUs. cpus_per_socket / sockets_per_node for "bind-to core"
!         print *, 'cudaistat ', cudaistat, mpprank,  MOD( mpprank /  ( 20 / 2 )  , 4 )
        cudaistat = cudaSetDevice( mpprank    )
        print *, 'cudaistat ', cudaistat, mpprank,   mpprank  


        !cudaistat = cudaSetDevice( 1 )
        !print *, 'cudaistat ', cudaistat, mpprank,  1
    END SUBROUTINE setdevice

!
!    SUBROUTINE alloc_pinned()
!        !from dom_oce.f90
!        allocate(mbkt(jpi,jpj)  , mbku(jpi,jpj)    , mbkv(jpi,jpj)                                     )
!        allocate(e1v(jpi,jpj)   , e2u(jpi,jpj)     , r1_e1v(jpi,jpj), r1_e2v(jpi,jpj) , r1_e1u(jpi,jpj))
!        allocate(e1e2t(jpi,jpj) , r1_e1e2t(jpi,jpj)                                                    )
!        allocate(e3t_n(jpi,jpj,jpk) ,   e3u_n(jpi,jpj,jpk) ,   e3v_n(jpi,jpj,jpk) ,  e3w_n(jpi,jpj,jpk))
!        !from oce.f90
!        allocate(tsn(jpi,jpj,jpk,2))
!        allocate(sshn(jpi,jpj))
!        allocate(un(jpi,jpj,jpk), wn(jpi,jpj,jpk), vn(jpi,jpj,jpk))
!        !from zdf_oce.f90
!        allocate(avm(jpi,jpj,jpk), avt(jpi,jpj,jpk), avs(jpi,jpj,jpk))
!        !from zdfdrg.f90
!        allocate(rCdU_bot(jpi,jpj))
!    END SUBROUTINE



END MODULE gpu_manager

