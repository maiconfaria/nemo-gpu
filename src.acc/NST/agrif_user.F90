#undef UPD_HIGH   /* MIX HIGH UPDATE */
#if defined key_agrif
   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/NST 4.0 , NEMO Consortium (2018)
   !! $Id$
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
   SUBROUTINE agrif_user
   END SUBROUTINE agrif_user

   
   SUBROUTINE agrif_before_regridding
   END SUBROUTINE agrif_before_regridding

   
   SUBROUTINE Agrif_InitWorkspace
   END SUBROUTINE Agrif_InitWorkspace

   
   SUBROUTINE Agrif_InitValues
      !!----------------------------------------------------------------------
      !!                 *** ROUTINE Agrif_InitValues ***
      !!----------------------------------------------------------------------
      USE nemogcm
      !!----------------------------------------------------------------------
      !
      CALL nemo_init       !* Initializations of each fine grid
      Kbb_a = Nbb; Kmm_a = Nnn; Krhs_a = Nrhs   ! agrif_oce module copies of time level indices
      !
      !                    !* Agrif initialization
      CALL Agrif_InitValues_cont
# if defined key_top
      CALL Agrif_InitValues_cont_top
# endif
# if defined key_si3
      CALL Agrif_InitValues_cont_ice
# endif
      !    
   END SUBROUTINE Agrif_initvalues

   
   SUBROUTINE Agrif_Istate( Kbb, Kmm, Kaa )
      !!----------------------------------------------------------------------
      !!                 *** ROUTINE agrif_istate ***
      !!----------------------------------------------------------------------
      USE domvvl
      USE domain
      USE par_oce
      USE agrif_oce
      USE agrif_oce_interp
      USE oce
      USE lib_mpp
      USE lbclnk
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(in)  :: Kbb, Kmm, Kaa
      INTEGER :: jn
      !!----------------------------------------------------------------------
      IF(lwp) WRITE(numout,*) ' '
      IF(lwp) WRITE(numout,*) 'AGRIF: interp child initial state from parent'
      IF(lwp) WRITE(numout,*) ' '

      l_ini_child           = .TRUE.
      Agrif_SpecialValue    = 0.0_wp
      Agrif_UseSpecialValue = .TRUE.
      uu(:,:,:,:) = 0.0_wp   ;   vv(:,:,:,:) = 0.0_wp   ;   ts(:,:,:,:,:) = 0.0_wp
       
      Krhs_a = Kbb   ;   Kmm_a = Kbb

      ! Brutal fix to pas 1x1 refinment. 
  !    IF(Agrif_Irhox() == 1) THEN
  !       CALL Agrif_Init_Variable(tsini_id, procname=agrif_initts)
  !    ELSE
      CALL Agrif_Init_Variable(tsini_id, procname=interptsn)

  !    ENDIF
! just for VORTEX because Parent velocities can actually be exactly zero
!      Agrif_UseSpecialValue = .FALSE.
      Agrif_UseSpecialValue = ln_spc_dyn
      use_sign_north = .TRUE.
      sign_north = -1.
      CALL Agrif_Init_Variable(uini_id , procname=interpun )
      CALL Agrif_Init_Variable(vini_id , procname=interpvn )
      use_sign_north = .FALSE.

      Agrif_UseSpecialValue = .FALSE.
      l_ini_child           = .FALSE.

      Krhs_a = Kaa   ;   Kmm_a = Kmm

      DO jn = 1, jpts
         ts(:,:,:,jn,Kbb) = ts(:,:,:,jn,Kbb)*tmask(:,:,:)
      END DO
      uu(:,:,:,Kbb) = uu(:,:,:,Kbb) * umask(:,:,:)     
      vv(:,:,:,Kbb) = vv(:,:,:,Kbb) * vmask(:,:,:) 


      CALL lbc_lnk_multi( 'agrif_istate', uu(:,:,:  ,Kbb), 'U', -1.0_wp , vv(:,:,:,Kbb), 'V', -1.0_wp )
      CALL lbc_lnk(       'agrif_istate', ts(:,:,:,:,Kbb), 'T',  1.0_wp )

   END SUBROUTINE Agrif_Istate

   
   SUBROUTINE agrif_declare_var_ini
      !!----------------------------------------------------------------------
      !!                 *** ROUTINE agrif_declare_var_ini ***
      !!----------------------------------------------------------------------
      USE agrif_util
      USE agrif_oce
      USE par_oce
      USE zdf_oce 
      USE oce
      USE dom_oce
      !
      IMPLICIT NONE
      !
      INTEGER :: ind1, ind2, ind3
      INTEGER :: its
      External :: nemo_mapping
      !!----------------------------------------------------------------------

! In case of East-West periodicity, prevent AGRIF interpolation at east and west boundaries
! The procnames will not be called at these boundaries
      IF (jperio == 1) THEN
         CALL Agrif_Set_NearCommonBorderX(.TRUE.)
         CALL Agrif_Set_DistantCommonBorderX(.TRUE.)
      ENDIF

      IF ( .NOT. lk_south ) THEN
         CALL Agrif_Set_NearCommonBorderY(.TRUE.)
      ENDIF

      ! 1. Declaration of the type of variable which have to be interpolated
      !---------------------------------------------------------------------
      ind1 =              nbghostcells
      ind2 = nn_hls + 2 + nbghostcells_x
      ind3 = nn_hls + 2 + nbghostcells_y_s

      CALL agrif_declare_variable((/2,2,0  /),(/ind2  ,ind3,0    /),(/'x','y','N'    /),(/1,1,1  /),(/jpi,jpj,jpk    /),   e3t_id)
      CALL agrif_declare_variable((/2,2    /),(/ind2  ,ind3      /),(/'x','y'        /),(/1,1    /),(/jpi,jpj        /),  mbkt_id)
      CALL agrif_declare_variable((/2,2    /),(/ind2  ,ind3      /),(/'x','y'        /),(/1,1    /),(/jpi,jpj        /),   ht0_id)

      CALL agrif_declare_variable((/1,2    /),(/ind2-1,ind3      /),(/'x','y'        /),(/1,1    /),(/jpi,jpj        /),   e1u_id)
      CALL agrif_declare_variable((/2,1    /),(/ind2  ,ind3-1    /),(/'x','y'        /),(/1,1    /),(/jpi,jpj        /),   e2v_id)
   
      ! Initial or restart velues
      its = jpts+1
      CALL agrif_declare_variable((/2,2,0,0/),(/ind2  ,ind3  ,0,0/),(/'x','y','N','N'/),(/1,1,1,1/),(/jpi,jpj,jpk,its/), tsini_id)
      CALL agrif_declare_variable((/1,2,0,0/),(/ind2-1,ind3  ,0,0/),(/'x','y','N','N'/),(/1,1,1,1/),(/jpi,jpj,jpk,2  /),  uini_id) 
      CALL agrif_declare_variable((/2,1,0,0/),(/ind2  ,ind3-1,0,0/),(/'x','y','N','N'/),(/1,1,1,1/),(/jpi,jpj,jpk,2  /),  vini_id)
      CALL agrif_declare_variable((/2,2    /),(/ind2  ,ind3      /),(/'x','y'        /),(/1,1    /),(/jpi,jpj        /),sshini_id)
      ! 
     
      ! 2. Type of interpolation
      !-------------------------
      CALL Agrif_Set_bcinterp(   e3t_id,interp =AGRIF_constant)

      CALL Agrif_Set_bcinterp(  mbkt_id,interp =AGRIF_constant)
      CALL Agrif_Set_interp  (  mbkt_id,interp =AGRIF_constant)
      CALL Agrif_Set_bcinterp(   ht0_id,interp =AGRIF_constant)
      CALL Agrif_Set_interp  (   ht0_id,interp =AGRIF_constant)

      CALL Agrif_Set_bcinterp(   e1u_id,interp1=Agrif_linear, interp2=AGRIF_ppm    )
      CALL Agrif_Set_bcinterp(   e2v_id,interp1=AGRIF_ppm   , interp2=Agrif_linear )

      ! Initial fields
      CALL Agrif_Set_bcinterp( tsini_id,interp =AGRIF_linear  )
      CALL Agrif_Set_interp  ( tsini_id,interp =AGRIF_linear  )
      CALL Agrif_Set_bcinterp(  uini_id,interp =AGRIF_linear  )
      CALL Agrif_Set_interp  (  uini_id,interp =AGRIF_linear  )
      CALL Agrif_Set_bcinterp(  vini_id,interp =AGRIF_linear  )
      CALL Agrif_Set_interp  (  vini_id,interp =AGRIF_linear  )
      CALL Agrif_Set_bcinterp(sshini_id,interp =AGRIF_linear  )
      CALL Agrif_Set_interp  (sshini_id,interp =AGRIF_linear  )

       ! 3. Location of interpolation
      !-----------------------------
!      CALL Agrif_Set_bc(  e3t_id, (/-nn_sponge_len*Agrif_irhox(),ind1-1/) )  
! JC: check near the boundary only until matching in sponge has been sorted out:
      CALL Agrif_Set_bc(    e3t_id, (/0,ind1-1/) )  

      ! extend the interpolation zone by 1 more point than necessary:
      ! RB check here
      CALL Agrif_Set_bc(   mbkt_id, (/-nn_sponge_len*Agrif_irhox()-2,ind1/) )
      CALL Agrif_Set_bc(    ht0_id, (/-nn_sponge_len*Agrif_irhox()-2,ind1/) )
      
      CALL Agrif_Set_bc(    e1u_id, (/0,ind1-1/) )
      CALL Agrif_Set_bc(    e2v_id, (/0,ind1-1/) )  

      CALL Agrif_Set_bc(  tsini_id, (/0,ind1-1/) ) ! if west,  rhox=3 and nbghost=3: columns 2 to 4
      CALL Agrif_Set_bc(   uini_id, (/0,ind1-1/) ) 
      CALL Agrif_Set_bc(   vini_id, (/0,ind1-1/) )
      CALL Agrif_Set_bc( sshini_id, (/0,ind1-1/) )

      ! 4. Update type
      !--------------- 
# if defined UPD_HIGH
      CALL Agrif_Set_Updatetype(e1u_id,update1 = Agrif_Update_Average       , update2=Agrif_Update_Full_Weighting)
      CALL Agrif_Set_Updatetype(e2v_id,update1 = Agrif_Update_Full_Weighting, update2=Agrif_Update_Average       )
#else
      CALL Agrif_Set_Updatetype(e1u_id,update1 = Agrif_Update_Copy          , update2=Agrif_Update_Average       )
      CALL Agrif_Set_Updatetype(e2v_id,update1 = Agrif_Update_Average       , update2=Agrif_Update_Copy          )
#endif
      
   !   CALL Agrif_Set_ExternalMapping(nemo_mapping)
      !
   END SUBROUTINE agrif_declare_var_ini


   SUBROUTINE Agrif_Init_Domain( Kbb, Kmm, Kaa ) 
      !!----------------------------------------------------------------------
      !!                 *** ROUTINE Agrif_Init_Domain ***
      !!----------------------------------------------------------------------
      USE agrif_oce_update
      USE agrif_oce_interp
      USE agrif_oce_sponge
      USE Agrif_Util
      USE oce 
      USE dom_oce
      USE zdf_oce
      USE nemogcm
      USE agrif_oce
      !
      USE lbclnk
      USE lib_mpp
      USE in_out_manager
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(in) ::  Kbb, Kmm, Kaa
      !
      LOGICAL :: check_namelist
      CHARACTER(len=15) :: cl_check1, cl_check2, cl_check3, cl_check4 
      REAL(wp), DIMENSION(jpi,jpj) ::   zk   ! workspace
      INTEGER :: ji, jj, jk
      !!----------------------------------------------------------------------
    
     ! CALL Agrif_Declare_Var_ini

      IF( agrif_oce_alloc()  > 0 )   CALL ctl_warn('agrif agrif_oce_alloc: allocation of arrays failed')

      ! Build consistent parent bathymetry and number of levels
      ! on the child grid 
      Agrif_UseSpecialValue = .FALSE.
      ht0_parent( :,:) = 0._wp
      mbkt_parent(:,:) = 0
      !
  !    CALL Agrif_Bc_variable(ht0_id ,calledweight=1.,procname=interpht0 )
  !    CALL Agrif_Bc_variable(mbkt_id,calledweight=1.,procname=interpmbkt)
      CALL Agrif_Init_Variable(ht0_id , procname=interpht0 )
      CALL Agrif_Init_Variable(mbkt_id, procname=interpmbkt)
      !
      ! Assume step wise change of bathymetry near interface
      ! TODO: Switch to linear interpolation of bathymetry in the s-coordinate case
      !       and no refinement
      DO_2D( 1, 0, 1, 0 )
         mbku_parent(ji,jj) = MIN( mbkt_parent(ji+1,jj  ), mbkt_parent(ji,jj) )
         mbkv_parent(ji,jj) = MIN( mbkt_parent(ji  ,jj+1), mbkt_parent(ji,jj) )
      END_2D
      IF ( ln_sco.AND.Agrif_Parent(ln_sco) ) THEN 
         DO_2D( 1, 0, 1, 0 )
            hu0_parent(ji,jj) = 0.5_wp * ( ht0_parent(ji,jj)+ht0_parent(ji+1,jj) )
            hv0_parent(ji,jj) = 0.5_wp * ( ht0_parent(ji,jj)+ht0_parent(ji,jj+1) )
         END_2D
      ELSE
         DO_2D( 1, 0, 1, 0 )
            hu0_parent(ji,jj) = MIN( ht0_parent(ji,jj), ht0_parent(ji+1,jj) )
            hv0_parent(ji,jj) = MIN( ht0_parent(ji,jj), ht0_parent(ji,jj+1) )
         END_2D
      ENDIF
      !
      CALL lbc_lnk_multi( 'Agrif_Init_Domain', hu0_parent, 'U', 1.0_wp, hv0_parent, 'V', 1.0_wp )
      DO_2D( 0, 0, 0, 0 )
         zk(ji,jj) = REAL( mbku_parent(ji,jj), wp )
      END_2D
      CALL lbc_lnk( 'Agrif_InitValues_cont', zk, 'U', 1.0_wp )
      mbku_parent(:,:) = MAX( NINT( zk(:,:) ), 1 )
      DO_2D( 0, 0, 0, 0 )
         zk(ji,jj) = REAL( mbkv_parent(ji,jj), wp )
      END_2D
      CALL lbc_lnk( 'Agrif_InitValues_cont', zk, 'V', 1.0_wp )
      mbkv_parent(:,:) = MAX( NINT( zk(:,:) ), 1 )   

      IF ( ln_init_chfrpar ) THEN 
         CALL Agrif_Init_Variable(sshini_id, procname=agrif_initssh)
         CALL lbc_lnk( 'Agrif_Init_Domain', ssh(:,:,Kbb), 'T', 1. )
         DO jk = 1, jpk
               e3t(:,:,jk,Kbb) =  e3t_0(:,:,jk) * ( ht_0(:,:) + ssh(:,:,Kbb)  ) &
                        &             / ( ht_0(:,:) + 1._wp - ssmask(:,:) ) * tmask(:,:,jk)   &
                        &              + e3t_0(:,:,jk) * ( 1._wp - tmask(:,:,jk) )
         END DO
      ENDIF

      ! check if masks and bathymetries match
      IF(ln_chk_bathy) THEN
         Agrif_UseSpecialValue = .FALSE.
         !
         IF(lwp) WRITE(numout,*) ' '
         IF(lwp) WRITE(numout,*) 'AGRIF: Check Bathymetry and masks near bdys. Level: ', Agrif_Level()
         !
         kindic_agr = 0
         IF( .NOT. l_vremap ) THEN
            !
            ! check if tmask and vertical scale factors agree with parent in sponge area:
            CALL Agrif_Bc_variable(e3t_id,calledweight=1.,procname=interpe3t)
            !
         ELSE
            !
            ! In case of vertical interpolation, check only that total depths agree between child and parent:
            DO ji = 1, jpi
               DO jj = 1, jpj
                  IF ((mbkt_parent(ji,jj)/=0).AND.(ABS(ht0_parent(ji,jj)-ht_0(ji,jj))>1.e-3)) kindic_agr = kindic_agr + 1
                  IF ((mbku_parent(ji,jj)/=0).AND.(ABS(hu0_parent(ji,jj)-hu_0(ji,jj))>1.e-3)) kindic_agr = kindic_agr + 1
                  IF ((mbkv_parent(ji,jj)/=0).AND.(ABS(hv0_parent(ji,jj)-hv_0(ji,jj))>1.e-3)) kindic_agr = kindic_agr + 1
               END DO
            END DO

            CALL mpp_sum( 'agrif_user', kindic_agr )
            IF( kindic_agr /= 0 ) THEN
               CALL ctl_stop('==> Child Bathymetry is NOT correct near boundaries.')
            ELSE
               IF(lwp) WRITE(numout,*) '==> Child Bathymetry is ok near boundaries.'
               IF(lwp) WRITE(numout,*) ' '
            ENDIF  
         ENDIF
      ENDIF

      IF( l_vremap ) THEN
      ! Additional constrain that should be removed someday:
         IF ( Agrif_Parent(jpk).GT.jpk ) THEN
            CALL ctl_stop( ' With l_vremap, child grids must have jpk greater or equal to the parent value' )
         ENDIF
      ENDIF
      !
   END SUBROUTINE Agrif_Init_Domain


   SUBROUTINE Agrif_InitValues_cont
      !!----------------------------------------------------------------------
      !!                 *** ROUTINE Agrif_InitValues_cont ***
      !!
      !! ** Purpose ::   Declaration of variables to be interpolated
      !!----------------------------------------------------------------------
      USE agrif_oce_update
      USE agrif_oce_interp
      USE agrif_oce_sponge
      USE Agrif_Util
      USE oce 
      USE dom_oce
      USE zdf_oce
      USE nemogcm
      USE agrif_oce
      !
      USE lbclnk
      USE lib_mpp
      USE in_out_manager
      !
      IMPLICIT NONE
      !
      LOGICAL :: check_namelist
      CHARACTER(len=15) :: cl_check1, cl_check2, cl_check3, cl_check4 
      REAL(wp), DIMENSION(jpi,jpj) ::   zk   ! workspace
      INTEGER :: ji, jj

      ! 1. Declaration of the type of variable which have to be interpolated
      !---------------------------------------------------------------------
      CALL agrif_declare_var

      ! 2. First interpolations of potentially non zero fields
      !-------------------------------------------------------
      Agrif_SpecialValue    = 0._wp
      Agrif_UseSpecialValue = .TRUE.
      CALL Agrif_Bc_variable(       tsn_id,calledweight=1.,procname=interptsn)
      CALL Agrif_Sponge
      tabspongedone_tsn = .FALSE.
      CALL Agrif_Bc_variable(tsn_sponge_id,calledweight=1.,procname=interptsn_sponge)
      ! reset tsa to zero
      ts(:,:,:,:,Krhs_a) = 0._wp

      Agrif_UseSpecialValue = ln_spc_dyn
      use_sign_north = .TRUE.
      sign_north = -1.
      CALL Agrif_Bc_variable(un_interp_id,calledweight=1.,procname=interpun)
      CALL Agrif_Bc_variable(vn_interp_id,calledweight=1.,procname=interpvn)
      tabspongedone_u = .FALSE.
      tabspongedone_v = .FALSE.
      CALL Agrif_Bc_variable(un_sponge_id,calledweight=1.,procname=interpun_sponge)
      tabspongedone_u = .FALSE.
      tabspongedone_v = .FALSE.
      CALL Agrif_Bc_variable(vn_sponge_id,calledweight=1.,procname=interpvn_sponge)
      use_sign_north = .FALSE.
      uu(:,:,:,Krhs_a) = 0._wp
      vv(:,:,:,Krhs_a) = 0._wp

      Agrif_UseSpecialValue = .TRUE.
      CALL Agrif_Bc_variable(sshn_id,calledweight=1., procname=interpsshn )
      hbdy(:,:) = 0._wp
      ssh(:,:,Krhs_a) = 0._wp

      IF ( ln_dynspg_ts ) THEN
         Agrif_UseSpecialValue = ln_spc_dyn
         use_sign_north = .TRUE.
         sign_north = -1.
         CALL Agrif_Bc_variable(        unb_id,calledweight=1.,procname=interpunb )
         CALL Agrif_Bc_variable(        vnb_id,calledweight=1.,procname=interpvnb )
         CALL Agrif_Bc_variable(ub2b_interp_id,calledweight=1.,procname=interpub2b)
         CALL Agrif_Bc_variable(vb2b_interp_id,calledweight=1.,procname=interpvb2b)
         use_sign_north = .FALSE.
         ubdy(:,:) = 0._wp
         vbdy(:,:) = 0._wp
      ENDIF
      Agrif_UseSpecialValue = .FALSE. 

      !-----------------
      check_namelist = .TRUE.

      IF( check_namelist ) THEN 
         ! Check free surface scheme
         IF ( ( Agrif_Parent(ln_dynspg_ts ).AND.ln_dynspg_exp ).OR.&
            & ( Agrif_Parent(ln_dynspg_exp).AND.ln_dynspg_ts ) ) THEN
            WRITE(cl_check1,*)  Agrif_Parent( ln_dynspg_ts )
            WRITE(cl_check2,*)  ln_dynspg_ts
            WRITE(cl_check3,*)  Agrif_Parent( ln_dynspg_exp )
            WRITE(cl_check4,*)  ln_dynspg_exp
            CALL ctl_stop( 'Incompatible free surface scheme between grids' ,  &
                  &               'parent grid ln_dynspg_ts  :'//cl_check1  ,  & 
                  &               'child  grid ln_dynspg_ts  :'//cl_check2  ,  &
                  &               'parent grid ln_dynspg_exp :'//cl_check3  ,  &
                  &               'child  grid ln_dynspg_exp :'//cl_check4  ,  &
                  &               'those logicals should be identical' )                 
            STOP
         ENDIF

         ! Check if identical linear free surface option
         IF ( ( Agrif_Parent(ln_linssh ).AND.(.NOT.ln_linssh )).OR.&
            & ( (.NOT.Agrif_Parent(ln_linssh)).AND.ln_linssh ) ) THEN
            WRITE(cl_check1,*)  Agrif_Parent(ln_linssh )
            WRITE(cl_check2,*)  ln_linssh
            CALL ctl_stop( 'Incompatible linearized fs option between grids',  &
                  &               'parent grid ln_linssh  :'//cl_check1     ,  &
                  &               'child  grid ln_linssh  :'//cl_check2     ,  &
                  &               'those logicals should be identical' )                  
            STOP
         ENDIF
      ENDIF

   END SUBROUTINE Agrif_InitValues_cont

   SUBROUTINE agrif_declare_var
      !!----------------------------------------------------------------------
      !!                 *** ROUTINE agrif_declare_var ***
      !!----------------------------------------------------------------------
      USE agrif_util
      USE agrif_oce
      USE par_oce
      USE zdf_oce 
      USE oce
      !
      IMPLICIT NONE
      !
      INTEGER :: ind1, ind2, ind3
      !!----------------------------------------------------------------------

      ! 1. Declaration of the type of variable which have to be interpolated
      !---------------------------------------------------------------------
      ind1 =              nbghostcells
      ind2 = nn_hls + 2 + nbghostcells_x
      ind3 = nn_hls + 2 + nbghostcells_y_s
# if defined key_vertical
      CALL agrif_declare_variable((/2,2,0,0/),(/ind2,ind3,0,0/),(/'x','y','N','N'/),(/1,1,1,1/),(/jpi,jpj,jpk,jpts+1/),tsn_id)
      CALL agrif_declare_variable((/2,2,0,0/),(/ind2,ind3,0,0/),(/'x','y','N','N'/),(/1,1,1,1/),(/jpi,jpj,jpk,jpts+1/),tsn_sponge_id)
      CALL agrif_declare_variable((/1,2,0,0/),(/ind2-1,ind3,0,0/),(/'x','y','N','N'/),(/1,1,1,1/),(/jpi,jpj,jpk,2/),un_interp_id)
      CALL agrif_declare_variable((/2,1,0,0/),(/ind2,ind3-1,0,0/),(/'x','y','N','N'/),(/1,1,1,1/),(/jpi,jpj,jpk,2/),vn_interp_id)
      CALL agrif_declare_variable((/1,2,0,0/),(/ind2-1,ind3,0,0/),(/'x','y','N','N'/),(/1,1,1,1/),(/jpi,jpj,jpk,2/),un_update_id)
      CALL agrif_declare_variable((/2,1,0,0/),(/ind2,ind3-1,0,0/),(/'x','y','N','N'/),(/1,1,1,1/),(/jpi,jpj,jpk,2/),vn_update_id)
      CALL agrif_declare_variable((/1,2,0,0/),(/ind2-1,ind3,0,0/),(/'x','y','N','N'/),(/1,1,1,1/),(/jpi,jpj,jpk,2/),un_sponge_id)
      CALL agrif_declare_variable((/2,1,0,0/),(/ind2,ind3-1,0,0/),(/'x','y','N','N'/),(/1,1,1,1/),(/jpi,jpj,jpk,2/),vn_sponge_id)
# else
      CALL agrif_declare_variable((/2,2,0,0/),(/ind2,ind3,0,0/),(/'x','y','N','N'/),(/1,1,1,1/),(/jpi,jpj,jpk,jpts/),tsn_id)
      CALL agrif_declare_variable((/2,2,0,0/),(/ind2,ind3,0,0/),(/'x','y','N','N'/),(/1,1,1,1/),(/jpi,jpj,jpk,jpts/),tsn_sponge_id)
      CALL agrif_declare_variable((/1,2,0,0/),(/ind2-1,ind3,0,0/),(/'x','y','N','N'/),(/1,1,1,1/),(/jpi,jpj,jpk,1/),un_interp_id)
      CALL agrif_declare_variable((/2,1,0,0/),(/ind2,ind3-1,0,0/),(/'x','y','N','N'/),(/1,1,1,1/),(/jpi,jpj,jpk,1/),vn_interp_id)
      CALL agrif_declare_variable((/1,2,0,0/),(/ind2-1,ind3,0,0/),(/'x','y','N','N'/),(/1,1,1,1/),(/jpi,jpj,jpk,1/),un_update_id)
      CALL agrif_declare_variable((/2,1,0,0/),(/ind2,ind3-1,0,0/),(/'x','y','N','N'/),(/1,1,1,1/),(/jpi,jpj,jpk,1/),vn_update_id)
      CALL agrif_declare_variable((/1,2,0,0/),(/ind2-1,ind3,0,0/),(/'x','y','N','N'/),(/1,1,1,1/),(/jpi,jpj,jpk,1/),un_sponge_id)
      CALL agrif_declare_variable((/2,1,0,0/),(/ind2,ind3-1,0,0/),(/'x','y','N','N'/),(/1,1,1,1/),(/jpi,jpj,jpk,1/),vn_sponge_id)
# endif
      CALL agrif_declare_variable((/1,2/),(/ind2-1,ind3/),(/'x','y'/),(/1,1/),(/jpi,jpj/),unb_id)
      CALL agrif_declare_variable((/2,1/),(/ind2,ind3-1/),(/'x','y'/),(/1,1/),(/jpi,jpj/),vnb_id)
      CALL agrif_declare_variable((/1,2/),(/ind2-1,ind3/),(/'x','y'/),(/1,1/),(/jpi,jpj/),ub2b_interp_id)
      CALL agrif_declare_variable((/2,1/),(/ind2,ind3-1/),(/'x','y'/),(/1,1/),(/jpi,jpj/),vb2b_interp_id)
      CALL agrif_declare_variable((/1,2/),(/ind2-1,ind3/),(/'x','y'/),(/1,1/),(/jpi,jpj/),ub2b_update_id)
      CALL agrif_declare_variable((/2,1/),(/ind2,ind3-1/),(/'x','y'/),(/1,1/),(/jpi,jpj/),vb2b_update_id)

!      CALL agrif_declare_variable((/2,2/),(/ind2,ind3/),(/'x','y'/),(/1,1/),(/jpi,jpj/),glamt_id)
!      CALL agrif_declare_variable((/2,2/),(/ind2,ind3/),(/'x','y'/),(/1,1/),(/jpi,jpj/),gphit_id)
      CALL agrif_declare_variable((/2,2/),(/ind2,ind3/),(/'x','y'/),(/1,1/),(/jpi,jpj/),sshn_id)


      IF( ln_zdftke.OR.ln_zdfgls ) THEN  ! logical not known at this point
!         CALL agrif_declare_variable((/2,2,0/),(/ind3,ind3,0/),(/'x','y','N'/),(/1,1,1/),(/jpi,jpj,jpk/), en_id)
!         CALL agrif_declare_variable((/2,2,0/),(/ind3,ind3,0/),(/'x','y','N'/),(/1,1,1/),(/jpi,jpj,jpk/),avt_id)
# if defined key_vertical
         CALL agrif_declare_variable((/2,2,0,0/),(/ind2,ind3,0,0/),(/'x','y','N','N'/),(/1,1,1,1/),(/jpi,jpj,jpk,2/),avm_id)
# else
         CALL agrif_declare_variable((/2,2,0,0/),(/ind2,ind3,0,0/),(/'x','y','N','N'/),(/1,1,1,1/),(/jpi,jpj,jpk,1/),avm_id)
# endif
      ENDIF
     
      ! 2. Type of interpolation
      !-------------------------
      CALL Agrif_Set_bcinterp(       tsn_id,interp =AGRIF_linear)
      CALL Agrif_Set_bcinterp( un_interp_id,interp1=Agrif_linear,interp2=AGRIF_ppm   )
      CALL Agrif_Set_bcinterp( vn_interp_id,interp1=AGRIF_ppm   ,interp2=Agrif_linear)

      CALL Agrif_Set_bcinterp( tsn_sponge_id,interp =AGRIF_linear)
      CALL Agrif_Set_bcinterp(  un_sponge_id,interp1=Agrif_linear,interp2=AGRIF_ppm   )
      CALL Agrif_Set_bcinterp(  vn_sponge_id,interp1=AGRIF_ppm   ,interp2=Agrif_linear)

      CALL Agrif_Set_bcinterp(       sshn_id,interp =AGRIF_linear)
      CALL Agrif_Set_bcinterp(        unb_id,interp1=Agrif_linear,interp2=AGRIF_ppm   )
      CALL Agrif_Set_bcinterp(        vnb_id,interp1=AGRIF_ppm   ,interp2=Agrif_linear)
      CALL Agrif_Set_bcinterp(ub2b_interp_id,interp1=Agrif_linear,interp2=AGRIF_ppm   )
      CALL Agrif_Set_bcinterp(vb2b_interp_id,interp1=AGRIF_ppm   ,interp2=Agrif_linear)
!
! > Divergence conserving alternative:
!      CALL Agrif_Set_bcinterp(sshn_id,interp=AGRIF_constant)
!      CALL Agrif_Set_bcinterp(unb_id,interp1=Agrif_linear,interp2=AGRIF_constant)
!      CALL Agrif_Set_bcinterp(vnb_id,interp1=AGRIF_constant,interp2=Agrif_linear)
!      CALL Agrif_Set_bcinterp(ub2b_interp_id,interp1=Agrif_linear,interp2=AGRIF_constant)
!      CALL Agrif_Set_bcinterp(vb2b_interp_id,interp1=AGRIF_constant,interp2=Agrif_linear)
!<

      IF( ln_zdftke.OR.ln_zdfgls )  CALL Agrif_Set_bcinterp( avm_id, interp=AGRIF_linear )
    

!      CALL Agrif_Set_bcinterp(gphit_id,interp=AGRIF_constant)
!      CALL Agrif_Set_bcinterp(glamt_id,interp=AGRIF_constant)

      ! 3. Location of interpolation
      !-----------------------------
      CALL Agrif_Set_bc(       tsn_id, (/0,ind1-1/) ) ! if west,  rhox=3 and nbghost=3: columns 2 to 4
      CALL Agrif_Set_bc( un_interp_id, (/0,ind1-1/) ) 
      CALL Agrif_Set_bc( vn_interp_id, (/0,ind1-1/) )

      CALL Agrif_Set_bc( tsn_sponge_id, (/-nn_sponge_len*Agrif_irhox()-1,0/) )  ! if west,  rhox=3, nn_sponge_len=2 
      CALL Agrif_Set_bc(  un_sponge_id, (/-nn_sponge_len*Agrif_irhox()-1,0/) )  ! and nbghost=3: 
      CALL Agrif_Set_bc(  vn_sponge_id, (/-nn_sponge_len*Agrif_irhox()-1,0/) )  ! columns 4 to 11

      CALL Agrif_Set_bc(        sshn_id, (/0,ind1-1/) )
      CALL Agrif_Set_bc(         unb_id, (/0,ind1-1/) )
      CALL Agrif_Set_bc(         vnb_id, (/0,ind1-1/) )
      CALL Agrif_Set_bc( ub2b_interp_id, (/0,ind1-1/) )
      CALL Agrif_Set_bc( vb2b_interp_id, (/0,ind1-1/) )

      IF( ln_zdftke.OR.ln_zdfgls ) CALL Agrif_Set_bc( avm_id, (/0,ind1/) )
!!$      CALL Agrif_Set_bc(glamt_id, (/0,ind1-1/) )  
!!$      CALL Agrif_Set_bc(gphit_id, (/0,ind1-1/) )  

      ! 4. Update type
      !--------------- 

# if defined UPD_HIGH
      CALL Agrif_Set_Updatetype(      tsn_id,update  = Agrif_Update_Full_Weighting)
      CALL Agrif_Set_Updatetype(un_update_id,update1 = Agrif_Update_Average       , update2 = Agrif_Update_Full_Weighting)
      CALL Agrif_Set_Updatetype(vn_update_id,update1 = Agrif_Update_Full_Weighting, update2 = Agrif_Update_Average       )

      CALL Agrif_Set_Updatetype(ub2b_update_id,update1 = Agrif_Update_Average       , update2 = Agrif_Update_Full_Weighting)
      CALL Agrif_Set_Updatetype(vb2b_update_id,update1 = Agrif_Update_Full_Weighting, update2 = Agrif_Update_Average       )
      CALL Agrif_Set_Updatetype(       sshn_id,update  = Agrif_Update_Full_Weighting)
      CALL Agrif_Set_Updatetype(        e3t_id,update  = Agrif_Update_Full_Weighting)

  !    IF( ln_zdftke.OR.ln_zdfgls ) THEN
!         CALL Agrif_Set_Updatetype( en_id, update = AGRIF_Update_Full_Weighting)
!         CALL Agrif_Set_Updatetype(avt_id, update = AGRIF_Update_Full_Weighting)
!         CALL Agrif_Set_Updatetype(avm_id, update = AGRIF_Update_Full_Weighting)
   !   ENDIF

#else
      CALL Agrif_Set_Updatetype(     tsn_id, update  = AGRIF_Update_Average)
      CALL Agrif_Set_Updatetype(un_update_id,update1 = Agrif_Update_Copy   , update2 = Agrif_Update_Average)
      CALL Agrif_Set_Updatetype(vn_update_id,update1 = Agrif_Update_Average, update2 = Agrif_Update_Copy   )

      CALL Agrif_Set_Updatetype(ub2b_update_id,update1 = Agrif_Update_Copy   , update2 = Agrif_Update_Average)
      CALL Agrif_Set_Updatetype(vb2b_update_id,update1 = Agrif_Update_Average, update2 = Agrif_Update_Copy   )
      CALL Agrif_Set_Updatetype(       sshn_id,update  = AGRIF_Update_Average)
      CALL Agrif_Set_Updatetype(        e3t_id,update  = AGRIF_Update_Average)

 !     IF( ln_zdftke.OR.ln_zdfgls ) THEN
!         CALL Agrif_Set_Updatetype( en_id, update = AGRIF_Update_Average)
!         CALL Agrif_Set_Updatetype(avt_id, update = AGRIF_Update_Average)
!         CALL Agrif_Set_Updatetype(avm_id, update = AGRIF_Update_Average)
 !     ENDIF

#endif
      !
   END SUBROUTINE agrif_declare_var

#if defined key_si3
   SUBROUTINE Agrif_InitValues_cont_ice
      !!----------------------------------------------------------------------
      !!                 *** ROUTINE Agrif_InitValues_cont_ice ***
      !!----------------------------------------------------------------------
      USE Agrif_Util
      USE sbc_oce, ONLY : nn_fsbc  ! clem: necessary otherwise Agrif_Parent(nn_fsbc) = nn_fsbc
      USE ice
      USE agrif_ice
      USE in_out_manager
      USE agrif_ice_interp
      USE lib_mpp
      !
      IMPLICIT NONE
      !
      !!----------------------------------------------------------------------
      ! Controls

      ! clem: For some reason, nn_fsbc(child)/=1 does not work properly (signal can be largely degraded by the agrif zoom)
      !          the run must satisfy CFL=Uice/(dx/dt) < 0.6/nn_fsbc(child)
      !          therefore, if nn_fsbc(child)>1 one must reduce the time-step in proportion to nn_fsbc(child), which is not acceptable
      !       If a solution is found, the following stop could be removed because the rest of the code take nn_fsbc(child) into account     
      IF( nn_fsbc > 1 )  CALL ctl_stop('nn_fsbc(child) must be set to 1 otherwise agrif and sea-ice may not work properly')

      ! stop if rhot * nn_fsbc(parent) /= N * nn_fsbc(child) with N being integer
      IF( MOD( Agrif_irhot() * Agrif_Parent(nn_fsbc), nn_fsbc ) /= 0 )  THEN
         CALL ctl_stop('rhot * nn_fsbc(parent) /= N * nn_fsbc(child), therefore nn_fsbc(child) should be set to 1 or nn_fsbc(parent)')
      ENDIF
      ! First Interpolations (using "after" ice subtime step => nbstep_ice=1)
      !----------------------------------------------------------------------
      nbstep_ice = ( Agrif_irhot() * Agrif_Parent(nn_fsbc) / nn_fsbc ) ! clem: to have calledweight=1 in interp (otherwise the western border of the zoom is wrong)
      CALL agrif_interp_ice('U') ! interpolation of ice velocities
      CALL agrif_interp_ice('V') ! interpolation of ice velocities
      CALL agrif_interp_ice('T') ! interpolation of ice tracers 
      nbstep_ice = 0   
      !
   END SUBROUTINE Agrif_InitValues_cont_ice

   
   SUBROUTINE agrif_declare_var_ice
      !!----------------------------------------------------------------------
      !!                 *** ROUTINE agrif_declare_var_ice ***
      !!----------------------------------------------------------------------
      USE Agrif_Util
      USE ice
      USE par_oce, ONLY : nbghostcells, nbghostcells_x, nbghostcells_y_s
      !
      IMPLICIT NONE
      !
      INTEGER :: ind1, ind2, ind3
      INTEGER :: ipl
      !!----------------------------------------------------------------------
      !
      ! 1. Declaration of the type of variable which have to be interpolated (parent=>child)
      !       agrif_declare_variable(position,1st point index,--,--,dimensions,name)
      !           ex.:  position=> 1,1 = not-centered (in i and j)
      !                            2,2 =     centered (    -     )
      !                 index   => 1,1 = one ghost line
      !                            2,2 = two ghost lines
      !-------------------------------------------------------------------------------------
      ind1 =              nbghostcells
      ind2 = nn_hls + 2 + nbghostcells_x
      ind3 = nn_hls + 2 + nbghostcells_y_s
      ipl = jpl*(9+nlay_s+nlay_i)
      CALL agrif_declare_variable((/2,2,0/),(/ind2,ind3,0/),(/'x','y','N'/),(/1,1,1/),(/jpi,jpj,ipl/),tra_ice_id)
      CALL agrif_declare_variable((/1,2/)  ,(/ind2-1,ind3/),(/'x','y'    /),(/1,1  /),(/jpi,jpj    /),  u_ice_id)
      CALL agrif_declare_variable((/2,1/)  ,(/ind2,ind3-1/),(/'x','y'    /),(/1,1  /),(/jpi,jpj    /),  v_ice_id)

      CALL agrif_declare_variable((/2,2,0/),(/ind3,ind3,0/),(/'x','y','N'/),(/1,1,1/),(/jpi,jpj,ipl/),tra_iceini_id)
      CALL agrif_declare_variable((/1,2/)  ,(/ind2-1,ind3/),(/'x','y'    /),(/1,1  /),(/jpi,jpj    /),  u_iceini_id)
      CALL agrif_declare_variable((/2,1/)  ,(/ind2,ind3-1/),(/'x','y'    /),(/1,1  /),(/jpi,jpj    /),  v_iceini_id)

      ! 2. Set interpolations (normal & tangent to the grid cell for velocities)
      !-----------------------------------
      CALL Agrif_Set_bcinterp(tra_ice_id, interp  = AGRIF_linear)
      CALL Agrif_Set_bcinterp(u_ice_id  , interp1 = Agrif_linear,interp2 = AGRIF_ppm   )
      CALL Agrif_Set_bcinterp(v_ice_id  , interp1 = AGRIF_ppm   ,interp2 = Agrif_linear)

      CALL Agrif_Set_bcinterp(tra_iceini_id, interp  = AGRIF_linear)
      CALL Agrif_Set_interp  (tra_iceini_id, interp  = AGRIF_linear)
      CALL Agrif_Set_bcinterp(u_iceini_id  , interp  = AGRIF_linear  )
      CALL Agrif_Set_interp  (u_iceini_id  , interp  = AGRIF_linear   )
      CALL Agrif_Set_bcinterp(v_iceini_id  , interp  = AGRIF_linear)
      CALL Agrif_Set_interp  (v_iceini_id  , interp  = AGRIF_linear)

      ! 3. Set location of interpolations
      !----------------------------------
      CALL Agrif_Set_bc(tra_ice_id,(/0,ind1/))
      CALL Agrif_Set_bc(u_ice_id  ,(/0,ind1/))
      CALL Agrif_Set_bc(v_ice_id  ,(/0,ind1/))

      CALL Agrif_Set_bc(tra_iceini_id,(/0,ind1/))
      CALL Agrif_Set_bc(u_iceini_id  ,(/0,ind1/))
      CALL Agrif_Set_bc(v_iceini_id  ,(/0,ind1/))

      ! 4. Set update type in case 2 ways (child=>parent) (normal & tangent to the grid cell for velocities)
      !--------------------------------------------------
# if defined UPD_HIGH
      CALL Agrif_Set_Updatetype(tra_ice_id, update  = Agrif_Update_Full_Weighting)
      CALL Agrif_Set_Updatetype(u_ice_id  , update1 = Agrif_Update_Average       , update2 = Agrif_Update_Full_Weighting)
      CALL Agrif_Set_Updatetype(v_ice_id  , update1 = Agrif_Update_Full_Weighting, update2 = Agrif_Update_Average       )
# else
      CALL Agrif_Set_Updatetype(tra_ice_id, update  = AGRIF_Update_Average)
      CALL Agrif_Set_Updatetype(u_ice_id  , update1 = Agrif_Update_Copy   , update2 = Agrif_Update_Average)
      CALL Agrif_Set_Updatetype(v_ice_id  , update1 = Agrif_Update_Average, update2 = Agrif_Update_Copy   )
# endif

   END SUBROUTINE agrif_declare_var_ice
#endif


# if defined key_top
   SUBROUTINE Agrif_InitValues_cont_top
      !!----------------------------------------------------------------------
      !!                 *** ROUTINE Agrif_InitValues_cont_top ***
      !!----------------------------------------------------------------------
      USE Agrif_Util
      USE oce 
      USE dom_oce
      USE nemogcm
      USE par_trc
      USE lib_mpp
      USE trc
      USE in_out_manager
      USE agrif_oce_sponge
      USE agrif_top_update
      USE agrif_top_interp
      USE agrif_top_sponge
      !
      IMPLICIT NONE
      !
      CHARACTER(len=10) :: cl_check1, cl_check2, cl_check3
      LOGICAL :: check_namelist
      !!----------------------------------------------------------------------

      ! 1. Declaration of the type of variable which have to be interpolated
      !---------------------------------------------------------------------
      CALL agrif_declare_var_top

      ! 2. First interpolations of potentially non zero fields
      !-------------------------------------------------------
      Agrif_SpecialValue=0._wp
      Agrif_UseSpecialValue = .TRUE.
      CALL Agrif_Bc_variable(trn_id,calledweight=1.,procname=interptrn)
      Agrif_UseSpecialValue = .FALSE.
      CALL Agrif_Sponge
      tabspongedone_trn = .FALSE.
      CALL Agrif_Bc_variable(trn_sponge_id,calledweight=1.,procname=interptrn_sponge)
      ! reset tsa to zero
      tra(:,:,:,:) = 0._wp

      ! 3. Some controls
      !-----------------
      check_namelist = .TRUE.

      IF( check_namelist ) THEN
         ! Check time steps
         IF( NINT(Agrif_Rhot()) * NINT(rdt) .NE. Agrif_Parent(rdt) ) THEN
            WRITE(cl_check1,*)  Agrif_Parent(rdt)
            WRITE(cl_check2,*)  rdt
            WRITE(cl_check3,*)  rdt*Agrif_Rhot()
            CALL ctl_stop( 'incompatible time step between grids',   &
               &               'parent grid value : '//cl_check1    ,   & 
               &               'child  grid value : '//cl_check2    ,   & 
               &               'value on child grid should be changed to  &
               &               :'//cl_check3  )
         ENDIF

         ! Check run length
         IF( Agrif_IRhot() * (Agrif_Parent(nitend)- &
            Agrif_Parent(nit000)+1) .NE. (nitend-nit000+1) ) THEN
            WRITE(cl_check1,*)  (Agrif_Parent(nit000)-1)*Agrif_IRhot() + 1
            WRITE(cl_check2,*)   Agrif_Parent(nitend)   *Agrif_IRhot()
            CALL ctl_warn( 'incompatible run length between grids'               ,   &
               &              ' nit000 on fine grid will be change to : '//cl_check1,   &
               &              ' nitend on fine grid will be change to : '//cl_check2    )
            nit000 = (Agrif_Parent(nit000)-1)*Agrif_IRhot() + 1
            nitend =  Agrif_Parent(nitend)   *Agrif_IRhot()
         ENDIF
      ENDIF
      !
   END SUBROUTINE Agrif_InitValues_cont_top


   SUBROUTINE agrif_declare_var_top
      !!----------------------------------------------------------------------
      !!                 *** ROUTINE agrif_declare_var_top ***
      !!----------------------------------------------------------------------
      USE agrif_util
      USE agrif_oce
      USE dom_oce
      USE trc
      !!
      IMPLICIT NONE
      !
      INTEGER :: ind1, ind2, ind3
      !!----------------------------------------------------------------------
!RB_CMEMS : declare here init for top      
      ! 1. Declaration of the type of variable which have to be interpolated
      !---------------------------------------------------------------------
      ind1 =              nbghostcells
      ind2 = nn_hls + 2 + nbghostcells_x
      ind3 = nn_hls + 2 + nbghostcells_y_s
# if defined key_vertical
      CALL agrif_declare_variable((/2,2,0,0/),(/ind2,ind3,0,0/),(/'x','y','N','N'/),(/1,1,1,1/),(/jpi,jpj,jpk,jptra+1/),trn_id)
      CALL agrif_declare_variable((/2,2,0,0/),(/ind2,ind3,0,0/),(/'x','y','N','N'/),(/1,1,1,1/),(/jpi,jpj,jpk,jptra+1/),trn_sponge_id)
# else
! LAURENT: STRANGE why (3,3) here ?
      CALL agrif_declare_variable((/2,2,0,0/),(/3,3,0,0/),(/'x','y','N','N'/),(/1,1,1,1/),(/jpi,jpj,jpk,jptra/),trn_id)
      CALL agrif_declare_variable((/2,2,0,0/),(/3,3,0,0/),(/'x','y','N','N'/),(/1,1,1,1/),(/jpi,jpj,jpk,jptra/),trn_sponge_id)
# endif

      ! 2. Type of interpolation
      !-------------------------
      CALL Agrif_Set_bcinterp(trn_id,interp=AGRIF_linear)
      CALL Agrif_Set_bcinterp(trn_sponge_id,interp=AGRIF_linear)

      ! 3. Location of interpolation
      !-----------------------------
      CALL Agrif_Set_bc(trn_id,(/0,ind1-1/))
      CALL Agrif_Set_bc(trn_sponge_id,(/-nn_sponge_len*Agrif_irhox()-1,0/))

      ! 4. Update type
      !--------------- 
# if defined UPD_HIGH
      CALL Agrif_Set_Updatetype(trn_id, update = Agrif_Update_Full_Weighting)
#else
      CALL Agrif_Set_Updatetype(trn_id, update = AGRIF_Update_Average)
#endif
   !
   END SUBROUTINE agrif_declare_var_top
# endif
   

   SUBROUTINE Agrif_detect( kg, ksizex )
      !!----------------------------------------------------------------------
      !!                      *** ROUTINE Agrif_detect ***
      !!----------------------------------------------------------------------
      INTEGER, DIMENSION(2) :: ksizex
      INTEGER, DIMENSION(ksizex(1),ksizex(2)) :: kg 
      !!----------------------------------------------------------------------
      !
      RETURN
      !
   END SUBROUTINE Agrif_detect

   
   SUBROUTINE agrif_nemo_init
      !!----------------------------------------------------------------------
      !!                     *** ROUTINE agrif_init ***
      !!----------------------------------------------------------------------
      USE agrif_oce 
      USE agrif_ice
      USE dom_oce
      USE in_out_manager
      USE lib_mpp
      !
      IMPLICIT NONE
      !
      INTEGER  ::   ios                 ! Local integer output status for namelist read
      NAMELIST/namagrif/ ln_agrif_2way, ln_init_chfrpar, rn_sponge_tra, rn_sponge_dyn, rn_trelax_tra, rn_trelax_dyn, &
                       & ln_spc_dyn, ln_chk_bathy
      !!--------------------------------------------------------------------------------------
      !
      READ  ( numnam_ref, namagrif, IOSTAT = ios, ERR = 901)
901 IF( ios /= 0 )   CALL ctl_nam ( ios , 'namagrif in reference namelist' )
      READ  ( numnam_cfg, namagrif, IOSTAT = ios, ERR = 902 )
902 IF( ios >  0 )   CALL ctl_nam ( ios , 'namagrif in configuration namelist' )
      IF(lwm) WRITE ( numond, namagrif )
      !
      IF(lwp) THEN                    ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'agrif_nemo_init : AGRIF parameters'
         WRITE(numout,*) '~~~~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namagrif : set AGRIF parameters'
         WRITE(numout,*) '      Two way nesting activated ln_agrif_2way         = ', ln_agrif_2way
         WRITE(numout,*) '      child initial state from parent ln_init_chfrpar = ', ln_init_chfrpar
         WRITE(numout,*) '      ad. sponge coeft for tracers      rn_sponge_tra = ', rn_sponge_tra
         WRITE(numout,*) '      ad. sponge coeft for dynamics     rn_sponge_tra = ', rn_sponge_dyn
         WRITE(numout,*) '      ad. time relaxation for tracers   rn_trelax_tra = ', rn_trelax_tra
         WRITE(numout,*) '      ad. time relaxation for dynamics  rn_trelax_dyn = ', rn_trelax_dyn
         WRITE(numout,*) '      use special values for dynamics   ln_spc_dyn    = ', ln_spc_dyn
         WRITE(numout,*) '      check bathymetry                  ln_chk_bathy  = ', ln_chk_bathy
      ENDIF

      lk_west  = .NOT. ( Agrif_Ix() == 1 )
      lk_east  = .NOT. ( Agrif_Ix() + nbcellsx/AGRIF_Irhox() == Agrif_Parent(jpiglo) -1 )
      lk_south = .NOT. ( Agrif_Iy() == 1 )
      lk_north = .NOT. ( Agrif_Iy() + nbcellsy/AGRIF_Irhoy() == Agrif_Parent(jpjglo) -1 )

      !
      ! Set the number of ghost cells according to periodicity
      nbghostcells_x   = nbghostcells
      nbghostcells_y_s = nbghostcells
      nbghostcells_y_n = nbghostcells
      !
      IF(   jperio == 1  )   nbghostcells_x   = 0
      IF( .NOT. lk_south )   nbghostcells_y_s = 0
      ! Some checks
      IF( jpiglo /= nbcellsx + 2 + 2*nn_hls + nbghostcells_x   + nbghostcells_x   )   CALL ctl_stop( 'STOP',    &
         &   'agrif_nemo_init: Agrif children requires jpiglo == nbcellsx + 2 + 2*nn_hls + 2*nbghostcells_x' )
      IF( jpjglo /= nbcellsy + 2 + 2*nn_hls + nbghostcells_y_s + nbghostcells_y_n )   CALL ctl_stop( 'STOP',    &
         &   'agrif_nemo_init: Agrif children requires jpjglo == nbcellsy + 2 + 2*nn_hls + nbghostcells_y_s + nbghostcells_y_n' )
      IF( ln_use_jattr )   CALL ctl_stop( 'STOP', 'agrif_nemo_init:Agrif children requires ln_use_jattr = .false. ' )
      !
   END SUBROUTINE agrif_nemo_init

   
# if defined key_mpp_mpi
   SUBROUTINE Agrif_InvLoc( indloc, nprocloc, i, indglob )
      !!----------------------------------------------------------------------
      !!                     *** ROUTINE Agrif_InvLoc ***
      !!----------------------------------------------------------------------
      USE dom_oce
      !!
      IMPLICIT NONE
      !
      INTEGER :: indglob, indloc, nprocloc, i
      !!----------------------------------------------------------------------
      !
      SELECT CASE( i )
      CASE(1)        ;   indglob = mig(indloc)
      CASE(2)        ;   indglob = mjg(indloc)
      CASE DEFAULT   ;   indglob = indloc
      END SELECT
      !
   END SUBROUTINE Agrif_InvLoc

   
   SUBROUTINE Agrif_get_proc_info( imin, imax, jmin, jmax )
      !!----------------------------------------------------------------------
      !!                 *** ROUTINE Agrif_get_proc_info ***
      !!----------------------------------------------------------------------
      USE par_oce
      !!
      IMPLICIT NONE
      !
      INTEGER, INTENT(out) :: imin, imax
      INTEGER, INTENT(out) :: jmin, jmax
      !!----------------------------------------------------------------------
      !
      imin = mig( 1 )
      jmin = mjg( 1 )
      imax = mig(jpi)
      jmax = mjg(jpj)
      ! 
   END SUBROUTINE Agrif_get_proc_info

   
   SUBROUTINE Agrif_estimate_parallel_cost(imin, imax,jmin, jmax, nbprocs, grid_cost)
      !!----------------------------------------------------------------------
      !!                 *** ROUTINE Agrif_estimate_parallel_cost ***
      !!----------------------------------------------------------------------
      USE par_oce
      !!
      IMPLICIT NONE
      !
      INTEGER,  INTENT(in)  :: imin, imax
      INTEGER,  INTENT(in)  :: jmin, jmax
      INTEGER,  INTENT(in)  :: nbprocs
      REAL(wp), INTENT(out) :: grid_cost
      !!----------------------------------------------------------------------
      !
      grid_cost = REAL(imax-imin+1,wp)*REAL(jmax-jmin+1,wp) / REAL(nbprocs,wp)
      !
   END SUBROUTINE Agrif_estimate_parallel_cost

# endif

   SUBROUTINE nemo_mapping(ndim,ptx,pty,bounds,bounds_chunks,correction_required,nb_chunks)
      !!----------------------------------------------------------------------
      !!                   *** ROUTINE Nemo_mapping ***
      !!----------------------------------------------------------------------
      USE dom_oce
      !!
      IMPLICIT NONE
      !
      INTEGER :: ndim
      INTEGER :: ptx, pty
      INTEGER, DIMENSION(ndim,2,2) :: bounds
      INTEGER, DIMENSION(:,:,:,:), ALLOCATABLE :: bounds_chunks
      LOGICAL, DIMENSION(:), ALLOCATABLE :: correction_required
      INTEGER :: nb_chunks
      !
      INTEGER :: i

      IF (agrif_debug_interp) THEN
         DO i=1,ndim
            WRITE(*,*) 'direction = ',i,bounds(i,1,2),bounds(i,2,2)
         ENDDO
      ENDIF

      IF( bounds(2,2,2) > jpjglo) THEN
         IF( bounds(2,1,2) <=jpjglo) THEN
            nb_chunks = 2
            ALLOCATE(bounds_chunks(nb_chunks,ndim,2,2))
            ALLOCATE(correction_required(nb_chunks))
            DO i = 1,nb_chunks
               bounds_chunks(i,:,:,:) = bounds
            END DO
        
      ! FIRST CHUNCK (for j<=jpjglo)   

      ! Original indices
            bounds_chunks(1,1,1,1) = bounds(1,1,2)
            bounds_chunks(1,1,2,1) = bounds(1,2,2)
            bounds_chunks(1,2,1,1) = bounds(2,1,2)
            bounds_chunks(1,2,2,1) = jpjglo

            bounds_chunks(1,1,1,2) = bounds(1,1,2)
            bounds_chunks(1,1,2,2) = bounds(1,2,2)
            bounds_chunks(1,2,1,2) = bounds(2,1,2)
            bounds_chunks(1,2,2,2) = jpjglo

      ! Correction required or not
            correction_required(1)=.FALSE.
       
      ! SECOND CHUNCK (for j>jpjglo)

      ! Original indices
            bounds_chunks(2,1,1,1) = bounds(1,1,2)
            bounds_chunks(2,1,2,1) = bounds(1,2,2)
            bounds_chunks(2,2,1,1) = jpjglo-2
            bounds_chunks(2,2,2,1) = bounds(2,2,2)

      ! Where to find them
      ! We use the relation TAB(ji,jj)=TAB(jpiglo-ji+2,jpjglo-2-(jj-jpjglo))

            IF( ptx == 2) THEN ! T, V points
               bounds_chunks(2,1,1,2) = jpiglo-bounds(1,2,2)+2
               bounds_chunks(2,1,2,2) = jpiglo-bounds(1,1,2)+2
            ELSE ! U, F points
               bounds_chunks(2,1,1,2) = jpiglo-bounds(1,2,2)+1
               bounds_chunks(2,1,2,2) = jpiglo-bounds(1,1,2)+1       
            ENDIF

            IF( pty == 2) THEN ! T, U points
               bounds_chunks(2,2,1,2) = jpjglo-2-(bounds(2,2,2) -jpjglo)
               bounds_chunks(2,2,2,2) = jpjglo-2-(jpjglo-2      -jpjglo)
            ELSE ! V, F points
               bounds_chunks(2,2,1,2) = jpjglo-3-(bounds(2,2,2) -jpjglo)
               bounds_chunks(2,2,2,2) = jpjglo-3-(jpjglo-2      -jpjglo)
            ENDIF
      ! Correction required or not
            correction_required(2)=.TRUE.

         ELSE
            nb_chunks = 1
            ALLOCATE(bounds_chunks(nb_chunks,ndim,2,2))
            ALLOCATE(correction_required(nb_chunks))
            DO i=1,nb_chunks
               bounds_chunks(i,:,:,:) = bounds
            END DO

            bounds_chunks(1,1,1,1) = bounds(1,1,2)
            bounds_chunks(1,1,2,1) = bounds(1,2,2)
            bounds_chunks(1,2,1,1) = bounds(2,1,2)
            bounds_chunks(1,2,2,1) = bounds(2,2,2)

            bounds_chunks(1,1,1,2) = jpiglo-bounds(1,2,2)+2
            bounds_chunks(1,1,2,2) = jpiglo-bounds(1,1,2)+2

            bounds_chunks(1,2,1,2) = jpjglo-2-(bounds(2,2,2)-jpjglo)
            bounds_chunks(1,2,2,2) = jpjglo-2-(bounds(2,1,2)-jpjglo)

            IF( ptx == 2) THEN ! T, V points
               bounds_chunks(1,1,1,2) = jpiglo-bounds(1,2,2)+2
               bounds_chunks(1,1,2,2) = jpiglo-bounds(1,1,2)+2
            ELSE ! U, F points
               bounds_chunks(1,1,1,2) = jpiglo-bounds(1,2,2)+1
               bounds_chunks(1,1,2,2) = jpiglo-bounds(1,1,2)+1    	
            ENDIF

            IF (pty == 2) THEN ! T, U points
               bounds_chunks(1,2,1,2) = jpjglo-2-(bounds(2,2,2) -jpjglo)
               bounds_chunks(1,2,2,2) = jpjglo-2-(bounds(2,1,2) -jpjglo)
            ELSE ! V, F points
               bounds_chunks(1,2,1,2) = jpjglo-3-(bounds(2,2,2) -jpjglo)
               bounds_chunks(1,2,2,2) = jpjglo-3-(bounds(2,1,2) -jpjglo)
            ENDIF

            correction_required(1)=.TRUE.          
         ENDIF

      ELSE IF (bounds(1,1,2) < 1) THEN
         IF (bounds(1,2,2) > 0) THEN
            nb_chunks = 2
            ALLOCATE(correction_required(nb_chunks))
            correction_required=.FALSE.
            ALLOCATE(bounds_chunks(nb_chunks,ndim,2,2))
            DO i=1,nb_chunks
               bounds_chunks(i,:,:,:) = bounds
            END DO
              
            bounds_chunks(1,1,1,2) = bounds(1,1,2)+jpiglo-2
            bounds_chunks(1,1,2,2) = 1+jpiglo-2
          
            bounds_chunks(1,1,1,1) = bounds(1,1,2)
            bounds_chunks(1,1,2,1) = 1
       
            bounds_chunks(2,1,1,2) = 2
            bounds_chunks(2,1,2,2) = bounds(1,2,2)
          
            bounds_chunks(2,1,1,1) = 2
            bounds_chunks(2,1,2,1) = bounds(1,2,2)

         ELSE
            nb_chunks = 1
            ALLOCATE(correction_required(nb_chunks))
            correction_required=.FALSE.
            ALLOCATE(bounds_chunks(nb_chunks,ndim,2,2))
            DO i=1,nb_chunks
               bounds_chunks(i,:,:,:) = bounds
            END DO    
            bounds_chunks(1,1,1,2) = bounds(1,1,2)+jpiglo-2
            bounds_chunks(1,1,2,2) = bounds(1,2,2)+jpiglo-2
          
            bounds_chunks(1,1,1,1) = bounds(1,1,2)
           bounds_chunks(1,1,2,1) = bounds(1,2,2)
         ENDIF
      ELSE
         nb_chunks=1  
         ALLOCATE(correction_required(nb_chunks))
         correction_required=.FALSE.
         ALLOCATE(bounds_chunks(nb_chunks,ndim,2,2))
         DO i=1,nb_chunks
            bounds_chunks(i,:,:,:) = bounds
         END DO
         bounds_chunks(1,1,1,2) = bounds(1,1,2)
         bounds_chunks(1,1,2,2) = bounds(1,2,2)
         bounds_chunks(1,2,1,2) = bounds(2,1,2)
         bounds_chunks(1,2,2,2) = bounds(2,2,2)
          
         bounds_chunks(1,1,1,1) = bounds(1,1,2)
         bounds_chunks(1,1,2,1) = bounds(1,2,2)
         bounds_chunks(1,2,1,1) = bounds(2,1,2)
         bounds_chunks(1,2,2,1) = bounds(2,2,2)              
      ENDIF
        
   END SUBROUTINE nemo_mapping

   FUNCTION agrif_external_switch_index(ptx,pty,i1,isens)

      USE dom_oce
      !
      IMPLICIT NONE

      INTEGER :: ptx, pty, i1, isens
      INTEGER :: agrif_external_switch_index
      !!----------------------------------------------------------------------

      IF( isens == 1 ) THEN
         IF( ptx == 2 ) THEN ! T, V points
            agrif_external_switch_index = jpiglo-i1+2
         ELSE ! U, F points
            agrif_external_switch_index = jpiglo-i1+1      
         ENDIF
      ELSE IF( isens ==2 ) THEN
         IF ( pty == 2 ) THEN ! T, U points
            agrif_external_switch_index = jpjglo-2-(i1 -jpjglo)
         ELSE ! V, F points
            agrif_external_switch_index = jpjglo-3-(i1 -jpjglo)
         ENDIF
      ENDIF

   END FUNCTION agrif_external_switch_index

   SUBROUTINE Correct_field(tab2d,i1,i2,j1,j2)
      !!----------------------------------------------------------------------
      !!                   *** ROUTINE Correct_field ***
      !!----------------------------------------------------------------------
      USE dom_oce
      USE agrif_oce
      !
      IMPLICIT NONE
      !
      INTEGER :: i1,i2,j1,j2
      REAL(wp), DIMENSION(i1:i2,j1:j2) :: tab2d
      !
      INTEGER :: i,j
      REAL(wp), DIMENSION(i1:i2,j1:j2) :: tab2dtemp
      !!----------------------------------------------------------------------

      tab2dtemp = tab2d

      IF( .NOT. use_sign_north ) THEN
         DO j=j1,j2
            DO i=i1,i2
               tab2d(i,j)=tab2dtemp(i2-(i-i1),j2-(j-j1))
            END DO
         END DO
      ELSE
         DO j=j1,j2
            DO i=i1,i2
               tab2d(i,j)=sign_north * tab2dtemp(i2-(i-i1),j2-(j-j1))
            END DO
         END DO
      ENDIF

   END SUBROUTINE Correct_field

#else
   SUBROUTINE Subcalledbyagrif
      !!----------------------------------------------------------------------
      !!                   *** ROUTINE Subcalledbyagrif ***
      !!----------------------------------------------------------------------
      WRITE(*,*) 'Impossible to be here'
   END SUBROUTINE Subcalledbyagrif
#endif
