
MODULE p3plugin_main_routines
  USE comin_plugin_interface,  ONLY : comin_parallel_get_host_mpi_rank, comin_descrdata_get_cell_indices

  USE p3plugin_types,          ONLY : t_dyn_vars_3dptr, t_mp_vars_3dptr, t_p3_vars_3dptr,                   &
    &                                 t_icon_tracer_3dptr, t_p3_tracer_3dptr
  USE p3plugin_global_vars,    ONLY : fastphystep, n_icecat, l3mom_ice, lliqfrac, dtime, p_global, p_patch, &
    &                                 dyn_vars, mp_vars, p3_vars, icon_tracer, p3_tracer
  USE p3plugin_utils,          ONLY : print_global_max

  USE microphy_p3,             ONLY : p3_main

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: p3_main_wrapper, set_reff_before_rad

  INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND(12,307)

CONTAINS

  ! --------------------------------------------------------------------
  ! ComIn constructor to run prepare and run P3 wrapper
  ! --------------------------------------------------------------------
  SUBROUTINE p3_main_wrapper()  BIND(C)

    INTEGER :: rank, i_icecat, jg, jk, jb
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end
    INTEGER :: n_diag_2d, n_diag_3d

    TYPE(t_dyn_vars_3dptr)    :: dyn_vars_3d
    TYPE(t_mp_vars_3dptr)     :: mp_vars_3d
    TYPE(t_p3_vars_3dptr)     :: p3_vars_3d(n_icecat)
    TYPE(t_icon_tracer_3dptr) :: icon_tracer_3d
    TYPE(t_p3_tracer_3dptr)   :: p3_tracer_3d(n_icecat)

    REAL :: dz_3d(p_global%nproma, p_patch%nlev, p_patch%cells%nblks)
    REAL :: temp_3d(p_global%nproma, p_patch%nlev, p_patch%cells%nblks)
    REAL :: temp_before_phys_3d(p_global%nproma, p_patch%nlev, p_patch%cells%nblks)
    REAL :: pres_3d(p_global%nproma, p_patch%nlev, p_patch%cells%nblks)
    REAL :: theta_3d(p_global%nproma, p_patch%nlev, p_patch%cells%nblks)
    REAL :: theta_old_3d(p_global%nproma, p_patch%nlev, p_patch%cells%nblks)
    REAL :: w_fl_3d(p_global%nproma, p_patch%nlev, p_patch%cells%nblks)

    REAL :: qv_3d(p_global%nproma, p_patch%nlev, p_patch%cells%nblks)
    REAL :: qv_old_3d(p_global%nproma, p_patch%nlev, p_patch%cells%nblks)
    REAL :: qc_3d(p_global%nproma, p_patch%nlev, p_patch%cells%nblks)
    REAL :: qr_3d(p_global%nproma, p_patch%nlev, p_patch%cells%nblks)
    REAL :: qnc_3d(p_global%nproma, p_patch%nlev, p_patch%cells%nblks)
    REAL :: qnr_3d(p_global%nproma, p_patch%nlev, p_patch%cells%nblks)

    REAL :: qitot_4d(p_global%nproma, p_patch%nlev, p_patch%cells%nblks, n_icecat)
    REAL :: qnitot_4d(p_global%nproma, p_patch%nlev, p_patch%cells%nblks, n_icecat)
    REAL :: qirim_4d(p_global%nproma, p_patch%nlev, p_patch%cells%nblks, n_icecat)
    REAL :: birim_4d(p_global%nproma, p_patch%nlev, p_patch%cells%nblks, n_icecat)
    REAL :: qzitot_4d(p_global%nproma, p_patch%nlev, p_patch%cells%nblks, n_icecat)
    REAL :: qiliq_4d(p_global%nproma, p_patch%nlev, p_patch%cells%nblks, n_icecat)

    REAL :: ssat(p_global%nproma, p_patch%nlev, p_patch%cells%nblks)
    REAL :: scf_3d(p_global%nproma, p_patch%nlev, p_patch%cells%nblks)
    REAL :: prt_liq_2d(p_global%nproma, p_patch%cells%nblks)
    REAL :: prt_sol_2d(p_global%nproma, p_patch%cells%nblks)
    REAL :: diag_ze_3d(p_global%nproma, p_patch%nlev, p_patch%cells%nblks)
    REAL :: diag_effc_3d(p_global%nproma, p_patch%nlev, p_patch%cells%nblks)
    REAL :: diag_effi_4d(p_global%nproma, p_patch%nlev, p_patch%cells%nblks, n_icecat)
    REAL :: diag_vmi_4d(p_global%nproma, p_patch%nlev, p_patch%cells%nblks, n_icecat)
    REAL :: diag_di_4d(p_global%nproma, p_patch%nlev, p_patch%cells%nblks, n_icecat)
    REAL :: diag_rhoi_4d(p_global%nproma, p_patch%nlev, p_patch%cells%nblks, n_icecat)
    REAL :: diag_dhmax_4d(p_global%nproma, p_patch%nlev, p_patch%cells%nblks, n_icecat)
    REAL :: dummysum(p_global%nproma, p_patch%nlev, p_patch%cells%nblks)
    REAL, ALLOCATABLE :: diag_2d(:,:,:)
    REAL, ALLOCATABLE :: diag_3d(:,:,:,:)

    REAL :: Rd, c_pd
    REAL :: clbfact_dep = 1.0
    REAL :: clbfact_sub = 1.0
    REAL :: scpf_pfrac = 1.0
    REAL :: scpf_resfact = 1.0
    REAL :: timer(20)
    LOGICAL  :: l2mom_clouddrops = .TRUE.
    LOGICAL  :: lscpf = .FALSE.
    LOGICAL  :: debug_on = .FALSE.
    CHARACTER(len=16) :: model = 'ICON'
    CHARACTER(len=20) :: timer_description(20)

    rank = comin_parallel_get_host_mpi_rank()
    IF (rank == 0) WRITE (0,*) 'run P3 microphysics'
    !IF (rank == 0) WRITE (0,'(a,i)') ' fastphystep', fastphystep
    !IF (rank == 0) WRITE (0,'(a,i)') ' ncells', p_patch%cells%ncells
    !IF (rank == 0) WRITE (0,'(a,i)') ' ncells_global', p_patch%cells%ncells_global


    CALL dyn_vars%dz%to_3d(dyn_vars_3d%dz)
    CALL dyn_vars%temp%to_3d(dyn_vars_3d%temp)
    CALL dyn_vars%rho%to_3d(dyn_vars_3d%rho)
    CALL dyn_vars%pres%to_3d(dyn_vars_3d%pres)
    CALL dyn_vars%exner%to_3d(dyn_vars_3d%exner)
    CALL dyn_vars%theta_old%to_3d(dyn_vars_3d%theta_old)
    CALL dyn_vars%ddt_temp_phys%to_3d(dyn_vars_3d%ddt_temp_phys)
    CALL dyn_vars%w_hl%to_3d(dyn_vars_3d%w_hl)

    CALL mp_vars%dmean_c%to_3d(mp_vars_3d%dmean_c)
    CALL mp_vars%dmean_r%to_3d(mp_vars_3d%dmean_r)
    CALL mp_vars%deff_c%to_3d(mp_vars_3d%deff_c)
    CALL mp_vars%deff_i%to_3d(mp_vars_3d%deff_i)
    CALL mp_vars%dhmax%to_3d(mp_vars_3d%dhmax)
    CALL mp_vars%dhmax_ground%to_3d(mp_vars_3d%dhmax_ground)
    CALL mp_vars%ze_p3%to_3d(mp_vars_3d%ze_p3)

    CALL mp_vars%prec_gsp_rate%to_3d(mp_vars_3d%prec_gsp_rate)
    CALL mp_vars%rain_gsp_rate%to_3d(mp_vars_3d%rain_gsp_rate)
    CALL mp_vars%snow_gsp_rate%to_3d(mp_vars_3d%snow_gsp_rate)
    CALL mp_vars%prec_gsp%to_3d(mp_vars_3d%prec_gsp)
    CALL mp_vars%prec_gsp_d%to_3d(mp_vars_3d%prec_gsp_d)
    CALL mp_vars%rain_gsp%to_3d(mp_vars_3d%rain_gsp)
    CALL mp_vars%snow_gsp%to_3d(mp_vars_3d%snow_gsp)

    CALL icon_tracer%qv%to_3d(icon_tracer_3d%qv)
    CALL icon_tracer%qv_old%to_3d(icon_tracer_3d%qv_old)
    CALL icon_tracer%qc%to_3d(icon_tracer_3d%qc)
    CALL icon_tracer%qr%to_3d(icon_tracer_3d%qr)
    CALL icon_tracer%qi%to_3d(icon_tracer_3d%qi)
    CALL icon_tracer%qnc%to_3d(icon_tracer_3d%qnc)
    CALL icon_tracer%qnr%to_3d(icon_tracer_3d%qnr)

    DO i_icecat = 1, n_icecat
      CALL p3_tracer(i_icecat)%qitot%to_3d(p3_tracer_3d(i_icecat)%qitot)
      CALL p3_tracer(i_icecat)%qnitot%to_3d(p3_tracer_3d(i_icecat)%qnitot)
      CALL p3_tracer(i_icecat)%qirim%to_3d(p3_tracer_3d(i_icecat)%qirim)
      CALL p3_tracer(i_icecat)%birim%to_3d(p3_tracer_3d(i_icecat)%birim)
      IF (l3mom_ice) THEN
        CALL p3_tracer(i_icecat)%qzitot%to_3d(p3_tracer_3d(i_icecat)%qzitot)
      ENDIF
      IF (lliqfrac) THEN
        CALL p3_tracer(i_icecat)%qiliq%to_3d(p3_tracer_3d(i_icecat)%qiliq)
      ENDIF

      CALL p3_vars(i_icecat)%dmean_i%to_3d(p3_vars_3d(i_icecat)%dmean_i)
      CALL p3_vars(i_icecat)%deff_i%to_3d(p3_vars_3d(i_icecat)%deff_i)
      CALL p3_vars(i_icecat)%rho_i%to_3d(p3_vars_3d(i_icecat)%rho_i)
      CALL p3_vars(i_icecat)%vm_i%to_3d(p3_vars_3d(i_icecat)%vm_i)
    END DO


    n_diag_2d = 1  ! not used
    n_diag_3d = 2  ! diag_3d contains dmean_c, dmean_r
    ALLOCATE(diag_2d(p_global%nproma, p_patch%cells%nblks, n_diag_2d))
    ALLOCATE(diag_3d(p_global%nproma, p_patch%nlev, p_patch%cells%nblks, n_diag_3d))

    Rd   = 287.04
    c_pd = Rd * 3.5
    ssat(:,:,:)       = 0.0
    prt_liq_2d(:,:)   = 0.0
    prt_sol_2d(:,:)   = 0.0
    diag_ze_3d(:,:,:) = 0.0
    diag_effc_3d(:,:,:) = 0.0
    diag_effi_4d(:,:,:,:) = 0.0
    diag_vmi_4d(:,:,:,:) = 0.0
    diag_di_4d(:,:,:,:) = 0.0
    diag_rhoi_4d(:,:,:,:) = 0.0
    diag_dhmax_4d(:,:,:,:) = 0.0
    diag_3d(:,:,:,:)  = 0.0

    ! calc theta, pres and w on full levels
    dyn_vars_3d%pres = dyn_vars_3d%rho * dyn_vars_3d%temp * Rd
    dyn_vars_3d%exner = (dyn_vars_3d%pres*1.0e-5) ** (Rd/c_pd)
    temp_before_phys_3d = dyn_vars_3d%temp
    theta_3d = dyn_vars_3d%temp / dyn_vars_3d%exner
    DO jk = 1, p_patch%nlev
        w_fl_3d(:,jk,:) = 0.5 * (dyn_vars_3d%w_hl(:,jk,:) + dyn_vars_3d%w_hl(:,jk+1,:))
    END DO

    ! initialize variables in first time step
    IF (fastphystep == 1) THEN
      dyn_vars_3d%theta_old = theta_3d
      icon_tracer_3d%qv_old = icon_tracer_3d%qv
    ENDIF

    IF ( mod((fastphystep-1)*dtime, 60.0) == 0.0 ) THEN
      mp_vars_3d%dhmax_ground = 0.0
    ENDIF


    dz_3d        = dyn_vars_3d%dz
    pres_3d      = dyn_vars_3d%pres
    theta_old_3d = dyn_vars_3d%theta_old

    qv_3d        = icon_tracer_3d%qv
    qv_old_3d    = icon_tracer_3d%qv_old
    qc_3d        = icon_tracer_3d%qc
    qr_3d        = icon_tracer_3d%qr
    qnc_3d       = icon_tracer_3d%qnc
    qnr_3d       = icon_tracer_3d%qnr

    DO i_icecat = 1, n_icecat
      qitot_4d(:,:,:,i_icecat) = p3_tracer_3d(i_icecat)%qitot
      qnitot_4d(:,:,:,i_icecat) = p3_tracer_3d(i_icecat)%qnitot
      qirim_4d(:,:,:,i_icecat) = p3_tracer_3d(i_icecat)%qirim
      birim_4d(:,:,:,i_icecat) = p3_tracer_3d(i_icecat)%birim
      IF (l3mom_ice) THEN
        qzitot_4d(:,:,:,i_icecat) = p3_tracer_3d(i_icecat)%qzitot
      ENDIF
      IF (lliqfrac) THEN
        qiliq_4d(:,:,:,i_icecat) = p3_tracer_3d(i_icecat)%qiliq
      ENDIF
    END DO


    ! block loop and call P3
    jg = 1
    rl_start = p_global%grf_bdywidth_c + 1
    rl_end   = p_global%min_rlcell_int

    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

    !IF (rank == 0) WRITE (0,*) 'i_startblk, i_endblk:', i_startblk, i_endblk

!!$OMP PARALLEL PRIVATE(jb,i_startidx,i_endidx) NUM_THREADS(1)
!!$OMP DO SCHEDULE(DYNAMIC, 1)

    DO jb = i_startblk, i_endblk
      !IF (rank == 0) WRITE (0,*) 'in do loop, jb=', jb
      CALL comin_descrdata_get_cell_indices(jg, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)

      CALL p3_main(qc             = qc_3d(i_startidx:i_endidx, 1:p_patch%nlev, jb),                       &
                   nc             = qnc_3d(i_startidx:i_endidx, 1:p_patch%nlev, jb),                      &
                   qr             = qr_3d(i_startidx:i_endidx, 1:p_patch%nlev, jb),                       &
                   nr             = qnr_3d(i_startidx:i_endidx, 1:p_patch%nlev, jb),                      &
                   th_old         = theta_old_3d(i_startidx:i_endidx, 1:p_patch%nlev, jb),                &
                   th             = theta_3d(i_startidx:i_endidx, 1:p_patch%nlev, jb),                    &
                   qv_old         = qv_old_3d(i_startidx:i_endidx, 1:p_patch%nlev, jb),                   &
                   qv             = qv_3d(i_startidx:i_endidx, 1:p_patch%nlev, jb),                       &
                   dt             = dtime,                                                                &
                   qitot          = qitot_4d(i_startidx:i_endidx, 1:p_patch%nlev, jb, 1:n_icecat),        &
                   qirim          = qirim_4d(i_startidx:i_endidx, 1:p_patch%nlev, jb, 1:n_icecat),        &
                   qiliq          = qiliq_4d(i_startidx:i_endidx, 1:p_patch%nlev, jb, 1:n_icecat),        &
                   nitot          = qnitot_4d(i_startidx:i_endidx, 1:p_patch%nlev, jb, 1:n_icecat),       &
                   birim          = birim_4d(i_startidx:i_endidx, 1:p_patch%nlev, jb, 1:n_icecat),        &
                   zitot          = qzitot_4d(i_startidx:i_endidx, 1:p_patch%nlev, jb, 1:n_icecat),       &
                   ssat           = ssat(i_startidx:i_endidx, 1:p_patch%nlev, jb),                        &
                   uzpl           = w_fl_3d(i_startidx:i_endidx, 1:p_patch%nlev, jb),                     &
                   pres           = pres_3d(i_startidx:i_endidx, 1:p_patch%nlev, jb),                     &
                   dzq            = dz_3d(i_startidx:i_endidx, 1:p_patch%nlev, jb),                       &
                   it             = fastphystep,                                                          &
                   prt_liq        = prt_liq_2d(i_startidx:i_endidx, jb),                                  &
                   prt_sol        = prt_sol_2d(i_startidx:i_endidx, jb),                                  &
                   its            = i_startidx,                                                           &
                   ite            = i_endidx,                                                             &
                   kts            = 1,                                                                    &
                   kte            = p_patch%nlev,                                                         &
                   nCat           = n_icecat,                                                             &
                   diag_ze        = diag_ze_3d(i_startidx:i_endidx, 1:p_patch%nlev, jb),                  &
                   diag_effc      = diag_effc_3d(i_startidx:i_endidx, 1:p_patch%nlev, jb),                &
                   diag_effi      = diag_effi_4d(i_startidx:i_endidx, 1:p_patch%nlev, jb, 1:n_icecat),    &
                   diag_vmi       = diag_vmi_4d(i_startidx:i_endidx, 1:p_patch%nlev, jb, 1:n_icecat),     &
                   diag_di        = diag_di_4d(i_startidx:i_endidx, 1:p_patch%nlev, jb, 1:n_icecat),      &
                   diag_rhoi      = diag_rhoi_4d(i_startidx:i_endidx, 1:p_patch%nlev, jb, 1:n_icecat),    &
                   n_diag_2d      = n_diag_2d,                                                            &
                   diag_2d        = diag_2d(i_startidx:i_endidx, jb, 1:n_diag_2d),                        &
                   n_diag_3d      = n_diag_3d,                                                            &
                   diag_3d        = diag_3d(i_startidx:i_endidx, 1:p_patch%nlev, jb, 1:n_diag_3d),        &
                   log_predictNc  = l2mom_clouddrops,                                                     &
                   model          = model,                                                                &
                   clbfact_dep    = clbfact_dep,                                                          &
                   clbfact_sub    = clbfact_sub,                                                          &
                   debug_on       = debug_on,                                                             &
                   scpf_on        = lscpf,                                                                &
                   scpf_pfrac     = scpf_pfrac,                                                           &
                   scpf_resfact   = scpf_resfact,                                                         &
                   SCF_out        = scf_3d(i_startidx:i_endidx, 1:p_patch%nlev, jb),                      &
                   log_3momentIce = l3mom_ice,                                                            &
                   log_LiquidFrac = lliqfrac,                                                             &
                   diag_dhmax     = diag_dhmax_4d(i_startidx:i_endidx, 1:p_patch%nlev, jb, 1:n_icecat),   &
                   timer          = timer,                                                                &
                   timer_description = timer_description                                                  &
                   )
    END DO

!!$OMP END DO
!!$OMP END PARALLEL

    !IF (rank == 0) WRITE (0,*) 'after do loop'

    ! update P3 vars and tracers
    icon_tracer_3d%qi = sum(qitot_4d, dim=4)
    DO i_icecat = 1, n_icecat
      p3_tracer_3d(i_icecat)%qitot = qitot_4d(:,:,:,i_icecat)
      p3_tracer_3d(i_icecat)%qnitot = qnitot_4d(:,:,:,i_icecat)
      p3_tracer_3d(i_icecat)%qirim = qirim_4d(:,:,:,i_icecat)
      p3_tracer_3d(i_icecat)%birim = birim_4d(:,:,:,i_icecat)
      IF (l3mom_ice) THEN
        p3_tracer_3d(i_icecat)%qzitot = qzitot_4d(:,:,:,i_icecat)
      ENDIF
      IF (lliqfrac) THEN
        p3_tracer_3d(i_icecat)%qiliq = qiliq_4d(:,:,:,i_icecat)
      ENDIF

      p3_vars_3d(i_icecat)%dmean_i = diag_di_4d(:,:,:,i_icecat)
      p3_vars_3d(i_icecat)%deff_i = diag_effi_4d(:,:,:,i_icecat) * 2.
      p3_vars_3d(i_icecat)%rho_i = diag_rhoi_4d(:,:,:,i_icecat)
      p3_vars_3d(i_icecat)%vm_i = diag_vmi_4d(:,:,:,i_icecat)
    END DO

    ! update warm-cloud tracers
    icon_tracer_3d%qc = qc_3d
    icon_tracer_3d%qr = qr_3d
    icon_tracer_3d%qnc = qnc_3d
    icon_tracer_3d%qnr = qnr_3d

    ! update thermodynamic vars
    ! updates to theta_v and temp_v are done in ICON after this call to p3plugin
    dyn_vars_3d%temp = theta_3d * dyn_vars_3d%exner
    dyn_vars_3d%theta_old = theta_old_3d
    icon_tracer_3d%qv = qv_3d
    icon_tracer_3d%qv_old = qv_old_3d

    ! calc heating rate (in K s-1)
    dyn_vars_3d%ddt_temp_phys = (dyn_vars_3d%temp - temp_before_phys_3d) / dtime

    ! update microphysical vars
    mp_vars_3d%dmean_c = diag_3d(:,:,:,1)
    mp_vars_3d%dmean_r = diag_3d(:,:,:,2)
    mp_vars_3d%deff_c = diag_effc_3d * 2.

    IF (n_icecat == 1) THEN
      mp_vars_3d%deff_i = p3_vars_3d(1)%deff_i
    ELSE
      dummysum = 0.0
      DO i_icecat = 1, n_icecat
        dummysum = dummysum + p3_tracer_3d(i_icecat)%qitot / p3_vars_3d(i_icecat)%deff_i
      END DO
      mp_vars_3d%deff_i = icon_tracer_3d%qi / dummysum
    ENDIF

    mp_vars_3d%dhmax = diag_dhmax_4d(:,:,:,0)
    mp_vars_3d%dhmax_ground(:,:,1) = max(mp_vars_3d%dhmax_ground(:,:,1), diag_dhmax_4d(:,p_patch%nlev,:,1))
    mp_vars_3d%ze_p3 = diag_ze_3d

    mp_vars_3d%prec_gsp_rate(:,:,1) = (prt_liq_2d + prt_sol_2d) * 1000.
    mp_vars_3d%rain_gsp_rate(:,:,1) = prt_liq_2d * 1000.
    mp_vars_3d%snow_gsp_rate(:,:,1) = prt_sol_2d * 1000.
    mp_vars_3d%prec_gsp(:,:,1) = mp_vars_3d%prec_gsp(:,:,1) + (prt_liq_2d + prt_sol_2d) * 1000.
    mp_vars_3d%prec_gsp_d(:,:,1) = mp_vars_3d%prec_gsp_d(:,:,1) + (prt_liq_2d + prt_sol_2d) * 1000.
    mp_vars_3d%rain_gsp(:,:,1) = mp_vars_3d%rain_gsp(:,:,1) + prt_liq_2d * 1000.
    mp_vars_3d%snow_gsp(:,:,1) = mp_vars_3d%snow_gsp(:,:,1) + prt_sol_2d * 1000.
    !IF (rank == 0) WRITE (0,*) 'shape(prt_sol_2d)', shape(prt_sol_2d)

    ! IF (rank == 0) THEN
    !   CALL print_global_max('w', dyn_vars_3d%w_hl)
    !   CALL print_global_max('qv', icon_tracer_3d%qv)
    !   CALL print_global_max('qc', icon_tracer_3d%qc)
    !   CALL print_global_max('qnc', icon_tracer_3d%qnc)
    !   CALL print_global_max('qr', icon_tracer_3d%qr)
    !   CALL print_global_max('qnr', icon_tracer_3d%qnr)
    !   CALL print_global_max('qi', icon_tracer_3d%qi)
    !   DO i_icecat = 1, n_icecat
    !     CALL print_global_max('qitot', p3_tracer_3d(i_icecat)%qitot)
    !     CALL print_global_max('dmean_i', p3_vars_3d(i_icecat)%dmean_i)
    !   END DO
    !   !CALL print_global_max('dhmax', mp_vars_3d%dhmax)
    !   !CALL print_global_max('dhmax_ground', mp_vars_3d%dhmax)

    !   !CALL print_global_max('rain_gsp_rate', mp_vars_3d%rain_gsp_rate, 3600.)
    !   !CALL print_global_max('snow_gsp_rate', mp_vars_3d%snow_gsp_rate, 3600.)
    !   !CALL print_global_max('ddt_temp_phys', dyn_vars_3d%ddt_temp_phys)
    ! ENDIF

    ! increment fastphysics step number
    fastphystep = fastphystep + 1

    ! clean up
    DEALLOCATE(diag_2d, diag_3d)

    CALL dyn_vars_3d%nullify()
    CALL mp_vars_3d%nullify()
    CALL icon_tracer_3d%nullify()
    DO i_icecat = 1, n_icecat
      CALL p3_vars_3d(i_icecat)%nullify()
      CALL p3_tracer_3d(i_icecat)%nullify()
    END DO

    !IF (rank == 0) WRITE (0,*) 'end of p3_main_wrapper'

  END SUBROUTINE p3_main_wrapper


  ! --------------------------------------------------------------------
  ! Set reff before the rad call to values from P3
  ! --------------------------------------------------------------------
  SUBROUTINE set_reff_before_rad()  BIND(C)
    INTEGER               :: rank
    TYPE(t_mp_vars_3dptr) :: mp_vars_3d

    rank = comin_parallel_get_host_mpi_rank()
    IF (rank == 0) WRITE (0,*) 'set reff_qc, reff_qi to P3 values'

    CALL mp_vars%reff_qc%to_3d(mp_vars_3d%reff_qc)
    CALL mp_vars%deff_c%to_3d(mp_vars_3d%deff_c)
    CALL mp_vars%reff_qi%to_3d(mp_vars_3d%reff_qi)
    CALL mp_vars%deff_i%to_3d(mp_vars_3d%deff_i)

    mp_vars_3d%reff_qc = mp_vars_3d%deff_c / 2.
    mp_vars_3d%reff_qi = mp_vars_3d%deff_i / 2.

    CALL mp_vars_3d%nullify()

  END SUBROUTINE set_reff_before_rad

END MODULE p3plugin_main_routines
