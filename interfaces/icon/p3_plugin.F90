
MODULE p3_plugin
  !USE omp_lib
  USE mpi
  USE netcdf,                  ONLY : nf90_open, nf90_close, nf90_inq_dimid, nf90_inquire,              &
    &                                 nf90_inq_varid, nf90_inquire_variable, nf90_inquire_dimension,    &
    &                                 nf90_get_var, NF90_FLOAT, NF90_DOUBLE,                            &
    &                                 NF90_NOWRITE, NF90_NOERR, NF90_MAX_VAR_DIMS

  USE comin_plugin_interface,  ONLY : comin_callback_register, comin_var_request_add,                   &
    &                                 comin_var_get, t_comin_var_descriptor, t_comin_var_handle,        &
    &                                 comin_parallel_get_host_mpi_rank,                                 &
    &                                 comin_parallel_get_host_mpi_comm,                                 &
    &                                 comin_parallel_get_plugin_mpi_comm,                               &
    &                                 t_comin_setup_version_info, comin_setup_get_version,              &
    &                                 comin_descrdata_get_domain, t_comin_descrdata_domain,             &
    &                                 comin_descrdata_get_global, t_comin_descrdata_global,             &
    &                                 EP_SECONDARY_CONSTRUCTOR, EP_ATM_INIT_FINALIZE, EP_DESTRUCTOR,    &
    &                                 EP_ATM_TURBULENCE_AFTER, EP_ATM_MICROPHYSICS_BEFORE,              &
    &                                 EP_ATM_RADIATION_BEFORE,                                          &
    &                                 EP_ATM_NUDGING_AFTER, EP_ATM_WRITE_OUTPUT_BEFORE,                 &
    &                                 COMIN_FLAG_READ, COMIN_FLAG_WRITE, COMIN_ZAXIS_3D, COMIN_ZAXIS_2D,&
    &                                 COMIN_VAR_DATATYPE_DOUBLE, COMIN_VAR_DATATYPE_FLOAT,              &
    &                                 t_comin_plugin_info, comin_current_get_plugin_info,               &
    &                                 comin_plugin_finish, comin_metadata_set,                          &
    &                                 comin_metadata_get, comin_descrdata_get_timesteplength,           &
    &                                 comin_descrdata_get_cell_indices

  USE p3_plugin_types,         ONLY : t_dyn_vars_handle, t_dyn_vars_3dptr,                              &
    &                                 t_mp_vars_handle, t_mp_vars_3dptr,                                &
    &                                 t_p3_vars_handle, t_p3_vars_3dptr,                                &
    &                                 t_icon_tracer_handle, t_icon_tracer_3dptr,                        &
    &                                 t_p3_tracer_handle, t_p3_tracer_3dptr

  USE microphy_p3,             ONLY : p3_init, p3_main, status_ok

  IMPLICIT NONE
  PRIVATE

  CHARACTER(*), PARAMETER :: icon_namelist_name = 'NAMELIST_NWP'
  CHARACTER(*), PARAMETER :: lookup_file_dir = '/home/hk-project-aci/nw5893/ICON/p3-microphysics-5.4-icon/lookup_tables'

  INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND(12,307)
  TYPE(t_comin_setup_version_info)        :: version
  TYPE(t_comin_descrdata_global), POINTER :: p_global
  TYPE(t_comin_descrdata_domain), POINTER :: p_patch

  TYPE(t_dyn_vars_handle)               :: dyn_vars
  TYPE(t_mp_vars_handle)                :: mp_vars
  TYPE(t_p3_vars_handle), ALLOCATABLE   :: p3_vars(:)
  TYPE(t_icon_tracer_handle)            :: icon_tracer, icon_tracer_ddt_turb
  TYPE(t_p3_tracer_handle), ALLOCATABLE :: p3_tracer(:), p3_tracer_ddt_turb(:)

  INTEGER        :: rank, fastphystep, i_icecat, n_icecat, ihydrometeor_ini
  REAL           :: dtime
  CHARACTER(20)  :: icecat_name, unit_name
  CHARACTER(999) :: tracer_ini_filename
  LOGICAL        :: l3mom_ice, lliqfrac, ltracer_turb

  NAMELIST /p3_nml/ n_icecat, l3mom_ice, lliqfrac, ihydrometeor_ini, tracer_ini_filename

CONTAINS

  ! --------------------------------------------------------------------
  ! ComIn primary constructor
  ! --------------------------------------------------------------------
  SUBROUTINE comin_main()  BIND(C)

    TYPE(t_comin_plugin_info) :: this_plugin
    INTEGER                   :: funit

    rank = comin_parallel_get_host_mpi_rank()

    version = comin_setup_get_version()
    IF (version%version_no_major > 1)  THEN
      CALL comin_plugin_finish('comin_main (p3_plugin)', 'incompatible version!')
    END IF

    ! print plugin id
    CALL comin_current_get_plugin_info(this_plugin)
    IF (rank == 0) WRITE (0,'(a,a,a,i4)') ' comin plugin ', this_plugin%name, ' has id: ', this_plugin%id

    p_global => comin_descrdata_get_global()
    p_patch  => comin_descrdata_get_domain(1)
    fastphystep = 1
    dtime = comin_descrdata_get_timesteplength(1)


    ! read p3 namelist
    OPEN(newunit=funit, file=icon_namelist_name, action='read', form='formatted')
    READ(funit, nml=p3_nml)
    CLOSE(funit)

    IF (rank == 0) WRITE (0,'(a)') ' read P3 settings:'
    IF (rank == 0) WRITE (0,'(a,i1)') ' n_icecat  = ', n_icecat
    IF (rank == 0) WRITE (0,'(a,l)') ' l3mom_ice =', l3mom_ice
    IF (rank == 0) WRITE (0,'(a,l)') ' lliqfrac  =', lliqfrac

    ALLOCATE(p3_vars(n_icecat))
    ALLOCATE(p3_tracer(n_icecat))
    ALLOCATE(p3_tracer_ddt_turb(n_icecat))


    ! create new vars in ICON

    CALL create_var('theta_old', 'K', '3d', 'dp')
    CALL create_var('qv_old', 'kg kg-1', '3d', 'dp')
    CALL create_var('ddt_temp_phys', 'K s-1', '3d', 'dp')
    CALL create_var('dmean_c', 'm', '3d', 'dp')
    CALL create_var('dmean_r', 'm', '3d', 'dp')
    CALL create_var('deff_c', 'm', '3d', 'dp')
    CALL create_var('deff_i', 'm', '3d', 'dp')
    CALL create_var('reff_qc', 'm', '3d', 'dp')
    CALL create_var('reff_qi', 'm', '3d', 'dp')
    CALL create_var('dhmax', 'm', '3d', 'dp')
    CALL create_var('dhmax_ground', 'm', '2d', 'dp')
    CALL create_var('ze_p3', 'dBZ', '3d', 'dp')

    DO i_icecat = 1, n_icecat
      WRITE(icecat_name, '(a,i0)') 'dmean_i', i_icecat
      CALL create_var(icecat_name, 'm', '3d', 'dp')
      WRITE(icecat_name, '(a,i0)') 'deff_i', i_icecat
      CALL create_var(icecat_name, 'm', '3d', 'dp')
      WRITE(icecat_name, '(a,i0)') 'rho_i', i_icecat
      CALL create_var(icecat_name, 'kg m-3', '3d', 'dp')
      WRITE(icecat_name, '(a,i0)') 'vm_i', i_icecat
      CALL create_var(icecat_name, 'm s-1', '3d', 'dp')
    END DO


    ! create new tracers in ICON

    ltracer_turb = .FALSE.
    !ltracer_turb = .TRUE.

    CALL create_tracer('qnc', 'kg-1', ltracer_turb)
    CALL create_tracer('qnr', 'kg-1', ltracer_turb)

    DO i_icecat = 1, n_icecat
      WRITE(icecat_name, '(a,i0)') 'qitot_', i_icecat
      WRITE(unit_name, '(a)')      'kg kg-1'
      CALL create_tracer(icecat_name, unit_name, ltracer_turb)
    END DO
    DO i_icecat = 1, n_icecat
      WRITE(icecat_name, '(a,i0)') 'qnitot_', i_icecat
      WRITE(unit_name, '(a)')      'kg-1'
      CALL create_tracer(icecat_name, unit_name, ltracer_turb)
    END DO
    DO i_icecat = 1, n_icecat
      WRITE(icecat_name, '(a,i0)') 'qirim_', i_icecat
      WRITE(unit_name, '(a)')      'kg kg-1'
      CALL create_tracer(icecat_name, unit_name, ltracer_turb)
    END DO
    DO i_icecat = 1, n_icecat
      WRITE(icecat_name, '(a,i0)') 'birim_', i_icecat
      WRITE(unit_name, '(a)')      'm3 kg-1'
      CALL create_tracer(icecat_name, unit_name, ltracer_turb)
    END DO
    IF (l3mom_ice) THEN
      DO i_icecat = 1, n_icecat
        WRITE(icecat_name, '(a,i0)') 'qzitot_', i_icecat
        WRITE(unit_name, '(a)')      'm6 kg-1'
        CALL create_tracer(icecat_name, unit_name, ltracer_turb)
      END DO
    ENDIF
    IF (lliqfrac) THEN
      DO i_icecat = 1, n_icecat
        WRITE(icecat_name, '(a,i0)') 'qiliq_', i_icecat
        WRITE(unit_name, '(a)')      'kg kg-1'
        CALL create_tracer(icecat_name, unit_name, ltracer_turb)
      END DO
    ENDIF


    ! register callbacks
    CALL comin_callback_register(EP_SECONDARY_CONSTRUCTOR, secondary_constructor)
    CALL comin_callback_register(EP_ATM_INIT_FINALIZE, call_p3_init)
    !CALL comin_callback_register(EP_ATM_TURBULENCE_AFTER, update_turb_tend)
    CALL comin_callback_register(EP_ATM_MICROPHYSICS_BEFORE, run_custom_microphysics)
    !CALL comin_callback_register(EP_ATM_RADIATION_BEFORE, set_reff_before_rad)
    !CALL comin_callback_register(EP_ATM_NUDGING_AFTER, update_ice_after_nudging)

  END SUBROUTINE comin_main


  ! --------------------------------------------------------------------
  ! ComIn secondary constructor to connect all vars and tracers
  ! --------------------------------------------------------------------
  SUBROUTINE secondary_constructor()  BIND(C)
    INTEGER :: ep_init = EP_ATM_INIT_FINALIZE
    INTEGER :: ep_turb = EP_ATM_TURBULENCE_AFTER
    INTEGER :: ep_mp   = EP_ATM_MICROPHYSICS_BEFORE
    INTEGER :: ep_reff = EP_ATM_RADIATION_BEFORE
    INTEGER :: ep_nudg = EP_ATM_NUDGING_AFTER
    INTEGER :: ep_out  = EP_ATM_WRITE_OUTPUT_BEFORE
    INTEGER :: FR = COMIN_FLAG_READ
    INTEGER :: FW = COMIN_FLAG_WRITE
    INTEGER :: id = 1

    IF (rank == 0) WRITE (0,*) 'run secondary constructor'

    CALL comin_var_get([ep_init], t_comin_var_descriptor('z_mc', id), FR, dyn_vars%hfl)
    CALL comin_var_get([ep_mp], t_comin_var_descriptor('ddqz_z_full', id), FR, dyn_vars%dz)
    CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor('temp', id), IOR(FR, FW), dyn_vars%temp)
    CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor('rho', id), FR, dyn_vars%rho)
    CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor('pres', id), FR, dyn_vars%pres)
    CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor('exner', id), FR, dyn_vars%exner)
    CALL comin_var_get([ep_mp], t_comin_var_descriptor('theta_old', id), IOR(FR, FW), dyn_vars%theta_old)
    CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor('ddt_temp_phys', id), IOR(FR, FW), dyn_vars%ddt_temp_phys)
    CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor('w', id), FR, dyn_vars%w_hl)

    CALL comin_var_get([ep_init, ep_mp, ep_out], t_comin_var_descriptor('qv', id), IOR(FR, FW), icon_tracer%qv)
    CALL comin_var_get([ep_mp], t_comin_var_descriptor('qv_old', id), IOR(FR, FW), icon_tracer%qv_old)
    CALL comin_var_get([ep_init, ep_mp, ep_out], t_comin_var_descriptor('qc', id), IOR(FR, FW), icon_tracer%qc)
    CALL comin_var_get([ep_init, ep_mp, ep_nudg, ep_out], t_comin_var_descriptor('qi', id), IOR(FR, FW), icon_tracer%qi)
    CALL comin_var_get([ep_init, ep_mp, ep_nudg, ep_out], t_comin_var_descriptor('qs', id), IOR(FR, FW), icon_tracer%qs)
    CALL comin_var_get([ep_init, ep_mp, ep_out], t_comin_var_descriptor('qr', id), IOR(FR, FW), icon_tracer%qr)
    CALL comin_var_get([ep_init, ep_turb, ep_mp, ep_out], t_comin_var_descriptor('qnc', id), IOR(FR, FW), icon_tracer%qnc)
    CALL comin_var_get([ep_init, ep_turb, ep_mp, ep_out], t_comin_var_descriptor('qnr', id), IOR(FR, FW), icon_tracer%qnr)
    !CALL comin_var_get([ep_turb], t_comin_var_descriptor('ddt_qnc_turb', id), FR, icon_tracer_ddt_turb%qnc)
    !CALL comin_var_get([ep_turb], t_comin_var_descriptor('ddt_qnr_turb', id), FR, icon_tracer_ddt_turb%qnr)

    CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor('dmean_c', id), IOR(FR, FW), mp_vars%dmean_c)
    CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor('dmean_r', id), IOR(FR, FW), mp_vars%dmean_r)
    CALL comin_var_get([ep_mp, ep_reff, ep_out], t_comin_var_descriptor('deff_c', id), IOR(FR, FW), mp_vars%deff_c)
    CALL comin_var_get([ep_mp, ep_reff, ep_out], t_comin_var_descriptor('deff_i', id), IOR(FR, FW), mp_vars%deff_i)
    CALL comin_var_get([ep_reff, ep_out], t_comin_var_descriptor('reff_qc', id), IOR(FR, FW), mp_vars%reff_qc)
    CALL comin_var_get([ep_reff, ep_out], t_comin_var_descriptor('reff_qi', id), IOR(FR, FW), mp_vars%reff_qi)
    CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor('dhmax', id), IOR(FR, FW), mp_vars%dhmax)
    CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor('dhmax_ground', id), IOR(FR, FW), mp_vars%dhmax_ground)
    CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor('ze_p3', id), IOR(FR, FW), mp_vars%ze_p3)

    CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor('prec_gsp_rate', id), IOR(FR, FW), mp_vars%prec_gsp_rate)
    CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor('rain_gsp_rate', id), IOR(FR, FW), mp_vars%rain_gsp_rate)
    CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor('snow_gsp_rate', id), IOR(FR, FW), mp_vars%snow_gsp_rate)
    CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor('prec_gsp', id), IOR(FR, FW), mp_vars%prec_gsp)
    CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor('prec_gsp_d', id), IOR(FR, FW), mp_vars%prec_gsp_d)
    CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor('rain_gsp', id), IOR(FR, FW), mp_vars%rain_gsp)
    CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor('snow_gsp', id), IOR(FR, FW), mp_vars%snow_gsp)
    !CALL comin_var_get([ep_mp], t_comin_var_descriptor('q_sedim', id), IOR(FR, FW), mp_vars%q_sedim)
    !CALL comin_var_get([ep_mp], t_comin_var_descriptor('twater', id), IOR(FR, FW), mp_vars%twater)

    DO i_icecat = 1, n_icecat
      WRITE(icecat_name, '(a,i0)') 'dmean_i', i_icecat
      CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor(TRIM(icecat_name), id), IOR(FR, FW), p3_vars(i_icecat)%dmean_i)
      WRITE(icecat_name, '(a,i0)') 'deff_i', i_icecat
      CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor(TRIM(icecat_name), id), IOR(FR, FW), p3_vars(i_icecat)%deff_i)
      WRITE(icecat_name, '(a,i0)') 'rho_i', i_icecat
      CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor(TRIM(icecat_name), id), IOR(FR, FW), p3_vars(i_icecat)%rho_i)
      WRITE(icecat_name, '(a,i0)') 'vm_i', i_icecat
      CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor(TRIM(icecat_name), id), IOR(FR, FW), p3_vars(i_icecat)%vm_i)

      WRITE(icecat_name, '(a,i0)') 'qitot_', i_icecat
      CALL comin_var_get([ep_init, ep_turb, ep_mp, ep_nudg, ep_out], t_comin_var_descriptor(TRIM(icecat_name), id), IOR(FR, FW), p3_tracer(i_icecat)%qitot)
      !CALL comin_var_get([ep_turb, ep_mp], t_comin_var_descriptor('ddt_'//TRIM(icecat_name)//'_turb', id), FR, p3_tracer_ddt_turb(i_icecat)%qitot)
      WRITE(icecat_name, '(a,i0)') 'qnitot_', i_icecat
      CALL comin_var_get([ep_init, ep_turb, ep_mp, ep_out], t_comin_var_descriptor(TRIM(icecat_name), id), IOR(FR, FW), p3_tracer(i_icecat)%qnitot)
      !CALL comin_var_get([ep_turb, ep_mp], t_comin_var_descriptor('ddt_'//TRIM(icecat_name)//'_turb', id), FR, p3_tracer_ddt_turb(i_icecat)%qnitot)
      WRITE(icecat_name, '(a,i0)') 'qirim_', i_icecat
      CALL comin_var_get([ep_init, ep_turb, ep_mp, ep_out], t_comin_var_descriptor(TRIM(icecat_name), id), IOR(FR, FW), p3_tracer(i_icecat)%qirim)
      !CALL comin_var_get([ep_turb, ep_mp], t_comin_var_descriptor('ddt_'//TRIM(icecat_name)//'_turb', id), FR, p3_tracer_ddt_turb(i_icecat)%qirim)
      WRITE(icecat_name, '(a,i0)') 'birim_', i_icecat
      CALL comin_var_get([ep_init, ep_turb, ep_mp, ep_out], t_comin_var_descriptor(TRIM(icecat_name), id), IOR(FR, FW), p3_tracer(i_icecat)%birim)
      !CALL comin_var_get([ep_turb, ep_mp], t_comin_var_descriptor('ddt_'//TRIM(icecat_name)//'_turb', id), FR, p3_tracer_ddt_turb(i_icecat)%birim)
      IF (l3mom_ice) THEN
        WRITE(icecat_name, '(a,i0)') 'qzitot_', i_icecat
        CALL comin_var_get([ep_init, ep_turb, ep_mp, ep_out], t_comin_var_descriptor(TRIM(icecat_name), id), IOR(FR, FW), p3_tracer(i_icecat)%qzitot)
        !CALL comin_var_get([ep_turb, ep_mp], t_comin_var_descriptor('ddt_'//TRIM(icecat_name)//'_turb', id), FR, p3_tracer_ddt_turb(i_icecat)%qzitot)
      ENDIF
      IF (lliqfrac) THEN
        WRITE(icecat_name, '(a,i0)') 'qiliq_', i_icecat
        CALL comin_var_get([ep_init, ep_turb, ep_mp, ep_out], t_comin_var_descriptor(TRIM(icecat_name), id), IOR(FR, FW), p3_tracer(i_icecat)%qiliq)
        !CALL comin_var_get([ep_turb, ep_mp], t_comin_var_descriptor('ddt_'//TRIM(icecat_name)//'_turb', id), FR, p3_tracer_ddt_turb(i_icecat)%qiliq)
      ENDIF

    END DO

  END SUBROUTINE secondary_constructor


  ! ------------------------------------------------------------------------
  ! Call p3_init to load lookup tables etc. and then initialize P3 tracers
  ! ------------------------------------------------------------------------
  SUBROUTINE call_p3_init()  BIND(C)

    CHARACTER(16) :: model        = 'ICON'
    CHARACTER(30) :: varname      = ''
    LOGICAL       :: abort_on_err = .TRUE.
    LOGICAL       :: dowr         = .FALSE.
    INTEGER       :: stat
    INTEGER       :: jg, jb, jk, jc, jglobal
    INTEGER       :: i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end
    REAL          :: hfl_3dpatch(p_global%nproma, p_patch%nlev, p_patch%cells%nblks)
    REAL          :: qs_ini_3dpatch(p_global%nproma, p_patch%nlev, p_patch%cells%nblks)
    REAL          :: qg_ini_3dpatch(p_global%nproma, p_patch%nlev, p_patch%cells%nblks)
    REAL          :: dmean_qc, dmean_qr, dmean_qi, dmean_qs, dmean_qg
    REAL          :: magicfac_qi, magicfac_qs, magicfac_qg
    REAL          :: rhow, rhop_qi, rhop_qs, rhop_qg, rhor_qs, rhor_qg
    REAL          :: frim_qs, frim_qg

    TYPE(t_dyn_vars_3dptr)    :: dyn_vars_3d
    TYPE(t_icon_tracer_3dptr) :: icon_tracer_3d
    TYPE(t_p3_tracer_3dptr)   :: p3_tracer_3d(n_icecat)

    IF (rank == 0) WRITE (0,*) 'call p3_init'

    CALL p3_init(lookup_file_dir, n_icecat, l3mom_ice, lliqfrac, model, stat, abort_on_err, dowr)
    IF (stat /= status_ok) CALL comin_plugin_finish('call_p3_init (p3_plugin)', 'failed!')


    IF (ihydrometeor_ini == 1) THEN
      IF (rank == 0) WRITE (0,*) 'initialize from 1M-scheme mass ice tracers qi, qs, qg'

      dmean_qc    = 15.e-6   ! number-mean diameter in m
      dmean_qr    = 600e-6   ! number-mean diameter in m (P3's internal upper limit is 1/inv_Drmax=2mm)
      rhow        = 1000.    ! density of liquid water in kg/m3

      dmean_qi    = 150.e-6  ! target number-mean diameter in m
      magicfac_qi = 17.      ! empirically estimated correction factor needed because ice is not spherical, not dense
      rhop_qi     = 917.     ! bulk density of whole cloud ice particle in kg/m3

      dmean_qs    = 2000.e-6  ! target number-mean diameter in m
      magicfac_qs = 208.      ! empirically estimated correction factor needed because ice is not spherical, not dense
      rhop_qs     = 200.      ! bulk density of whole snow particle in kg/m3
      rhor_qs     = 500.      ! bulk density of rimed ice part in kg/m3
      frim_qs     = 0.1       ! bulk rime mass fraction

      dmean_qg    = 1200.e-6  ! target number-mean diameter in m
      magicfac_qg = 15.3      ! empirically estimated correction factor needed because ice is not spherical, not dense
      rhop_qg     = 500.      ! bulk density of whole graupel particle in kg/m3
      rhor_qg     = 500.      ! bulk density of rimed ice part in kg/m3
      frim_qg     = 0.9       ! bulk rime mass fraction

      CALL dyn_vars%hfl%to_3d(dyn_vars_3d%hfl)
      hfl_3dpatch = dyn_vars_3d%hfl
      CALL read_vinterp_ini_var(tracer_ini_filename, 'QS', hfl_3dpatch, qs_ini_3dpatch)
      CALL read_vinterp_ini_var(tracer_ini_filename, 'QG', hfl_3dpatch, qg_ini_3dpatch)
    ENDIF

    !IF (rank == 0) WRITE (0,'(a)') 'successfully read values of ' // TRIM(varname)
    !IF (rank == 0) WRITE (0,'(a,F10.5,F10.5)') 'min max:', minval(qc_ini_3d(:,120)), maxval(qc_ini_3d(:,120))


    !IF (rank == 0) WRITE (0,*) 'set initial qnc & qnr to zero'
    !IF (rank == 0) WRITE (0,*) 'set initial qitot_1 to read in qi values, other qitot_x to zero'
    !IF (rank == 0) WRITE (0,*) 'set initial qnitot_x, qirim_x, birim_x, qzitot_x, qiliq_x all to zero'

    CALL icon_tracer%qv%to_3d(icon_tracer_3d%qv)
    CALL icon_tracer%qc%to_3d(icon_tracer_3d%qc)
    CALL icon_tracer%qnc%to_3d(icon_tracer_3d%qnc)
    CALL icon_tracer%qr%to_3d(icon_tracer_3d%qr)
    CALL icon_tracer%qnr%to_3d(icon_tracer_3d%qnr)
    CALL icon_tracer%qi%to_3d(icon_tracer_3d%qi)
    CALL icon_tracer%qs%to_3d(icon_tracer_3d%qs)

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
    END DO

    jg = 1
    rl_start = p_global%grf_bdywidth_c + 1
    rl_end   = p_global%min_rlcell_int

    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

!!$OMP PARALLEL PRIVATE(jb,i_startidx,i_endidx,jk,jc,i_icecat)
!!$OMP DO SCHEDULE(STATIC, 1)

    DO jb = i_startblk, i_endblk
      CALL comin_descrdata_get_cell_indices(jg, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)
      DO jk = 1, p_patch%nlev
        DO jc = i_startidx, i_endidx

          SELECT CASE (ihydrometeor_ini)
          CASE (0)
            ! initialization of warm phase:
            ! qv, qc, qr were already loaded and set in initicon if in ini file, therefore reset qc, qr to zero

            icon_tracer_3d%qc(jc,jk,jb) = 0.0
            icon_tracer_3d%qr(jc,jk,jb) = 0.0
            icon_tracer_3d%qnc(jc,jk,jb) = 0.0
            icon_tracer_3d%qnr(jc,jk,jb) = 0.0

            ! initialization of cold phase:

            icon_tracer_3d%qi(jc,jk,jb) = 0.0
            icon_tracer_3d%qs(jc,jk,jb) = 0.0

            DO i_icecat = 1, n_icecat
              p3_tracer_3d(i_icecat)%qitot(jc,jk,jb) = 0.0
              p3_tracer_3d(i_icecat)%qnitot(jc,jk,jb) = 0.0
              p3_tracer_3d(i_icecat)%qirim(jc,jk,jb) = 0.0
              p3_tracer_3d(i_icecat)%birim(jc,jk,jb) = 0.0
              IF (l3mom_ice) THEN
                p3_tracer_3d(i_icecat)%qzitot(jc,jk,jb) = 0.0
              ENDIF
              IF (lliqfrac) THEN
                p3_tracer_3d(i_icecat)%qiliq(jc,jk,jb) = 0.0
              ENDIF
            END DO

          CASE (1)
            ! initialization of warm phase:
            ! qv, qc, qr were already loaded and set in initicon

            icon_tracer_3d%qnc(jc,jk,jb) = icon_tracer_3d%qc(jc,jk,jb) / (rhow*3.14) * dmean_qc**-3
            icon_tracer_3d%qnr(jc,jk,jb) = icon_tracer_3d%qr(jc,jk,jb) / (rhow*3.14) * dmean_qr**-3
            ! dmean_qc = (qc / (qnc*rhow*3.14))**(1./3.)

            ! initialization of cold phase:
            ! 1: initialize from 1M-scheme mass tracers (qi, qs, qg)
            ! if n_icecat == 1: init only cloud ice qi and ignore precipitating types qs, qg
            ! if n_icecat == 2: use the first icecat for qi and the second icecat for merged qs + qg
            ! if n_icecat >= 3: use one icecat for qi, qs, qg each
            ! all other icecats (if more available) are kept empty

            IF (n_icecat == 1) THEN
              p3_tracer_3d(1)%qitot(jc,jk,jb) = icon_tracer_3d%qi(jc,jk,jb)
              p3_tracer_3d(1)%qnitot(jc,jk,jb) = icon_tracer_3d%qi(jc,jk,jb) / (rhop_qi*3.14) &
                &                                               * magicfac_qi * dmean_qi**-3
              p3_tracer_3d(1)%qirim(jc,jk,jb) = 0.0
              p3_tracer_3d(1)%birim(jc,jk,jb) = 0.0
              IF (l3mom_ice) THEN
                p3_tracer_3d(1)%qzitot(jc,jk,jb) = 0.0
              ENDIF
              IF (lliqfrac) THEN
                p3_tracer_3d(1)%qiliq(jc,jk,jb) = 0.0
              ENDIF

            ELSE
              ! cloud ice, assume no riming part present
              p3_tracer_3d(1)%qitot(jc,jk,jb) = icon_tracer_3d%qi(jc,jk,jb)
              p3_tracer_3d(1)%qnitot(jc,jk,jb) = icon_tracer_3d%qi(jc,jk,jb) / (rhop_qi*3.14) &
                &                                               * magicfac_qi * dmean_qi**-3
              p3_tracer_3d(1)%qirim(jc,jk,jb) = 0.0
              p3_tracer_3d(1)%birim(jc,jk,jb) = 0.0
              IF (l3mom_ice) THEN
                p3_tracer_3d(1)%qzitot(jc,jk,jb) = 0.0
              ENDIF
              IF (lliqfrac) THEN
                p3_tracer_3d(1)%qiliq(jc,jk,jb) = 0.0
              ENDIF

              ! precipitating ice, use two separate icecats if available, merge into one if not
              IF (n_icecat == 2) THEN
                p3_tracer_3d(2)%qitot(jc,jk,jb) = qs_ini_3dpatch(jc,jk,jb) + qg_ini_3dpatch(jc,jk,jb)
                p3_tracer_3d(2)%qnitot(jc,jk,jb) = qs_ini_3dpatch(jc,jk,jb) / (rhop_qs*3.14)  &
                  &                                * magicfac_qs * dmean_qs**-3               &
                  &                              + qg_ini_3dpatch(jc,jk,jb) / (rhop_qg*3.14)  &
                  &                                * magicfac_qg * dmean_qg**-3
                p3_tracer_3d(2)%qirim(jc,jk,jb) = qs_ini_3dpatch(jc,jk,jb) * frim_qs  &
                  &                             + qg_ini_3dpatch(jc,jk,jb) * frim_qg
                p3_tracer_3d(2)%birim(jc,jk,jb) = qs_ini_3dpatch(jc,jk,jb) * frim_qs / rhor_qs  &
                  &                             + qg_ini_3dpatch(jc,jk,jb) * frim_qg / rhor_qg
                IF (l3mom_ice) THEN
                  p3_tracer_3d(2)%qzitot(jc,jk,jb) = 0.0
                ENDIF
                IF (lliqfrac) THEN
                  p3_tracer_3d(2)%qiliq(jc,jk,jb) = 0.0
                ENDIF

              ELSE
                p3_tracer_3d(2)%qitot(jc,jk,jb) = qs_ini_3dpatch(jc,jk,jb)
                p3_tracer_3d(2)%qnitot(jc,jk,jb) = qs_ini_3dpatch(jc,jk,jb) / (rhop_qs*3.14) &
                  &                                * magicfac_qs * dmean_qs**-3
                p3_tracer_3d(2)%qirim(jc,jk,jb) = qs_ini_3dpatch(jc,jk,jb) * frim_qs
                p3_tracer_3d(2)%birim(jc,jk,jb) = qs_ini_3dpatch(jc,jk,jb) * frim_qs / rhor_qs
                IF (l3mom_ice) THEN
                  p3_tracer_3d(2)%qzitot(jc,jk,jb) = 0.0
                ENDIF
                IF (lliqfrac) THEN
                  p3_tracer_3d(2)%qiliq(jc,jk,jb) = 0.0
                ENDIF

                p3_tracer_3d(3)%qitot(jc,jk,jb) = qg_ini_3dpatch(jc,jk,jb)
                p3_tracer_3d(3)%qnitot(jc,jk,jb) = qg_ini_3dpatch(jc,jk,jb) / (rhop_qg*3.14) &
                  &                                * magicfac_qg * dmean_qg**-3
                p3_tracer_3d(3)%qirim(jc,jk,jb) = qg_ini_3dpatch(jc,jk,jb) * frim_qg
                p3_tracer_3d(3)%birim(jc,jk,jb) = qg_ini_3dpatch(jc,jk,jb) * frim_qg / rhor_qg
                IF (l3mom_ice) THEN
                  p3_tracer_3d(3)%qzitot(jc,jk,jb) = 0.0
                ENDIF
                IF (lliqfrac) THEN
                  p3_tracer_3d(3)%qiliq(jc,jk,jb) = 0.0
                ENDIF
              ENDIF

              IF (n_icecat > 3) THEN
                DO i_icecat = 4, n_icecat
                  p3_tracer_3d(i_icecat)%qitot(jc,jk,jb) = 0.0
                  p3_tracer_3d(i_icecat)%qnitot(jc,jk,jb) = 0.0
                  p3_tracer_3d(i_icecat)%qirim(jc,jk,jb) = 0.0
                  p3_tracer_3d(i_icecat)%birim(jc,jk,jb) = 0.0
                  IF (l3mom_ice) THEN
                    p3_tracer_3d(i_icecat)%qzitot(jc,jk,jb) = 0.0
                  ENDIF
                  IF (lliqfrac) THEN
                    p3_tracer_3d(i_icecat)%qiliq(jc,jk,jb) = 0.0
                  ENDIF
                END DO
              ENDIF

            ENDIF

          END SELECT

          ! set qs to zero and re-sum all qitot into qi
          icon_tracer_3d%qs(jc,jk,jb) = 0.0
          icon_tracer_3d%qi(jc,jk,jb) = 0.0
          DO i_icecat = 1, n_icecat
            icon_tracer_3d%qi(jc,jk,jb) = icon_tracer_3d%qi(jc,jk,jb) + p3_tracer_3d(i_icecat)%qitot(jc,jk,jb)
          END DO

        END DO
      END DO
    END DO

!!$OMP END DO
!!$OMP END PARALLEL

    CALL dyn_vars_3d%nullify()
    CALL icon_tracer_3d%nullify()
    DO i_icecat = 1, n_icecat
      CALL p3_tracer_3d(i_icecat)%nullify()
    END DO

  END SUBROUTINE call_p3_init


  ! --------------------------------------------------------------------
  ! Update the turbulence tendencies to all tracers
  ! --------------------------------------------------------------------
  SUBROUTINE update_turb_tend()  BIND(C)

    INTEGER :: jg, jb, jk, jc
    INTEGER :: i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end

    TYPE(t_icon_tracer_3dptr) :: icon_tracer_3d, icon_tracer_ddt_turb_3d
    TYPE(t_p3_tracer_3dptr)   :: p3_tracer_3d(n_icecat), p3_tracer_ddt_turb_3d(n_icecat)

    IF (rank == 0) WRITE (0,*) 'update turbulence tendencies of tracers'

    CALL icon_tracer%qnc%to_3d(icon_tracer_3d%qnc)
    CALL icon_tracer%qnr%to_3d(icon_tracer_3d%qnr)

    CALL icon_tracer_ddt_turb%qnc%to_3d(icon_tracer_ddt_turb_3d%qnc)
    CALL icon_tracer_ddt_turb%qnr%to_3d(icon_tracer_ddt_turb_3d%qnr)

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

      CALL p3_tracer_ddt_turb(i_icecat)%qitot%to_3d(p3_tracer_ddt_turb_3d(i_icecat)%qitot)
      CALL p3_tracer_ddt_turb(i_icecat)%qnitot%to_3d(p3_tracer_ddt_turb_3d(i_icecat)%qnitot)
      CALL p3_tracer_ddt_turb(i_icecat)%qirim%to_3d(p3_tracer_ddt_turb_3d(i_icecat)%qirim)
      CALL p3_tracer_ddt_turb(i_icecat)%birim%to_3d(p3_tracer_ddt_turb_3d(i_icecat)%birim)
      CALL p3_tracer_ddt_turb(i_icecat)%qzitot%to_3d(p3_tracer_ddt_turb_3d(i_icecat)%qzitot)
      CALL p3_tracer_ddt_turb(i_icecat)%qiliq%to_3d(p3_tracer_ddt_turb_3d(i_icecat)%qiliq)
    END DO

    !IF (rank == 0) WRITE(0,'(a,e12.5)') ' max(p3_tracer_3d(1)%qitot)', maxval(p3_tracer_3d(1)%qitot)

    jg = 1
    rl_start = p_global%grf_bdywidth_c + 1
    rl_end   = p_global%min_rlcell_int

    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

!!$OMP PARALLEL PRIVATE(jb,i_startidx,i_endidx,jk,jc,i_icecat)
!!$OMP DO SCHEDULE(STATIC, 1)

    DO jb = i_startblk, i_endblk

      CALL comin_descrdata_get_cell_indices(jg, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)
      DO jk = 1, p_patch%nlev
        DO jc = i_startidx, i_endidx
          icon_tracer_3d%qnc(jc,jk,jb) = max(0.0, icon_tracer_3d%qnc(jc,jk,jb) + dtime * icon_tracer_ddt_turb_3d%qnc(jc,jk,jb))
          icon_tracer_3d%qnr(jc,jk,jb) = max(0.0, icon_tracer_3d%qnr(jc,jk,jb) + dtime * icon_tracer_ddt_turb_3d%qnr(jc,jk,jb))

          DO i_icecat = 1, n_icecat
            p3_tracer_3d(i_icecat)%qitot(jc,jk,jb) = max(0.0, p3_tracer_3d(i_icecat)%qitot(jc,jk,jb) + dtime * p3_tracer_ddt_turb_3d(i_icecat)%qitot(jc,jk,jb))
            p3_tracer_3d(i_icecat)%qnitot(jc,jk,jb) = max(0.0, p3_tracer_3d(i_icecat)%qnitot(jc,jk,jb) + dtime * p3_tracer_ddt_turb_3d(i_icecat)%qnitot(jc,jk,jb))
            p3_tracer_3d(i_icecat)%qirim(jc,jk,jb) = max(0.0, p3_tracer_3d(i_icecat)%qirim(jc,jk,jb) + dtime * p3_tracer_ddt_turb_3d(i_icecat)%qirim(jc,jk,jb))
            p3_tracer_3d(i_icecat)%birim(jc,jk,jb) = max(0.0, p3_tracer_3d(i_icecat)%birim(jc,jk,jb) + dtime * p3_tracer_ddt_turb_3d(i_icecat)%birim(jc,jk,jb))
            IF (l3mom_ice) THEN
              p3_tracer_3d(i_icecat)%qzitot(jc,jk,jb) = max(0.0, p3_tracer_3d(i_icecat)%qzitot(jc,jk,jb) + dtime * p3_tracer_ddt_turb_3d(i_icecat)%qzitot(jc,jk,jb))
            ENDIF
            IF (lliqfrac) THEN
              p3_tracer_3d(i_icecat)%qiliq(jc,jk,jb) = max(0.0, p3_tracer_3d(i_icecat)%qiliq(jc,jk,jb) + dtime * p3_tracer_ddt_turb_3d(i_icecat)%qiliq(jc,jk,jb))
            ENDIF
          END DO

        END DO
      END DO
    END DO

!!$OMP END DO
!!$OMP END PARALLEL

    CALL icon_tracer_3d%nullify()
    CALL icon_tracer_ddt_turb_3d%nullify()
    DO i_icecat = 1, n_icecat
      CALL p3_tracer_3d(i_icecat)%nullify()
      CALL p3_tracer_ddt_turb_3d(i_icecat)%nullify()
    END DO

  END SUBROUTINE update_turb_tend


  ! --------------------------------------------------------------------
  ! Set reff before the rad call to values from P3
  ! --------------------------------------------------------------------
  SUBROUTINE set_reff_before_rad()  BIND(C)

    TYPE(t_mp_vars_3dptr) :: mp_vars_3d

    IF (rank == 0) WRITE (0,*) 'set reff_qc, reff_qi to P3 values'

    CALL mp_vars%reff_qc%to_3d(mp_vars_3d%reff_qc)
    CALL mp_vars%deff_c%to_3d(mp_vars_3d%deff_c)
    CALL mp_vars%reff_qi%to_3d(mp_vars_3d%reff_qi)
    CALL mp_vars%deff_i%to_3d(mp_vars_3d%deff_i)

    mp_vars_3d%reff_qc = mp_vars_3d%deff_c / 2.
    mp_vars_3d%reff_qi = mp_vars_3d%deff_i / 2.

    CALL mp_vars_3d%nullify()

  END SUBROUTINE set_reff_before_rad


  ! -----------------------------------------------------------------------------
  ! Update P3 ice to the changes that come from the changes in the nudging zone
  ! -----------------------------------------------------------------------------
  SUBROUTINE update_ice_after_nudging()  BIND(C)

    INTEGER  :: jg, jb, jk, jc
    INTEGER  :: i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end
    REAL(wp) :: qi_sum
    REAL(wp) :: qsmall = 1.0e-12    ! minimum threshold in kg/kg

    TYPE(t_icon_tracer_3dptr) :: icon_tracer_3d
    TYPE(t_p3_tracer_3dptr)   :: p3_tracer_3d(n_icecat)

    IF (rank == 0) WRITE (0,*) 'update qitot_x values by changes of nudging'

    CALL icon_tracer%qi%to_3d(icon_tracer_3d%qi)

    DO i_icecat = 1, n_icecat
      CALL p3_tracer(i_icecat)%qitot%to_3d(p3_tracer_3d(i_icecat)%qitot)
    END DO

    jg = 1
    rl_start = p_global%grf_bdywidth_c + 1
    rl_end   = p_global%min_rlcell_int

    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

!!$OMP PARALLEL PRIVATE(jb,i_startidx,i_endidx,jk,jc,i_icecat)
!!$OMP DO SCHEDULE(STATIC, 1)

    DO jb = i_startblk, i_endblk
      CALL comin_descrdata_get_cell_indices(jg, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)
      DO jk = 1, p_patch%nlev
        DO jc = i_startidx, i_endidx

          qi_sum = 0.0
          DO i_icecat = 1, n_icecat
            qi_sum = qi_sum + p3_tracer_3d(i_icecat)%qitot(jc,jk,jb)
          END DO

          IF (qi_sum > qsmall) THEN
            DO i_icecat = 1, n_icecat
              p3_tracer_3d(i_icecat)%qitot(jc,jk,jb) = icon_tracer_3d%qi(jc,jk,jb) * p3_tracer_3d(i_icecat)%qitot(jc,jk,jb) / qi_sum
            END DO
          ELSE
            p3_tracer_3d(1)%qitot(jc,jk,jb) = icon_tracer_3d%qi(jc,jk,jb)
          ENDIF

        END DO
      END DO
    END DO

!!$OMP END DO
!!$OMP END PARALLEL

    CALL icon_tracer_3d%nullify()
    DO i_icecat = 1, n_icecat
      CALL p3_tracer_3d(i_icecat)%nullify()
    END DO

  END SUBROUTINE update_ice_after_nudging


  ! --------------------------------------------------------------------
  ! ComIn constructor to run prepare and run P3 wrapper
  ! --------------------------------------------------------------------
  SUBROUTINE run_custom_microphysics()  BIND(C)

    INTEGER :: jg, jk, jb
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
    n_diag_3d = 3  ! diag_3d contains dmean_c, dmean_r
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
    ! updates to theta_v and temp_v are done in ICON after this call to p3_plugin
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

    CALL print_global_max('w', dyn_vars_3d%w_hl)
    CALL print_global_max('qv', icon_tracer_3d%qv)
    CALL print_global_max('qc', icon_tracer_3d%qc)
    CALL print_global_max('qnc', icon_tracer_3d%qnc)
    CALL print_global_max('qr', icon_tracer_3d%qr)
    CALL print_global_max('qnr', icon_tracer_3d%qnr)
    CALL print_global_max('qi', icon_tracer_3d%qi)
    DO i_icecat = 1, n_icecat
      CALL print_global_max('qitot', p3_tracer_3d(i_icecat)%qitot)
      CALL print_global_max('dmean_i', p3_vars_3d(i_icecat)%dmean_i)
    END DO
    !CALL print_global_max('dhmax', mp_vars_3d%dhmax)
    !CALL print_global_max('dhmax_ground', mp_vars_3d%dhmax)

    !CALL print_global_max('rain_gsp_rate', mp_vars_3d%rain_gsp_rate, 3600.)
    !CALL print_global_max('snow_gsp_rate', mp_vars_3d%snow_gsp_rate, 3600.)
    !CALL print_global_max('ddt_temp_phys', dyn_vars_3d%ddt_temp_phys)

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

    !IF (rank == 0) WRITE (0,*) 'end of run_custom_microphysics'

  END SUBROUTINE run_custom_microphysics

  ! --------------------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------------------

  SUBROUTINE print_global_max(var_name, var_3dptr, factor)
    CHARACTER(*), INTENT(IN)      :: var_name
    !TYPE(t_comin_var_handle), INTENT(IN) :: var
    REAL(wp), POINTER, INTENT(IN) :: var_3dptr(:,:,:)
    REAL, INTENT(IN), OPTIONAL    :: factor

    REAL                          :: local_max, global_max
    INTEGER                       :: comm, root, ierr

    !CALL var%get_ptr(var_ptr)

    IF (PRESENT(factor)) THEN
      local_max = maxval(var_3dptr * factor)
    ELSE
      local_max = maxval(var_3dptr)
    END IF

    root = 0
    comm = comin_parallel_get_plugin_mpi_comm()

    CALL MPI_REDUCE(local_max, global_max, 1, MPI_REAL, MPI_MAX, root, comm, ierr)
    IF (ierr /= 0) CALL comin_plugin_finish('print_global_max (p3_plugin)', 'failed!')

    !IF (rank == root) WRITE(0,*) 'local_max('//var_name//')', local_max
    IF (rank == root) WRITE(0,*) 'global_max('//var_name//')', global_max
  END SUBROUTINE print_global_max


  SUBROUTINE create_var(var_name, unit_name, axis_type, datatype_precision)
    CHARACTER(*), INTENT(IN) :: var_name, unit_name, axis_type, datatype_precision
    INTEGER                  :: zaxis_id, datatype

    IF (axis_type == '3d') THEN
      zaxis_id = COMIN_ZAXIS_3D
    ELSE IF (axis_type == '2d') THEN
      zaxis_id = COMIN_ZAXIS_2D
    END IF

    IF (datatype_precision == 'dp') THEN
      datatype = COMIN_VAR_DATATYPE_DOUBLE
    ELSE IF (axis_type == 'sp') THEN
      datatype = COMIN_VAR_DATATYPE_FLOAT
    END IF

    CALL comin_var_request_add_wrapper(descriptor=t_comin_var_descriptor(name=TRIM(var_name), id=1), &
      &                                units=TRIM(unit_name), lmode_exclusive=.TRUE., zaxis_id=zaxis_id, &
      &                                ltracer=.FALSE., lrestart=.FALSE., datatype=datatype)
  END SUBROUTINE create_var


  SUBROUTINE create_tracer(var_name, unit_name, ltracer_turb)
    CHARACTER(*), INTENT(IN) :: var_name, unit_name
    LOGICAL, INTENT(IN)      :: ltracer_turb

    CALL comin_var_request_add_wrapper(descriptor=t_comin_var_descriptor(name=TRIM(var_name), id=-1), &
      &                                units=TRIM(unit_name), lmode_exclusive=.FALSE., zaxis_id=COMIN_ZAXIS_3D, &
                                       ltracer=.TRUE., lrestart=.FALSE., ltracer_turb=ltracer_turb)
  END SUBROUTINE create_tracer


  SUBROUTINE comin_var_request_add_wrapper(descriptor, units, lmode_exclusive, zaxis_id, &
                                           ltracer, lrestart, datatype, ltracer_turb)
    TYPE(t_comin_var_descriptor), INTENT(IN)  :: descriptor
    LOGICAL,            OPTIONAL, INTENT(IN)  :: lmode_exclusive
    LOGICAL,            OPTIONAL, INTENT(IN)  :: ltracer_turb
    INTEGER,            OPTIONAL, INTENT(IN)  :: zaxis_id
    INTEGER,            OPTIONAL, INTENT(IN)  :: datatype
    LOGICAL,            OPTIONAL, INTENT(IN)  :: ltracer
    LOGICAL,            OPTIONAL, INTENT(IN)  :: lrestart
    CHARACTER(LEN=*),   OPTIONAL, INTENT(IN)  :: units
    LOGICAL                                   :: lexclusive

    IF (PRESENT(lmode_exclusive)) THEN
      lexclusive = lmode_exclusive
    ELSE
      lexclusive = .FALSE.
    END IF

    CALL comin_var_request_add(descriptor, lexclusive)

    IF (PRESENT(ltracer_turb)) THEN
      CALL comin_metadata_set(descriptor, "tracer_turb", ltracer_turb)
    END IF
    IF (PRESENT(zaxis_id)) THEN
      CALL comin_metadata_set(descriptor, "zaxis_id", zaxis_id)
    END IF
    IF (PRESENT(ltracer)) THEN
      CALL comin_metadata_set(descriptor, "tracer", ltracer)

      IF (ltracer) THEN
        CALL comin_metadata_set(descriptor, "tracer_hlimit", 4)
        CALL comin_metadata_set(descriptor, "tracer_hadv", 2)
        CALL comin_metadata_set(descriptor, "tracer_vadv", 3)
      END IF
    END IF
    IF (PRESENT(lrestart)) THEN
      CALL comin_metadata_set(descriptor, "restart", lrestart)
    END IF
    IF (PRESENT(units)) THEN
      CALL comin_metadata_set(descriptor, "units", TRIM(units))
    END IF
  END SUBROUTINE comin_var_request_add_wrapper


  SUBROUTINE read_vinterp_ini_var(filename, varname, hfl_3dpatch_outlevs, var_3dpatch_outlevs)
    CHARACTER(*), INTENT(IN) :: filename, varname
    REAL, INTENT(IN)         :: hfl_3dpatch_outlevs(:, :, :)
    REAL, INTENT(OUT)        :: var_3dpatch_outlevs(:, :, :)

    CHARACTER(30)            :: varname_dummy, dimname
    INTEGER                  :: nc_status, ncid, ncells_global, nlev_in !, nDimensions, nVariables, nAttributes
    INTEGER                  :: i, varid, hhlid, xtype, ndims, natts, dimpos_time, dimpos_ncells, dimid_ncells
    INTEGER                  :: dimids(NF90_MAX_VAR_DIMS), dimlen(NF90_MAX_VAR_DIMS)
    INTEGER                  :: jg, jb, jk, jc, jglobal
    INTEGER                  :: i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end
    REAL, ALLOCATABLE        :: varvalues_file_2d(:, :), hhlvalues_file_2d(:, :), varvalues_file_3d(:, :, :)
    REAL, ALLOCATABLE        :: var_global_inlevs(:, :), hhl_global_inlevs(:, :), hfl_1d_inlevs(:)

    rank = comin_parallel_get_host_mpi_rank()
    !IF (rank == 0) WRITE (0,'(a)') 'read from filename: ', TRIM(filename)

    ncid = -99
    nc_status = nf90_open(TRIM(filename), NF90_NOWRITE, ncid)
    IF (nc_status /= NF90_NOERR) &
      & CALL comin_plugin_finish('read_vinterp_ini_var (p3_plugin)', 'Could not read ini file: ' // TRIM(filename))

    ! read ncells of ini file and check if equals to icon model's ncells_global
    nc_status = nf90_inq_dimid(ncid, 'ncells', dimid_ncells)
    IF (nc_status /= NF90_NOERR) &
      & CALL comin_plugin_finish('read_vinterp_ini_var (p3_plugin)', &
                               & 'Could not find dimension ncells in ini file!')
    nc_status = nf90_inquire_dimension(ncid, dimid_ncells, dimname, ncells_global)
    IF (ncells_global /= p_patch%cells%ncells_global) &
      & CALL comin_plugin_finish('read_vinterp_ini_var (p3_plugin)', &
                               & 'ncells number of ini file does not match with model!')

    ! read in variable varname and handle dimensions
    !nc_status = nf90_inquire(ncid, nDimensions, nVariables, nAttributes)
    nc_status = nf90_inq_varid(ncid, varname, varid)
    IF (nc_status /= NF90_NOERR) &
      & CALL comin_plugin_finish('read_vinterp_ini_var (p3_plugin)', &
                               & 'Could not find variable "' // TRIM(varname) // '" in ini file: ' // TRIM(filename))

    nc_status = nf90_inquire_variable(ncid, varid, varname_dummy, xtype, ndims, dimids, natts)
    IF (nc_status /= NF90_NOERR) &
      & CALL comin_plugin_finish('read_vinterp_ini_var (p3_plugin)', &
                               & 'Could not inquire variable: ' // TRIM(varname))
    IF (xtype /= NF90_DOUBLE .and. xtype /= NF90_FLOAT) &
      & CALL comin_plugin_finish('read_vinterp_ini_var (p3_plugin)', &
                               & 'Variable type not double or float: ' // TRIM(varname))

    SELECT CASE (ndims)
    CASE (1)
      CALL comin_plugin_finish('read_vinterp_ini_var (p3_plugin)', &
                             & 'Variable has only one dimension: ' // TRIM(varname))
    CASE (2)
      dimpos_ncells = 0
      DO i = 1, 2
        nc_status = nf90_inquire_dimension(ncid, dimids(i), dimname, dimlen(i))
        !IF (rank == 0) WRITE (0,'(a,i2,a,i9)') "dim len of dimid", dimids(i), ", " // TRIM(dimname) // ":", dimlen(i)
        IF (TRIM(dimname) == 'time') THEN
          CALL comin_plugin_finish('read_vinterp_ini_var (p3_plugin)', &
                                 & 'Found dimension "time" in 2d array: ' // TRIM(varname))
        ELSE IF (TRIM(dimname) == 'ncells') THEN
          dimpos_ncells = i
        ELSE
          nlev_in = dimlen(i)
        ENDIF
      END DO

      IF (dimpos_ncells == 0) &
        & CALL comin_plugin_finish('read_vinterp_ini_var (p3_plugin)', &
                                 & 'Could not find dimension "ncells" in:' // TRIM(varname))

      ALLOCATE(var_global_inlevs(ncells_global, nlev_in))
      ALLOCATE(varvalues_file_2d(dimlen(1), dimlen(2)))
      nc_status = nf90_get_var(ncid, varid, varvalues_file_2d)
      !IF (rank == 0) WRITE (0,'(a)') 'successfully read values of ' // TRIM(varname)
      !IF (rank == 0) WRITE (0,'(a,F10.5,F10.5)') 'min max:', minval(varvalues_file_2d), maxval(varvalues_file_2d)

      SELECT CASE (dimpos_ncells)
      CASE (1)
        var_global_inlevs(:, :) = varvalues_file_2d(:, :)
      CASE (2)
        var_global_inlevs(:, :) = transpose(varvalues_file_2d(:, :))
      END SELECT

    CASE (3)
      dimpos_time = 0
      dimpos_ncells = 0
      DO i = 1, 3
        nc_status = nf90_inquire_dimension(ncid, dimids(i), dimname, dimlen(i))
        !IF (rank == 0) WRITE (0,'(a,i2,a,i9)') "dim len of dimid", dimids(i), ", " // TRIM(dimname) // ":", dimlen(i)
        IF (TRIM(dimname) == 'time') THEN
          dimpos_time = i
          IF (dimlen(i) > 1) THEN
            IF (rank == 0) WRITE (0,'(a)') 'Dimension "time" has more than one time step, choosing the first'
          ENDIF
        ELSE IF (TRIM(dimname) == 'ncells') THEN
          dimpos_ncells = i
        ELSE
          nlev_in = dimlen(i)
        ENDIF
      END DO

      IF (dimpos_time == 0) &
        & CALL comin_plugin_finish('read_vinterp_ini_var (p3_plugin)', &
                                 & 'Could not find dimension "time" in:' // TRIM(varname))
      IF (dimpos_ncells == 0) &
        & CALL comin_plugin_finish('read_vinterp_ini_var (p3_plugin)', &
                                 & 'Could not find dimension "ncells" in:' // TRIM(varname))

      ALLOCATE(var_global_inlevs(ncells_global, nlev_in))
      ALLOCATE(varvalues_file_3d(dimlen(1), dimlen(2), dimlen(3)))
      nc_status = nf90_get_var(ncid, varid, varvalues_file_3d)
      !IF (rank == 0) WRITE (0,'(a)') 'successfully read values of ' // TRIM(varname)
      !IF (rank == 0) WRITE (0,'(a,F10.5,F10.5)') 'min max:', minval(varvalues_file_3d), maxval(varvalues_file_3d)

      SELECT CASE (dimpos_time)
      CASE (1)
        IF (dimpos_ncells == 2) THEN
          var_global_inlevs(:, :) = varvalues_file_3d(1, :, :)
        ELSE
          var_global_inlevs(:, :) = transpose(varvalues_file_3d(1, :, :))
        ENDIF
      CASE (2)
        IF (dimpos_ncells == 1) THEN
          var_global_inlevs(:, :) = varvalues_file_3d(:, 1, :)
        ELSE
          var_global_inlevs(:, :) = transpose(varvalues_file_3d(:, 1, :))
        ENDIF
      CASE (3)
        IF (dimpos_ncells == 1) THEN
          var_global_inlevs(:, :) = varvalues_file_3d(:, :, 1)
        ELSE
          var_global_inlevs(:, :) = transpose(varvalues_file_3d(:, :, 1))
        ENDIF
      END SELECT

    CASE (4:)
      CALL comin_plugin_finish('read_vinterp_ini_var (p3_plugin)', &
                             & 'Variable has more than 3 dimensions: ' // TRIM(varname))
    END SELECT
    ! end of reading variable


    ! read in hhl of ini file
    nc_status = nf90_inq_varid(ncid, 'HHL', hhlid)
    IF (nc_status /= NF90_NOERR) &
      & nc_status = nf90_inq_varid(ncid, 'hhl', hhlid)
    IF (nc_status /= NF90_NOERR) &
      & CALL comin_plugin_finish('read_vinterp_ini_var (p3_plugin)', &
                               & 'Could not find hhl in ini file: ' // TRIM(filename))

    nc_status = nf90_inquire_variable(ncid, hhlid, varname_dummy, xtype, ndims, dimids, natts)
    IF (nc_status /= NF90_NOERR) &
      & CALL comin_plugin_finish('read_vinterp_ini_var (p3_plugin)', 'Could not inquire hhl!')
    IF (xtype /= NF90_DOUBLE .and. xtype /= NF90_FLOAT) &
      & CALL comin_plugin_finish('read_vinterp_ini_var (p3_plugin)', 'hhl type not double or float!')
    IF (ndims /= 2) &
      & CALL comin_plugin_finish('read_vinterp_ini_var (p3_plugin)', 'hhl does not have 2 dimensions!')
    dimpos_ncells = 0
    DO i = 1, 2
      nc_status = nf90_inquire_dimension(ncid, dimids(i), dimname, dimlen(i))
      !IF (rank == 0) WRITE (0,'(a,i2,a,i9)') "dim len of dimid", dimids(i), ", " // TRIM(dimname) // ":", dimlen(i)
      IF (TRIM(dimname) == 'ncells') dimpos_ncells = i
    END DO
    IF (dimpos_ncells == 0) &
      & CALL comin_plugin_finish('read_vinterp_ini_var (p3_plugin)', 'Could not find dimension "ncells" in hhl!')

    ALLOCATE(hhl_global_inlevs(ncells_global, nlev_in+1))
    ALLOCATE(hhlvalues_file_2d(dimlen(1), dimlen(2)))
    nc_status = nf90_get_var(ncid, hhlid, hhlvalues_file_2d)
    !IF (rank == 0) WRITE (0,'(a)') 'successfully read values of hhl'

    SELECT CASE (dimpos_ncells)
    CASE (1)
      hhl_global_inlevs(:, :) = hhlvalues_file_2d(:, :)
    CASE (2)
      hhl_global_inlevs(:, :) = transpose(hhlvalues_file_2d(:, :))
    END SELECT

    nc_status = nf90_close(ncid)
    IF (nc_status /= NF90_NOERR) &
      & CALL comin_plugin_finish('read_vinterp_ini_var (p3_plugin)', 'File closing not successful!')

    ! vertical interpolate from inlevs to model outlevs
    ALLOCATE(hfl_1d_inlevs(nlev_in))
    jg = 1
    rl_start = p_global%grf_bdywidth_c + 1
    rl_end   = p_global%min_rlcell_int

    i_startblk = p_patch%cells%start_block(rl_start)
    i_endblk   = p_patch%cells%end_block(rl_end)

    DO jb = i_startblk, i_endblk
      CALL comin_descrdata_get_cell_indices(jg, jb, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end)
      DO jc = i_startidx, i_endidx
        jglobal = p_patch%cells%glb_index((jb-1)*p_global%nproma+jc)
        hfl_1d_inlevs(:) = (hhl_global_inlevs(jglobal, 1:nlev_in) + hhl_global_inlevs(jglobal, 2:nlev_in+1)) / 2.

        CALL vert_intp_linear_1d(hfl_1d_inlevs(:),                var_global_inlevs(jglobal, :), &
                               & hfl_3dpatch_outlevs(jc, :, jb),  var_3dpatch_outlevs(jc, :, jb))
      END DO
    END DO

  END SUBROUTINE read_vinterp_ini_var


  ! linear vertical interpolation: height levels za -> zb
  ! extrapolates if no data given for zb > za
  ! taken from mo_nh_vert_interp_les:vert_intp_linear_1d (there taken from UCLA-LES)
  SUBROUTINE vert_intp_linear_1d(za, xa, zb, xb)
     REAL, INTENT(IN)  :: za(:), zb(:), xa(:)
     REAL, INTENT(OUT) :: xb(:)
     REAL              :: wt
     INTEGER           :: l, k

     l = SIZE(za)
     DO k = SIZE(zb), 1, -1
       IF (zb(k) <= za(1)) THEN
          DO WHILE ( zb(k) > za(l-1) .AND. l > 1)
             l = l-1
          END DO
          wt=(zb(k)-za(l))/(za(l-1)-za(l))
          xb(k)=xa(l)+(xa(l-1)-xa(l))*wt
       ELSE
          wt=(zb(k)-za(1))/(za(2)-za(1))
          xb(k)=xa(1)+(xa(2)-xa(1))*wt
       END IF
    END DO
  END SUBROUTINE vert_intp_linear_1d

END MODULE p3_plugin
