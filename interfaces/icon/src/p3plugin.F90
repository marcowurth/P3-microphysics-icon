
MODULE p3plugin
  USE mpi,                     ONLY : MPI_Comm_split_type, MPI_Comm_rank, MPI_Comm_size, MPI_INFO_NULL,   &
    &                                 MPI_COMM_TYPE_SHARED, MPI_Allreduce, MPI_Allgather, MPI_INT, MPI_MAX

  USE comin_plugin_interface,  ONLY : comin_callback_register, comin_var_get, t_comin_var_descriptor,     &
    &                                 comin_parallel_get_host_mpi_comm, comin_parallel_get_host_mpi_rank, &
    &                                 t_comin_setup_version_info, comin_setup_get_version,                &
    &                                 comin_descrdata_get_domain, comin_descrdata_get_global,             &
    &                                 EP_SECONDARY_CONSTRUCTOR, EP_ATM_INIT_FINALIZE, EP_DESTRUCTOR,      &
    &                                 EP_ATM_MICROPHYSICS_BEFORE, EP_ATM_TURBULENCE_AFTER,                &
    &                                 EP_ATM_NUDGING_AFTER, EP_ATM_RADIATION_BEFORE,                      &
    &                                 EP_ATM_WRITE_OUTPUT_BEFORE, COMIN_FLAG_READ, COMIN_FLAG_WRITE,      &
    &                                 t_comin_plugin_info, comin_current_get_plugin_info,                 &
    &                                 comin_plugin_finish

  USE p3plugin_global_vars,    ONLY : comm_world, comm_insidenode, rank_world, rank_insidenode,           &
    &                                 numprocs_insidenode, max_patch_size,                                &
    &                                 node_patches_sizes, node_patches_idx,                               &
    &                                 n_icecat, itracer_ini, l3mom_ice, lliqfrac, autoAccr_param_in,      &
    &                                 tracer_ini_filename, lookup_tables_path, p_global, p_patch,         &
    &                                 dyn_vars, mp_vars, p3_vars,                                         &
    &                                 icon_tracer, icon_tracer_ddt_turb, p3_tracer, p3_tracer_ddt_turb
  USE p3plugin_utils,          ONLY : create_var, create_tracer
  USE p3plugin_tracer_init,    ONLY : init_p3_and_tracer
  USE p3plugin_main_routines,  ONLY : p3_main_wrapper, set_reff_before_rad
!! JM_20260331 >> using new comin handle for ice-phase diagnostics
  USE p3plugin_global_vars,    ONLY : p3_ice_diag
!! << JM_20260331

  IMPLICIT NONE
  PRIVATE

  INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND(12,307)

  NAMELIST /p3_nml/ n_icecat, l3mom_ice, lliqfrac, itracer_ini, tracer_ini_filename, lookup_tables_path, autoAccr_param_in

CONTAINS

  ! --------------------------------------------------------------------
  ! ComIn primary constructor
  ! --------------------------------------------------------------------
  SUBROUTINE comin_main()  BIND(C)
    TYPE(t_comin_setup_version_info) :: version
    TYPE(t_comin_plugin_info)        :: this_plugin
    INTEGER                          :: k, ierr, loc_patch_size, funit, i_icecat
    CHARACTER(20)                    :: icecat_name, unit_name
    CHARACTER(99)                    :: icon_namelist_name
    LOGICAL                          :: ltracer_turb
    INTEGER, ALLOCATABLE             :: idxmap(:), idxvec_local(:), idxvec_global_padded(:)

    comm_world = comin_parallel_get_host_mpi_comm()
    rank_world = comin_parallel_get_host_mpi_rank()

    CALL MPI_Comm_split_type(comm_world, MPI_COMM_TYPE_SHARED, rank_world, MPI_INFO_NULL, comm_insidenode, ierr)
    CALL MPI_Comm_rank(comm_insidenode, rank_insidenode, ierr)
    CALL MPI_Comm_size(comm_insidenode, numprocs_insidenode, ierr)

    version = comin_setup_get_version()
    IF (version%version_no_major > 1)  THEN
      CALL comin_plugin_finish('comin_main (p3plugin)', 'incompatible version!')
    END IF

    ! print plugin id
    CALL comin_current_get_plugin_info(this_plugin)
    IF (rank_world == 0) WRITE (0,'(a,a,a,i4)') ' comin plugin ', this_plugin%name, ' has id: ', this_plugin%id

    p_global => comin_descrdata_get_global()
    p_patch  => comin_descrdata_get_domain(1)

    idxvec_local = [(k, k=1, p_patch%cells%ncells)]
    loc_patch_size = size(idxvec_local)

    CALL MPI_Allreduce(loc_patch_size, max_patch_size, 1, MPI_INT, MPI_MAX, comm_insidenode, ierr)

    ALLOCATE(idxvec_global_padded(max_patch_size))
    ALLOCATE(node_patches_sizes(numprocs_insidenode))
    ALLOCATE(node_patches_idx(max_patch_size, numprocs_insidenode))

    idxvec_global_padded(1:loc_patch_size) = p_patch%cells%glb_index(idxvec_local)
    idxvec_global_padded(loc_patch_size+1:max_patch_size) = -1    ! pad index array with negative numbers

    CALL MPI_Allgather(loc_patch_size,     1, MPI_INT, &
      &                node_patches_sizes, 1, MPI_INT, comm_insidenode, ierr)
    IF (ierr /= 0) CALL comin_plugin_finish('gathering node patches sizes (p3plugin)', 'failed!')

    CALL MPI_Allgather(idxvec_global_padded, max_patch_size, MPI_INT, &
      &                node_patches_idx,     max_patch_size, MPI_INT, comm_insidenode, ierr)
    IF (ierr /= 0) CALL comin_plugin_finish('gathering global indices of patches (p3plugin)', 'failed!')


    ! change name of main nwp namelist here (&master_model_nml:model_namelist_filename):
    icon_namelist_name = 'NAMELIST_NWP'

    ! read p3 namelist
    OPEN(newunit=funit, file=TRIM(icon_namelist_name), action='read', form='formatted')
    READ(funit, nml=p3_nml)
    CLOSE(funit)

    IF (rank_world == 0) WRITE (0,'(a)')    ' read P3 settings:'
    IF (rank_world == 0) WRITE (0,'(a,i1)') ' n_icecat  =', n_icecat
    IF (rank_world == 0) WRITE (0,'(a,l)')  ' l3mom_ice =', l3mom_ice
    IF (rank_world == 0) WRITE (0,'(a,l)')  ' lliqfrac  =', lliqfrac
    IF (rank_world == 0) WRITE (0,'(a,i1)') ' autoAccr_param_in  =', autoAccr_param_in

    ALLOCATE(p3_vars(n_icecat))
    ALLOCATE(p3_tracer(n_icecat))
    ALLOCATE(p3_tracer_ddt_turb(n_icecat))
!! JM_20260331 >> allocate memory for new comin handle for ice-phase diagnostics
    ALLOCATE(p3_ice_diag(n_icecat))
!! << JM_20260331


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
!! JM_20260323 >> adding new diagnostic
    ! --- warm-rain process rates ---
    ! nucleation
    CALL create_var('d_qcnuc',       'kg kg-1', '3d', 'dp') ! activation of cloud droplets from CCN
    CALL create_var('d_ncnuc',       '#  kg-1', '3d', 'dp') ! change in cloud droplet number from activation of CCN
    ! condensation
    CALL create_var('d_qccon',       'kg kg-1', '3d', 'dp') ! cloud droplet condensation
    CALL create_var('d_qrcon',       'kg kg-1', '3d', 'dp') ! rain drop condensation
    ! evaporation
    CALL create_var('d_qcevp',       'kg kg-1', '3d', 'dp') ! cloud droplet evaporation
    CALL create_var('d_qrevp',       'kg kg-1', '3d', 'dp') ! rain evaporation
    CALL create_var('d_nrevp',       '#  kg-1', '3d', 'dp') ! change in rain drop number from evaporation
    ! collection (warm)
    CALL create_var('d_qcacc',       'kg kg-1', '3d', 'dp') ! cloud droplet accretion (acc) by rain
    CALL create_var('d_ncacc',       '#  kg-1', '3d', 'dp') ! change in cloud droplet number from accretion by rain
    CALL create_var('d_qcaut',       'kg kg-1', '3d', 'dp') ! cloud droplet autoconversion (aut) to rain
    CALL create_var('d_ncautc',      '#  kg-1', '3d', 'dp') ! change in cloud droplet number from autoconversion
    CALL create_var('d_ncautr',      '#  kg-1', '3d', 'dp') ! change in rain drop number from autoconversion of cloud water
    ! self-collection
    CALL create_var('d_ncslf',       '#  kg-1', '3d', 'dp') ! change in cloud droplet number from self-collection
    CALL create_var('d_nrslf',       '#  kg-1', '3d', 'dp') ! change in rain drop number from self-collection
!! << JM_20260323

    DO i_icecat = 1, n_icecat
      WRITE(icecat_name, '(a,i0)') 'dmean_i', i_icecat
      CALL create_var(icecat_name, 'm', '3d', 'dp')
      WRITE(icecat_name, '(a,i0)') 'deff_i', i_icecat
      CALL create_var(icecat_name, 'm', '3d', 'dp')
      WRITE(icecat_name, '(a,i0)') 'rho_i', i_icecat
      CALL create_var(icecat_name, 'kg m-3', '3d', 'dp')
      WRITE(icecat_name, '(a,i0)') 'vm_i', i_icecat
      CALL create_var(icecat_name, 'm s-1', '3d', 'dp')
!! JM_20260331 >> create new comin variabe (need to be adjusted here when adding more diagnostics)
      WRITE(icecat_name, '(a,i0)') 'd_qidep_', i_icecat
      CALL create_var(icecat_name, 'kg kg-1', '3d', 'dp')
      WRITE(icecat_name, '(a,i0)') 'd_qisub_', i_icecat
      CALL create_var(icecat_name, 'kg kg-1', '3d', 'dp')
      WRITE(icecat_name, '(a,i0)') 'd_qinuc_', i_icecat
      CALL create_var(icecat_name, 'kg kg-1', '3d', 'dp')
      WRITE(icecat_name, '(a,i0)') 'd_qchetc_', i_icecat
      CALL create_var(icecat_name, 'kg kg-1', '3d', 'dp')
      WRITE(icecat_name, '(a,i0)') 'd_qcheti_', i_icecat
      CALL create_var(icecat_name, 'kg kg-1', '3d', 'dp')
      WRITE(icecat_name, '(a,i0)') 'd_qrhetc_', i_icecat
      CALL create_var(icecat_name, 'kg kg-1', '3d', 'dp')
      WRITE(icecat_name, '(a,i0)') 'd_qrheti_', i_icecat
      CALL create_var(icecat_name, 'kg kg-1', '3d', 'dp')
      WRITE(icecat_name, '(a,i0)') 'd_qimlt_', i_icecat
      CALL create_var(icecat_name, 'kg kg-1', '3d', 'dp')
      WRITE(icecat_name, '(a,i0)') 'd_qccol_', i_icecat
      CALL create_var(icecat_name, 'kg kg-1', '3d', 'dp')
      WRITE(icecat_name, '(a,i0)') 'd_qrcol_', i_icecat
      CALL create_var(icecat_name, 'kg kg-1', '3d', 'dp')
      WRITE(icecat_name, '(a,i0)') 'd_qwgrth_', i_icecat
      CALL create_var(icecat_name, 'kg kg-1', '3d', 'dp')
      WRITE(icecat_name, '(a,i0)') 'd_qcshd_', i_icecat
      CALL create_var(icecat_name, 'kg kg-1', '3d', 'dp')
      WRITE(icecat_name, '(a,i0)') 'd_qcmul_', i_icecat
      CALL create_var(icecat_name, 'kg kg-1', '3d', 'dp')
      WRITE(icecat_name, '(a,i0)') 'd_qrmul_', i_icecat
      CALL create_var(icecat_name, 'kg kg-1', '3d', 'dp')
      WRITE(icecat_name, '(a,i0)') 'd_nimul_', i_icecat
      CALL create_var(icecat_name, '# kg-1', '3d', 'dp')
!! << JM_20260331
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
    CALL comin_callback_register(EP_ATM_INIT_FINALIZE, init_p3_and_tracer)
    CALL comin_callback_register(EP_ATM_MICROPHYSICS_BEFORE, p3_main_wrapper)
    !CALL comin_callback_register(EP_ATM_RADIATION_BEFORE, set_reff_before_rad)

  END SUBROUTINE comin_main


  ! --------------------------------------------------------------------
  ! ComIn secondary constructor to connect all vars and tracers
  ! --------------------------------------------------------------------
  SUBROUTINE secondary_constructor()  BIND(C)
    INTEGER       :: ep_init = EP_ATM_INIT_FINALIZE
    INTEGER       :: ep_turb = EP_ATM_TURBULENCE_AFTER
    INTEGER       :: ep_mp   = EP_ATM_MICROPHYSICS_BEFORE
    INTEGER       :: ep_reff = EP_ATM_RADIATION_BEFORE
    INTEGER       :: ep_nudg = EP_ATM_NUDGING_AFTER
    INTEGER       :: ep_out  = EP_ATM_WRITE_OUTPUT_BEFORE
    INTEGER       :: FR = COMIN_FLAG_READ
    INTEGER       :: FW = COMIN_FLAG_WRITE
    INTEGER       :: id = 1
    INTEGER       :: i_icecat
    CHARACTER(20) :: icecat_name

    IF (rank_world == 0) WRITE (0,*) 'run secondary constructor'

    CALL comin_var_get([ep_init], t_comin_var_descriptor('z_mc', id), FR, dyn_vars%hfl)
    CALL comin_var_get([ep_mp], t_comin_var_descriptor('ddqz_z_full', id), FR, dyn_vars%dz)
    CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor('temp', id), IOR(FR, FW), dyn_vars%temp)
    CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor('rho', id), FR, dyn_vars%rho)
    CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor('pres', id), FR, dyn_vars%pres)
    CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor('exner', id), FR, dyn_vars%exner)
    CALL comin_var_get([ep_mp], t_comin_var_descriptor('theta_old', id), IOR(FR, FW), dyn_vars%theta_old)
    CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor('ddt_temp_phys', id), IOR(FR, FW), dyn_vars%ddt_temp_phys)
    CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor('w', id), FR, dyn_vars%w_hl)

    CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor('qv', id), IOR(FR, FW), icon_tracer%qv)
    CALL comin_var_get([ep_mp], t_comin_var_descriptor('qv_old', id), IOR(FR, FW), icon_tracer%qv_old)
    CALL comin_var_get([ep_init, ep_mp, ep_out], t_comin_var_descriptor('qc', id), IOR(FR, FW), icon_tracer%qc)
    CALL comin_var_get([ep_init, ep_mp, ep_out], t_comin_var_descriptor('qr', id), IOR(FR, FW), icon_tracer%qr)
    CALL comin_var_get([ep_init, ep_turb, ep_mp, ep_out], t_comin_var_descriptor('qnc', id), IOR(FR, FW), icon_tracer%qnc)
    CALL comin_var_get([ep_init, ep_turb, ep_mp, ep_out], t_comin_var_descriptor('qnr', id), IOR(FR, FW), icon_tracer%qnr)
    CALL comin_var_get([ep_init, ep_mp, ep_nudg, ep_out], t_comin_var_descriptor('qi', id), IOR(FR, FW), icon_tracer%qi)
    CALL comin_var_get([ep_init, ep_mp, ep_nudg, ep_out], t_comin_var_descriptor('qs', id), IOR(FR, FW), icon_tracer%qs)
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
    !! JM_20260323 >> adding new diagnostics
    ! --- warm-rain process rates ---
    ! nucleation
    CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor('d_qcnuc', id), IOR(FR, FW), mp_vars%d_qcnuc)
    CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor('d_ncnuc', id), IOR(FR, FW), mp_vars%d_ncnuc)
    ! condensation
    CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor('d_qccon', id), IOR(FR, FW), mp_vars%d_qccon)
    CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor('d_qrcon', id), IOR(FR, FW), mp_vars%d_qrcon)
    ! evaporation
    CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor('d_qcevp', id), IOR(FR, FW), mp_vars%d_qcevp)
    CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor('d_qrevp', id), IOR(FR, FW), mp_vars%d_qrevp)
    CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor('d_nrevp', id), IOR(FR, FW), mp_vars%d_nrevp)
    ! collection (warm)
    CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor('d_qcacc',  id), IOR(FR, FW), mp_vars%d_qcacc)
    CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor('d_ncacc',  id), IOR(FR, FW), mp_vars%d_ncacc)
    CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor('d_qcaut', id), IOR(FR, FW), mp_vars%d_qcaut)
    CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor('d_ncautc', id), IOR(FR, FW), mp_vars%d_ncautc)
    CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor('d_ncautr', id), IOR(FR, FW), mp_vars%d_ncautr)
    ! self-collection
    CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor('d_ncslf', id), IOR(FR, FW), mp_vars%d_ncslf)
    CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor('d_nrslf', id), IOR(FR, FW), mp_vars%d_nrslf)
!! << JM_20260323

    CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor('prec_gsp_rate', id), IOR(FR, FW), mp_vars%prec_gsp_rate)
    CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor('rain_gsp_rate', id), IOR(FR, FW), mp_vars%rain_gsp_rate)
    CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor('snow_gsp_rate', id), IOR(FR, FW), mp_vars%snow_gsp_rate)
    CALL comin_var_get([ep_init, ep_out], t_comin_var_descriptor('ice_gsp_rate', id), IOR(FR, FW), mp_vars%ice_gsp_rate)
    CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor('prec_gsp', id), IOR(FR, FW), mp_vars%prec_gsp)
    CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor('prec_gsp_d', id), IOR(FR, FW), mp_vars%prec_gsp_d)
    CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor('rain_gsp', id), IOR(FR, FW), mp_vars%rain_gsp)
    CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor('snow_gsp', id), IOR(FR, FW), mp_vars%snow_gsp)
    CALL comin_var_get([ep_init, ep_out], t_comin_var_descriptor('ice_gsp', id), IOR(FR, FW), mp_vars%ice_gsp)
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
!! JM_20260331 >> registering new ice-phase diagnostics (need to be adjusted here when adding more diagnostics)
      WRITE(icecat_name, '(a,i0)') 'd_qidep_', i_icecat
      CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor(TRIM(icecat_name), id), IOR(FR, FW), p3_ice_diag(i_icecat)%d_qidep_)
      WRITE(icecat_name, '(a,i0)') 'd_qisub_', i_icecat
      CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor(TRIM(icecat_name), id), IOR(FR, FW), p3_ice_diag(i_icecat)%d_qisub_)
      WRITE(icecat_name, '(a,i0)') 'd_qinuc_', i_icecat
      CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor(TRIM(icecat_name), id), IOR(FR, FW), p3_ice_diag(i_icecat)%d_qinuc_)
      WRITE(icecat_name, '(a,i0)') 'd_qchetc_', i_icecat
      CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor(TRIM(icecat_name), id), IOR(FR, FW), p3_ice_diag(i_icecat)%d_qchetc_)
      WRITE(icecat_name, '(a,i0)') 'd_qcheti_', i_icecat
      CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor(TRIM(icecat_name), id), IOR(FR, FW), p3_ice_diag(i_icecat)%d_qcheti_)
      WRITE(icecat_name, '(a,i0)') 'd_qrhetc_', i_icecat
      CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor(TRIM(icecat_name), id), IOR(FR, FW), p3_ice_diag(i_icecat)%d_qrhetc_)
      WRITE(icecat_name, '(a,i0)') 'd_qrheti_', i_icecat
      CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor(TRIM(icecat_name), id), IOR(FR, FW), p3_ice_diag(i_icecat)%d_qrheti_)
      WRITE(icecat_name, '(a,i0)') 'd_qimlt_', i_icecat
      CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor(TRIM(icecat_name), id), IOR(FR, FW), p3_ice_diag(i_icecat)%d_qimlt_)
      WRITE(icecat_name, '(a,i0)') 'd_qccol_', i_icecat
      CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor(TRIM(icecat_name), id), IOR(FR, FW), p3_ice_diag(i_icecat)%d_qccol_)
      WRITE(icecat_name, '(a,i0)') 'd_qrcol_', i_icecat
      CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor(TRIM(icecat_name), id), IOR(FR, FW), p3_ice_diag(i_icecat)%d_qrcol_)
      WRITE(icecat_name, '(a,i0)') 'd_qwgrth_', i_icecat
      CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor(TRIM(icecat_name), id), IOR(FR, FW), p3_ice_diag(i_icecat)%d_qwgrth_)
      WRITE(icecat_name, '(a,i0)') 'd_qcshd_', i_icecat
      CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor(TRIM(icecat_name), id), IOR(FR, FW), p3_ice_diag(i_icecat)%d_qcshd_)
      WRITE(icecat_name, '(a,i0)') 'd_qcmul_', i_icecat
      CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor(TRIM(icecat_name), id), IOR(FR, FW), p3_ice_diag(i_icecat)%d_qcmul_)
      WRITE(icecat_name, '(a,i0)') 'd_qrmul_', i_icecat
      CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor(TRIM(icecat_name), id), IOR(FR, FW), p3_ice_diag(i_icecat)%d_qrmul_)
      WRITE(icecat_name, '(a,i0)') 'd_nimul_', i_icecat
      CALL comin_var_get([ep_mp, ep_out], t_comin_var_descriptor(TRIM(icecat_name), id), IOR(FR, FW), p3_ice_diag(i_icecat)%d_nimul_)
!! << JM_20260331

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

END MODULE p3plugin
