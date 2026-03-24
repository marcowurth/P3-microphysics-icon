
MODULE p3plugin_tracer_init
  USE mpi,                     ONLY : MPI_Wtime, MPI_Bcast, MPI_Scatter, MPI_INT, MPI_REAL, MPI_LOGICAL

  USE netcdf,                  ONLY : nf90_open, nf90_close, nf90_inq_dimid, nf90_inquire,                  &
    &                                 nf90_inq_varid, nf90_inquire_variable, nf90_inquire_dimension,        &
    &                                 nf90_get_var, NF90_FLOAT, NF90_DOUBLE,                                &
    &                                 NF90_NOWRITE, NF90_NOERR, NF90_MAX_VAR_DIMS

  USE comin_plugin_interface,  ONLY : comin_descrdata_get_timesteplength, comin_descrdata_get_cell_indices, &
    &                                 comin_plugin_finish

  USE p3plugin_utils,          ONLY : uppercase, lowercase
  USE p3plugin_types,          ONLY : t_dyn_vars_3dptr, t_icon_tracer_3dptr, t_mp_vars_3dptr,               &
    &                                 t_p3_tracer_3dptr
  USE p3plugin_global_vars,    ONLY : comm_world, comm_insidenode, rank_world, rank_insidenode,             &
    &                                 numprocs_insidenode, max_patch_size,                                  &
    &                                 node_patches_sizes, node_patches_idx,                                 &
    &                                 dtime, fastphystep, n_icecat, itracer_ini, l3mom_ice, lliqfrac,       &
    &                                 tracer_ini_filename, lookup_tables_path, p_global, p_patch,           &
    &                                 dyn_vars, icon_tracer, mp_vars, p3_tracer

  USE microphy_p3,             ONLY : p3_init, status_ok

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: init_p3_and_tracer

CONTAINS

  ! ------------------------------------------------------------------------
  ! Call p3_init to load lookup tables etc. and then initialize P3 tracers
  ! ------------------------------------------------------------------------
  SUBROUTINE init_p3_and_tracer()  BIND(C)

    CHARACTER(16)     :: model        = 'ICON'
    CHARACTER(30)     :: varname      = ''
    CHARACTER(20)     :: icecat_name
    LOGICAL           :: abort_on_err = .TRUE.
    LOGICAL           :: dowr         = .FALSE.
    LOGICAL           :: lvarfound, l3mom_ice_ini, lliqfrac_ini
    INTEGER           :: stat, i_icecat, n_icecat_ini
    INTEGER           :: jg, jb, jk, jc
    INTEGER           :: i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end
    DOUBLE PRECISION  :: start, finish
    REAL              :: dmean_qc, dmean_qr, dmean_qi, dmean_qs, dmean_qg
    REAL              :: magicfac_qi, magicfac_qs, magicfac_qg
    REAL              :: rhow, rhop_qi, rhop_qs, rhop_qg, rhor_qs, rhor_qg
    REAL              :: frim_qs, frim_qg
    REAL, ALLOCATABLE :: hfl_3dpatch(:, :, :)
    REAL, ALLOCATABLE :: qnc_ini_3dpatch(:, :, :)
    REAL, ALLOCATABLE :: qr_ini_3dpatch(:, :, :)
    REAL, ALLOCATABLE :: qnr_ini_3dpatch(:, :, :)
    REAL, ALLOCATABLE :: qs_ini_3dpatch(:, :, :)
    REAL, ALLOCATABLE :: qg_ini_3dpatch(:, :, :)
    REAL, ALLOCATABLE :: qitot_ini_3dpatch(:, :, :, :)
    REAL, ALLOCATABLE :: qnitot_ini_3dpatch(:, :, :, :)
    REAL, ALLOCATABLE :: qirim_ini_3dpatch(:, :, :, :)
    REAL, ALLOCATABLE :: birim_ini_3dpatch(:, :, :, :)
    REAL, ALLOCATABLE :: qzitot_ini_3dpatch(:, :, :, :)
    REAL, ALLOCATABLE :: qiliq_ini_3dpatch(:, :, :, :)

    TYPE(t_icon_tracer_3dptr) :: icon_tracer_3d
    TYPE(t_p3_tracer_3dptr)   :: p3_tracer_3d(n_icecat)
    TYPE(t_dyn_vars_3dptr)    :: dyn_vars_3d
    TYPE(t_mp_vars_3dptr)     :: mp_vars_3d

    IF (rank_world == 0) WRITE (0,*) 'call p3_init'

    CALL p3_init(TRIM(lookup_tables_path), n_icecat, l3mom_ice, lliqfrac, model, stat, abort_on_err, dowr, autoAccr_param_in)
    IF (stat /= status_ok) CALL comin_plugin_finish('init_p3_and_tracer (p3plugin)', 'calling failed!')

    dtime = comin_descrdata_get_timesteplength(1)
    fastphystep = 1

    CALL mp_vars%ice_gsp_rate%to_3d(mp_vars_3d%ice_gsp_rate)
    CALL mp_vars%ice_gsp%to_3d(mp_vars_3d%ice_gsp)

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

    start = MPI_Wtime()

    SELECT CASE (itracer_ini)
    CASE (0)
      IF (rank_world == 0) WRITE (0,*) 'initialize without clouds and precipitation, "dry"'

    CASE (1)
      IF (rank_world == 0) WRITE (0,*) 'initialize from 1M-scheme mass tracers qc, qr, qi, qs and if available qg'

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

      ALLOCATE(hfl_3dpatch(p_global%nproma, p_patch%nlev, p_patch%cells%nblks))
      ALLOCATE(qr_ini_3dpatch(p_global%nproma, p_patch%nlev, p_patch%cells%nblks))
      ALLOCATE(qs_ini_3dpatch(p_global%nproma, p_patch%nlev, p_patch%cells%nblks))
      ALLOCATE(qg_ini_3dpatch(p_global%nproma, p_patch%nlev, p_patch%cells%nblks))

      CALL dyn_vars%hfl%to_3d(dyn_vars_3d%hfl)
      hfl_3dpatch = dyn_vars_3d%hfl

      CALL read_vinterp_ini_var(tracer_ini_filename, 'qr', hfl_3dpatch, qr_ini_3dpatch)
      CALL read_vinterp_ini_var(tracer_ini_filename, 'qs', hfl_3dpatch, qs_ini_3dpatch)

      CALL query_file_for_var(tracer_ini_filename, 'qg', lvarfound)

      IF (lvarfound) THEN
        IF (rank_world == 0) WRITE (0,*) 'found tracer qg in tracer_ini_filename, reading it in...'
        CALL read_vinterp_ini_var(tracer_ini_filename, 'qg', hfl_3dpatch, qg_ini_3dpatch)
      ELSE
        IF (rank_world == 0) WRITE (0,*) 'did not find tracer qg in tracer_ini_filename, initialize P3 ice with qi & qs only'
      ENDIF

    CASE (3)
      IF (rank_world == 0) WRITE (0,*) 'initialize from P3 tracers'

      ALLOCATE(hfl_3dpatch(p_global%nproma, p_patch%nlev, p_patch%cells%nblks))
      ALLOCATE(qnc_ini_3dpatch(p_global%nproma, p_patch%nlev, p_patch%cells%nblks))
      ALLOCATE(qr_ini_3dpatch(p_global%nproma, p_patch%nlev, p_patch%cells%nblks))
      ALLOCATE(qnr_ini_3dpatch(p_global%nproma, p_patch%nlev, p_patch%cells%nblks))

      CALL dyn_vars%hfl%to_3d(dyn_vars_3d%hfl)
      hfl_3dpatch = dyn_vars_3d%hfl

      CALL read_vinterp_ini_var(tracer_ini_filename, 'qnc', hfl_3dpatch, qnc_ini_3dpatch)
      CALL read_vinterp_ini_var(tracer_ini_filename, 'qr', hfl_3dpatch, qr_ini_3dpatch)
      CALL read_vinterp_ini_var(tracer_ini_filename, 'qnr', hfl_3dpatch, qnr_ini_3dpatch)

      n_icecat_ini = 0
      DO i_icecat = 1, n_icecat
        WRITE(icecat_name, '(a,i0)') 'qitot_', i_icecat
        CALL query_file_for_var(tracer_ini_filename, TRIM(icecat_name), lvarfound)
        IF (lvarfound) n_icecat_ini = i_icecat
      END DO

      IF (rank_world == 0) WRITE (0,'(a,i0)') 'found number of P3 ice categories in tracer_ini_filename: ', n_icecat_ini
      IF (n_icecat < n_icecat_ini) THEN
        IF (rank_world == 0) WRITE (0,*) 'found more ini P3 ice categories than used in this model run, will ignore additional fields'
      ELSE IF (n_icecat > n_icecat_ini) THEN
        IF (n_icecat_ini == 0) CALL comin_plugin_finish('init_p3_and_tracer (p3plugin)', 'found no P3 ice categories!')
        IF (rank_world == 0) WRITE (0,*) 'found less P3 ice categories than used in this model run, will initialize the rest empty'
      ENDIF

      l3mom_ice_ini = .FALSE.
      IF (l3mom_ice) THEN
        CALL query_file_for_var(tracer_ini_filename, 'qzitot_1', lvarfound)
        IF (lvarfound) THEN
          l3mom_ice_ini = .TRUE.
          IF (rank_world == 0) WRITE (0,*) 'found qzitot, reading it in...'
        ELSE
          IF (rank_world == 0) WRITE (0,*) 'did not find qzitot, leaving it zero...'
        ENDIF
      ENDIF

      lliqfrac_ini = .FALSE.
      IF (lliqfrac) THEN
        CALL query_file_for_var(tracer_ini_filename, 'qiliq_1', lvarfound)
        IF (lvarfound) THEN
          lliqfrac_ini = .TRUE.
          IF (rank_world == 0) WRITE (0,*) 'found qiliq, reading it in...'
        ELSE
          IF (rank_world == 0) WRITE (0,*) 'did not find qiliq, leaving it zero...'
        ENDIF
      ENDIF

      ALLOCATE(qitot_ini_3dpatch(p_global%nproma, p_patch%nlev, p_patch%cells%nblks, n_icecat_ini))
      ALLOCATE(qnitot_ini_3dpatch(p_global%nproma, p_patch%nlev, p_patch%cells%nblks, n_icecat_ini))
      ALLOCATE(qirim_ini_3dpatch(p_global%nproma, p_patch%nlev, p_patch%cells%nblks, n_icecat_ini))
      ALLOCATE(birim_ini_3dpatch(p_global%nproma, p_patch%nlev, p_patch%cells%nblks, n_icecat_ini))
      IF (l3mom_ice_ini) THEN
        ALLOCATE(qzitot_ini_3dpatch(p_global%nproma, p_patch%nlev, p_patch%cells%nblks, n_icecat_ini))
      ENDIF
      IF (l3mom_ice_ini) THEN
        ALLOCATE(qiliq_ini_3dpatch(p_global%nproma, p_patch%nlev, p_patch%cells%nblks, n_icecat_ini))
      ENDIF

      DO i_icecat = 1, n_icecat_ini
        WRITE(icecat_name, '(a,i0)') 'qitot_', i_icecat
        CALL read_vinterp_ini_var(tracer_ini_filename, icecat_name, hfl_3dpatch, qitot_ini_3dpatch(:, :, :, i_icecat))

        WRITE(icecat_name, '(a,i0)') 'qnitot_', i_icecat
        CALL read_vinterp_ini_var(tracer_ini_filename, icecat_name, hfl_3dpatch, qnitot_ini_3dpatch(:, :, :, i_icecat))

        WRITE(icecat_name, '(a,i0)') 'qirim_', i_icecat
        CALL read_vinterp_ini_var(tracer_ini_filename, icecat_name, hfl_3dpatch, qirim_ini_3dpatch(:, :, :, i_icecat))

        WRITE(icecat_name, '(a,i0)') 'birim_', i_icecat
        CALL read_vinterp_ini_var(tracer_ini_filename, icecat_name, hfl_3dpatch, birim_ini_3dpatch(:, :, :, i_icecat))

        IF (l3mom_ice_ini) THEN
          WRITE(icecat_name, '(a,i0)') 'qzitot_', i_icecat
          CALL read_vinterp_ini_var(tracer_ini_filename, icecat_name, hfl_3dpatch, qzitot_ini_3dpatch(:, :, :, i_icecat))
        ENDIF

        IF (l3mom_ice_ini) THEN
          WRITE(icecat_name, '(a,i0)') 'qiliq_', i_icecat
          CALL read_vinterp_ini_var(tracer_ini_filename, icecat_name, hfl_3dpatch, qiliq_ini_3dpatch(:, :, :, i_icecat))
        ENDIF
      END DO

    END SELECT

    finish = MPI_Wtime()
    IF (rank_insidenode == 0) &
     & WRITE (0,'(a,i3,a,F7.3,a)') 'whole file reading with rank_world ', rank_world, ' completed in ', finish - start, 'sec'

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

          SELECT CASE (itracer_ini)
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
            ! qv, qc were already read from file and set in initicon

            icon_tracer_3d%qr(jc,jk,jb) = qr_ini_3dpatch(jc,jk,jb)

            icon_tracer_3d%qnc(jc,jk,jb) = icon_tracer_3d%qc(jc,jk,jb) / (rhow*3.14) * dmean_qc**-3
            icon_tracer_3d%qnr(jc,jk,jb) = icon_tracer_3d%qr(jc,jk,jb) / (rhow*3.14) * dmean_qr**-3
            ! dmean_qc = (qc / (qnc*rhow*3.14))**(1./3.)

            ! initialization of cold phase:
            ! qi was already read from file and set in initicon, use it as ice category 1
            ! use the read in 1M-scheme mass tracer qs [and qg]
            ! if n_icecat == 1: init only cloud ice qi and ignore precipitating types qs, qg
            ! if n_icecat == 2: use the first icecat for qi and the second icecat for merged qs + qg
            ! if n_icecat >= 3: use one icecat for qi, qs, qg each
            ! all other icecats (if more available) are kept empty

            ! cloud ice, assume no riming part present: Fr = 0
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

            ! precipitating ice, use two separate icecats if available (n_icecat > 2), if not merge into one.
            ! in case QG was not in the ini file, qg_ini_3dpatch equals zero and ice category 2 is only snow
            IF (lvarfound == .FALSE.) THEN
              qg_ini_3dpatch(jc,jk,jb) = 0.0
            ENDIF

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

            ELSE IF (n_icecat > 2) THEN
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

          CASE (3)
            ! initialization of warm phase:
            ! qv, qc were already read from file and set in initicon

            icon_tracer_3d%qnc(jc,jk,jb) = qnc_ini_3dpatch(jc,jk,jb)
            icon_tracer_3d%qr(jc,jk,jb) = qr_ini_3dpatch(jc,jk,jb)
            icon_tracer_3d%qnr(jc,jk,jb) = qnr_ini_3dpatch(jc,jk,jb)

            ! initialization of cold phase:
            ! use read in P3 tracers, leave additional categories empty
            DO i_icecat = 1, n_icecat_ini
              p3_tracer_3d(i_icecat)%qitot(jc,jk,jb) = qitot_ini_3dpatch(jc,jk,jb,i_icecat)
              p3_tracer_3d(i_icecat)%qnitot(jc,jk,jb) = qnitot_ini_3dpatch(jc,jk,jb,i_icecat)
              p3_tracer_3d(i_icecat)%qirim(jc,jk,jb) = qirim_ini_3dpatch(jc,jk,jb,i_icecat)
              p3_tracer_3d(i_icecat)%birim(jc,jk,jb) = birim_ini_3dpatch(jc,jk,jb,i_icecat)
              IF (l3mom_ice) THEN
                IF (l3mom_ice_ini) THEN
                  p3_tracer_3d(i_icecat)%qzitot(jc,jk,jb) = qzitot_ini_3dpatch(jc,jk,jb,i_icecat)
                ELSE
                  p3_tracer_3d(i_icecat)%qzitot(jc,jk,jb) = 0.0
                ENDIF
              ENDIF
              IF (lliqfrac) THEN
                IF (lliqfrac_ini) THEN
                  p3_tracer_3d(i_icecat)%qiliq(jc,jk,jb) = qiliq_ini_3dpatch(jc,jk,jb,i_icecat)
                ELSE
                  p3_tracer_3d(i_icecat)%qiliq(jc,jk,jb) = 0.0
                ENDIF
              ENDIF
            END DO

            IF (n_icecat > n_icecat_ini) THEN
              DO i_icecat = n_icecat_ini + 1, n_icecat
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


    ! set ice_gsp_rate to zero, it is needed in the sfc interface for terra (blowing snow for snow tiles)
    mp_vars_3d%ice_gsp_rate(:,:,1) = 0.0
    mp_vars_3d%ice_gsp(:,:,1) = 0.0


    CALL dyn_vars_3d%nullify()
    CALL mp_vars_3d%nullify()
    CALL icon_tracer_3d%nullify()
    DO i_icecat = 1, n_icecat
      CALL p3_tracer_3d(i_icecat)%nullify()
    END DO

  END SUBROUTINE init_p3_and_tracer


  SUBROUTINE query_file_for_var(filename, varname, lvarfound)
    CHARACTER(*), INTENT(IN) :: filename, varname
    LOGICAL, INTENT(OUT)     :: lvarfound

    INTEGER                  :: nc_status, ncid, varid, ierr

    IF (rank_world == 0) THEN
      ncid = -99
      nc_status = nf90_open(TRIM(filename), NF90_NOWRITE, ncid)
      IF (nc_status /= NF90_NOERR) &
        & CALL comin_plugin_finish('query_file_for_var (p3plugin)', 'Could not read file: ' // TRIM(filename))

      lvarfound = .FALSE.
      nc_status = nf90_inq_varid(ncid, varname, varid)
      IF (nc_status /= NF90_NOERR) &
        & nc_status = nf90_inq_varid(ncid, uppercase(TRIM(varname)), varid)
      IF (nc_status /= NF90_NOERR) &
        & nc_status = nf90_inq_varid(ncid, lowercase(TRIM(varname)), varid)
      IF (nc_status == NF90_NOERR) lvarfound = .TRUE.
    END IF

    CALL MPI_Bcast(lvarfound, 1, MPI_LOGICAL, 0, comm_world, ierr)
    IF (ierr /= 0) CALL comin_plugin_finish('MPI_Bcast(lvarfound, 1, MPI_LOGICAL, 0, comm_world) (p3plugin)', 'failed!')

  END SUBROUTINE


  SUBROUTINE read_vinterp_ini_var(filename, varname, hfl_3dpatch_outlevs, var_3dpatch_outlevs)
    CHARACTER(*), INTENT(IN) :: filename, varname
    REAL, INTENT(IN)         :: hfl_3dpatch_outlevs(:, :, :)
    REAL, INTENT(OUT)        :: var_3dpatch_outlevs(:, :, :)

    DOUBLE PRECISION         :: start, finish
    CHARACTER(30)            :: dimname
    INTEGER                  :: nc_status, ncid, dimid_ncells, ncells_global, nlev_in, nlev_in_hhl, ierr
    INTEGER                  :: i, jg, jb, jk, jc, jlocal
    INTEGER                  :: i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end
    REAL, ALLOCATABLE        :: var_patches_padded(:, :, :), hhl_patches_padded(:, :, :)
    REAL, ALLOCATABLE        :: var_patch_padded(:, :), hhl_patch_padded(:, :)
    REAL, ALLOCATABLE        :: hfl_1d_inlevs(:)

    ! file reading done by node head process only (therefore once per node)
    IF (rank_insidenode == 0) THEN
      IF (rank_world == 0) WRITE (0,'(a,a)') 'read from file, varname: ', TRIM(varname)
      start = MPI_Wtime()

      ncid = -99
      nc_status = nf90_open(TRIM(filename), NF90_NOWRITE, ncid)
      IF (nc_status /= NF90_NOERR) &
        & CALL comin_plugin_finish('read_vinterp_ini_var (p3plugin)', 'Could not read ini file: ' // TRIM(filename))

      ! read ncells of ini file and check if equals to icon model's ncells_global
      nc_status = nf90_inq_dimid(ncid, 'ncells', dimid_ncells)
      IF (nc_status /= NF90_NOERR) &
        & CALL comin_plugin_finish('read_vinterp_ini_var (p3plugin)', &
                                 & 'Could not find dimension ncells in ini file!')
      nc_status = nf90_inquire_dimension(ncid, dimid_ncells, dimname, ncells_global)
      IF (ncells_global /= p_patch%cells%ncells_global) &
        & CALL comin_plugin_finish('read_vinterp_ini_var (p3plugin)', &
                                 & 'ncells number of ini file does not match with model!')

      CALL read_netcdf_nlev(varname, ncid, nlev_in)
      CALL read_netcdf_nlev('hhl', ncid, nlev_in_hhl)

      ALLOCATE(var_patches_padded(max_patch_size, nlev_in, numprocs_insidenode))
      ALLOCATE(hhl_patches_padded(max_patch_size, nlev_in_hhl, numprocs_insidenode))

      CALL read_netcdf_var(varname, ncid, var_patches_padded)
      CALL read_netcdf_var('hhl', ncid, hhl_patches_padded)

      nc_status = nf90_close(ncid)
      IF (nc_status /= NF90_NOERR) &
        & CALL comin_plugin_finish('read_vinterp_ini_var (p3plugin)', 'File closing not successful!')

      finish = MPI_Wtime()
      WRITE (0,'(a,i3,a,F7.3,a)') 'file reading with rank_world ', rank_world, ' completed in ', finish - start, 'sec'
    ENDIF
    ! end of node head reading work

    ! communicate the number of levels, so the other MPI processes can iterate over it
    CALL MPI_Bcast(nlev_in, 1, MPI_INT, 0, comm_insidenode, ierr)
    IF (ierr /= 0) CALL comin_plugin_finish('MPI_Bcast(nlev_in, ...) (p3plugin)', 'failed!')
    CALL MPI_Bcast(nlev_in_hhl, 1, MPI_INT, 0, comm_insidenode, ierr)
    IF (ierr /= 0) CALL comin_plugin_finish('MPI_Bcast(nlev_in_hhl, ...) (p3plugin)', 'failed!')

    ! ! allocate patch fields in all MPI processes inside node
    ALLOCATE(var_patch_padded(max_patch_size, nlev_in))
    ALLOCATE(hhl_patch_padded(max_patch_size, nlev_in_hhl))

    ! start = MPI_Wtime()

    DO jk = 1, nlev_in
      CALL MPI_Scatter(var_patches_padded(:, jk, :), max_patch_size, MPI_REAL, &
        &              var_patch_padded(:, jk),      max_patch_size, MPI_REAL, 0, comm_insidenode, ierr)
        IF (ierr /= 0) CALL comin_plugin_finish('MPI_Scatter(var_patches_padded, ...) (p3plugin)', 'failed!')
    END DO
    DO jk = 1, nlev_in_hhl
      CALL MPI_Scatter(hhl_patches_padded(:, jk, :), max_patch_size, MPI_REAL, &
        &              hhl_patch_padded(:, jk),      max_patch_size, MPI_REAL, 0, comm_insidenode, ierr)
      IF (ierr /= 0) CALL comin_plugin_finish('MPI_Scatter(hhl_patches_padded, ...) (p3plugin)', 'failed!')
    END DO

    finish = MPI_Wtime()
    IF (rank_insidenode == 0) &
     & WRITE (0,'(a,i3,a,F5.3,a)') 'communicating fields with rank_world ', rank_world, ' completed in ', finish - start, 'sec'

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
        jlocal = (jb-1)*p_global%nproma+jc
        hfl_1d_inlevs(:) = (hhl_patch_padded(jlocal, 1:nlev_in) + hhl_patch_padded(jlocal, 2:nlev_in+1)) / 2.

        CALL vert_intp_linear_1d(hfl_1d_inlevs(:),               var_patch_padded(jlocal, :),    &
          &                      hfl_3dpatch_outlevs(jc, :, jb), var_3dpatch_outlevs(jc, :, jb))
      END DO
    END DO

    IF (rank_insidenode == 0) DEALLOCATE(var_patches_padded, hhl_patches_padded)
    DEALLOCATE(var_patch_padded, hhl_patch_padded, hfl_1d_inlevs)

  END SUBROUTINE read_vinterp_ini_var


  SUBROUTINE read_netcdf_var(varname, ncid, var_patches_padded)
    CHARACTER(*), INTENT(IN) :: varname
    INTEGER, INTENT(IN)      :: ncid
    REAL, INTENT(INOUT)      :: var_patches_padded(:, :, :)

    INTEGER                  :: i, nc_status, varid, xtype, ndims, natts, ncells_global, dimpos_time, dimpos_ncells, nlev
    INTEGER                  :: read_start(3), read_count(3)
    INTEGER                  :: dimids(NF90_MAX_VAR_DIMS), dimlen(NF90_MAX_VAR_DIMS)
    INTEGER, ALLOCATABLE     :: varglobal_read_first(:, :), varglobal_read_last(:, :)
    REAL                     :: fillval
    CHARACTER(30)            :: varname_dummy, dimname
    LOGICAL                  :: ncells_first

    IF (rank_insidenode /= 0) &
      & CALL comin_plugin_finish('read_netcdf_var (p3plugin)', &
                               & 'error: calling routine not with node head process!')

    ! read in variable varname and handle dimensions
    nc_status = nf90_inq_varid(ncid, TRIM(varname), varid)
    IF (nc_status /= NF90_NOERR) &
      & nc_status = nf90_inq_varid(ncid, uppercase(TRIM(varname)), varid)
    IF (nc_status /= NF90_NOERR) &
      & nc_status = nf90_inq_varid(ncid, lowercase(TRIM(varname)), varid)
    IF (nc_status /= NF90_NOERR) &
      & CALL comin_plugin_finish('read_netcdf_var (p3plugin)', &
                               & 'Could not find variable "' // TRIM(varname))

    nc_status = nf90_inquire_variable(ncid, varid, varname_dummy, xtype, ndims, dimids, natts)
    IF (nc_status /= NF90_NOERR) &
      & CALL comin_plugin_finish('read_netcdf_var (p3plugin)', &
                               & 'Could not inquire variable: ' // TRIM(varname))
    IF (xtype /= NF90_DOUBLE .and. xtype /= NF90_FLOAT) &
      & CALL comin_plugin_finish('read_netcdf_var (p3plugin)', &
                               & 'Variable type not double or float: ' // TRIM(varname))

    fillval = -999.

    SELECT CASE (ndims)
    CASE (1)
      CALL comin_plugin_finish('read_netcdf_var (p3plugin)', &
                             & 'Variable has only one dimension: ' // TRIM(varname))
    CASE (2)
      dimpos_ncells = 0
      DO i = 1, 2
        nc_status = nf90_inquire_dimension(ncid, dimids(i), dimname, dimlen(i))
        IF (TRIM(dimname) == 'time') THEN
          CALL comin_plugin_finish('read_netcdf_var (p3plugin)', &
                                 & 'Found dimension "time" in 2d array: ' // TRIM(varname))
        ELSE IF (TRIM(dimname) == 'ncells') THEN
          dimpos_ncells = i
          ncells_global = dimlen(i)
        ELSE
          nlev = dimlen(i)
        ENDIF
      END DO

      IF (dimpos_ncells == 0) &
        & CALL comin_plugin_finish('read_netcdf_var (p3plugin)', &
                                 & 'Could not find dimension "ncells" in:' // TRIM(varname))
      IF (ncells_global /= p_patch%cells%ncells_global) &
        & CALL comin_plugin_finish('read_netcdf_var (p3plugin)', &
                                 & 'ncells number of ini file does not match with model!')

      SELECT CASE (dimpos_ncells)
      CASE (1)
        ALLOCATE(varglobal_read_first(ncells_global, nlev))
        nc_status = nf90_get_var(ncid, varid, varglobal_read_first)
        DO i = 1, numprocs_insidenode
          var_patches_padded(1:node_patches_sizes(i), :, i) &
            &  = varglobal_read_first(node_patches_idx(1:node_patches_sizes(i), i), :)
          var_patches_padded(node_patches_sizes(i)+1:max_patch_size, :, i) = fillval
        END DO
        DEALLOCATE(varglobal_read_first)
      CASE (2)
        ALLOCATE(varglobal_read_last(nlev, ncells_global))
        nc_status = nf90_get_var(ncid, varid, varglobal_read_last)
        DO i = 1, numprocs_insidenode
          var_patches_padded(1:node_patches_sizes(i), :, i) &
            &  = varglobal_read_last(:, node_patches_idx(1:node_patches_sizes(i), i))
          var_patches_padded(node_patches_sizes(i)+1:max_patch_size, :, i) = fillval
        END DO
        DEALLOCATE(varglobal_read_last)
      END SELECT

    CASE (3)
      dimpos_time = 0
      dimpos_ncells = 0
      DO i = 1, 3
        nc_status = nf90_inquire_dimension(ncid, dimids(i), dimname, dimlen(i))
        IF (TRIM(dimname) == 'time') THEN
          dimpos_time = i
          IF (dimlen(i) > 1) THEN
            WRITE (0,'(a)') 'Reading of 3D var: dimension "time" has more than one time step, choosing the first'
          ENDIF
        ELSE IF (TRIM(dimname) == 'ncells') THEN
          dimpos_ncells = i
          ncells_global = dimlen(i)
        ELSE
          nlev = dimlen(i)
        ENDIF
      END DO

      IF (dimpos_time == 0) &
        & CALL comin_plugin_finish('read_netcdf_var (p3plugin)', &
                                 & 'Could not find dimension "time" in:' // TRIM(varname))
      IF (dimpos_ncells == 0) &
        & CALL comin_plugin_finish('read_netcdf_var (p3plugin)', &
                                 & 'Could not find dimension "ncells" in:' // TRIM(varname))
      IF (ncells_global /= p_patch%cells%ncells_global) &
        & CALL comin_plugin_finish('read_netcdf_var (p3plugin)', &
                                 & 'ncells number of ini file does not match with model!')

      read_start = (/ 1, 1, 1 /)
      SELECT CASE (dimpos_ncells)
      CASE (1)
        SELECT CASE (dimpos_time)
        CASE (2)
          read_count = (/ ncells_global, 1, nlev /)
        CASE (3)
          read_count = (/ ncells_global, nlev, 1 /)
        END SELECT
        ncells_first = .TRUE.
      CASE (2)
        SELECT CASE (dimpos_time)
        CASE (1)
          read_count = (/ 1, ncells_global, nlev /)
          ncells_first = .TRUE.
        CASE (3)
          read_count = (/ nlev, ncells_global, 1 /)
          ncells_first = .FALSE.
        END SELECT
      CASE (3)
        SELECT CASE (dimpos_time)
        CASE (1)
          read_count = (/ 1, nlev, ncells_global /)
        CASE (2)
          read_count = (/ nlev, 1, ncells_global /)
        END SELECT
        ncells_first = .FALSE.
      END SELECT

      IF (ncells_first) THEN
        ALLOCATE(varglobal_read_first(ncells_global, nlev))
        nc_status = nf90_get_var(ncid, varid, varglobal_read_first, start=read_start, count=read_count)
        DO i = 1, numprocs_insidenode
          var_patches_padded(1:node_patches_sizes(i), :, i) &
            &  = varglobal_read_first(node_patches_idx(1:node_patches_sizes(i), i), :)
          var_patches_padded(node_patches_sizes(i)+1:max_patch_size, :, i) = fillval
        END DO
        DEALLOCATE(varglobal_read_first)
      ELSE
        ALLOCATE(varglobal_read_last(nlev, ncells_global))
        nc_status = nf90_get_var(ncid, varid, varglobal_read_last, start=read_start, count=read_count)
        DO i = 1, numprocs_insidenode
          var_patches_padded(1:node_patches_sizes(i), :, i) &
            &  = varglobal_read_last(:, node_patches_idx(1:node_patches_sizes(i), i))
          var_patches_padded(node_patches_sizes(i)+1:max_patch_size, :, i) = fillval
        END DO
        DEALLOCATE(varglobal_read_last)
      END IF

    CASE (4:)
      CALL comin_plugin_finish('read_netcdf_var (p3plugin)', &
                             & 'Variable has more than 3 dimensions: ' // TRIM(varname))
    END SELECT
    ! end of reading variable

  END SUBROUTINE read_netcdf_var


  SUBROUTINE read_netcdf_nlev(varname, ncid, nlev)
    CHARACTER(*), INTENT(IN) :: varname
    INTEGER, INTENT(IN)      :: ncid
    INTEGER, INTENT(OUT)     :: nlev

    INTEGER                  :: nc_status, varid, xtype, ndims, natts, i
    INTEGER                  :: dimids(NF90_MAX_VAR_DIMS), dimlen(NF90_MAX_VAR_DIMS)
    CHARACTER(30)            :: varname_dummy, dimname

    IF (rank_insidenode /= 0) &
      & CALL comin_plugin_finish('read_netcdf_nlev (p3plugin)', &
                               & 'error: calling routine not with node head process!')

    ! read in variable varname and handle dimensions
    nc_status = nf90_inq_varid(ncid, TRIM(varname), varid)
    IF (nc_status /= NF90_NOERR) &
      & nc_status = nf90_inq_varid(ncid, uppercase(TRIM(varname)), varid)
    IF (nc_status /= NF90_NOERR) &
      & nc_status = nf90_inq_varid(ncid, lowercase(TRIM(varname)), varid)
    IF (nc_status /= NF90_NOERR) &
      & CALL comin_plugin_finish('read_netcdf_nlev (p3plugin)', &
                               & 'Could not find variable "' // TRIM(varname))

    nc_status = nf90_inquire_variable(ncid, varid, varname_dummy, xtype, ndims, dimids, natts)
    IF (nc_status /= NF90_NOERR) &
      & CALL comin_plugin_finish('read_netcdf_nlev (p3plugin)', &
                               & 'Could not inquire variable: ' // TRIM(varname))
    IF (xtype /= NF90_DOUBLE .and. xtype /= NF90_FLOAT) &
      & CALL comin_plugin_finish('read_netcdf_nlev (p3plugin)', &
                               & 'Variable type not double or float: ' // TRIM(varname))

    SELECT CASE (ndims)
    CASE (1)
      CALL comin_plugin_finish('read_netcdf_nlev (p3plugin)', &
                             & 'Variable has only one dimension: ' // TRIM(varname))
    CASE (2)
      DO i = 1, 2
        nc_status = nf90_inquire_dimension(ncid, dimids(i), dimname, dimlen(i))
        IF (TRIM(dimname) /= 'time' .and. TRIM(dimname) /= 'ncells') THEN
          nlev = dimlen(i)
        ENDIF
      END DO

    CASE (3)
      DO i = 1, 3
        nc_status = nf90_inquire_dimension(ncid, dimids(i), dimname, dimlen(i))
        IF (TRIM(dimname) /= 'time' .and. TRIM(dimname) /= 'ncells') THEN
          nlev = dimlen(i)
        ENDIF
      END DO

    CASE (4:)
      CALL comin_plugin_finish('read_netcdf_nlev (p3plugin)', &
                             & 'Variable has more than 3 dimensions: ' // TRIM(varname))
    END SELECT

  END SUBROUTINE read_netcdf_nlev


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

END MODULE p3plugin_tracer_init
