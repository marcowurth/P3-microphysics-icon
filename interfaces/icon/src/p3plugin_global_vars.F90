
MODULE p3plugin_global_vars
  USE comin_plugin_interface,  ONLY : t_comin_descrdata_global, t_comin_descrdata_domain

  USE p3plugin_types,          ONLY : t_dyn_vars_handle, t_mp_vars_handle, t_p3_vars_handle,   &
    &                                 t_icon_tracer_handle, t_p3_tracer_handle

!! JM_20260407 >> adding new comin handle type for 2-moment warm-phase, 2-moment ice-phase, 2-moment ice-ice collision, 2-moment ice-liquid and 3-moment ice-phase diagnostics
  USE p3plugin_types,          ONLY : t_p3_diag_wrm_2mom_handle, t_p3_diag_ice_2mom_handle, &
                                      t_p3_diag_ice_2mom_coll_handle, t_p3_diag_ice_2mom_liqfrac_handle, &
                                      t_p3_diag_ice_3mom_handle
!! << JM_20260407

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: comm_world, comm_insidenode, rank_world, rank_insidenode
  PUBLIC :: numprocs_insidenode, max_patch_size, node_patches_sizes, node_patches_idx
  PUBLIC :: n_icecat, itracer_ini, fastphystep, dtime
  PUBLIC :: l3mom_ice, lliqfrac
  PUBLIC :: tracer_ini_filename, lookup_tables_path
  PUBLIC :: p_global, p_patch
  PUBLIC :: dyn_vars, mp_vars, p3_vars
  PUBLIC :: icon_tracer, p3_tracer
  PUBLIC :: icon_tracer_ddt_turb, p3_tracer_ddt_turb
  PUBLIC :: autoAccr_param_in
!! JM_20260629 >> adding n_inact as argument for depletion of INPs
  PUBLIC :: n_inact
!! << JM_20260629
!! JM_20260407 >> defining new comin handle for 2-moment warm-phase, 2-moment ice-phase, 2-moment ice-ice collision, 2-moment ice-liquid and 3-moment ice-phase diagnostics
  PUBLIC :: p3_diag_wrm_2mom, p3_diag_ice_2mom, p3_diag_ice_2mom_coll, p3_diag_ice_2mom_liqfrac, p3_diag_ice_3mom
!! << JM_20260407
!! JM_20260723 >> new integer for warm-/ice-phase diagnostics
  PUBLIC :: n_diag_wrm_2mom, n_diag_ice_2mom, n_diag_ice_2mom_coll, n_diag_ice_2mom_liqfrac, n_diag_ice_3mom
  PUBLIC :: diag_wrm_2mom, diag_ice_2mom, diag_ice_2mom_coll, diag_ice_2mom_liqfrac, diag_ice_3mom
!! << JM_20260723

  INTEGER        :: comm_world, comm_insidenode, rank_world, rank_insidenode
  INTEGER        :: numprocs_insidenode, max_patch_size
  INTEGER        :: n_icecat, itracer_ini, fastphystep
  REAL           :: dtime
  LOGICAL        :: l3mom_ice, lliqfrac
  CHARACTER(999) :: tracer_ini_filename, lookup_tables_path
  INTEGER        :: autoAccr_param_in = 2

!! JM_20260629 >> adding n_inact as argument for depletion of INPs
  REAL, ALLOCATABLE :: n_inact(:,:,:)
!! << JM_20260629
!! JM_20260723 >> trying to make warm-phase and ice-phase diagnostics accumulated
  INTEGER           :: n_diag_wrm_2mom = 7
  INTEGER           :: n_diag_ice_2mom = 39
  INTEGER           :: n_diag_ice_2mom_coll = 1
  INTEGER           :: n_diag_ice_2mom_liqfrac = 15
  INTEGER           :: n_diag_ice_3mom = 7
  REAL, ALLOCATABLE :: diag_wrm_2mom(:,:,:,:), diag_ice_2mom(:,:,:,:,:), diag_ice_2mom_coll(:,:,:,:,:,:), diag_ice_2mom_liqfrac(:,:,:,:,:), diag_ice_3mom(:,:,:,:,:)
!! << JM_20260723

  INTEGER, ALLOCATABLE :: node_patches_sizes(:)
  INTEGER, ALLOCATABLE :: node_patches_idx(:, :)

  TYPE(t_comin_descrdata_global), POINTER :: p_global
  TYPE(t_comin_descrdata_domain), POINTER :: p_patch

  TYPE(t_dyn_vars_handle)               :: dyn_vars
  TYPE(t_mp_vars_handle)                :: mp_vars
  TYPE(t_p3_vars_handle), ALLOCATABLE   :: p3_vars(:)
  TYPE(t_icon_tracer_handle)            :: icon_tracer, icon_tracer_ddt_turb
  TYPE(t_p3_tracer_handle), ALLOCATABLE :: p3_tracer(:), p3_tracer_ddt_turb(:)

!! JM_20260407 >> defining new comin handle for 2-moment warm-phase, 2-moment ice-phase, 2-moment ice-ice collision, 2-moment ice-liquid and 3-moment ice-phase diagnostics
  TYPE(t_p3_diag_wrm_2mom_handle)                      :: p3_diag_wrm_2mom
  TYPE(t_p3_diag_ice_2mom_handle), ALLOCATABLE         :: p3_diag_ice_2mom(:)
  TYPE(t_p3_diag_ice_2mom_coll_handle), ALLOCATABLE    :: p3_diag_ice_2mom_coll(:,:)
  TYPE(t_p3_diag_ice_2mom_liqfrac_handle), ALLOCATABLE :: p3_diag_ice_2mom_liqfrac(:)
  TYPE(t_p3_diag_ice_3mom_handle), ALLOCATABLE         :: p3_diag_ice_3mom(:)
!! << JM_20260407

END MODULE p3plugin_global_vars
