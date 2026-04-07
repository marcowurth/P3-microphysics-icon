
MODULE p3plugin_global_vars
  USE comin_plugin_interface,  ONLY : t_comin_descrdata_global, t_comin_descrdata_domain

  USE p3plugin_types,          ONLY : t_dyn_vars_handle, t_mp_vars_handle, t_p3_vars_handle,   &
    &                                 t_icon_tracer_handle, t_p3_tracer_handle

!! JM_20260407 >> adding new comin handle type for 2moment ice-phase, 2mom ice-ice collision, 2mom ice-liquid and 3mom ice-phase diagnostics
  USE p3plugin_types,          ONLY : t_p3_ice_diag_2mom_handle, t_p3_ice_diag_2mom_coll_handle, t_p3_ice_diag_2mom_liqfrac_handle, &
                                      t_p3_ice_diag_3mom_handle
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
!! JM_20260407 >> defining new comin handle for 2moment ice-phase, 2mom ice-ice collision, 2mom ice-liquid and 3mom ice-phase diagnostics
  PUBLIC :: p3_ice_diag_2mom, p3_ice_diag_2mom_coll, p3_ice_diag_2mom_liqfrac, p3_ice_diag_3mom
!! << JM_20260407

  INTEGER        :: comm_world, comm_insidenode, rank_world, rank_insidenode
  INTEGER        :: numprocs_insidenode, max_patch_size
  INTEGER        :: n_icecat, itracer_ini, fastphystep
  REAL           :: dtime
  LOGICAL        :: l3mom_ice, lliqfrac
  CHARACTER(999) :: tracer_ini_filename, lookup_tables_path
!! JM_20260331 >> adding default value here, as it will be initialized with 0 when no value is given in the p3-nml block
  INTEGER        :: autoAccr_param_in = 2
!! << JM_20260331

  INTEGER, ALLOCATABLE :: node_patches_sizes(:)
  INTEGER, ALLOCATABLE :: node_patches_idx(:, :)

  TYPE(t_comin_descrdata_global), POINTER :: p_global
  TYPE(t_comin_descrdata_domain), POINTER :: p_patch

  TYPE(t_dyn_vars_handle)               :: dyn_vars
  TYPE(t_mp_vars_handle)                :: mp_vars
  TYPE(t_p3_vars_handle), ALLOCATABLE   :: p3_vars(:)
  TYPE(t_icon_tracer_handle)            :: icon_tracer, icon_tracer_ddt_turb
  TYPE(t_p3_tracer_handle), ALLOCATABLE :: p3_tracer(:), p3_tracer_ddt_turb(:)

!! JM_20260407 >> defining new comin handle for 2moment ice-phase, 2mom ice-ice collision, 2mom ice-liquid and 3mom ice-phase diagnostics
  TYPE(t_p3_ice_diag_2mom_handle), ALLOCATABLE         :: p3_ice_diag_2mom(:)
  TYPE(t_p3_ice_diag_2mom_coll_handle), ALLOCATABLE    :: p3_ice_diag_2mom_coll(:,:)
  TYPE(t_p3_ice_diag_2mom_liqfrac_handle), ALLOCATABLE :: p3_ice_diag_2mom_liqfrac(:)
  TYPE(t_p3_ice_diag_3mom_handle), ALLOCATABLE         :: p3_ice_diag_3mom(:)
!! << JM_20260407

END MODULE p3plugin_global_vars
