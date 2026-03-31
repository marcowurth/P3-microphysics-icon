
MODULE p3plugin_types
  USE comin_plugin_interface,  ONLY : t_comin_var_handle

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_dyn_vars_handle, t_dyn_vars_3dptr
  PUBLIC :: t_mp_vars_handle, t_mp_vars_3dptr
  PUBLIC :: t_p3_vars_handle, t_p3_vars_3dptr
  PUBLIC :: t_icon_tracer_handle, t_icon_tracer_3dptr
  PUBLIC :: t_p3_tracer_handle, t_p3_tracer_3dptr

  INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND(12,307)

  TYPE :: t_dyn_vars_handle
    TYPE(t_comin_var_handle)            :: hfl, dz, temp, rho, pres, exner, theta_old, ddt_temp_phys, w_hl
  END type t_dyn_vars_handle

  TYPE :: t_dyn_vars_3dptr
    REAL(wp), POINTER, DIMENSION(:,:,:) :: hfl, dz, temp, rho, pres, exner, theta_old, ddt_temp_phys, w_hl
  CONTAINS
    PROCEDURE :: nullify => nullify_dyn_vars_3dptr
  END type t_dyn_vars_3dptr

  TYPE :: t_mp_vars_handle
    TYPE(t_comin_var_handle)            :: dmean_c, dmean_r, deff_c, deff_i, reff_qc, reff_qi
    TYPE(t_comin_var_handle)            :: dhmax, dhmax_ground, ze_p3
    TYPE(t_comin_var_handle)            :: prec_gsp_rate, rain_gsp_rate, snow_gsp_rate, ice_gsp_rate
    TYPE(t_comin_var_handle)            :: prec_gsp, prec_gsp_d, rain_gsp, snow_gsp, ice_gsp
    TYPE(t_comin_var_handle)            :: q_sedim, twater
!! JM_20260323 >> adding new diagnostics
    ! --- warm-rain process rates ---
    TYPE(t_comin_var_handle)            :: d_qcnuc, d_ncnuc  ! nucleation
    TYPE(t_comin_var_handle)            :: d_qccon, d_qrcon  ! condensation
    TYPE(t_comin_var_handle)            :: d_qcevp, d_qrevp, d_nrevp ! evaporation
    TYPE(t_comin_var_handle)            :: d_qcacc, d_ncacc, d_qcaut, d_ncautc, d_ncautr ! collection (warm)
    TYPE(t_comin_var_handle)            :: d_ncslf, d_nrslf ! self-collection
!! << JM_20260323
  END type t_mp_vars_handle

  TYPE :: t_mp_vars_3dptr
    REAL(wp), POINTER, DIMENSION(:,:,:) :: dmean_c, dmean_r, deff_c, deff_i, reff_qc, reff_qi
    REAL(wp), POINTER, DIMENSION(:,:,:) :: dhmax, dhmax_ground, ze_p3
    REAL(wp), POINTER, DIMENSION(:,:,:) :: prec_gsp_rate, rain_gsp_rate, snow_gsp_rate, ice_gsp_rate
    REAL(wp), POINTER, DIMENSION(:,:,:) :: prec_gsp, prec_gsp_d, rain_gsp, snow_gsp, ice_gsp
    REAL(wp), POINTER, DIMENSION(:,:,:) :: q_sedim, twater
!! JM_20260323 >> adding new diagnostics
    ! --- warm-rain process rates ---
    REAL(wp), POINTER, DIMENSION(:,:,:) :: d_qcnuc, d_ncnuc  ! nucleation
    REAL(wp), POINTER, DIMENSION(:,:,:) :: d_qccon, d_qrcon  ! condensation
    REAL(wp), POINTER, DIMENSION(:,:,:) :: d_qcevp, d_qrevp, d_nrevp ! evaporation
    REAL(wp), POINTER, DIMENSION(:,:,:) :: d_qcacc, d_ncacc, d_qcaut, d_ncautc, d_ncautr ! collection (warm)
    REAL(wp), POINTER, DIMENSION(:,:,:) :: d_ncslf, d_nrslf ! self-collection
!! << JM_20260323
  CONTAINS
    PROCEDURE :: nullify => nullify_mp_vars_3dptr
  END type t_mp_vars_3dptr

  TYPE :: t_p3_vars_handle
    TYPE(t_comin_var_handle)            :: dmean_i, deff_i, rho_i, vm_i
  END type t_p3_vars_handle

  TYPE :: t_p3_vars_3dptr
    REAL(wp), POINTER, DIMENSION(:,:,:) :: dmean_i, deff_i, rho_i, vm_i
  CONTAINS
    PROCEDURE :: nullify => nullify_p3_vars_3dptr
  END type t_p3_vars_3dptr

  TYPE :: t_icon_tracer_handle
    TYPE(t_comin_var_handle)            :: qv, qv_old, qc, qi, qr, qs, qnc, qnr
  END type t_icon_tracer_handle

  TYPE :: t_icon_tracer_3dptr
    REAL(wp), POINTER, DIMENSION(:,:,:) :: qv, qv_old, qc, qi, qr, qs, qnc, qnr
  CONTAINS
    PROCEDURE :: nullify => nullify_icon_tracer_3dptr
  END type t_icon_tracer_3dptr

  TYPE :: t_p3_tracer_handle
    TYPE(t_comin_var_handle)            :: qitot, qnitot, qirim, birim, qzitot, qiliq
  END type t_p3_tracer_handle

  TYPE :: t_p3_tracer_3dptr
    REAL(wp), POINTER, DIMENSION(:,:,:) :: qitot, qnitot, qirim, birim, qzitot, qiliq
  CONTAINS
    PROCEDURE :: nullify => nullify_p3_tracer_3dptr
  END type t_p3_tracer_3dptr


CONTAINS

  SUBROUTINE nullify_dyn_vars_3dptr(this)
    CLASS(t_dyn_vars_3dptr), INTENT(inout) :: this
    NULLIFY(this%hfl, this%dz, this%temp, this%rho, this%pres)
    NULLIFY(this%exner, this%theta_old, this%ddt_temp_phys, this%w_hl)
  END SUBROUTINE nullify_dyn_vars_3dptr

  SUBROUTINE nullify_mp_vars_3dptr(this)
    CLASS(t_mp_vars_3dptr), INTENT(inout) :: this
    NULLIFY(this%dmean_c, this%dmean_r, this%deff_c, this%deff_i, this%reff_qc, this%reff_qi)
    NULLIFY(this%dhmax, this%dhmax_ground, this%ze_p3)
    NULLIFY(this%prec_gsp_rate, this%rain_gsp_rate, this%snow_gsp_rate, this%ice_gsp_rate)
    NULLIFY(this%prec_gsp, this%prec_gsp_d, this%rain_gsp, this%snow_gsp, this%ice_gsp)
    NULLIFY(this%q_sedim, this%twater)
!! JM_20260323 >> adding new diagnostics
    ! --- warm-rain process rates ---
    NULLIFY(this%d_qcnuc, this%d_ncnuc) ! nucleation
    NULLIFY(this%d_qccon, this%d_qrcon) ! condensation
    NULLIFY(this%d_qcevp, this%d_qrevp, this%d_nrevp) ! evaporation
    NULLIFY(this%d_qcacc, this%d_ncacc, this%d_qcaut, this%d_ncautc, this%d_ncautr) ! collection (warm)
    NULLIFY(this%d_ncslf, this%d_nrslf) ! self-collection
!! << JM_20260323
  END SUBROUTINE nullify_mp_vars_3dptr

  SUBROUTINE nullify_p3_vars_3dptr(this)
    CLASS(t_p3_vars_3dptr), INTENT(inout) :: this
    NULLIFY(this%dmean_i, this%deff_i, this%rho_i, this%vm_i)
  END SUBROUTINE nullify_p3_vars_3dptr

  SUBROUTINE nullify_icon_tracer_3dptr(this)
    CLASS(t_icon_tracer_3dptr), INTENT(inout) :: this
    NULLIFY(this%qv, this%qv_old, this%qc, this%qi, this%qr, this%qs, this%qnc, this%qnr)
  END SUBROUTINE nullify_icon_tracer_3dptr

  SUBROUTINE nullify_p3_tracer_3dptr(this)
    CLASS(t_p3_tracer_3dptr), INTENT(inout) :: this
    NULLIFY(this%qitot, this%qnitot, this%qirim, this%birim, this%qzitot, this%qiliq)
  END SUBROUTINE nullify_p3_tracer_3dptr

END MODULE p3plugin_types
