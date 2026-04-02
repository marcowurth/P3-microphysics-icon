
MODULE p3plugin_types
  USE comin_plugin_interface,  ONLY : t_comin_var_handle

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: t_dyn_vars_handle, t_dyn_vars_3dptr
  PUBLIC :: t_mp_vars_handle, t_mp_vars_3dptr
  PUBLIC :: t_p3_vars_handle, t_p3_vars_3dptr
  PUBLIC :: t_icon_tracer_handle, t_icon_tracer_3dptr
  PUBLIC :: t_p3_tracer_handle, t_p3_tracer_3dptr
!! JM_:20260331 >> defining new comin handle and pointer for ice-phase diagnostics
  PUBLIC :: t_p3_ice_diag_handle, t_p3_ice_diag_3dptr
!! << JM_20260331
!! JM_20260402 >> defining new comin handle and pointer for ice-phase 3moment diagnostics
  PUBLIC :: t_p3_ice_3mom_diag_handle, t_p3_ice_3mom_diag_3dptr
!! << JM_20260402
!! JM_20260401 >> defining new comin handle and pointer for ice-ice collisions diagnostics
  PUBLIC :: t_p3_ice_coll_handle, t_p3_ice_coll_3dptr
!! << JM_20260401
!! JM_20260402 >> defining new comin handle and pointer for ice-liquid diagnostics
  PUBLIC :: t_p3_ice_liqfrac_handle, t_p3_ice_liqfrac_3dptr
!! << JM_20260402

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
    TYPE(t_comin_var_handle)            :: d_qcnuc, d_qccon, d_qrcon, d_qcevp, d_qrevp, d_qcacc, d_qcaut
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
    REAL(wp), POINTER, DIMENSION(:,:,:) :: d_qcnuc, d_qccon, d_qrcon, d_qcevp, d_qrevp, d_qcacc, d_qcaut
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

!! JM_20260331 >> defining new comin handle and pointer for ice-phase diagnostics (need to be adjusted here when adding more diagnostics)
  TYPE :: t_p3_ice_diag_handle
    TYPE(t_comin_var_handle)            :: d_qidep_, d_qisub_, d_qinuc_, d_qchetc_, d_qcheti_, d_qrhetc_, d_qrheti_, d_qrmlt_, d_qccol_, d_qrcol_, d_qwgrth_, d_qcshd_
    TYPE(t_comin_var_handle)            :: d_qcmul_, d_qrmul_, d_nimul_
  END type t_p3_ice_diag_handle

  TYPE :: t_p3_ice_diag_3dptr
    REAL(wp), POINTER, DIMENSION(:,:,:) :: d_qidep_, d_qisub_, d_qinuc_, d_qchetc_, d_qcheti_, d_qrhetc_, d_qrheti_, d_qrmlt_, d_qccol_, d_qrcol_, d_qwgrth_, d_qcshd_
    REAL(wp), POINTER, DIMENSION(:,:,:) :: d_qcmul_, d_qrmul_, d_nimul_
  CONTAINS
    PROCEDURE :: nullify => nullify_p3_ice_diag_3dptr
  END type t_p3_ice_diag_3dptr
!! JM_20260402 >> defining new comin handle and pointer for ice-phase 3moment diagnostics (need to be adjusted here when adding more diagnostics)
  TYPE :: t_p3_ice_3mom_diag_handle
    TYPE(t_comin_var_handle)            :: d_zidep_, d_zisub_, d_zimlt_, d_zislf_, d_zishd_, d_zqccol_, d_zqrcol_
  END type t_p3_ice_3mom_diag_handle

  TYPE :: t_p3_ice_3mom_diag_3dptr
    REAL(wp), POINTER, DIMENSION(:,:,:) :: d_zidep_, d_zisub_, d_zimlt_, d_zislf_, d_zishd_, d_zqccol_, d_zqrcol_
  CONTAINS
    PROCEDURE :: nullify => nullify_p3_ice_3mom_diag_3dptr
  END type t_p3_ice_3mom_diag_3dptr
!! << JM_20260402
!! JM_20260331 >> defining new comin handle and pointer for ice-phase diagnostics (need to be adjusted here when adding more diagnostics)
  TYPE :: t_p3_ice_coll_handle
    TYPE(t_comin_var_handle)            :: d_qicol_
  END type t_p3_ice_coll_handle

  TYPE :: t_p3_ice_coll_3dptr
    REAL(wp), POINTER, DIMENSION(:,:,:) :: d_qicol_
  CONTAINS
    PROCEDURE :: nullify => nullify_p3_ice_coll_3dptr
  END type t_p3_ice_coll_3dptr
!! << JM_20260401
!! JM_20260402 >> defining new comin handle and pointer for ice-liquid diagnostics (need to be adjusted here when adding more diagnostics)
  TYPE :: t_p3_ice_liqfrac_handle
    TYPE(t_comin_var_handle)            :: d_qimlt_, d_qwgrth1_, d_qwgrth1c_, d_qwgrth1r_, d_qlshd_, d_qlcon_, d_qlevp_, d_qifrz_, d_qccoll_, d_qrcoll_
  END type t_p3_ice_liqfrac_handle

  TYPE :: t_p3_ice_liqfrac_3dptr
    REAL(wp), POINTER, DIMENSION(:,:,:) :: d_qimlt_, d_qwgrth1_, d_qwgrth1c_, d_qwgrth1r_, d_qlshd_, d_qlcon_, d_qlevp_, d_qifrz_, d_qccoll_, d_qrcoll_
  CONTAINS
    PROCEDURE :: nullify => nullify_p3_ice_liqfrac_3dptr
  END type t_p3_ice_liqfrac_3dptr
!! << JM_20260402

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
!! JM_20260323 >> adding warm-rain diagnostics
    NULLIFY(this%d_qcnuc, this%d_qccon, this%d_qrcon, this%d_qcevp, this%d_qrevp, this%d_qcacc, this%d_qcaut)
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

!! JM_20260331 >> cleanup 4D ice-phase diagnostics (need to be adjusted here when adding more diagnostics)
  SUBROUTINE nullify_p3_ice_diag_3dptr(this)
    CLASS(t_p3_ice_diag_3dptr), INTENT(inout) :: this
    NULLIFY(this%d_qidep_, this%d_qisub_, this%d_qinuc_, this%d_qchetc_, this%d_qcheti_, this%d_qrhetc_, this%d_qrheti_, this%d_qrmlt_, this%d_qccol_, this%d_qrcol_, this%d_qwgrth_, this%d_qcshd_, this%d_qcmul_, this%d_qrmul_, this%d_nimul_)
  END SUBROUTINE nullify_p3_ice_diag_3dptr
!! << JM_20260331
!! JM_20260402 >> cleanup 4D ice-phase diagnostics (need to be adjusted here when adding more diagnostics)
  SUBROUTINE nullify_p3_ice_3mom_diag_3dptr(this)
    CLASS(t_p3_ice_3mom_diag_3dptr), INTENT(inout) :: this
    NULLIFY(this%d_zidep_, this%d_zisub_, this%d_zimlt_, this%d_zislf_, this%d_zishd_, this%d_zqccol_, this%d_zqrcol_)
  END SUBROUTINE nullify_p3_ice_3mom_diag_3dptr
!! << JM_20260402
!! JM_20260402 >> cleanup 4D ice-liquid diagnostics (need to be adjusted here when adding more diagnostics)
  SUBROUTINE nullify_p3_ice_liqfrac_3dptr(this)
    CLASS(t_p3_ice_liqfrac_3dptr), INTENT(inout) :: this
    NULLIFY(this%d_qimlt_, this%d_qwgrth1_, this%d_qwgrth1c_, this%d_qwgrth1r_, this%d_qlshd_, this%d_qlcon_, this%d_qlevp_, this%d_qifrz_, this%d_qccoll_, this%d_qrcoll_)
  END SUBROUTINE nullify_p3_ice_liqfrac_3dptr
!! << JM_20260402
!! JM_20260401 >> cleanup 5D ice-ice collisions diagnostics (need to be adjusted here when adding more diagnostics)
  SUBROUTINE nullify_p3_ice_coll_3dptr(this)
   CLASS(t_p3_ice_coll_3dptr), INTENT(inout) :: this
   NULLIFY(this%d_qicol_)
  END SUBROUTINE nullify_p3_ice_coll_3dptr
!! << JM_20260401
END MODULE p3plugin_types
