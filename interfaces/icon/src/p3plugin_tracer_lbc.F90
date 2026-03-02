
MODULE p3plugin_tracer_lbc
  USE comin_plugin_interface,  ONLY : comin_parallel_get_host_mpi_rank, comin_descrdata_get_cell_indices

  USE p3plugin_types,          ONLY : t_icon_tracer_3dptr, t_p3_tracer_3dptr
  USE p3plugin_global_vars,    ONLY : n_icecat, p_global, p_patch, icon_tracer, p3_tracer

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: update_lbc_ice_after_nudging

  INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND(12,307)

CONTAINS

  ! -----------------------------------------------------------------------------
  ! Update P3 ice to the changes that come from the changes in the nudging zone
  ! -----------------------------------------------------------------------------
  ! this subroutine is under development and currently not used.
  ! it can be included with the following callback line:
  !CALL comin_callback_register(EP_ATM_NUDGING_AFTER, update_lbc_ice_after_nudging)

  SUBROUTINE update_lbc_ice_after_nudging()  BIND(C)

    INTEGER  :: jg, jb, jk, jc
    INTEGER  :: rank, i_icecat, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end
    REAL(wp) :: qi_sum
    REAL(wp) :: qsmall = 1.0e-12    ! minimum threshold in kg/kg

    TYPE(t_icon_tracer_3dptr) :: icon_tracer_3d
    TYPE(t_p3_tracer_3dptr)   :: p3_tracer_3d(n_icecat)

    rank = comin_parallel_get_host_mpi_rank()
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

  END SUBROUTINE update_lbc_ice_after_nudging

END MODULE p3plugin_tracer_lbc
