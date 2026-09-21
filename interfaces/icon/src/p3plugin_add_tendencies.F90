
MODULE p3plugin_add_tendencies
  USE comin_plugin_interface,  ONLY : comin_descrdata_get_cell_indices

  USE p3plugin_types,          ONLY : t_icon_tracer_3dptr, t_p3_tracer_3dptr
  USE p3plugin_global_vars,    ONLY : rank_world, n_icecat, l3mom_ice, lliqfrac, dtime, p_global, p_patch, &
    &                                 icon_tracer, p3_tracer,                                              &
    &                                 icon_tracer_ddt_turb, p3_tracer_ddt_turb

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: add_turb_tendencies

  INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND(12,307)

CONTAINS

  ! --------------------------------------------------------------------
  ! Add the turbulence tendencies to all tracers
  ! --------------------------------------------------------------------
  ! this subroutine is under development and currently not used.
  ! it can be included with the following callback line:
  !CALL comin_callback_register(EP_ATM_TURBULENCE_AFTER, add_turb_tendencies)
  SUBROUTINE add_turb_tendencies()  BIND(C)

    INTEGER :: jg, jb, jk, jc
    INTEGER :: i_icecat, i_startblk, i_endblk, i_startidx, i_endidx, rl_start, rl_end

    TYPE(t_icon_tracer_3dptr) :: icon_tracer_3d, icon_tracer_ddt_turb_3d
    TYPE(t_p3_tracer_3dptr)   :: p3_tracer_3d(n_icecat), p3_tracer_ddt_turb_3d(n_icecat)

    IF (rank_world == 0) WRITE (0,*) 'update turbulence tendencies of tracers'

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
      IF (l3mom_ice) THEN
        CALL p3_tracer_ddt_turb(i_icecat)%qzitot%to_3d(p3_tracer_ddt_turb_3d(i_icecat)%qzitot)
      ENDIF
      IF (lliqfrac) THEN
        CALL p3_tracer_ddt_turb(i_icecat)%qiliq%to_3d(p3_tracer_ddt_turb_3d(i_icecat)%qiliq)
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

  END SUBROUTINE add_turb_tendencies

END MODULE p3plugin_add_tendencies
