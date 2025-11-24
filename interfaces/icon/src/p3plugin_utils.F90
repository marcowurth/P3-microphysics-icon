
MODULE p3plugin_utils
  USE mpi

  USE comin_plugin_interface,  ONLY : comin_parallel_get_host_mpi_rank, comin_parallel_get_plugin_mpi_comm, &
    &                                 comin_plugin_finish, t_comin_var_descriptor,                          &
    &                                 COMIN_ZAXIS_3D, COMIN_ZAXIS_2D,                                       &
    &                                 COMIN_VAR_DATATYPE_DOUBLE, COMIN_VAR_DATATYPE_FLOAT,                  &
    &                                 comin_var_request_add, comin_metadata_set

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: print_global_max, create_var, create_tracer

  INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND(12,307)

CONTAINS

  SUBROUTINE print_global_max(var_name, var_3dptr, factor)
    CHARACTER(*), INTENT(IN)      :: var_name
    REAL(wp), POINTER, INTENT(IN) :: var_3dptr(:,:,:)
    REAL, INTENT(IN), OPTIONAL    :: factor
    REAL                          :: local_max, global_max
    INTEGER                       :: comm, root, rank, ierr

    IF (PRESENT(factor)) THEN
      local_max = maxval(var_3dptr * factor)
    ELSE
      local_max = maxval(var_3dptr)
    END IF

    root = 0
    comm = comin_parallel_get_plugin_mpi_comm()
    rank = comin_parallel_get_host_mpi_rank()

    CALL MPI_REDUCE(local_max, global_max, 1, MPI_REAL, MPI_MAX, root, comm, ierr)
    IF (ierr /= 0) CALL comin_plugin_finish('print_global_max (p3plugin)', 'failed!')

    IF (rank == root) WRITE(0,*) 'global_max(' // var_name // ')', global_max
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

END MODULE p3plugin_utils
