
MODULE p3plugin_utils
  USE mpi,                     ONLY : MPI_Reduce, MPI_REAL, MPI_MAX

  USE comin_plugin_interface,  ONLY : comin_plugin_finish, t_comin_var_descriptor,                          &
    &                                 COMIN_ZAXIS_3D, COMIN_ZAXIS_2D,                                       &
    &                                 COMIN_VAR_DATATYPE_DOUBLE, COMIN_VAR_DATATYPE_FLOAT,                  &
    &                                 comin_var_request_add, comin_metadata_set
  USE p3plugin_global_vars,    ONLY : rank_world, comm_world

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: print_global_max, create_var, create_tracer
  PUBLIC :: uppercase, lowercase

  INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND(12,307)

CONTAINS

  SUBROUTINE print_global_max(var_name, var_3dptr, factor)
    CHARACTER(*), INTENT(IN)      :: var_name
    REAL(wp), POINTER, INTENT(IN) :: var_3dptr(:,:,:)
    REAL, INTENT(IN), OPTIONAL    :: factor
    REAL                          :: local_max, global_max
    INTEGER                       :: root, ierr

    IF (PRESENT(factor)) THEN
      local_max = maxval(var_3dptr * factor)
    ELSE
      local_max = maxval(var_3dptr)
    END IF

    root = 0
    CALL MPI_Reduce(local_max, global_max, 1, MPI_REAL, MPI_MAX, root, comm_world, ierr)
    IF (ierr /= 0) CALL comin_plugin_finish('print_global_max (p3plugin)', 'failed!')

    IF (rank_world == root) WRITE(0,*) 'global_max(' // var_name // ')', global_max
  END SUBROUTINE print_global_max


  SUBROUTINE create_var(var_name, unit_name, axis_type, datatype_precision, long_name)
    CHARACTER(*), INTENT(IN) :: var_name, unit_name, axis_type, datatype_precision
    CHARACTER(*), INTENT(IN), OPTIONAL :: long_name
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
      &                                ltracer=.FALSE., lrestart=.FALSE., datatype=datatype, &
      &                                long_name=long_name)
  END SUBROUTINE create_var


  SUBROUTINE create_tracer(var_name, unit_name, ltracer_turb, long_name)
    CHARACTER(*), INTENT(IN) :: var_name, unit_name
    CHARACTER(*), INTENT(IN), OPTIONAL :: long_name
    LOGICAL, INTENT(IN)      :: ltracer_turb

    CALL comin_var_request_add_wrapper(descriptor=t_comin_var_descriptor(name=TRIM(var_name), id=-1), &
      &                                units=TRIM(unit_name), lmode_exclusive=.FALSE., zaxis_id=COMIN_ZAXIS_3D, &
                                       ltracer=.TRUE., lrestart=.FALSE., ltracer_turb=ltracer_turb, &
                                       long_name=long_name)
  END SUBROUTINE create_tracer


  SUBROUTINE comin_var_request_add_wrapper(descriptor, units, lmode_exclusive, zaxis_id, &
                                           ltracer, lrestart, datatype, ltracer_turb, &
                                           long_name)
    TYPE(t_comin_var_descriptor), INTENT(IN)  :: descriptor
    LOGICAL,            OPTIONAL, INTENT(IN)  :: lmode_exclusive
    LOGICAL,            OPTIONAL, INTENT(IN)  :: ltracer_turb
    INTEGER,            OPTIONAL, INTENT(IN)  :: zaxis_id
    INTEGER,            OPTIONAL, INTENT(IN)  :: datatype
    LOGICAL,            OPTIONAL, INTENT(IN)  :: ltracer
    LOGICAL,            OPTIONAL, INTENT(IN)  :: lrestart
    CHARACTER(LEN=*),   OPTIONAL, INTENT(IN)  :: units
    CHARACTER(LEN=*),   OPTIONAL, INTENT(IN)  :: long_name
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
!! JM_20260415 >> adding long name
    IF (PRESENT(long_name)) THEN
      CALL comin_metadata_set(descriptor, "long_name", TRIM(long_name))
    END IF
!! << JM_20260415

  END SUBROUTINE comin_var_request_add_wrapper

  FUNCTION uppercase(word_in) RESULT(word_out)
    CHARACTER(len=*), INTENT(IN) :: word_in
    CHARACTER(:), ALLOCATABLE    :: word_out
    CHARACTER(len=1)   :: c
    INTEGER            :: i

    ALLOCATE(CHARACTER(len_trim(word_in)) :: word_out)
    DO i = 1, len_trim(word_in)
      c = word_in(i:i)
      IF (c >= 'a' .and. c <= 'z') THEN
        word_out(i:i) = char(ichar(c) - 32)
      ELSE
        word_out(i:i) = c
      END IF
    END DO
  END FUNCTION uppercase

  FUNCTION lowercase(word_in) RESULT(word_out)
    CHARACTER(len=*), INTENT(IN) :: word_in
    CHARACTER(:), ALLOCATABLE    :: word_out
    CHARACTER(len=1)   :: c
    INTEGER            :: i

    ALLOCATE(CHARACTER(len_trim(word_in)) :: word_out)
    DO i = 1, len_trim(word_in)
      c = word_in(i:i)
      IF (c >= 'A' .and. c <= 'Z') THEN
        word_out(i:i) = char(ichar(c) + 32)
      ELSE
        word_out(i:i) = c
      END IF
    END DO
  END FUNCTION lowercase

END MODULE p3plugin_utils
