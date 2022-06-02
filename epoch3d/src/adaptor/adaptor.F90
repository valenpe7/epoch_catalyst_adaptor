MODULE coprocessor
  
  USE, INTRINSIC :: iso_c_binding
  USE fields
  USE calc_df
 
  IMPLICIT NONE

INTERFACE

  SUBROUTINE add_input_description(grid_name) bind(c)

    USE, INTRINSIC :: iso_c_binding
    IMPLICIT NONE

    CHARACTER(kind=c_char), DIMENSION(*), INTENT(IN) :: grid_name

  END SUBROUTINE add_input_description

  SUBROUTINE build_uniform_grid(rank, size, extent, zero_extent, whole_extent, &
    origin, spacing, grid_name) bind(c)

    USE, INTRINSIC :: iso_c_binding
    IMPLICIT NONE

    INTEGER(c_int), INTENT(IN) :: rank, size
    INTEGER(c_int), DIMENSION(*), INTENT(IN) :: extent, zero_extent, whole_extent
    REAL(c_double), DIMENSION(*), INTENT(IN) :: origin, spacing
    CHARACTER(kind=c_char), DIMENSION(*), INTENT(IN) :: grid_name

  END SUBROUTINE build_uniform_grid

  SUBROUTINE build_rectilinear_grid(rank, size, extent, &
    x_coords, y_coords, z_coords, ghost_levels, grid_name) bind(c)

    USE, INTRINSIC :: iso_c_binding
    IMPLICIT NONE

    INTEGER(c_int), INTENT(IN) :: rank, size, ghost_levels
    INTEGER(c_int), DIMENSION(*), INTENT(IN) :: extent
    REAL(c_double), DIMENSION(*), INTENT(IN) :: x_coords, y_coords, z_coords
    CHARACTER(kind=c_char), DIMENSION(*), INTENT(IN) :: grid_name

  END SUBROUTINE build_rectilinear_grid

  SUBROUTINE build_unstructured_grid(count, grid_name) bind(c)

    USE, INTRINSIC :: iso_c_binding
    IMPLICIT NONE

    INTEGER(c_long_long), INTENT(IN) :: count
    CHARACTER(kind=c_char), DIMENSION(*), INTENT(IN) :: grid_name

  END SUBROUTINE build_unstructured_grid

  SUBROUTINE set_time_data(step, time) bind(c)

    USE, INTRINSIC :: iso_c_binding
    IMPLICIT NONE

    INTEGER(c_int), INTENT(IN) :: step
    REAL(c_double), INTENT(IN) :: time

  END SUBROUTINE set_time_data

  SUBROUTINE set_scalar_field(rank, scalar, name, grid_name) bind(c)

    USE, INTRINSIC :: iso_c_binding
    IMPLICIT NONE

    INTEGER(c_int), INTENT(IN) :: rank
    REAL(c_double), DIMENSION(*), INTENT(IN) :: scalar
    CHARACTER(kind=c_char), DIMENSION(*), INTENT(IN) :: name, grid_name

  END SUBROUTINE set_scalar_field

  SUBROUTINE set_vector_field(rank, x_comp, y_comp, z_comp, name, grid_name) bind(c)
  
    USE, INTRINSIC :: iso_c_binding
    IMPLICIT NONE

    INTEGER(c_int), INTENT(IN) :: rank
    REAL(c_double), DIMENSION(*), INTENT(IN) :: x_comp, y_comp, z_comp
    CHARACTER(kind=c_char), DIMENSION(*), INTENT(IN) :: name, grid_name

  END SUBROUTINE set_vector_field

  SUBROUTINE set_particle_position(id, x, y, z, grid_name) bind(c)

    USE, INTRINSIC :: iso_c_binding
    IMPLICIT NONE

    INTEGER(c_long_long), INTENT(IN) :: id
    REAL(c_double), INTENT(IN) :: x, y, z
    CHARACTER(kind=c_char), DIMENSION(*), INTENT(IN) :: grid_name

  END SUBROUTINE set_particle_position

  SUBROUTINE set_scalar_particle(id, scalar, name, grid_name) bind(c)

    USE, INTRINSIC :: iso_c_binding
    IMPLICIT NONE

    REAL(c_double), INTENT(IN) :: scalar
    INTEGER(c_long_long), INTENT(IN) :: id
    CHARACTER(kind=c_char), DIMENSION(*), INTENT(IN) :: name, grid_name

  END SUBROUTINE set_scalar_particle

  SUBROUTINE set_vector_particle(id, x_comp, y_comp, z_comp, name, grid_name) bind(c)

    USE, INTRINSIC :: iso_c_binding
    IMPLICIT NONE

    REAL(c_double), INTENT(IN) :: x_comp, y_comp, z_comp
    INTEGER(c_long_long), INTENT(IN) :: id
    CHARACTER(kind=c_char), DIMENSION(*), INTENT(IN) :: name, grid_name

  END SUBROUTINE set_vector_particle

END INTERFACE

CONTAINS

  SUBROUTINE init_coproc(step, time)  
    INTEGER, INTENT(in) :: step
    REAL(num), INTENT(in) :: time
    INTEGER :: ilen, i
    CHARACTER(len=200) :: arg
    CALL coprocessorinitialize()
    DO i = 1, iargc()
       CALL getarg(i, arg)
       ilen = len_trim(arg)
       arg(ilen+1:) = C_NULL_CHAR
       CALL coprocessoraddpythonscript(arg, ilen)
    ENDDO
    CALL set_time_data(step, time)
    CALL add_input_description("Grid"//C_NULL_CHAR)
   ! DO i = 1, n_species
   !   CALL add_input_description("Particles - "//TRIM(species_list(i)%name)//C_NULL_CHAR)
   ! ENDDO
    
    ALLOCATE(dens(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng, n_species))
  
  END SUBROUTINE init_coproc

  SUBROUTINE run_coproc(step, time)
    INTEGER, INTENT(in) :: step
    REAL(num), INTENT(in) :: time
    INTEGER :: flag, ispecies
    INTEGER(i8) :: ipart
    INTEGER, DIMENSION(6) :: extent, zero_extent, whole_extent
    REAL(num), DIMENSION(3) :: origin, spacings
    !TYPE(particle), POINTER :: current

#ifdef NO_IO
    RETURN
#endif

    extent = (/ nx_global_min - ng, nx_global_max + ng, ny_global_min - ng, ny_global_max + ng, & 
                nz_global_min - ng, nz_global_max + ng /)
    zero_extent = (/ nx_global_min, nx_global_max, ny_global_min, ny_global_max, nz_global_min, nz_global_max /)
    whole_extent = (/ 1 - ng, nx_global + ng, 1 - ng, ny_global + ng, 1 - ng, nz_global + ng /)
    origin = (/ x_grid_min, y_grid_min, z_grid_min /)
    spacings = (/ dx, dy, dz /)

    CALL set_time_data(step, time)

    CALL requestdatadescription(step, time, flag)
    IF (flag /= 0) THEN
      
      CALL build_uniform_grid(rank, nproc, extent, zero_extent, whole_extent, origin, spacings, "Grid"//C_NULL_CHAR)
     ! CALL build_rectilinear_grid(rank, nproc, extent, x(0:nx), y(0:ny), z(0:nz), ng, "EM Field"//C_NULL_CHAR)

     ! DO ispecies = 1, n_species
     !   CALL build_unstructured_grid(species_list(ispecies)%attached_list%count, &
     !     "Particles - "//TRIM(species_list(ispecies)%name)//C_NULL_CHAR)
     ! ENDDO

     ! CALL set_vector_field(rank, ex, ey, ez, "Electric Field (V/m)"//C_NULL_CHAR, "EM Field"//C_NULL_CHAR)
     ! CALL set_vector_field(rank, bx, by, bz, "Magnetic Field (T)"//C_NULL_CHAR, "EM Field"//C_NULL_CHAR)
     ! CALL set_vector_field(rank, jx, jy, jz, "Current Density (A/m^2)"//C_NULL_CHAR, "EM Field"//C_NULL_CHAR)
      CALL set_scalar_field(rank, ex, "Electric Field - E_x (V/m)"//C_NULL_CHAR, "Grid"//C_NULL_CHAR)
      CALL set_scalar_field(rank, ez, "Electric Field - E_z (V/m)"//C_NULL_CHAR, "Grid"//C_NULL_CHAR)
      
      DO ispecies = 1, n_species
        CALL calc_number_density(dens(:, :, :, ispecies), ispecies)
        CALL do_field_mpi_with_lengths(dens(:, :, :, ispecies), ng, nx, ny, nz)
        CALL set_scalar_field(rank, dens(:, :, :, ispecies), &
          "Number Density - "//TRIM(species_list(ispecies)%name)//" (m^-3)"//C_NULL_CHAR, "Grid"//C_NULL_CHAR)
      ENDDO

     ! DO ispecies = 1, n_species
     !   current => species_list(ispecies)%attached_list%head
     !   DO ipart = 0, species_list(ispecies)%attached_list%count - 1
          
     !     CALL set_particle_position(ipart, current%part_pos(1), current%part_pos(2), &
     !       current%part_pos(3), "Particles - "//TRIM(species_list(ispecies)%name)//C_NULL_CHAR)

     !     CALL set_vector_particle(ipart, current%part_p(1), current%part_p(2), current%part_p(3), & 
     !       "Momentum (kg*m/s)"//C_NULL_CHAR, "Particles - "//TRIM(species_list(ispecies)%name)//C_NULL_CHAR)

     !     CALL set_scalar_particle(ipart, current%weight, "Weight (-)"//C_NULL_CHAR, &
     !       "Particles - "//TRIM(species_list(ispecies)%name)//C_NULL_CHAR)
          
     !     current => current%next
     !   ENDDO
     ! ENDDO
      
      CALL coprocess()

    ENDIF
  END SUBROUTINE run_coproc

  SUBROUTINE finalise_coproc()
    DEALLOCATE(dens)
    CALL coprocessorfinalize()
  END SUBROUTINE finalise_coproc

END MODULE coprocessor
