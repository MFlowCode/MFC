!>
!! @file m_boundary_conditions.fpp
!! @brief Contains module m_boundary_conditions

#:include 'macros.fpp'

!> @brief The purpose of the module is to apply noncharacteristic and processor
!! boundary condiitons
module m_boundary_conditions

    ! Dependencies =============================================================
    use m_derived_types        !< Definitions of the derived types

    use m_global_parameters    !< Definitions of the global parameters

    use m_mpi_proxy            !< Definitions of the MPI proxy

    use m_constants

    use m_boundary_conditions_common

    ! ==========================================================================

    implicit none

    private; 
    public :: s_populate_prim_buffers, &
              s_populate_capillary_buffers, &
              s_initialize_boundary_conditions_module

contains

    subroutine s_populate_capillary_buffers(c_divs)

        type(scalar_field), dimension(num_dims + 1), intent(inout) :: c_divs

        @:BOUNDARY_CONDITION_INTEGER_DECLARATIONS()

        #:for dir, loc in itertools.product([1, 2, 3], [-1, 1])
            #:block ITERATE_OVER_BUFFER_REGION_SIDED(dir=dir, loc=loc, outer_loops=[("i", 1, "num_dims + 1")])
                select case (bc_id_sfs(${dir}$, ${loc}$)%sf(exlhs, eylhs, ezlhs)%type)
                case (-13:-3); 
                    c_divs(i)%sf(x, y, z) = c_divs(i)%sf(ex, ey, ez)
                case (-2); !< slip wall or reflective
                    if (i == 1) then
                        c_divs(i)%sf(x, y, z) = -c_divs(i)%sf(sx, sy, sz)
                    else
                        c_divs(i)%sf(x, y, z) = +c_divs(i)%sf(sx, sy, sz)
                    end if
                case (-1); 
                    c_divs(i)%sf(x, y, z) = c_divs(i)%sf(px, py, pz)
                end select
            #:endblock

            if (${dir}$ <= num_dims) then
                call s_mpi_sendrecv_capilary_variables_buffers(c_divs, bc_id_sfs, ${dir}$, ${loc}$)
            end if
        #:endfor

    end subroutine s_populate_capillary_buffers

end module m_boundary_conditions
