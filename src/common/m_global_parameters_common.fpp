module m_global_parameters_common

    use m_constants
    use m_derived_types

    implicit none

contains

    subroutine s_bc_assign_default_values_to_user_inputs(num_bc_patches, patch_bc)

        integer, intent(inout) :: num_bc_patches
        type(bc_patch_parameters), intent(inout) :: patch_bc(num_bc_patches_max)

        integer :: i

        num_bc_patches = 0
        do i = 1, num_bc_patches_max
            patch_bc(i)%type = dflt_int
            patch_bc(i)%geometry = dflt_int
            patch_bc(i)%dir = dflt_int
            patch_bc(i)%loc = dflt_int
        end do

    end subroutine s_bc_assign_default_values_to_user_inputs

end module m_global_parameters_common
