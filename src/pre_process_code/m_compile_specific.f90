!!       __  _______________
!!      /  |/  / ____/ ____/
!!     / /|_/ / /_  / /
!!    / /  / / __/ / /___
!!   /_/  /_/_/    \____/
!!
!!  This file is part of MFC.
!!
!!  MFC is the legal property of its developers, whose names
!!  are listed in the copyright file included with this source
!!  distribution.
!!
!!  MFC is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU General Public License as published
!!  by the Free Software Foundation, either version 3 of the license
!!  or any later version.
!!
!!  MFC is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!!  GNU General Public License for more details.
!!
!!  You should have received a copy of the GNU General Public License
!!  along with MFC (LICENSE).
!!  If not, see <http://www.gnu.org/licenses/>.

!>
!! @file m_compile_specific.f90
!! @brief Contains module m_compile_specific
!! @author S. Bryngelson, K. Schimdmayer, V. Coralic, J. Meng, K. Maeda, T. Colonius
!! @version 1.0
!! @date JUNE 06 2019

!> @brief This module contains subroutines that are compiler specific
module m_compile_specific

    implicit none

contains

    !>  Inquires on the existence of a directory
        !!  @param fileloc File directory location
        !!  @param dircheck Switch that indicates if directory exists
    subroutine my_inquire(fileloc, dircheck)
        character(LEN=*), intent(IN) :: fileloc
        logical, intent(INOUT) :: dircheck

! #ifdef __INTEL_COMPILER
        ! inquire (DIRECTORY=trim(fileloc), EXIST=dircheck)   !Intel
! #else
        inquire (FILE=trim(fileloc), EXIST=dircheck)        !GCC
! #endif

    end subroutine my_inquire

end module m_compile_specific
