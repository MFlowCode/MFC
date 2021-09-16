!!       __  _______________
!!      /  |/  / ____/ ____/
!!     / /|_/ / /_  / /     
!!    / /  / / __/ / /___   
!!   /_/  /_/_/    \____/   
!!                       
!!  This file is part of MFC.
!!
!! Copyright 2021
!!
!! Permission is hereby granted, free of charge, to any person 
!! obtaining a copy of this software and associated documentation 
!! files (the "Software"), to deal in the Software without 
!! restriction, including without limitation the rights to use, 
!! copy, modify, merge, publish, distribute, sublicense, 
!! and/or sell copies of the Software, and to permit persons 
!! to whom the Software is furnished to do so, subject to the 
!! following conditions:
!!
!! The above copyright notice and this permission notice shall 
!! be included in all copies or substantial portions of the Software.
!!
!! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
!! EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
!! MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
!! IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
!! FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION 
!! OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
!! CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN 
!! THE SOFTWARE.

!>
!! @file m_compile_specific.f90
!! @brief Contains module m_compile_specific
!! @author S. Bryngelson, K. Schimdmayer, V. Coralic, J. Meng, K. Maeda, T. Colonius
!! @version 1.0
!! @date JUNE 06 2019

!> @brief This module contains subroutines that are compiler specific
MODULE m_compile_specific

    IMPLICIT NONE

    CONTAINS

        !>  Inquires on the existence of a directory
        !!  @param fileloc File directory location
        !!  @param dircheck Switch that indicates if directory exists
        SUBROUTINE my_inquire(fileloc,dircheck)
            CHARACTER(LEN=*), INTENT(IN) :: fileloc
            LOGICAL, INTENT(INOUT) :: dircheck

#ifdef __INTEL_COMPILER
    INQUIRE(DIRECTORY=TRIM(fileloc),EXIST=dircheck)   !Intel
#else
    INQUIRE(FILE=TRIM(fileloc),EXIST=dircheck)        !GCC
#endif

        END SUBROUTINE my_inquire

END MODULE m_compile_specific
