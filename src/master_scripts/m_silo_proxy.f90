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
!! @file m_mpi_proxy.f90
!! @brief Contains module m_silo_proxy.f90 
!! @author S. Bryngelson, K. Schimdmayer, V. Coralic, J. Meng, K. Maeda, T. Colonius
!! @version 1.0
!! @date JUNE 06 2019

!> @brief The purpose of this module is to serve as a replacement framework
!!             for the Silo library and header file when those are not available
!!             on the platform on which the raw simulation data is to be post-
!!             processed. Note that this module merely allows for the code for
!!             the post-process to be compiled and executed and does not in
!!             actuality provide support for the creation of Silo-HDF5 database
!!             files. This means that when using this module during the post-
!!             process, the user must select another database output format.
!!  INSTRUCTIONS: To utilize this module, first place a copy of it in the
!!  post-process code directory. Next, modify the module m_data_output.f90 to
!!  use m_silo_proxy.f90 and erase the include command referencing the header
!!  file silo.inc. Finally, modify the makefile so that it includes in the
!!  compilation the Silo proxy module and remove any linker flags referencing
!!  the Silo library. No further changes should be necessary to compile and
!!  execute the post-process code.
MODULE m_silo_proxy
    
    
    IMPLICIT NONE
    
    
    !> @name Refer to Silo's user guide (10/2007, v4.6) for the variables' definitions
    !! and the header file silo.inc for their choice of values
    !> @{
        INTEGER, PARAMETER :: DB_CLOBBER         = 0
    INTEGER, PARAMETER :: DB_COLLINEAR       = 130
    INTEGER, PARAMETER :: DB_FLOAT           = 19
    INTEGER, PARAMETER :: DB_DOUBLE          = 20
    INTEGER, PARAMETER :: DB_F77NULL         = -99
    INTEGER, PARAMETER :: DB_HDF5            = 7
    INTEGER, PARAMETER :: DB_LOCAL           = 0
    INTEGER, PARAMETER :: DB_QUAD_RECT       = 130
    INTEGER, PARAMETER :: DB_QUADVAR         = 501
    INTEGER, PARAMETER :: DB_ZONECENT        = 111
    INTEGER, PARAMETER :: DBOPT_EXTENTS      = 300
    INTEGER, PARAMETER :: DBOPT_EXTENTS_SIZE = 299
    INTEGER, PARAMETER :: DBOPT_HI_OFFSET    = 265
    INTEGER, PARAMETER :: DBOPT_LO_OFFSET    = 266
    !> @}
    
    CONTAINS
        
        
        !> @brief  Refer to page 28 of Silo's user guide (10/2007, v4.6) for
        !!             information about this subroutine
        FUNCTION DBCREATE( pathname, lpathname,   mode  , target, & !! ----------
                           fileinfo, lfileinfo, filetype, status  )
            
            
            INTEGER                        :: DBCREATE
            CHARACTER(LEN = *), INTENT(IN) :: pathname
            INTEGER           , INTENT(IN) :: lpathname
            INTEGER           , INTENT(IN) :: mode
            INTEGER           , INTENT(IN) :: target
            CHARACTER(LEN = *), INTENT(IN) :: fileinfo
            INTEGER           , INTENT(IN) :: lfileinfo
            INTEGER           , INTENT(IN) :: filetype
            INTEGER           , INTENT(IN) :: status
            
            
            PRINT '(A)', 'Creation of Silo-HDF5 database file is ' // &
                         'inconsistent with use of Silo proxy module. ' // &
                         'Exiting ...'
            STOP
            
            
        END FUNCTION DBCREATE !! ------------------------------------------------
        
        
        
        
        !> @brief  Refer to page 235 of Silo's user guide (10/2007, v4.6)
        !!             for information about this subroutine        
        FUNCTION DBGET2DSTRLEN() !! ---------------------------------------------

            
            
            INTEGER :: DBGET2DSTRLEN
            
            
            PRINT '(A)', 'Attempt to query 2D string length is ' // &
                         'inconsistent with use of Silo proxy module. ' // &
                         'Exiting ...'
            STOP
            
            
        END FUNCTION DBGET2DSTRLEN !! -------------------------------------------
        
        
        
        
        !> @brief  Refer to page 234 of Silo's user guide (10/2007, v4.6)
        !!             for information about this subroutine        
        FUNCTION DBSET2DSTRLEN(len) !! ------------------------------------------

            
            
            INTEGER             :: DBSET2DSTRLEN
            INTEGER, INTENT(IN) :: len
            
            
            PRINT '(A)', 'Attempt to set 2D string length is ' // &
                         'inconsistent with use of Silo proxy module. ' // &
                         'Exiting ...'
            STOP
            
            
        END FUNCTION DBSET2DSTRLEN !! -------------------------------------------
        
        
        
        
        !> @brief  Refer to page 185 of Silo's user guide (10/2007, v4.6)
        !!             for information about this subroutine        
        FUNCTION DBMKOPTLIST(maxopts, optlist_id) !! ----------------------------

            
            
            INTEGER             :: DBMKOPTLIST
            INTEGER, INTENT(IN) :: maxopts
            INTEGER, INTENT(IN) :: optlist_id
            
            
            PRINT '(A)', 'Allocation of an options list is inconsistent ' // &
                         'with use of Silo proxy module. Exiting ...'
            STOP
            
            
        END FUNCTION DBMKOPTLIST !! ---------------------------------------------
        
        
        
        
        !> @brief  Refer to page 186 of Silo's user guide (10/2007, v4.6)
        !!             for information about this subroutine        
        FUNCTION DBADDIOPT(optlist_id, option, ivalue) !! -----------------------

            
            
            INTEGER             :: DBADDIOPT
            INTEGER, INTENT(IN) :: optlist_id
            INTEGER, INTENT(IN) :: option
            INTEGER, INTENT(IN) :: ivalue
            
            
            PRINT '(A)', 'Adding an option to an options list is ' // &
                         'inconsistent with use of Silo proxy module. ' // &
                         'Exiting ...'
            STOP
            
            
        END FUNCTION DBADDIOPT !! -----------------------------------------------
        
        
        
        
        !> @brief  Refer to page 186 of Silo's user guide (10/2007, v4.6)
        !!             for information about this subroutine        
        FUNCTION DBADDDOPT(optlist_id, option, dvalue) !! -----------------------

            
            
            INTEGER                             :: DBADDDOPT
            INTEGER                , INTENT(IN) :: optlist_id
            INTEGER                , INTENT(IN) :: option
            INTEGER, DIMENSION(:,:), INTENT(IN) :: dvalue
            
            
            PRINT '(A)', 'Adding an option to an options list is ' // &
                         'inconsistent with use of Silo proxy module. ' // &
                         'Exiting ...'
            STOP
            
            
        END FUNCTION DBADDDOPT !! -----------------------------------------------
        
        
        
        
        !> @brief  Refer to page 121 of Silo's user guide (10/2007, v4.6)
        !!             for information about this subroutine        
        FUNCTION DBPUTMMESH(    dbid, name, lname, nmesh, meshnames,   & !! -----
                             lmeshnames, meshtypes, optlist_id, status )

            
            
            INTEGER                                      :: DBPUTMMESH
            INTEGER                         , INTENT(IN) :: dbid
            CHARACTER(LEN = *)              , INTENT(IN) :: name
            INTEGER                         , INTENT(IN) :: lname
            INTEGER                         , INTENT(IN) :: nmesh
            CHARACTER(LEN = *), DIMENSION(:), INTENT(IN) :: meshnames
            INTEGER           , DIMENSION(:), INTENT(IN) :: lmeshnames
            INTEGER           , DIMENSION(:), INTENT(IN) :: meshtypes
            INTEGER                         , INTENT(IN) :: optlist_id
            INTEGER                         , INTENT(IN) :: status
            
            
            PRINT '(A)', 'Writing of multi-block mesh object into ' // &
                         'Silo-HDF5 database file is inconsistent with ' // &
                         'use of Silo proxy module. Exiting ...'
            STOP
            
            
        END FUNCTION DBPUTMMESH !! ----------------------------------------------
        
        
        
        
        !> @brief  Refer to page 189 of Silo's user guide (10/2007, v4.6)
        !!             for information about this subroutine    
        FUNCTION DBFREEOPTLIST(optlist_id) !! -----------------------------------

            
            
            INTEGER             :: DBFREEOPTLIST
            INTEGER, INTENT(IN) :: optlist_id
            
            
            PRINT '(A)', 'Deallocation of an options list is inconsistent ' // &
                         'with use of Silo proxy module. Exiting ...'
            STOP
            
            
        END FUNCTION DBFREEOPTLIST !! -------------------------------------------
        
        
        
        
        !> @brief  Refer to page 57 of Silo's user guide (10/2007, v4.6) for
        !!             information about this subroutine
        FUNCTION DBPUTQM( dbid, name, lname, xname, lxname, yname, lyname, & !! -
                           zname, lzname, x, y, z, dims, ndims, datatype,  &
                                    coordtype, optlist_id, status          )
            
            
            INTEGER                                      :: DBPUTQM
            INTEGER           ,               INTENT(IN) :: dbid
            CHARACTER(LEN = *),               INTENT(IN) :: name
            INTEGER           ,               INTENT(IN) :: lname
            CHARACTER(LEN = *),               INTENT(IN) :: xname
            INTEGER           ,               INTENT(IN) :: lxname
            CHARACTER(LEN = *),               INTENT(IN) :: yname
            INTEGER           ,               INTENT(IN) :: lyname
            CHARACTER(LEN = *),               INTENT(IN) :: zname
            INTEGER                         , INTENT(IN) :: lzname
            REAL(KIND(0d0))   , DIMENSION(:), INTENT(IN) :: x
            REAL(KIND(0d0))   , DIMENSION(:), INTENT(IN) :: y
            REAL(KIND(0d0))   , DIMENSION(:), INTENT(IN) :: z
            INTEGER           , DIMENSION(:), INTENT(IN) :: dims
            INTEGER           ,               INTENT(IN) :: ndims
            INTEGER           ,               INTENT(IN) :: datatype
            INTEGER           ,               INTENT(IN) :: coordtype
            INTEGER           ,               INTENT(IN) :: optlist_id
            INTEGER           ,               INTENT(IN) :: status
            
            
            PRINT '(A)', 'Writing of quad mesh object into Silo-HDF5 ' // &
                         'database file is inconsistent with use of Silo ' // &
                         'proxy module. Exiting ...'
            STOP
            
            
        END FUNCTION DBPUTQM !! -------------------------------------------------
        
        
        
        
        !> @brief  Refer to page 46 of Silo's user guide (10/2007, v4.6) for
        !!             information about this subroutine
        FUNCTION DBPUTCURVE( dbid, curvename, lcurvename, xvals, yvals, & !! ----
                               datatype, npoints, optlist_id, status    )

            
            
            INTEGER                                      :: DBPUTCURVE
            INTEGER           ,               INTENT(IN) :: dbid
            CHARACTER(LEN = *),               INTENT(IN) :: curvename
            INTEGER           ,               INTENT(IN) :: lcurvename
            REAL(KIND(0d0))   , DIMENSION(:), INTENT(IN) :: xvals
            REAL(KIND(0d0))   , DIMENSION(:), INTENT(IN) :: yvals
            INTEGER           ,               INTENT(IN) :: datatype
            INTEGER           ,               INTENT(IN) :: npoints
            INTEGER           ,               INTENT(IN) :: optlist_id
            INTEGER           ,               INTENT(IN) :: status
            
            
            PRINT '(A)', 'Writing of curve object into Silo-HDF5 database ' // &
                         'file is inconsistent with use of Silo proxy ' // &
                         'module. Exiting ...'
            STOP
            
            
        END FUNCTION DBPUTCURVE !! ----------------------------------------------
        
        
        
        
        !> @brief  Refer to page 130 of Silo's user guide (10/2007, v4.6)
        !!             for information about this subroutine
        FUNCTION DBPUTMVAR( dbid, name, lname, nvar, varnames, lvarnames, & !! --
                                    vartypes, optlist_id, status          )

            
            INTEGER                                      :: DBPUTMVAR
            INTEGER           ,               INTENT(IN) :: dbid
            CHARACTER(LEN = *),               INTENT(IN) :: name
            INTEGER           ,               INTENT(IN) :: lname
            INTEGER           ,               INTENT(IN) :: nvar
            CHARACTER(LEN = *), DIMENSION(:), INTENT(IN) :: varnames
            INTEGER           , DIMENSION(:), INTENT(IN) :: lvarnames
            INTEGER           , DIMENSION(:), INTENT(IN) :: vartypes
            INTEGER           ,               INTENT(IN) :: optlist_id
            INTEGER           ,               INTENT(IN) :: status
            
            
            PRINT '(A)', 'Writing of multi-block variable into Silo-HDF5 ' // &
                         'database file is inconsistent with use of Silo ' // &
                         'proxy module. Exiting ...'
            STOP
            
            
        END FUNCTION DBPUTMVAR !! -----------------------------------------------
        
        
        
        
        !> @brief  Refer to page 64 of Silo's user guide (10/2007, v4.6) for
        !!             information about this subroutine    
        FUNCTION DBPUTQV1( dbid, name, lname, meshname, lmeshname, var, & !! ----
                              dims, ndims, mixvar, mixlen, datatype,    &
                                  centering, optlist_id, status         )

            
            
            INTEGER                                          :: DBPUTQV1
            INTEGER           ,                   INTENT(IN) :: dbid
            CHARACTER(LEN = *),                   INTENT(IN) :: name
            INTEGER           ,                   INTENT(IN) :: lname
            CHARACTER(LEN = *),                   INTENT(IN) :: meshname
            INTEGER           ,                   INTENT(IN) :: lmeshname
            REAL(KIND(0d0))   , DIMENSION(:,:,:), INTENT(IN) :: var
            INTEGER           , DIMENSION(:)    , INTENT(IN) :: dims
            INTEGER           ,                   INTENT(IN) :: ndims
            INTEGER           ,                   INTENT(IN) :: mixvar
            INTEGER           ,                   INTENT(IN) :: mixlen
            INTEGER           ,                   INTENT(IN) :: datatype
            INTEGER           ,                   INTENT(IN) :: centering
            INTEGER           ,                   INTENT(IN) :: optlist_id
            INTEGER           ,                   INTENT(IN) :: status
            
            
            PRINT '(A)', 'Writing a scalar quad variable into a Silo-HDF5 ' // &
                         'database file is inconsistent with use of Silo ' // &
                         'proxy module. Exiting ...'
            STOP
            
            
        END FUNCTION DBPUTQV1 !! ------------------------------------------------
        
        
        
        
        !> @brief  Refer to page 31 of Silo's user guide (10/2007, v4.6) for
        !!             information about this subroutine
        FUNCTION DBCLOSE(dbid) !! -----------------------------------------------

            
            
            INTEGER             :: DBCLOSE
            INTEGER, INTENT(IN) :: dbid
            
            
            PRINT '(A)', 'Attempt to close Silo-HDF5 database file is ' // &
                         'inconsistent with use of Silo proxy module. ' // &
                         'Exiting ...'
            STOP
            
            
        END FUNCTION DBCLOSE !! -------------------------------------------------
        
        
        
        
        
END MODULE m_silo_proxy
