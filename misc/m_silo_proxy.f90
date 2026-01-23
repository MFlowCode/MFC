!>
!! @file m_mpi_proxy.f90
!! @brief Contains module m_silo_proxy.f90

!> @brief The purpose of this module is to serve as a replacement framework
!!                           for the Silo library and header file when those are not available
!!                           on the platform on which the raw simulation data is to be post-
!!                           processed. Note that this module merely allows for the code for
!!                           the post-process to be compiled and executed and does not in
!!                           actuality provide support for the creation of Silo-HDF5 database
!!                           files. This means that when using this module during the post-
!!                           process, the user must select another database output format.
!!  INSTRUCTIONS: To utilize this module, first place a copy of it in the
!!        post-process code directory. Next, modify the module m_data_output.f90 to
!!        use m_silo_proxy.f90 and erase the include command referencing the header
!!        file silo.inc. Finally, modify the makefile so that it includes in the
!!        compilation the Silo proxy module and remove any linker flags referencing
!!        the Silo library. No further changes should be necessary to compile and
!!        execute the post-process code.
module m_silo_proxy

    implicit none

    !> @name Refer to Silo's user guide (10/2007, v4.6) for the variables' definitions
        !! and the header file silo.inc for their choice of values
    !> @{
    integer, parameter :: DB_CLOBBER = 0
    integer, parameter :: DB_COLLINEAR = 130
    integer, parameter :: DB_FLOAT = 19
    integer, parameter :: DB_DOUBLE = 20
    integer, parameter :: DB_F77NULL = -99
    integer, parameter :: DB_HDF5 = 7
    integer, parameter :: DB_LOCAL = 0
    integer, parameter :: DB_QUAD_RECT = 130
    integer, parameter :: DB_QUADVAR = 501
    integer, parameter :: DB_ZONECENT = 111
    integer, parameter :: DBOPT_EXTENTS = 300
    integer, parameter :: DBOPT_EXTENTS_SIZE = 299
    integer, parameter :: DBOPT_HI_OFFSET = 265
    integer, parameter :: DBOPT_LO_OFFSET = 266
    !> @}

contains

    !> @brief  Refer to page 28 of Silo's user guide (10/2007, v4.6) for
                !!                           information about this subroutine
   impure  function DBCREATE(pathname, lpathname, mode, target, &
                      fileinfo, lfileinfo, filetype, status)

        integer                                                   :: DBCREATE
        character(LEN=*), intent(IN) :: pathname
        integer, intent(IN) :: lpathname
        integer, intent(IN) :: mode
        integer, intent(IN) :: target
        character(LEN=*), intent(IN) :: fileinfo
        integer, intent(IN) :: lfileinfo
        integer, intent(IN) :: filetype
        integer, intent(IN) :: status

        print '(A)', 'Creation of Silo-HDF5 database file is '// &
            'inconsistent with use of Silo proxy module. '// &
            'Exiting.'
        stop

    end function DBCREATE

    !> @brief  Refer to page 235 of Silo's user guide (10/2007, v4.6)
                !!                           for information about this subroutine
    impure function DBGET2DSTRLEN()

        integer :: DBGET2DSTRLEN

        print '(A)', 'Attempt to query 2D string length is '// &
            'inconsistent with use of Silo proxy module. '// &
            'Exiting.'
        stop

    end function DBGET2DSTRLEN

    !> @brief  Refer to page 234 of Silo's user guide (10/2007, v4.6)
                !!                           for information about this subroutine
    impure function DBSET2DSTRLEN(len)

        integer                                :: DBSET2DSTRLEN
        integer, intent(IN) :: len

        print '(A)', 'Attempt to set 2D string length is '// &
            'inconsistent with use of Silo proxy module. '// &
            'Exiting.'
        stop

    end function DBSET2DSTRLEN

    !> @brief  Refer to page 185 of Silo's user guide (10/2007, v4.6)
                !!                           for information about this subroutine
    impure function DBMKOPTLIST(maxopts, optlist_id)

        integer                                :: DBMKOPTLIST
        integer, intent(IN) :: maxopts
        integer, intent(IN) :: optlist_id

        print '(A)', 'Allocation of an options list is inconsistent '// &
            'with use of Silo proxy module. Exiting.'
        stop

    end function DBMKOPTLIST

    !> @brief  Refer to page 186 of Silo's user guide (10/2007, v4.6)
                !!                           for information about this subroutine
    impure function DBADDIOPT(optlist_id, option, ivalue)

        integer                                :: DBADDIOPT
        integer, intent(IN) :: optlist_id
        integer, intent(IN) :: option
        integer, intent(IN) :: ivalue

        print '(A)', 'Adding an option to an options list is '// &
            'inconsistent with use of Silo proxy module. '// &
            'Exiting.'
        stop

    end function DBADDIOPT

    !> @brief  Refer to page 186 of Silo's user guide (10/2007, v4.6)
                !!                           for information about this subroutine
    impure function DBADDDOPT(optlist_id, option, dvalue)

        integer                                                                :: DBADDDOPT
        integer, intent(IN) :: optlist_id
        integer, intent(IN) :: option
        integer, dimension(:, :), intent(IN) :: dvalue

        print '(A)', 'Adding an option to an options list is '// &
            'inconsistent with use of Silo proxy module. '// &
            'Exiting.'
        stop

    end function DBADDDOPT

    !> @brief  Refer to page 121 of Silo's user guide (10/2007, v4.6)
                !!                           for information about this subroutine
    impure function DBPUTMMESH(dbid, name, lname, nmesh, meshnames, &
                        lmeshnames, meshtypes, optlist_id, status)

        integer                                                                                 :: DBPUTMMESH
        integer, intent(IN) :: dbid
        character(LEN=*), intent(IN) :: name
        integer, intent(IN) :: lname
        integer, intent(IN) :: nmesh
        character(LEN=*), dimension(:), intent(IN) :: meshnames
        integer, dimension(:), intent(IN) :: lmeshnames
        integer, dimension(:), intent(IN) :: meshtypes
        integer, intent(IN) :: optlist_id
        integer, intent(IN) :: status

        print '(A)', 'Writing of multi-block mesh object into '// &
            'Silo-HDF5 database file is inconsistent with '// &
            'use of Silo proxy module. Exiting.'
        stop

    end function DBPUTMMESH

    !> @brief  Refer to page 189 of Silo's user guide (10/2007, v4.6)
                !!                           for information about this subroutine
    impure function DBFREEOPTLIST(optlist_id)

        integer                                :: DBFREEOPTLIST
        integer, intent(IN) :: optlist_id

        print '(A)', 'Deallocation of an options list is inconsistent '// &
            'with use of Silo proxy module. Exiting.'
        stop

    end function DBFREEOPTLIST

    !> @brief  Refer to page 57 of Silo's user guide (10/2007, v4.6) for
                !!                           information about this subroutine
    impure function DBPUTQM(dbid, name, lname, xname, lxname, yname, lyname, &
                     zname, lzname, x, y, z, dims, ndims, datatype, &
                     coordtype, optlist_id, status)

        integer                                                                                 :: DBPUTQM
        integer, intent(IN) :: dbid
        character(LEN=*), intent(IN) :: name
        integer, intent(IN) :: lname
        character(LEN=*), intent(IN) :: xname
        integer, intent(IN) :: lxname
        character(LEN=*), intent(IN) :: yname
        integer, intent(IN) :: lyname
        character(LEN=*), intent(IN) :: zname
        integer, intent(IN) :: lzname
        real(wp), dimension(:), intent(IN) :: x
        real(wp), dimension(:), intent(IN) :: y
        real(wp), dimension(:), intent(IN) :: z
        integer, dimension(:), intent(IN) :: dims
        integer, intent(IN) :: ndims
        integer, intent(IN) :: datatype
        integer, intent(IN) :: coordtype
        integer, intent(IN) :: optlist_id
        integer, intent(IN) :: status

        print '(A)', 'Writing of quad mesh object into Silo-HDF5 '// &
            'database file is inconsistent with use of Silo '// &
            'proxy module. Exiting.'
        stop

    end function DBPUTQM

    !> @brief  Refer to page 46 of Silo's user guide (10/2007, v4.6) for
                !!                           information about this subroutine
    impure function DBPUTCURVE(dbid, curvename, lcurvename, xvals, yvals, &
                        datatype, npoints, optlist_id, status)

        integer                                                                                 :: DBPUTCURVE
        integer, intent(IN) :: dbid
        character(LEN=*), intent(IN) :: curvename
        integer, intent(IN) :: lcurvename
        real(wp), dimension(:), intent(IN) :: xvals
        real(wp), dimension(:), intent(IN) :: yvals
        integer, intent(IN) :: datatype
        integer, intent(IN) :: npoints
        integer, intent(IN) :: optlist_id
        integer, intent(IN) :: status

        print '(A)', 'Writing of curve object into Silo-HDF5 database '// &
            'file is inconsistent with use of Silo proxy '// &
            'module. Exiting.'
        stop

    end function DBPUTCURVE

    !> @brief  Refer to page 130 of Silo's user guide (10/2007, v4.6)
                !!                           for information about this subroutine
    impure function DBPUTMVAR(dbid, name, lname, nvar, varnames, lvarnames, &
                       vartypes, optlist_id, status)

        integer                                                                                 :: DBPUTMVAR
        integer, intent(IN) :: dbid
        character(LEN=*), intent(IN) :: name
        integer, intent(IN) :: lname
        integer, intent(IN) :: nvar
        character(LEN=*), dimension(:), intent(IN) :: varnames
        integer, dimension(:), intent(IN) :: lvarnames
        integer, dimension(:), intent(IN) :: vartypes
        integer, intent(IN) :: optlist_id
        integer, intent(IN) :: status

        print '(A)', 'Writing of multi-block variable into Silo-HDF5 '// &
            'database file is inconsistent with use of Silo '// &
            'proxy module. Exiting.'
        stop

    end function DBPUTMVAR

    !> @brief  Refer to page 64 of Silo's user guide (10/2007, v4.6) for
                !!                           information about this subroutine
    impure function DBPUTQV1(dbid, name, lname, meshname, lmeshname, var, &
                      dims, ndims, mixvar, mixlen, datatype, &
                      centering, optlist_id, status)

        integer                                                                                         :: DBPUTQV1
        integer, intent(IN) :: dbid
        character(LEN=*), intent(IN) :: name
        integer, intent(IN) :: lname
        character(LEN=*), intent(IN) :: meshname
        integer, intent(IN) :: lmeshname
        real(wp), dimension(:, :, :), intent(IN) :: var
        integer, dimension(:), intent(IN) :: dims
        integer, intent(IN) :: ndims
        integer, intent(IN) :: mixvar
        integer, intent(IN) :: mixlen
        integer, intent(IN) :: datatype
        integer, intent(IN) :: centering
        integer, intent(IN) :: optlist_id
        integer, intent(IN) :: status

        print '(A)', 'Writing a scalar quad variable into a Silo-HDF5 '// &
            'database file is inconsistent with use of Silo '// &
            'proxy module. Exiting.'
        stop

    end function DBPUTQV1

    !> @brief  Refer to page 31 of Silo's user guide (10/2007, v4.6) for
                !!                           information about this subroutine
    impure function DBCLOSE(dbid)

        integer                                :: DBCLOSE
        integer, intent(IN) :: dbid

        print '(A)', 'Attempt to close Silo-HDF5 database file is '// &
            'inconsistent with use of Silo proxy module. '// &
            'Exiting.'
        stop

    end function DBCLOSE

end module m_silo_proxy
