! MFC v3.0 - Master Scripts: m_silo_proxy.f90
! Description: The purpose of this module is to serve as a replacement framework
!			   for the Silo library and header file when those are not available
!			   on the platform on which the raw simulation data is to be post-
!			   processed. Note that this module merely allows for the code for
!			   the post-process to be compiled and executed and does not in
!			   actuality provide support for the creation of Silo-HDF5 database
!			   files. This means that when using this module during the post-
!			   process, the user must select another database output format.
! Author: Vedran Coralic
! Date: 06/08/12


MODULE m_silo_proxy
	
	
	IMPLICIT NONE
	
	
	! INSTRUCTIONS: To utilize this module, first place a copy of it in the
	! post-process code directory. Next, modify the module m_data_output.f90 to
	! use m_silo_proxy.f90 and erase the include command referencing the header
	! file silo.inc. Finally, modify the makefile so that it includes in the
	! compilation the Silo proxy module and remove any linker flags referencing
	! the Silo library. No further changes should be necessary to compile and
	! execute the post-process code.
	
	
	! Refer to Silo's user guide (10/2007, v4.6) for the variables' definitions
	! and the header file silo.inc for their choice of values
	INTEGER, PARAMETER :: DB_CLOBBER		 = 0
	INTEGER, PARAMETER :: DB_COLLINEAR		 = 130
	INTEGER, PARAMETER :: DB_FLOAT			 = 19
	INTEGER, PARAMETER :: DB_DOUBLE			 = 20
	INTEGER, PARAMETER :: DB_F77NULL		 = -99
	INTEGER, PARAMETER :: DB_HDF5			 = 7
	INTEGER, PARAMETER :: DB_LOCAL			 = 0
	INTEGER, PARAMETER :: DB_QUAD_RECT		 = 130
	INTEGER, PARAMETER :: DB_QUADVAR		 = 501
	INTEGER, PARAMETER :: DB_ZONECENT		 = 111
	INTEGER, PARAMETER :: DBOPT_EXTENTS		 = 300
	INTEGER, PARAMETER :: DBOPT_EXTENTS_SIZE = 299
	INTEGER, PARAMETER :: DBOPT_HI_OFFSET	 = 265
	INTEGER, PARAMETER :: DBOPT_LO_OFFSET	 = 266
	
	
	CONTAINS
		
		
		
		
		
		FUNCTION DBCREATE( pathname, lpathname,   mode  , target, & ! ----------
						   fileinfo, lfileinfo, filetype, status  )
		! Description: Refer to page 28 of Silo's user guide (10/2007, v4.6) for
		!			   information about this subroutine
			
			
			INTEGER						   :: DBCREATE
			CHARACTER(LEN = *), INTENT(IN) :: pathname
			INTEGER			  , INTENT(IN) :: lpathname
			INTEGER			  , INTENT(IN) :: mode
			INTEGER			  , INTENT(IN) :: target
			CHARACTER(LEN = *), INTENT(IN) :: fileinfo
			INTEGER			  , INTENT(IN) :: lfileinfo
			INTEGER			  , INTENT(IN) :: filetype
			INTEGER			  , INTENT(IN) :: status
			
			
			PRINT '(A)', 'Creation of Silo-HDF5 database file is ' // &
						 'inconsistent with use of Silo proxy module. ' // &
						 'Exiting ...'
			STOP
			
			
		END FUNCTION DBCREATE ! ------------------------------------------------
		
		
		
		
		
		FUNCTION DBGET2DSTRLEN() ! ---------------------------------------------
		! Description: Refer to page 235 of Silo's user guide (10/2007, v4.6)
		!			   for information about this subroutine
			
			
			INTEGER :: DBGET2DSTRLEN
			
			
			PRINT '(A)', 'Attempt to query 2D string length is ' // &
						 'inconsistent with use of Silo proxy module. ' // &
						 'Exiting ...'
			STOP
			
			
		END FUNCTION DBGET2DSTRLEN ! -------------------------------------------
		
		
		
		
		
		FUNCTION DBSET2DSTRLEN(len) ! ------------------------------------------
		! Description: Refer to page 234 of Silo's user guide (10/2007, v4.6)
		!			   for information about this subroutine
			
			
			INTEGER				:: DBSET2DSTRLEN
			INTEGER, INTENT(IN) :: len
			
			
			PRINT '(A)', 'Attempt to set 2D string length is ' // &
						 'inconsistent with use of Silo proxy module. ' // &
						 'Exiting ...'
			STOP
			
			
		END FUNCTION DBSET2DSTRLEN ! -------------------------------------------
		
		
		
		
		
		FUNCTION DBMKOPTLIST(maxopts, optlist_id) ! ----------------------------
		! Description: Refer to page 185 of Silo's user guide (10/2007, v4.6)
		!			   for information about this subroutine
			
			
			INTEGER				:: DBMKOPTLIST
			INTEGER, INTENT(IN) :: maxopts
			INTEGER, INTENT(IN) :: optlist_id
			
			
			PRINT '(A)', 'Allocation of an options list is inconsistent ' // &
						 'with use of Silo proxy module. Exiting ...'
			STOP
			
			
		END FUNCTION DBMKOPTLIST ! ---------------------------------------------
		
		
		
		
		
		FUNCTION DBADDIOPT(optlist_id, option, ivalue) ! -----------------------
		! Description: Refer to page 186 of Silo's user guide (10/2007, v4.6)
		!			   for information about this subroutine
			
			
			INTEGER				:: DBADDIOPT
			INTEGER, INTENT(IN) :: optlist_id
			INTEGER, INTENT(IN) :: option
			INTEGER, INTENT(IN) :: ivalue
			
			
			PRINT '(A)', 'Adding an option to an options list is ' // &
						 'inconsistent with use of Silo proxy module. ' // &
						 'Exiting ...'
			STOP
			
			
		END FUNCTION DBADDIOPT ! -----------------------------------------------
		
		
		
		
		
		FUNCTION DBADDDOPT(optlist_id, option, dvalue) ! -----------------------
		! Description: Refer to page 186 of Silo's user guide (10/2007, v4.6)
		!			   for information about this subroutine
			
			
			INTEGER								:: DBADDDOPT
			INTEGER				   , INTENT(IN) :: optlist_id
			INTEGER				   , INTENT(IN) :: option
			INTEGER, DIMENSION(:,:), INTENT(IN) :: dvalue
			
			
			PRINT '(A)', 'Adding an option to an options list is ' // &
						 'inconsistent with use of Silo proxy module. ' // &
						 'Exiting ...'
			STOP
			
			
		END FUNCTION DBADDDOPT ! -----------------------------------------------
		
		
		
		
		
		FUNCTION DBPUTMMESH(    dbid, name, lname, nmesh, meshnames,   & ! -----
							 lmeshnames, meshtypes, optlist_id, status )
		! Description: Refer to page 121 of Silo's user guide (10/2007, v4.6)
		!			   for information about this subroutine
			
			
			INTEGER										 :: DBPUTMMESH
			INTEGER							, INTENT(IN) :: dbid
			CHARACTER(LEN = *)				, INTENT(IN) :: name
			INTEGER							, INTENT(IN) :: lname
			INTEGER							, INTENT(IN) :: nmesh
			CHARACTER(LEN = *), DIMENSION(:), INTENT(IN) :: meshnames
			INTEGER			  , DIMENSION(:), INTENT(IN) :: lmeshnames
			INTEGER			  , DIMENSION(:), INTENT(IN) :: meshtypes
			INTEGER							, INTENT(IN) :: optlist_id
			INTEGER							, INTENT(IN) :: status
			
			
			PRINT '(A)', 'Writing of multi-block mesh object into ' // &
						 'Silo-HDF5 database file is inconsistent with ' // &
						 'use of Silo proxy module. Exiting ...'
			STOP
			
			
		END FUNCTION DBPUTMMESH ! ----------------------------------------------
		
		
		
		
		
		FUNCTION DBFREEOPTLIST(optlist_id) ! -----------------------------------
		! Description: Refer to page 189 of Silo's user guide (10/2007, v4.6)
		!			   for information about this subroutine
			
			
			INTEGER				:: DBFREEOPTLIST
			INTEGER, INTENT(IN) :: optlist_id
			
			
			PRINT '(A)', 'Deallocation of an options list is inconsistent ' // &
						 'with use of Silo proxy module. Exiting ...'
			STOP
			
			
		END FUNCTION DBFREEOPTLIST ! -------------------------------------------
		
		
		
		
		
		FUNCTION DBPUTQM( dbid, name, lname, xname, lxname, yname, lyname, & ! -
						   zname, lzname, x, y, z, dims, ndims, datatype,  &
						  			coordtype, optlist_id, status		   )
		! Description: Refer to page 57 of Silo's user guide (10/2007, v4.6) for
		!			   information about this subroutine
			
			
			INTEGER										 :: DBPUTQM
			INTEGER			  ,				  INTENT(IN) :: dbid
			CHARACTER(LEN = *),				  INTENT(IN) :: name
			INTEGER			  ,				  INTENT(IN) :: lname
			CHARACTER(LEN = *),				  INTENT(IN) :: xname
			INTEGER			  ,				  INTENT(IN) :: lxname
			CHARACTER(LEN = *),				  INTENT(IN) :: yname
			INTEGER			  ,				  INTENT(IN) :: lyname
			CHARACTER(LEN = *),				  INTENT(IN) :: zname
			INTEGER							, INTENT(IN) :: lzname
			REAL(KIND(0d0))	  , DIMENSION(:), INTENT(IN) :: x
			REAL(KIND(0d0))	  , DIMENSION(:), INTENT(IN) :: y
			REAL(KIND(0d0))	  , DIMENSION(:), INTENT(IN) :: z
			INTEGER			  , DIMENSION(:), INTENT(IN) :: dims
			INTEGER			  ,				  INTENT(IN) :: ndims
			INTEGER			  ,				  INTENT(IN) :: datatype
			INTEGER			  ,				  INTENT(IN) :: coordtype
			INTEGER			  ,				  INTENT(IN) :: optlist_id
			INTEGER			  ,				  INTENT(IN) :: status
			
			
			PRINT '(A)', 'Writing of quad mesh object into Silo-HDF5 ' // &
						 'database file is inconsistent with use of Silo ' // &
						 'proxy module. Exiting ...'
			STOP
			
			
		END FUNCTION DBPUTQM ! -------------------------------------------------
		
		
		
		
		
		FUNCTION DBPUTCURVE( dbid, curvename, lcurvename, xvals, yvals, & ! ----
							   datatype, npoints, optlist_id, status	)
		! Description: Refer to page 46 of Silo's user guide (10/2007, v4.6) for
		!			   information about this subroutine
			
			
			INTEGER										 :: DBPUTCURVE
			INTEGER			  ,				  INTENT(IN) :: dbid
			CHARACTER(LEN = *),				  INTENT(IN) :: curvename
			INTEGER			  ,				  INTENT(IN) :: lcurvename
			REAL(KIND(0d0))	  , DIMENSION(:), INTENT(IN) :: xvals
			REAL(KIND(0d0))	  , DIMENSION(:), INTENT(IN) :: yvals
			INTEGER			  ,				  INTENT(IN) :: datatype
			INTEGER			  ,				  INTENT(IN) :: npoints
			INTEGER			  ,				  INTENT(IN) :: optlist_id
			INTEGER			  ,				  INTENT(IN) :: status
			
			
			PRINT '(A)', 'Writing of curve object into Silo-HDF5 database ' // &
						 'file is inconsistent with use of Silo proxy ' // &
						 'module. Exiting ...'
			STOP
			
			
		END FUNCTION DBPUTCURVE ! ----------------------------------------------
		
		
		
		
		
		FUNCTION DBPUTMVAR( dbid, name, lname, nvar, varnames, lvarnames, & ! --
									vartypes, optlist_id, status		  )
		! Description: Refer to page 130 of Silo's user guide (10/2007, v4.6)
		!			   for information about this subroutine
			
			
			INTEGER										 :: DBPUTMVAR
			INTEGER			  ,				  INTENT(IN) :: dbid
			CHARACTER(LEN = *),				  INTENT(IN) :: name
			INTEGER			  ,				  INTENT(IN) :: lname
			INTEGER			  ,				  INTENT(IN) :: nvar
			CHARACTER(LEN = *), DIMENSION(:), INTENT(IN) :: varnames
			INTEGER			  , DIMENSION(:), INTENT(IN) :: lvarnames
			INTEGER			  , DIMENSION(:), INTENT(IN) :: vartypes
			INTEGER			  ,				  INTENT(IN) :: optlist_id
			INTEGER			  ,				  INTENT(IN) :: status
			
			
			PRINT '(A)', 'Writing of multi-block variable into Silo-HDF5 ' // &
						 'database file is inconsistent with use of Silo ' // &
						 'proxy module. Exiting ...'
			STOP
			
			
		END FUNCTION DBPUTMVAR ! -----------------------------------------------
		
		
		
		
		
		FUNCTION DBPUTQV1( dbid, name, lname, meshname, lmeshname, var, & ! ----
						      dims, ndims, mixvar, mixlen, datatype,	&
						          centering, optlist_id, status			)
		! Description: Refer to page 64 of Silo's user guide (10/2007, v4.6) for
		!			   information about this subroutine
			
			
			INTEGER											 :: DBPUTQV1
			INTEGER			  ,					  INTENT(IN) :: dbid
			CHARACTER(LEN = *),					  INTENT(IN) :: name
			INTEGER			  ,					  INTENT(IN) :: lname
			CHARACTER(LEN = *),					  INTENT(IN) :: meshname
			INTEGER			  ,					  INTENT(IN) :: lmeshname
			REAL(KIND(0d0))	  , DIMENSION(:,:,:), INTENT(IN) :: var
			INTEGER			  , DIMENSION(:)	, INTENT(IN) :: dims
			INTEGER			  ,					  INTENT(IN) :: ndims
			INTEGER			  ,					  INTENT(IN) :: mixvar
			INTEGER			  ,					  INTENT(IN) :: mixlen
			INTEGER			  ,					  INTENT(IN) :: datatype
			INTEGER			  ,					  INTENT(IN) :: centering
			INTEGER			  ,					  INTENT(IN) :: optlist_id
			INTEGER			  ,					  INTENT(IN) :: status
			
			
			PRINT '(A)', 'Writing a scalar quad variable into a Silo-HDF5 ' // &
						 'database file is inconsistent with use of Silo ' // &
						 'proxy module. Exiting ...'
			STOP
			
			
		END FUNCTION DBPUTQV1 ! ------------------------------------------------
		
		
		
		
		
		FUNCTION DBCLOSE(dbid) ! -----------------------------------------------
		! Description: Refer to page 31 of Silo's user guide (10/2007, v4.6) for
		!			   information about this subroutine
			
			
			INTEGER				:: DBCLOSE
			INTEGER, INTENT(IN) :: dbid
			
			
			PRINT '(A)', 'Attempt to close Silo-HDF5 database file is ' // &
						 'inconsistent with use of Silo proxy module. ' // &
						 'Exiting ...'
			STOP
			
			
		END FUNCTION DBCLOSE ! -------------------------------------------------
		
		
		
		
		
END MODULE m_silo_proxy
