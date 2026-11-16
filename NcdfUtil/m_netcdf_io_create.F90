!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: m_netcdf_io_create.F90
!
! !INTERFACE:
!
module m_netcdf_io_create
!
! !USES:
!
  use iso_c_binding, only: c_int, c_size_t, c_ptr, c_char, c_null_char, &
                          c_null_ptr, c_associated

  implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
  public  Nccr_Wr
  public  Ncdo_Sync
  public  Ncdo_Memory_Close
!
! !PRIVATE MEMBER FUNCTIONS:
!
  private :: memory_create, memory_close, memory_free
  private :: disk_open, disk_write, disk_close
!
! !PRIVATE C INTEROPERABILITY SYMBOLS:
!
  private :: c_int, c_size_t, c_ptr, c_char, c_null_char, c_null_ptr
  private :: c_associated
!
! !DESCRIPTION: Routines for creating and synchronizing netCDF files,
!  including in-memory files saved to disk when closed.
!\\
!\\
! !AUTHOR:
!  Jules Kouatchou
!
! !REMARKS:
!  This file is based on code from NASA/GSFC, SIVO, Code 610.3
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/ncdfutil for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !PRIVATE PARAMETERS:
!
  ! Maximum number of in-memory files that can be open at the same time
  integer, parameter, private :: max_memory_files = 256
!
! !PRIVATE VARIABLES:
!
  ! Map each open memory-file ID to its destination on disk
  integer, save, private :: memory_ids(max_memory_files) = -1
  character(len=1024), save, private :: memory_names(max_memory_files) = ''
!
! !PRIVATE TYPES:
!
  ! C-compatible netCDF memory-image descriptor (NC_memio)
  type, bind(C), private :: memory_image
     integer(c_size_t) :: size
     type(c_ptr)      :: memory
     integer(c_int)   :: flags
  end type memory_image
!
! !PRIVATE INTERFACES:
!
  ! netCDF memory-image operations and C file I/O used for persistence
  interface
     function memory_create(path, mode, initial_size, ncid) &
          bind(C, name='nc_create_mem') result(status)
       import :: c_char, c_int, c_size_t
       character(kind=c_char), intent(in) :: path(*)
       integer(c_int), value              :: mode
       integer(c_size_t), value           :: initial_size
       integer(c_int), intent(out)        :: ncid
       integer(c_int)                     :: status
     end function memory_create

     function memory_close(ncid, image) &
          bind(C, name='nc_close_memio') result(status)
       import :: c_int, memory_image
       integer(c_int), value             :: ncid
       type(memory_image), intent(inout) :: image
       integer(c_int)                    :: status
     end function memory_close

     function disk_open(path, mode) bind(C, name='fopen') result(stream)
       import :: c_char, c_ptr
       character(kind=c_char), intent(in) :: path(*), mode(*)
       type(c_ptr)                       :: stream
     end function disk_open

     function disk_write(buffer, element_size, count, stream) &
          bind(C, name='fwrite') result(written)
       import :: c_ptr, c_size_t
       type(c_ptr), value       :: buffer, stream
       integer(c_size_t), value :: element_size, count
       integer(c_size_t)        :: written
     end function disk_write

     function disk_close(stream) bind(C, name='fclose') result(status)
       import :: c_ptr, c_int
       type(c_ptr), value :: stream
       integer(c_int)    :: status
     end function disk_close

     subroutine memory_free(buffer) bind(C, name='free')
       import :: c_ptr
       type(c_ptr), value :: buffer
     end subroutine memory_free
  end interface
!
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Nccr_Wr
!
! !INTERFACE:
!
  subroutine Nccr_Wr (ncid, filname, WRITE_NC4, MemoryFile)
!
! !USES:
!
    use netCDF
    use m_do_err_out
!
! !INPUT PARAMETERS:
!   ncid       : opened netCDF file id
!   filname    : name of netCDF file to open for writing
!   WRITE_NC4  : optional flag to create a netCDF-4 file
!   MemoryFile : optional flag to keep the file in RAM until it is closed
    integer          , intent(INOUT) :: ncid
    character (len=*), intent(IN)    :: filname
    LOGICAL, OPTIONAL, INTENT(IN)    :: WRITE_NC4
    LOGICAL, OPTIONAL, INTENT(IN)    :: MemoryFile
!
! !DESCRIPTION: Creates a netCDF file for writing and does some error checking.
!\\
!\\
! !AUTHOR:
!  John Tannahill (LLNL) and Jules Kouatchou
!
! !REMARKS:
!  If the netCDF4 library is used, then the NF90_CLOBBER flag will write
!  a classic (i.e. netCDF3) file.  Use OR(NF_NETCDF4,NF_CLASSIC_MODEL) to
!  create netCDF-4 file that supports compression and uses "classic"
!  netcdf data model (no groups, no user-defined types)
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/ncdfutil for complete history
!EOP
!-------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
    character (len=128) :: err_msg
    integer             :: ierr
    INTEGER             :: mode
    LOGICAL             :: TMP_NC4
    LOGICAL             :: keep_in_memory
    INTEGER             :: slot
!
    ! Save the value of the optional WRITE_NC4 variable in
    ! a local shadow variable (bmy, 11/7/11)
    IF ( PRESENT( WRITE_NC4 ) ) THEN
       TMP_NC4 = WRITE_NC4
    ELSE
       TMP_NC4 = .FALSE.
    ENDIF

    IF ( TMP_NC4 ) THEN
#ifdef NC_HAS_COMPRESSION
       mode = IOR( NF90_NETCDF4, NF90_CLASSIC_MODEL )
#else
       mode = NF90_64BIT_OFFSET
#endif
    ELSE
       mode = NF90_CLOBBER
    ENDIF

    keep_in_memory = .FALSE.
    IF ( PRESENT( MemoryFile ) ) keep_in_memory = MemoryFile
    IF ( keep_in_memory ) THEN
       slot = 1
       DO WHILE ( slot <= max_memory_files )
          IF ( memory_ids(slot) == -1 ) EXIT
          slot = slot + 1
       ENDDO
       IF ( slot > max_memory_files ) THEN
          err_msg = 'In Nccr_Wr: too many in-memory HISTORY files'
          call Do_Err_Out (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0)
       ENDIF
       ierr = memory_create( TRIM( filname ) // c_null_char, mode, &
                             0_c_size_t, ncid )
       IF ( ierr == NF90_NOERR ) THEN
          memory_ids(slot) = ncid
          memory_names(slot) = filname
       ENDIF
    ELSE
       ierr = NF90_Create(filname, mode, ncid)
    ENDIF

    if (ierr /= NF90_NOERR) then
       err_msg = 'In Nccr_Wr, cannot create:  ' // Trim (filname)
       call Do_Err_Out (err_msg, .true., 0, 0, 0, 0 , 0.0d0, 0.0d0)
    end if

    return

  end subroutine Nccr_Wr
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncdo_Sync
!
! !INTERFACE:
!
  subroutine Ncdo_Sync (ncid)
!
! !USES:
!
    use netCDF
    use m_do_err_out
!
! !INPUT PARAMETERS:
!!  ncid : netCDF file id
    integer, intent(in)   :: ncid
!
! !DESCRIPTION: Synchronizes a netCDF file.
!\\
!\\
! !AUTHOR:
!  John Tannahill (LLNL) and Jules Kouatchou
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/ncdfutil for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
    character (len=128) :: err_msg
    integer             :: ierr
!
    ierr = NF90_Sync (ncid)

    if (ierr /= NF90_NOERR) then
       err_msg = 'In Ncdo_Sync:  ' // NF90_Strerror (ierr)
       call Do_Err_Out (err_msg, .true., 1, ncid, 0, 0, 0.0d0, 0.0d0)
    end if

  end subroutine Ncdo_Sync
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Ncdo_Memory_Close
!
! !INTERFACE:
!
  subroutine Ncdo_Memory_Close (ncid, Handled, ierr)
!
! !USES:
!
    use netCDF, only: NF90_NOERR, NF90_EIO
!
! !INPUT PARAMETERS:
!
    integer, intent(in)  :: ncid     ! netCDF file ID
!
! !OUTPUT PARAMETERS:
!
    logical, intent(out) :: Handled  ! True if this is a registered memory file
    integer, intent(out) :: ierr     ! netCDF status code
!
! !DESCRIPTION: Closes an in-memory netCDF file, saves its image to disk,
!  and releases the image buffer. Ordinary files return Handled = false
!  so that the caller can use its existing netCDF close operation.
!
! !REMARKS:
!  The file must have been registered by Nccr_Wr with MemoryFile = true.
!  File-open, short-write, and stream-close failures return NF90_EIO.
!
! !REVISION HISTORY:
!  See https://github.com/geoschem/ncdfutil for complete history
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:
!
    ! Scalars
    integer           :: slot
    integer(c_int)    :: close_status
    integer(c_size_t) :: written

    ! File image and output stream
    type(memory_image) :: image
    type(c_ptr)        :: stream

    ! Initialize status and search for a registered memory file
    Handled = .false.
    ierr    = NF90_NOERR
    do slot = 1, max_memory_files
       if ( memory_ids(slot) /= ncid ) cycle
       Handled = .true.

       ! Close the netCDF image and obtain its buffer
       image%size   = 0_c_size_t
       image%memory = c_null_ptr
       image%flags  = 0_c_int
       ierr = memory_close (ncid, image)
       if ( ierr /= NF90_NOERR ) return

       ! Save the completed image to its destination on disk
       stream = disk_open( trim( memory_names(slot) ) // c_null_char, &
                           'wb' // c_null_char )
       if ( c_associated( stream ) ) then
          written = disk_write (image%memory, 1_c_size_t, image%size, stream)
          close_status = disk_close (stream)
          if ( written /= image%size .or. close_status /= 0 ) ierr = NF90_EIO
       else
          ierr = NF90_EIO
       end if

       ! Release the buffer and clear the file's registry entry
       if ( c_associated( image%memory ) ) call memory_free (image%memory)
       memory_ids(slot)   = -1
       memory_names(slot) = ''
       return
    end do

  end subroutine Ncdo_Memory_Close
!EOC
end module m_netcdf_io_create
