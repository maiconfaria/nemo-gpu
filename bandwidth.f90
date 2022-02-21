program BandwidthTest

  use cudafor
  implicit none

  integer :: nElements

  ! host arrays
  real :: a_pageable(4*1024*1024), b_pageable(4*1024*1024)
  real, allocatable :: a_pinned(:), b_pinned(:)

  ! device arrays
  real, device :: a_d(4*1024*1024)

  ! events for timing
  type (cudaEvent) :: startEvent, stopEvent

  ! misc
  type (cudaDeviceProp) :: prop
  real :: time
  integer :: istat, i
  logical :: pinnedFlag

  INTEGER, parameter                ::  nstreams = 3
  INTEGER(kind=cuda_stream_kind)    ::  stream(nstreams), str               ! Stream ID

  nElements = 4*1024*1024
  ! allocate and initialize
  do i = 1, nElements
    a_pageable(i) = i
  end do
  b_pageable = 0.0

  allocate(a_pinned(nElements), b_pinned(nElements), &
           STAT=istat)
  DO i = 1, nstreams
          istat = cudaStreamCreate(stream(i))
          IF(istat /= 0) print *, 'Error in Stream creation dia_hsb_init', i
  END DO
! test reallocation to pin
  istat = cudaHostRegister(C_LOC(a_pinned ),sizeof(a_pinned),cudaHostRegisterMapped)
  istat = cudaHostRegister(C_LOC(b_pinned ),sizeof(b_pinned),cudaHostRegisterMapped)
  pinnedFlag =  .true.
  !

if (istat /= 0) then
    write(*,*) 'Allocation of a_pinned/b_pinned failed'
    pinnedFlag = .false.
  else
    if (.not. pinnedFlag) write(*,*) 'Pinned allocation failed'
  end if

  if (pinnedFlag) then
    a_pinned = a_pageable
    b_pinned = 0.0
  endif

  istat = cudaEventCreate(startEvent)
  istat = cudaEventCreate(stopEvent)

  ! output device info and transfer size
  istat = cudaGetDeviceProperties(prop, 0)

  write(*,*)
  write(*,*) 'Device: ', trim(prop%name)
  write(*,*) 'Transfer size (MB): ', 4*nElements/1024./1024.

  ! pageable data transfers
  write(*,*)
  write(*,*) 'Pageable transfers'

  istat = cudaEventRecord(startEvent, 0)
  a_d = a_pageable
  istat = cudaEventRecord(stopEvent, 0)
  istat = cudaEventSynchronize(stopEvent)

  istat = cudaEventElapsedTime(time, startEvent, stopEvent)
  write(*,*) '  Host to Device bandwidth (GB/s): ', &
    nElements*4*1e-6/time

  istat = cudaEventRecord(startEvent, 0)
  b_pageable = a_d
  istat = cudaEventRecord(stopEvent, 0)
  istat = cudaEventSynchronize(stopEvent)

  istat = cudaEventElapsedTime(time, startEvent, stopEvent)
  write(*,*) '  Device to Host bandwidth (GB/s): ', &
    nElements*4*1e-6/time

  if (any(a_pageable /= b_pageable)) &
    write(*,*) '*** Pageable transfers failed ***'

  ! pinned data transfers
  if (pinnedFlag) then
    write(*,*)
    write(*,*) 'Pinned transfers'

    istat = cudaEventRecord(startEvent, 0)
    a_d = a_pinned
    istat = cudaEventRecord(stopEvent, 0)
    istat = cudaEventSynchronize(stopEvent)

    istat = cudaEventElapsedTime(time, startEvent, stopEvent)
    write(*,*) '  Host to Device bandwidth (GB/s): ', &
      nElements*4*1e-6/time

    istat = cudaEventRecord(startEvent, 0)
    b_pinned = a_d
    istat = cudaEventRecord(stopEvent, 0)
    istat = cudaEventSynchronize(stopEvent)

    istat = cudaEventElapsedTime(time, startEvent, stopEvent)
    write(*,*) '  Device to Host bandwidth (GB/s): ', &
      nElements*4*1e-6/time

    if (any(a_pinned /= b_pinned)) &
      write(*,*) '*** Pinned transfers failed ***'
  end if

  write(*,*)

  ! cleanup
  if (allocated(a_pinned)) deallocate(a_pinned)
  if (allocated(b_pinned)) deallocate(b_pinned)
  istat = cudaEventDestroy(startEvent)
  istat = cudaEventDestroy(stopEvent)

end program BandwidthTest
