! x6op7.h
!     Strings, arrays, and such for the menu-style ("D") input
!     format for version 7.
!     Options (iopt, iopg, iopr, and iodb):
!       nop6pa = maximum number of options (the sum of all iopt,
!                  iopg, iopr, and iodb options)
!       nod6pa = maximum number of choices for all options
!                  in the above group
!       uopt6  = array of strings describing the above options
!       udesc6 = array of strings describing option choices
!       uvar6  = array of names of option arrays (iopt, iopg, iopr, or
!                  iodb) corresponding to the options
!       iopti6 = array of indices of option array elements corresponding
!                  to strings describing option choices
!       ivalu6 = array of values of option array elements corresponding
!                  to strings describing option choices
!       index6 = pointer array giving the index of the option
!                  corresponding to a string describing an option choice
!     Development options (iopt, iopg, iopr, and iodb):
!       ndv6pa = maximum number of development options
!       udevl6 = array of strings describing development options
!       udv6vr = array of names of option arrays (iopt, iopg, iopr, or
!                  iodb) corresponding to the development options
!       idev6i = array of values of option array elements corresponding
!                  to strings describing development option choices
!     Tolerances:
!       nto6pa = maximum number of tolerance variables
!       utol6  = array of tolerance descriptor strings
integer :: nop6pa
integer :: nod6pa
integer :: ndv6pa
integer :: nto6pa
parameter (nop6pa = 34,nod6pa = 85,ndv6pa = 2,nto6pa = 23)

integer :: iopti6
integer :: ivalu6
integer :: index6
integer :: idev6i

common /x6op7/ iopti6(nod6pa),ivalu6(nod6pa),index6(nop6pa),idev6i(ndv6pa)

character(len=80) :: udesc6
character(len=80) :: uopt6
character(len=80) :: udevl6
character(len=32) :: utol6
character(len=8) :: uvar6
character(len=8) :: udv6vr

common /x6op7c/ udesc6(nod6pa),udevl6(ndv6pa),udv6vr(ndv6pa),utol6(nto6pa),uopt6(nop6pa),uvar6(nop6pa)

