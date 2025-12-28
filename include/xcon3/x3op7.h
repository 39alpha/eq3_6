! x3op7.h
!     Strings, arrays, and such for the menu-style ("D") input
!     format for version 7.
!     Jflag strings:
!       njf7pa = the largest jflag value, including "hidden" values
!                  (> 30) that map to standard ones (0-30)
!       ujflg7 = array of strings corresponding to jflag values
!     Options (iopt, iopg, and iopr):
!       nop3pa = maximum number of options (the sum of all iopt,
!                  iopg, and iopr options)
!       nod3pa = maximum number of choices for all options in
!                  the above group
!       uopt3  = array of strings describing the above options
!       udesc3 = array of strings describing option choices
!       uvar3  = array of names of option arrays (iopt, iopg, or
!                  iopr) corresponding to the options
!       iopti3 = index of the option array element corresponding to
!                  a string describing an option choice
!       ivalu3 = array of values of option array elements corresponding
!                  to strings describing option choices
!       index3 = pointer array giving the index of the option
!                  corresponding to a string describing an option choice
!     Debug options (iodb):
!       ndb3pa = maximum number of debug options
!       udebug = array of strings describing debug options
!       idbugi = array of indices of the iodb array corresponding to
!                  strings describing debug options
!     Development options (iopt, iopg, iopr, and iodb):
!       None
!     Tolerances:
!       nto3pa = maximum number of tolerance strings
!       utol3  = array of strings describing tolerance variables
integer :: njf7pa
integer :: nop3pa
integer :: nod3pa
integer :: ndb3pa
integer :: nto3pa
parameter (njf7pa = 32,nop3pa = 15,nod3pa = 40,ndb3pa = 8,nto3pa = 4)

integer :: iopti3
integer :: ivalu3
integer :: index3
integer :: idbugi

common /x3op7/ iopti3(nod3pa),ivalu3(nod3pa),index3(nop3pa),idbugi(ndb3pa)

character(len=80) :: uopt3
character(len=80) :: udesc3
character(len=80) :: udebug
character(len=40) :: ujflg7
character(len=32) :: utol3
character(len=8) :: uvar3

common /x3op7c/ uopt3(nop3pa),udesc3(nod3pa),udebug(ndb3pa),ujflg7(-1:njf7pa),utol3(nto3pa),uvar3(nop3pa)

