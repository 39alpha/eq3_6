! eql1s.h
!     The presssure on the 1.013-bar/steam-saturation curve. The
!     variable is presh. The coefficients for computing it are in the
!     apresh array, which has fixed dimensions (apresh(5,2)). The values
!     of the coefficients are stored in a DATA statement in the EQLIB
!     INCLUDE file eql1sd.h.
!     Note that presg and apresg (which are in EQLIB/eqlxf.h) are,
!     respectively, the pressure and the corresponding coefficients) for
!     the data file reference pressure curve. This may or may not be the
!     1.013-bar/steam-saturation curve. Even if it is, the two
!     representations may not exactly coincide.
integer :: narxth

common /eql1sa/ narxth(2)

real(kind=8) :: apresh

common /eql1sb/ apresh(5,2)

real(kind=8) :: presh

common /eql1sc/ presh

