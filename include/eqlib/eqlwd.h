! eqlwd.h
!     This constants needed for the WIPP brine density model:
!       density (g/L) = a + b x TDS (g/L)
!     This is based on data for pure NaCl solutions at 20C.
!     The respective variables corresponding to a and b are adwipp
!     and bdwipp. The values of the coefficients are stored in a
!     DATA statement in the EQLIB INCLUDE file eqlwdd.h.
real(kind=8) :: adwipp
real(kind=8) :: bdwipp

common /eqlwda/ adwipp,bdwipp

