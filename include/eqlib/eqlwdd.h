! eqlwdd.h
!     This is a pseudo-data file for constants needed for the
!     WIPP brine density model:
!       density (g/L) = a + b x TDS (g/L)
!     This is based on data for pure NaCl solutions at 20C.
!     The adwipp (a) and bdwipp (b) variables are declared in
!     EQLIB/eqlwd.h.
!     This pseudo-data file is referenced by:
!       EQ3NR/bkdeq3.f
!       EQ6/bkdeq6.f
data adwipp /1000.96/, bdwipp /0.639963/

