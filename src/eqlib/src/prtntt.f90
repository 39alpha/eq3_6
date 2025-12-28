subroutine prtntt(nat,nata,natmax,nbt,nbta,nbtmax,nct,ncta,nctmax,ngt,ngta,ngtmax,nlt,nlta,nltmax,nmt,nmta,nmtmax,noutpt,npt,npta,nptmax,nst,nsta,nstmax,nxt,nxta,nxtmax)
    !! This subroutine prints a table of statistics for species, phases,
    !! and groups thereof showing for each entity the number on the data
    !! file, the number the software is dimensioned for, and the number
    !! appearing in the current problem.
    !! This subroutine is called by:
    !!   EQ3NR/echox.f
    !!   EQ6/echoz.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: natmax
    integer :: nbtmax
    integer :: nctmax
    integer :: ngtmax
    integer :: nltmax
    integer :: nmtmax
    integer :: nptmax
    integer :: nstmax
    integer :: nxtmax

    integer :: noutpt

    integer :: nat
    integer :: nata
    integer :: nbt
    integer :: nbta
    integer :: nct
    integer :: ncta
    integer :: ngt
    integer :: ngta
    integer :: nlt
    integer :: nlta
    integer :: nmt
    integer :: nmta
    integer :: npt
    integer :: npta
    integer :: nst
    integer :: nsta
    integer :: nxt
    integer :: nxta

    ! Local variable declarations.
    character(len=24) :: ux

    write (noutpt,1000)
1000 format(/7x,'--- Numbers of Phases, Species, and Groups Thereof','---',//7x,'Entity',15x,'Date Base',4x,'Dimension',3x,'Current Problem',/)

    ux = 'Chemical Elements       '
    write (noutpt,1010) ux,ncta,nctmax,nct
1010 format(3x,a24,3x,i5,7x,i5,7x,i5)

    ux = 'Basis Species'
    write (noutpt,1010) ux,nbta,nbtmax,nbt

    ux = 'Phases'
    write (noutpt,1010) ux,npta,nptmax,npt

    ux = 'Species'
    write (noutpt,1010) ux,nsta,nstmax,nst

    ux = 'Aqueous Species'
    write (noutpt,1010) ux,nata,natmax,nat

    ux = 'Pure Minerals'
    write (noutpt,1010) ux,nmta,nmtmax,nmt

    ux = 'Pure Liquids'
    write (noutpt,1010) ux,nlta,nltmax,nlt

    ux = 'Gas Species'
    write (noutpt,1010) ux,ngta,ngtmax,ngt

    ux = 'Solid Soutions'
    write (noutpt,1020) ux,nxta,nxtmax,nxt
1020 format(3x,a24,3x,i5,7x,i5,7x,i5,/)
end subroutine prtntt