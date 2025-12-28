subroutine prreac(cdrs,ndrs,ndrsmx,ndrsr,nf,ns,nstmax,uspec)
    !! This subroutine writes the n-th reaction in a set on the file
    !! whose unit number isf nf.
    !! This subroutine is called by:
    !!   EQLIB/alters.f
    !!   EQLIB/echolk.f
    !!   EQLIB/switch.f
    !!   EQLIB/swtchb.f
    !!   EQLIB/swtchk.f
    !!   EQ3NR/arrsim.f
    !!   EQ3NR/dawfix.f
    !!   EQ3NR/echox.f
    !!   EQ3NR/eq3nr.f
    !! Principal input:
    !!   cdrs   = the array of reaction coefficients
    !!   ndrs   = the array of corresponding species indices
    !!   ndrsmx = dimension of the cdrs and ndrs arrays
    !!   ndrsr  = aray giving the range in the cdrs and ndrs arrays
    !!              corresponding to the reaction for a given species
    !!   nf     = the unit number of the file to write on
    !!   ns     = the index of the species whose reaction is to
    !!              be written
    !!   uspec  = array of species names
    !! Principal output:
    !!   None
    implicit none

    ! Calling sequence variable declarations.
    integer :: ndrsmx
    integer :: nstmax

    integer :: ndrs(ndrsmx)
    integer :: ndrsr(2,nstmax)
    integer :: nf
    integer :: ns

    character(len=48) :: uspec(nstmax)

    real(kind=8) :: cdrs(ndrsmx)

    ! Local variable declarations.
    integer :: jlen
    integer :: n
    integer :: nr1
    integer :: nr2
    integer :: nss
    integer :: nt

    logical :: qdtach
    logical :: qfirst

    character(len=56) :: uspn56

    real(kind=8) :: cx

    write (nf,1000)
1000 format(1x)

    nr1 = ndrsr(1,ns)
    nr2 = ndrsr(2,ns)

    nt = nr2 - nr1 + 1

    if (nt .lt. 2) then
        ! Calling sequence substitutions:
        !   uspec(ns) for unam48
        call fmspnx(jlen,uspec(ns),uspn56)

        write (nf,1010) uspn56(1:jlen)
1010 format(3x,a,' is a strict basis species and has no reaction.',/)

        go to 999
    end if

    qdtach = .false.

    do n = nr1 + 1,nr2
        nss = ndrs(n)

        if (nss .eq. 0) then
            qdtach = .true.
            go to 100
        end if
    end do

100 continue

    if (qdtach) then
        ! Calling sequence substitutions:
        !   uspec(ns) for unam48
        call fmspnx(jlen,uspec(ns),uspn56)

        write (nf,1020) uspn56(1:jlen)
1020 format(3x,a,' is a detached auxiliary basis species and',' in effect',/3x,'has no reaction.',/)

        go to 999
    end if

    ! Write the coefficients and names of the reactants.
    qfirst = .true.

    do n = nr1,nr2
        cx = -cdrs(n)
        nss = ndrs(n)

        if (cx .gt. 0.) then
            ! Calling sequence substitutions:
            !   uspec(nss) for unam48
            call fmspnx(jlen,uspec(nss),uspn56)

            if (qfirst) then
                write (nf,1030) cx,uspn56(1:jlen)
1030 format(4x,f7.3,2x,a)

                qfirst = .false.
            else
                write (nf,1050) cx,uspn56(1:jlen)
1050 format(2x,'+ ',f7.3,2x,a)
            end if
        end if
    end do

    write (nf,1070)
1070 format(10x,'==')

    ! Write the coefficients and names of the products.
    qfirst = .true.

    do n = nr1,nr2
        cx = cdrs(n)
        nss = ndrs(n)

        if (cx .gt. 0.) then
            ! Calling sequence substitutions:
            !   uspec(nss) for unam48
            call fmspnx(jlen,uspec(nss),uspn56)

            if (qfirst) then
                write (nf,1030) cx,uspn56(1:jlen)
                qfirst = .false.
            else
                write (nf,1050) cx,uspn56(1:jlen)
            end if
        end if
    end do

999 continue
end subroutine prreac