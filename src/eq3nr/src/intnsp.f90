subroutine intnsp(coval,covali,ier,jflag,narn1,narn2,nbasp,nbt,nbti,nbtmax,nchlor,ncosp,ndecsp,nhydr,noutpt,nst,nstmax,nttyo,ucospi,uspec)
    !! This subroutine finds the indices of species required to evaluate
    !! certain types of constraints, such as a mean activity or a
    !! heterogeneous equilibrium.
    !! This subroutine is called by:
    !!   EQ3NR/eq3nr.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: nbtmax
    integer :: nstmax

    integer :: noutpt
    integer :: nttyo

    integer :: jflag(nstmax)
    integer :: nbasp(nbtmax)
    integer :: ncosp(nbtmax)
    integer :: ndecsp(nbtmax)

    integer :: ier
    integer :: narn1
    integer :: narn2
    integer :: nbt
    integer :: nbti
    integer :: nchlor
    integer :: nhydr
    integer :: nst

    character(len=48) :: ucospi(nbtmax)
    character(len=48) :: uspec(nstmax)

    real(kind=8) :: coval(nbtmax)
    real(kind=8) :: covali(nbtmax)

    ! Local variable declarations.
    integer :: jfl
    integer :: jlen
    integer :: j2
    integer :: n
    integer :: nbi
    integer :: nb1
    integer :: nb2
    integer :: ns1
    integer :: ns2

    character(len=56) :: uspn56
    character(len=24) :: ublk24

    data ublk24 /'                        '/

    ier = 0

    do nbi = 1,nbti
        nb1 = ndecsp(nbi)

        if (nb1 .le. 0) then
            go to 300
        end if

        ns1 = nbasp(nb1)
        jfl = jflag(ns1)
        coval(nb1) = covali(nbi)

        if (jfl.eq.17 .or. jfl.eq.18) then
            ! Find the species index of the counterion for a mean log
            ! activity or neutral log activity combination.
            do nb2 = 1,nbt
                ns2 = nbasp(nb2)

                if (ucospi(nbi)(1:24) .eq. uspec(ns2)(1:24)) then
                    ncosp(nb1) = ns2
                    go to 300
                end if
            end do

            ier = ier + 1

            ! Calling sequence substitutions:
            !   ucospi(nbi) for unam48
            call fmspnx(jlen,ucospi(nbi),uspn56)

            if (jfl .eq. 17) then
                write (noutpt,1000) uspn56(1:jlen)
                write (nttyo,1000) uspn56(1:jlen)
1000 format(/' * Error - (EQ3NR/intnsp) The species ',a,/7x,'is required for an abs(zj)*log ai + abs(zi)*log aj',/7x,"constraint (jflag= 17), but it isn't otherwise",' present in the',/7x,'current system.')
            else if (jfl .eq. 18) then
                write (noutpt,1010) uspn56(1:jlen)
                write (nttyo,1010) uspn56(1:jlen)
1010 format(/' * Error - (EQ3NR/intnsp) The species ',a,/7x,'is required for a log a(+/-,ij) constraint',' (jflag= 18),',/7x,"but it isn't otherwise present",' in the current system.')
            end if
        else if (jfl .eq. 21) then
            ! Find the counterion for the pHCl constraint. Here the
            ! constraint may be applied to either H+ or Cl-.
            if (ns1 .eq. nhydr) then
                ucospi(nbi) = 'Cl-'
                j2 = 3

                if (nchlor .gt. 0) then
                    ncosp(nb1) = nchlor
                else
                    ier = ier + 1
                    write (noutpt,1020) ucospi(nbi)(1:j2)
                    write (nttyo,1020) ucospi(nbi)(1:j2)
1020 format(/' * Error - (EQ3NR/intnsp) The species ',a,/7x,'is required for a pHCl constraint (jflag= 21), but',/7x,"it isn't otherwise present in the current system.")
                end if
            else if (ns1 .eq. nchlor) then
                ucospi(nbi) = 'H+'
                j2 = 2

                if (nhydr .gt. 0) then
                    ncosp(nb1) = nhydr
                else
                    ier = ier + 1
                    write (noutpt,1020) ucospi(nbi)(1:j2)
                    write (nttyo,1020) ucospi(nbi)(1:j2)
                end if
            else
                ier = ier + 1

                ! Calling sequence substitutions:
                !   uspec(ns1) for unam48
                call fmspnx(jlen,uspec(ns1),uspn56)
                write (noutpt,1040) uspn56(1:jlen)
                write (nttyo,1040) uspn56(1:jlen)
1040 format(/' * Error - (EQ3NR/intnsp) The pHCl constraint',' (jflag= 21)',/7x,"can't be applied to the species ",a,'.',/7x,'It can only be applied to H+ or Cl-.')
            end if
        else if (jfl .eq. 25) then
            ! Find the species index of the species whose heterogeneous
            ! reaction is used to define an equilibrium constraint.
            n = narn1 - 1

            if (n .ge. 1) then
                do ns2 = 1,n
                    if ( ucospi(nbi)(1:24) .eq. uspec(ns2)(1:24) .and.           ( ucospi(nbi)(25:48) .eq. uspec(ns2)(25:48) .or.             ucospi(nbi)(25:48) .eq. ublk24(1:24) ) ) then
                        ncosp(nb1) = ns2
                        go to 300
                    end if
                end do
            end if

            n = narn2 + 1

            if (nst .ge. n) then
                do ns2 = n,nst
                    if ( ucospi(nbi)(1:24) .eq. uspec(ns2)(1:24) .and.           ( ucospi(nbi)(25:48) .eq. uspec(ns2)(25:48) .or.             ucospi(nbi)(25:48) .eq. ublk24(1:24) ) ) then
                        ncosp(nb1) = ns2
                        go to 300
                    end if
                end do
            end if

            ier = ier + 1

            ! Calling sequence substitutions:
            !   ucospi(nbi) for unam48
            call fmspnx(jlen,ucospi(nbi),uspn56)
            write (noutpt,1050) uspn56(1:jlen)
            write (nttyo,1050) uspn56(1:jlen)
1050 format(/' * Error - (EQ3NR/intnsp) The species ',a,/7x,'is required for a heterogeneous equilibrium',' constraint',/7x,"(jflag= 25), but it isn't otherwise",' present in the current system.')
        end if

300 continue
    end do
end subroutine intnsp