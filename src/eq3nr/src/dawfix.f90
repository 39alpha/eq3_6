subroutine dawfix(aamatr,cdrs,eps100,gmmatr,iindx1,iodb,irdxc3,jflag,jjndex,kbt,kkndex,kmax,narn1,nbasp,nbtmax,ncosp,ndrs,ndrsmx,ndrsr,nelect,nhydr,nodbmx,no2gaq,noutpt,nstmax,qawfix,uspec)
    !! This subroutine determines if the activity of water is directly
    !! or indirectly fixed. For example, the user might directly specify
    !! the activity of water. Alternatively, a set of solubility
    !! constraints (e.g., gysum + anhydrite) might indirectly fix it.
    !! If the activity of water is fixed directly or indirectly, the
    !! the flag 'qawfix' is returned with a value of '.true.'.
    !! This subroutine is called by:
    !!   EQ3NR/arrset.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: kmax
    integer :: nbtmax
    integer :: ndrsmx
    integer :: nodbmx
    integer :: nstmax

    integer :: noutpt

    integer :: iindx1(kmax)
    integer :: iodb(nodbmx)
    integer :: jjndex(nbtmax)
    integer :: jflag(nstmax)
    integer :: kkndex(nbtmax)
    integer :: ncosp(nbtmax)
    integer :: nbasp(nbtmax)
    integer :: ndrs(ndrsmx)
    integer :: ndrsr(2,nstmax)

    integer :: irdxc3
    integer :: kbt
    integer :: narn1
    integer :: nelect
    integer :: nhydr
    integer :: no2gaq

    logical :: qawfix

    character(len=48) :: uspec(nstmax)

    real(kind=8) :: aamatr(kmax,kmax)
    real(kind=8) :: cdrs(ndrsmx)
    real(kind=8) :: gmmatr(kmax,kmax)

    real(kind=8) :: eps100

    ! Local variable declarations.
    integer :: ibt
    integer :: icol
    integer :: ihydr
    integer :: io2gaq
    integer :: irow
    integer :: irow1
    integer :: irow2
    integer :: iwater
    integer :: j2
    integer :: jcol
    integer :: jcol1
    integer :: jcol2
    integer :: krow
    integer :: nb
    integer :: nb1
    integer :: ns
    integer :: ns1
    integer :: ns2

    integer :: ilnobl

    logical :: qldep
    logical :: qx

    real(kind=8) :: coefdr

    qawfix = .false.

    ! Check for directly specified water activity. Note that the
    ! following assumes that water is the first aqueous species
    ! (e.g., it has the species index narn1).
    if (jflag(narn1) .eq. 16) then
        qawfix = .true.
        go to 999
    end if

    ! Build a matrix to test for water activity fixed by simultaneous
    ! equilibria.
    ibt = 0
    iwater = 0
    ihydr = 0
    io2gaq = 0

    do krow = 1,kbt
        nb = iindx1(krow)
        ns = nbasp(nb)
        kkndex(nb) = 0
        qx = .false.

        if (ns .eq. narn1) then
            qx = .true.
        else if (ns.eq.nelect .or. ns.eq.no2gaq) then
            qx = irdxc3 .ne. 0
        else
            qx = jflag(ns).eq.25 .or. jflag(ns).eq.27
        end if

        if (qx) then
            ibt = ibt + 1
            jjndex(ibt) = nb
            kkndex(nb) = 1

            if (ns .eq. narn1) then
                iwater = ibt
            end if

            if (ns .eq. nhydr) then
                ihydr = ibt
            end if

            if (ns .eq. no2gaq) then
                io2gaq = ibt
            end if
        end if
    end do

    ! Test the matrix dimension. Unless it is greater than or equal to
    ! two, the activity of water is not fixed indirectly.
    if (ibt .lt. 2) then
        go to 999
    end if

    ! Build the matrix.
    do icol = 1,ibt
        do irow = 1,ibt
            aamatr(irow,icol) = 0.
        end do
    end do

    do irow = 1,ibt
        nb = jjndex(irow)
        ns = nbasp(nb)

        if (ns .eq. narn1) then
            ! Water.
            aamatr(irow,irow) = 1.0
        else if (jflag(ns) .eq. 27) then
            ! Aqueous homogenous equilibrium.
            do icol = 1,ibt
                nb1 = jjndex(icol)
                ns1 = nbasp(nb1)

                ! Calling sequence substitutions:
                !   ns1 for nse
                aamatr(irow,icol)       = coefdr(cdrs,ndrs,ndrsmx,ndrsr,ns1,ns,nstmax)
            end do
        else if (jflag(ns) .eq. 25) then
            ! Heterogeneous equilibrium.
            ns2 = ncosp(nb)

            do icol = 1,ibt
                ns1 = jjndex(icol)

                ! Calling sequence substitutions:
                !   ns1 for nse
                !   ns2 for ns
                aamatr(irow,icol)       = coefdr(cdrs,ndrs,ndrsmx,ndrsr,ns1,ns2,nstmax)
            end do
        else if (irdxc3 .lt. 0) then
            ! Eh (an input pe- has been converted to Eh).
            aamatr(irow,iwater) = -2.
            aamatr(irow,io2gaq) = 1.

            if (ihydr .gt. 0) then
                aamatr(irow,ihydr) = 4.
            end if
        else
            ! Aqueous redox couple.
            ns2 = irdxc3

            do icol = 1,ibt
                ns1 = jjndex(icol)

                if (kkndex(ns1) .ge. 1) then
                    ! Calling sequence substitutions:
                    !   ns1 for nse
                    !   ns2 for ns
                    aamatr(irow,icol)        = coefdr(cdrs,ndrs,ndrsmx,ndrsr,ns1,ns2,nstmax)
                end if
            end do
        end if
    end do

    if (iodb(3) .ge. 3) then
        write (noutpt,1000)
1000 format(/10x,'--- Matrix to Test for Fixed Activity of Water',' ---',/)

        do irow = 1,ibt
            write (noutpt,1010) (aamatr(irow,icol), icol = 1,ibt)
1010 format(2x,10(f7.2,2x))
        end do

        write (noutpt,1020)
1020 format(/1x)

        do irow = 1,ibt
            nb = jjndex(irow)
            ns = nbasp(nb)
            j2 = ilnobl(uspec(ns)(1:24))
            write (noutpt,1030) irow,uspec(ns)(1:j2)
1030 format(2x,i3,2x,a)

            if (ns .eq. narn1) then
                write (noutpt,1040)
1040 format(/10x,'Mole fraction definition',/)
            else if (jflag(ns) .eq. 27) then
                ! Calling sequence substitutions:
                !   noutpt for nf
                call prreac(cdrs,ndrs,ndrsmx,ndrsr,noutpt,ns,nstmax,uspec)
            else if (jflag(ns) .eq. 25) then
                ns2 = ncosp(ns)

                ! Calling sequence substitutions:
                !   noutpt for nf
                !   ns2 for ns
                call prreac(cdrs,ndrs,ndrsmx,ndrsr,noutpt,ns2,nstmax,uspec)
            else if (irdxc3 .lt. 0) then
                write (noutpt,1050)
1050 format(/10x,'2 H2O(l) = 4 H+ + 4 e- + O2(g)',/)
            else
                ns2 = irdxc3

                ! Calling sequence substitutions:
                !   noutpt for nf
                !   ns2 for ns
                call prreac(cdrs,ndrs,ndrsmx,ndrsr,noutpt,ns2,nstmax,uspec)
            end if

            write (noutpt,1020)
        end do
    end if

    ! Make a copy of the matrix.
    do irow = 1,ibt
        do jcol = 1,ibt
            gmmatr(irow,jcol) = aamatr(irow,jcol)
        end do
    end do

    ! Test for linear dependence omitting the water row and column.
    ! The original matrix is overwritten in the process.
    irow1 = 2
    irow2 = ibt
    jcol1 = 2
    jcol2 = ibt
    call lindep(aamatr,eps100,irow1,irow2,jcol1,jcol2,kmax,qldep)

    if (qldep) then
        ! Have linear dependence omitting the water column. Test for
        ! linear dependence with the water column (still omitting the
        ! water row).
        ! Restore the matrix from the copy.
        do irow = 1,ibt
            do jcol = 1,ibt
                aamatr(irow,jcol) = gmmatr(irow,jcol)
            end do
        end do

        jcol1 = 1
        call lindep(aamatr,eps100,irow1,irow2,jcol1,jcol2,kmax,qldep)
        qawfix = .not.qldep
    end if

999 continue
end subroutine dawfix