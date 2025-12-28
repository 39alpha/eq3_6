subroutine alters(afcnst,apresg,axlks,cdrs,kxmod,narxmx,narxt,ndrs,ndrsmx,ndrsr,noutpt,npt,nptmax,nst,nstmax,ntpr,ntprmx,nttyo,nxmdmx,nxmod,tempc,uphase,uspec,uxmod,xlkmod)
    !! This subroutine alters the log K values of reactions for the
    !! destruction of the specified species after they are read from
    !! the data file. This is done before the reactions are rewritten
    !! by any basis switching. The reactions are identified on the
    !! input file by the names (uxmod) of the associated species.
    !! The log K values (actually here the corresponding interpolating
    !! polynomial coefficients) are altered according the following
    !! options:
    !!   kxmod = 0   Replace the log K value by xlkmod
    !!   kxmod = 1   Augment the log K value by xlkmod
    !!   kxmod = 2   Destablize the associated species by xlkmod kcal
    !!                (at the temperature at the start of the run)
    !! This subroutine maps any kxmod = 0 or 2 options into the
    !! equivalent kxmod = 1 option. An unaltered pickup file (written
    !! by EQ3NR or EQ6) therefore may have only the kxmod = 1 alter
    !! option in addition to the kxmod = -1 suppression option (which
    !! is handled by EQLIB/supprs.f). The effect of this is that if
    !! temperature in EQ6 is a function of reaction progress, the
    !! log K at any temperature is augmented by a constant amount.
    !! Note that uxmod is a 48-character name with a 24-character
    !! species part followed by a 24-character phase part. If the phase
    !! part is blank, which would normally be the case, name matching
    !! is done only on the species part. That way, for example,
    !! directing the code to change the log K of 'albite' means that
    !! the code makes the change for 'albite'  in 'albite,' 'albite' in
    !! 'plagioclase,' and 'albite' in (Na,K)-feldspar.'
    !! This subroutine is called by:
    !!   EQ3NR/eq3nr.f
    !!   EQ6/eq6.f
    !! Principal input:
    !!   axlks  = array of polynomial coefficients for computing log K
    !!              values (altered by the options carried out by this
    !!              subroutine)
    !!   kxmod  = array of flag variables defining the type of nxmod
    !!              options
    !!   nxmod  = the number of 'nxmod' alter/suppress options
    !!   uxmod  = array of names of species whose thermodynamic data
    !!              are to be altered
    !!   xlkmod = array of values used by the 'nxmod' options to replace
    !!              or augment the the log K values
    !! Principal output:
    !!   axlks  = array of polynomial coefficients for computing log K
    !!              values
    implicit none

    ! Calling sequence variable declarations.
    integer :: narxmx
    integer :: ndrsmx
    integer :: nptmax
    integer :: nstmax
    integer :: ntprmx
    integer :: nxmdmx

    integer :: noutpt
    integer :: nttyo

    integer :: kxmod(nxmdmx)
    integer :: narxt(ntprmx)
    integer :: ndrs(ndrsmx)
    integer :: ndrsr(2,nstmax)
    integer :: npt
    integer :: nst
    integer :: ntpr
    integer :: nxmod

    character(len=48) :: uxmod(nxmdmx)
    character(len=48) :: uspec(nstmax)
    character(len=24) :: uphase(nptmax)

    real(kind=8) :: axlks(narxmx,ntprmx,nstmax)
    real(kind=8) :: apresg(narxmx,ntprmx)
    real(kind=8) :: cdrs(ndrsmx)
    real(kind=8) :: xlkmod(nxmdmx)
    real(kind=8) :: afcnst
    real(kind=8) :: tempc

    ! Local variable declarations.
    integer :: j
    integer :: jlen
    integer :: j2
    integer :: n
    integer :: nchar
    integer :: nerr
    integer :: nhits
    integer :: nn
    integer :: np
    integer :: nr1
    integer :: nr2
    integer :: ns
    integer :: nt

    integer :: ilnobl

    character(len=56) :: uspn56
    character(len=48) :: unam48
    character(len=24) :: ublk24
    character(len=8) :: ufix

    real(kind=8) :: cds
    real(kind=8) :: pgrid
    real(kind=8) :: xlkdel
    real(kind=8) :: xlknew
    real(kind=8) :: xlkold

    data ublk24 /'                        '/
    data ufix   /'fix_    '/

    nerr = 0

    if (nxmod .le. 0) then
        go to 999
    end if

    ! Compute the grid pressure.
    ! Calling sequence substitutions:
    !   apresg for arr
    !   pgrid for prop
    call evdat2(apresg,narxmx,narxt,ntpr,ntprmx,pgrid,tempc)

    do n = 1,nxmod
        if (kxmod(n) .lt. 0) then
            go to 130
        end if

        unam48 = uxmod(n)
        xlkdel = xlkmod(n)

        nchar = 24

        if (unam48(25:48) .ne. ublk24(1:24)) then
            nchar = 48
        end if

        nhits = 0

        do ns = 1,nst
            if (unam48(1:nchar) .ne. uspec(ns)(1:nchar)) then
                go to 120
            end if

            nhits = nhits + 1
            nr1 = ndrsr(1,ns)
            nr2 = ndrsr(2,ns)
            nt = nr2 - nr1 + 1

            if (nt .lt. 2) then
                call fmspnm(jlen,unam48,uspn56)
                write (noutpt,1005) uspn56(1:jlen)
                write (nttyo,1005) uspn56(1:jlen)
1005 format(/' * Error - (EQLIB/alters) The species ',a,/7x,'is in the strict basis set, so it can not be affected',/7x,'by the  specifed nxmod alter option.')

                nerr = nerr + 1

                if (nchar .eq. 48) then
                    go to 130
                end if

                go to 120
            end if

            ! Calling sequence substitutions:
            !   axlks for arr
            !   ns for k
            !   nstmax for nmax
            !   xlkold for prop
            call evdat3(axlks,ns,nstmax,narxmx,narxt,ntpr,ntprmx,xlkold,tempc)

            ! Map the kxmod = 0 and kxmod = 2 options to the kxmod = 1
            ! option.
            if (kxmod(n) .eq. 0) then
                kxmod(n) = 1
                xlkdel = xlkdel - xlkold
                xlkmod(n) = xlkdel
            else if (kxmod(n) .eq. 2) then
                do nn = nr1,nr2
                    if (ndrs(nn) .eq. ns) then
                        cds = cdrs(nn)
                        go to 110
                    end if
                end do

110 continue
                if (cds .le. 0) then
                    call fmspnm(jlen,unam48,uspn56)
                    write (noutpt,1010) uspn56(1:jlen)
                    write (nttyo,1010) uspn56(1:jlen)
1010 format(/" * Error - (EQLIB/alters) Couldn't find a",' non-zero reaction coefficient',/7x,'for the species ',a,', which is',/7x,'specified in an nxmod alter option.'," Therefore can't alter the",/7x,'corresponding',' equilibrium constant.')

                    nerr = nerr + 1

                    if (nchar .eq. 48) then
                        go to 130
                    end if

                    go to 120
                end if

                xlkdel = -cds*xlkdel/afcnst
                kxmod(n) = 1
                xlkmod(n) = xlkdel
            end if

            ! Change the polynomial coefficients.
            do j = 1,ntprmx
                axlks(1,j,ns) = axlks(1,j,ns) + xlkdel
            end do

            ! Calling sequence substitutions:
            !   noutpt for nf
            call prreac(cdrs,ndrs,ndrsmx,ndrsr,noutpt,ns,nstmax,uspec)

            ! Calling sequence substitutions:
            !   nttyo for nf
            call prreac(cdrs,ndrs,ndrsmx,ndrsr,nttyo,ns,nstmax,uspec)
            xlknew = xlkold + xlkdel
            write (noutpt,1020) tempc,pgrid,xlkold,xlknew
            write (nttyo,1020) tempc,pgrid,xlkold,xlknew
1020 format(/7x,'The log K of the above reaction at ',f10.3,'C',/7x,'and ',f10.3,' bars was changed from',f10.4,' to ',f10.4,'.',/)

            if (nchar .eq. 48) then
                go to 130
            end if

120 continue
        end do

        if (nhits .le. 0) then
            call fmspnm(jlen,unam48,uspn56)
            write (noutpt,1025) uspn56(1:jlen)
            write (nttyo,1025) uspn56(1:jlen)
1025 format(/" * Error- (EQLIB/alters) Can't find the species",/7x,a,', which is specified in an nxmod alter option.',' Check to make sure',/7x,'that a species of this name',' appears on the supporting data file.')

            j2 = ilnobl(unam48(1:24))

            do np = 1,npt
                if (unam48(1:24) .eq. uphase(np)(1:24)) then
                    write (noutpt,1030) unam48(1:j2)
                    write (nttyo,1030) unam48(1:j2)
1030 format(/' * Note - (EQLIB/alters) The entity "',a,'"',' specified',/7x,'in an nxmod alter option is a phase,',' not a species. Such an option applies only to species.')

                    go to 130
                end if
            end do
        end if

        ! Check to see that the species is not a fictive fugacity
        ! fixing species.
        if (unam48(1:4) .eq. ufix(1:4)) then
            call fmspnm(jlen,unam48,uspn56)
            write (noutpt,1035) uspn56(1:jlen)
            write (nttyo,1035) uspn56(1:jlen)
1035 format(/' * Note - (EQLIB/alters) The species "',a,'" specified',/7x,'in an nxmod alter option comprises a',' fictive fugacity-fixing.',/7x,' phase. Such an option'," can't be applied to this kind of species.")
        end if

130 continue
    end do

    if (nerr .gt. 0) then
        stop
    end if

999 continue
end subroutine alters