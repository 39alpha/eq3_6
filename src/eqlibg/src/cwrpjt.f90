subroutine cwrpjt(noutpt)
    !! This subroutine calculates and writes tables of the higher-order
    !! electrostatic term functions J(x) and J'(x). These functions may
    !! be computed from two different approximate formulations. The
    !! first is one given by Pitzer (1975, eq. 47). The second is given
    !! by Harvie (1981). Almost all modern application of the Pitzer
    !! equations for aqueous electrolyte thermodynamics uses the more
    !! recent Harvie formulation. These formulations are all
    !! approximate in nature.
    !! It is important to note that the tabulation in Table II of
    !! Pitzer (1975) is not for the formulation represented by
    !! eq. 47. Therefore, the results given here for the eq. 47
    !! approximation will not precisely match the results given in
    !! that table. Pitzer (1975) does not provide such a table for
    !! the eq. 47 formulation.
    !!                        References
    !! Harvie, C. E. 1981. Theoretical Investigations in Geochemistry
    !!   and Atom Surface Scattering. Ph.D. dissertation, University
    !!   of California, San Diego (Available as #8203026 from
    !!   University Microfilms International, Ann Arbor, Michigan).
    !! Pitzer, K.S. 1975. Thermodynamics of electrolytes. V. Effects
    !!   of higher-order electrostatic terms. Journal of Solution
    !!   Chemistry, v. 4, p. 249-265.
    !! This subroutine is called by:
    !!   EQ6/eq6.f
    !!   EQ3NR/eq3nr.f
    !! Principal input:
    !!   noutpt = unit number of the output file
    !! Principal output:
    !!   None returned
    implicit none

    ! Calling sequence variable declarations.
    integer :: noutpt

    ! Local variable declarations.
    integer :: ktable
    integer :: n

    logical :: qpit75

    real(kind=8) :: dhj0
    real(kind=8) :: d2hj0
    real(kind=8) :: hj0
    real(kind=8) :: hj1
    real(kind=8) :: hj2
    real(kind=8) :: x

    write (noutpt,1000)
1000 format(/1x,"Tables of J(x) and J'(x) (used with Pitzer's",' equations)')

    ! Tables using the Pitzer (1975) table format.
    write (noutpt,1010)
1010 format(//3x,'Tables in the format of Pitzer (1975, Table II)',/3x,'(Note: that table was for a different approximation than',/3x,'than Pitzer, 1975, eq. 47, so results for eq. 47 will',/3x,'not precisely match those in Table II)')

    ktable = 1
    qpit75 = .true.
100 continue

    if (qpit75) then
        write (noutpt,1020)
1020 format(//5x,'Pitzer (1975, eq. 47) formulation')
    else
        write (noutpt,1030)
1030 format(//5x,'Harvie (1981, Appendix B) formulation')
    end if

    write (noutpt,1040)
1040 format(/7x,"x           J(x)        J'(x)",/)

    x = 0.

    do n = 1,10
        x = x + 0.01

        if (.not.qpit75) then
            call ghj0(dhj0,d2hj0,hj0,hj1,hj2,x)
        else
            call gpj0(dhj0,d2hj0,hj0,x)
        end if

        write (noutpt,1050) x,hj0,dhj0
1050 format(1x,f9.3,3x,f13.8,3x,f7.4)
    end do

    ! x should now be 0.1
    do n = 1,5
        x = x + 0.02

        if (.not.qpit75) then
            call ghj0(dhj0,d2hj0,hj0,hj1,hj2,x)
        else
            call gpj0(dhj0,d2hj0,hj0,x)
        end if

        write (noutpt,1050) x,hj0,dhj0
    end do

    ! x should now be 0.2
    do n = 1,10
        x = x + 0.04

        if (.not.qpit75) then
            call ghj0(dhj0,d2hj0,hj0,hj1,hj2,x)
        else
            call gpj0(dhj0,d2hj0,hj0,x)
        end if

        write (noutpt,1050) x,hj0,dhj0
    end do

    ! x should now be 0.6
    do n = 1,7
        x = x + 0.2

        if (.not.qpit75) then
            call ghj0(dhj0,d2hj0,hj0,hj1,hj2,x)
        else
            call gpj0(dhj0,d2hj0,hj0,x)
        end if

        write (noutpt,1050) x,hj0,dhj0
    end do

    ! x should now be 2.0
    do n = 1,8
        x = x + 1.0

        if (.not.qpit75) then
            call ghj0(dhj0,d2hj0,hj0,hj1,hj2,x)
        else
            call gpj0(dhj0,d2hj0,hj0,x)
        end if

        write (noutpt,1050) x,hj0,dhj0
    end do

    ! x should now be 10.0
    x = x + 2.0

    if (.not.qpit75) then
        call ghj0(dhj0,d2hj0,hj0,hj1,hj2,x)
    else
        call gpj0(dhj0,d2hj0,hj0,x)
    end if

    write (noutpt,1050) x,hj0,dhj0

    ! x should now be 12.0
    do n = 1,7
        x = x + 4.0

        if (.not.qpit75) then
            call ghj0(dhj0,d2hj0,hj0,hj1,hj2,x)
        else
            call gpj0(dhj0,d2hj0,hj0,x)
        end if

        write (noutpt,1050) x,hj0,dhj0
    end do

    ! x should now be 40.0
    do n = 1,6
        x = x + 10.0

        if (.not.qpit75) then
            call ghj0(dhj0,d2hj0,hj0,hj1,hj2,x)
        else
            call gpj0(dhj0,d2hj0,hj0,x)
        end if

        write (noutpt,1050) x,hj0,dhj0
    end do

    ! x should now be 100.0
    x = x + 100.0

    if (.not.qpit75) then
        call ghj0(dhj0,d2hj0,hj0,hj1,hj2,x)
    else
        call gpj0(dhj0,d2hj0,hj0,x)
    end if

    write (noutpt,1050) x,hj0,dhj0

    ! x should now be 200.0
    do n = 1,4
        x = x + 200.0

        if (.not.qpit75) then
            call ghj0(dhj0,d2hj0,hj0,hj1,hj2,x)
        else
            call gpj0(dhj0,d2hj0,hj0,x)
        end if

        write (noutpt,1050) x,hj0,dhj0
    end do

    ! x should now be 1000.0
    x = x + 1000.0

    if (.not.qpit75) then
        call ghj0(dhj0,d2hj0,hj0,hj1,hj2,x)
    else
        call gpj0(dhj0,d2hj0,hj0,x)
    end if

    write (noutpt,1050) x,hj0,dhj0

    ! x should now be 2000.0
    do n = 1,4
        x = x + 2000.0

        if (.not.qpit75) then
            call ghj0(dhj0,d2hj0,hj0,hj1,hj2,x)
        else
            call gpj0(dhj0,d2hj0,hj0,x)
        end if

        write (noutpt,1050) x,hj0,dhj0
    end do

    ! x should now be 10000.0
    if (ktable .le. 1) then
        qpit75 = .false.
        ktable = ktable + 1
        go to 100
    end if

    ! Tables using the Harvie (1981) table format (extended).
    write (noutpt,1060)
1060 format(//3x,'Tables in the format of Harvie (1981, Appendix B)',/3x,'(extended format)')

    ktable = 1
    qpit75 = .true.
120 continue

    if (qpit75) then
        write (noutpt,1020)
    else
        write (noutpt,1030)
    end if

    write (noutpt,1040)

    x = 0.

    do n = 1,10
        x = x + 0.001

        if (.not.qpit75) then
            call ghj0(dhj0,d2hj0,hj0,hj1,hj2,x)
        else
            call gpj0(dhj0,d2hj0,hj0,x)
        end if

        write (noutpt,1070) x,hj0,dhj0
1070 format(1x,f9.3,3x,f12.7,4x,f10.7)
    end do

    ! x should now be 0.01
    do n = 1,9
        x = x + 0.01

        if (.not.qpit75) then
            call ghj0(dhj0,d2hj0,hj0,hj1,hj2,x)
        else
            call gpj0(dhj0,d2hj0,hj0,x)
        end if

        write (noutpt,1070) x,hj0,dhj0
    end do

    ! x should now be 0.1
    do n = 1,9
        x = x + 0.1

        if (.not.qpit75) then
            call ghj0(dhj0,d2hj0,hj0,hj1,hj2,x)
        else
            call gpj0(dhj0,d2hj0,hj0,x)
        end if

        write (noutpt,1070) x,hj0,dhj0
    end do

    ! x should now be 1.0
    do n = 1,9
        x = x + 1.0

        if (.not.qpit75) then
            call ghj0(dhj0,d2hj0,hj0,hj1,hj2,x)
        else
            call gpj0(dhj0,d2hj0,hj0,x)
        end if

        write (noutpt,1070) x,hj0,dhj0
    end do

    ! x should now be 10.0
    do n = 1,9
        x = x + 10.0

        if (.not.qpit75) then
            call ghj0(dhj0,d2hj0,hj0,hj1,hj2,x)
        else
            call gpj0(dhj0,d2hj0,hj0,x)
        end if

        write (noutpt,1070) x,hj0,dhj0
    end do

    ! x should now be 100.0
    do n = 1,9
        x = x + 100.0

        if (.not.qpit75) then
            call ghj0(dhj0,d2hj0,hj0,hj1,hj2,x)
        else
            call gpj0(dhj0,d2hj0,hj0,x)
        end if

        write (noutpt,1070) x,hj0,dhj0
    end do

    ! x should now be 1000.0
    do n = 1,9
        x = x + 1000.0

        if (.not.qpit75) then
            call ghj0(dhj0,d2hj0,hj0,hj1,hj2,x)
        else
            call gpj0(dhj0,d2hj0,hj0,x)
        end if

        write (noutpt,1070) x,hj0,dhj0
    end do

    ! x should now be 10000.0
    if (ktable .le. 1) then
        qpit75 = .false.
        ktable = ktable + 1
        go to 120
    end if

    write (noutpt,2000)
2000 format(/' - - - - - - - - - - - - - - - - - - - - - - - - ','- - - - - - ',/)
end subroutine cwrpjt