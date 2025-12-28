subroutine artrip(iz1,iz2,iz3,na,nc,nn,n1,n2,n3,u1,u2,u3,z1,z2,z3)
    !! This suboutine arranges the members of a species triplet
    !! according to the following rules:
    !!   If exactly one neutral is present:
    !!     neutral, cation, anion
    !!   If exactly two neutrals are present:
    !!     neutral 1, neutral 1, neutral2
    !!   If two cations are present:
    !!     cation 1, cation 2, anion
    !!   If two anions are present:
    !!     anion 1, anion 2, cation
    !! This suboutine is called by:
    !!   EQPT/rdpz3.f
    !! Principal input:
    !!   iz1    = (integer) charge of the first species in a triplet
    !!   iz2    = (integer) charge of the second species in a triplet
    !!   iz3    = (integer) charge of the third species in a triplet
    !!   na     = number of anions in a species triplet
    !!   nc     = number of cations in a species triplet
    !!   nn     = number of neutral species in a species triplet
    !!   u1     = first name in a species triplet
    !!   u2     = second name in a species triplet
    !!   u3     = third name in a species triplet
    !!   z1     = charge of the first species in a triplet
    !!   z2     = charge of the second species in a triplet
    !!   z3     = charge of the third species in a triplet
    !! Principal output:
    !!   none     (u1, u2, u3 are rearranged)
    implicit none

    ! Calling sequence variable declarations.
    integer :: iz1
    integer :: iz2
    integer :: iz3
    integer :: na
    integer :: nc
    integer :: nn
    integer :: n1
    integer :: n2
    integer :: n3

    character(len=24) :: u1
    character(len=24) :: u2
    character(len=24) :: u3

    real(kind=8) :: z1
    real(kind=8) :: z2
    real(kind=8) :: z3

    ! Local variable declarations.
    integer :: nu
    integer :: nx
    integer :: izx

    character(len=24) :: ux24

    real(kind=8) :: zx

    if (nn .eq. 3) then
        if (n1.ne.n2 .or. n1.ne.n3) then
            ! Have the combination nnn'.
            ! Make sure that the unique neutral is in third place.
            ! Find the unique neutral.
            if (n2 .eq. n3) then
                nu = 1
            end if

            if (n1 .eq. n3) then
                nu = 2
            end if

            if (n1 .eq. n2) then
                nu = 3
            end if

            if (nu .eq. 1) then
                ! The unique neutral is the first species. Switch it
                ! with the third one.
                ux24 = u1
                u1 = u3
                u3 = ux24

                nx = n1
                n1 = n3
                n3 = nx

                zx = z1
                z1 = z3
                z3 = zx

                izx = iz1
                iz1 = iz3
                iz3 = izx
            else if (nu .eq. 2) then
                ! The unique neutral is the second species. Switch it
                ! with the third one.
                ux24 = u2
                u2 = u3
                u3 = ux24

                nx = n2
                n2 = n3
                n3 = nx

                zx = z2
                z2 = z3
                z3 = zx

                izx = iz2
                iz2 = iz3
                iz3 = izx
            end if
        end if
    end if

    if (nn.eq.1 .and. nc.eq.1 .and. na.eq.1) then
        ! Have the combination nca.
        ! Make sure that the neutral is in first place.
        if (iz2 .eq. 0) then
            ux24 = u1
            u1 = u2
            u2 = ux24

            nx = n1
            n1 = n2
            n2 = nx

            zx = z1
            z1 = z2
            z2 = zx

            izx = iz1
            iz1 = iz2
            iz2 = izx
        else if (iz3 .eq. 0) then
            ux24 = u1
            u1 = u3
            u3 = ux24

            nx = n1
            n1 = n3
            n3 = nx

            zx = z1
            z1 = z3
            z3 = zx

            izx = iz1
            iz1 = iz3
            iz3 = izx
        end if

        ! Make sure that the cation is in second place. It
        ! must presently be in second or third place.
        if (iz3 .gt. 0) then
            ux24 = u2
            u2 = u3
            u3 = ux24

            nx = n2
            n2 = n3
            n3 = nx

            zx = z2
            z2 = z3
            z3 = zx

            izx = iz2
            iz2 = iz3
            iz3 = izx
        end if
    end if

    if (nc.eq.2 .and. na.eq.1) then
        ! Have the combination cc'a.
        ! Make sure that the anion is in third place.
        if (iz1 .lt. 0) then
            ux24 = u1
            u1 = u3
            u3 = ux24

            nx = n1
            n1 = n3
            n3 = nx

            zx = z1
            z1 = z3
            z3 = zx

            izx = iz1
            iz1 = iz3
            iz3 = izx
        else if (iz2 .lt. 0) then
            ux24 = u2
            u2 = u3
            u3 = ux24

            nx = n2
            n2 = n3
            n3 = nx

            zx = z2
            z2 = z3
            z3 = zx

            izx = iz2
            iz2 = iz3
            iz3 = izx
        end if

        ! Put the two cations in alphabetical order.
        if (u1(1:24) .gt. u2(1:24)) then
            ux24 = u1
            u1 = u2
            u2 = ux24

            nx = n1
            n1 = n2
            n2 = nx

            zx = z1
            z1 = z2
            z2 = zx

            izx = iz1
            iz1 = iz2
            iz2 = izx
        end if
    end if

    if (na.eq.2 .and. nc.eq.1) then
        ! Have the combination aa'c.
        ! Make sure that the cation is in third place.
        if (iz1 .gt. 0) then
            ux24 = u1
            u1 = u3
            u3 = ux24

            nx = n1
            n1 = n3
            n3 = nx

            zx = z1
            z1 = z3
            z3 = zx

            izx = iz1
            iz1 = iz3
            iz3 = izx
        else if (iz2 .gt. 0) then
            ux24 = u2
            u2 = u3
            u3 = ux24

            nx = n2
            n2 = n3
            n3 = nx

            zx = z2
            z2 = z3
            z3 = zx

            izx = iz2
            iz2 = iz3
            iz3 = izx
        end if

        ! Put the two anions in alphabetical order.
        if (u1(1:24) .gt. u2(1:24)) then
            ux24 = u1
            u1 = u2
            u2 = ux24

            nx = n1
            n1 = n2
            n2 = nx

            zx = z1
            z1 = z2
            z2 = zx

            izx = iz1
            iz1 = iz2
            iz2 = izx
        end if
    end if
end subroutine artrip