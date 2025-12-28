subroutine cnvndx(nlim1,nlim1a,nlim2,nlim2a,nsmap,nstmax,ntot)
    !! This subroutine finds the first and last species in a reduced
    !! range corresponding to an original set defined by the species
    !! index limits nlim1a, nlim2a. Note that this is not a
    !! straightforward mapping of such as "nlim1 = nsmap(nlim1a),"
    !! because the species whose original index is nlim1a may have been
    !! eliminated in the reduction.
    !! This subroutine is called by:
    !!   EQLIB/cmpdat.f
    !! Principal input:
    !!   nlim1a = starting index of the range for the original set
    !!   nlim2a = ending index of the range for the original set
    !!   nsmap  = pointer array mapping species indices from the
    !!              original set to the reduced set
    !! Principal output:
    !!   nlim1  = starting index of the range for the reduced set
    !!   nlim2  = ending index of the range for the reduced set
    !!   ntot   = number of species in the range for the reduced set
    implicit none

    ! Calling sequence variable declarations.
    integer :: nstmax

    integer :: nsmap(nstmax)
    integer :: nlim1
    integer :: nlim1a
    integer :: nlim2
    integer :: nlim2a
    integer :: ntot

    ! Local variable declarations.
    integer :: j

    ! Search for the first species that is in the list.
    do j = nlim1a,nlim2a,1
        if (nsmap(j) .ne. 0) then
            go to 120
        end if
    end do

    ! None found, error.
    nlim1 = 0
    nlim2 = -1
    ntot = 0
    go to 999

120 continue

    ! Save the first species.
    nlim1 = nsmap(j)

    ! Search for last species that is in the list.
    do j = nlim2a,nlim1a,-1
        if (nsmap(j) .ne. 0) then
            go to 140
        end if
    end do

140 continue

    ! Save the last species.
    nlim2 = nsmap(j)
    ntot = nlim2 - nlim1 + 1

999 continue
end subroutine cnvndx