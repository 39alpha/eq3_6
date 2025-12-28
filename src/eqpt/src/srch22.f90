subroutine srch22(jpair,unam1,unam2,upair,npx2mx,npx2t)
    !! This subroutine searches for the species pair corresponding
    !! to unam1, unam2 or unam2, unam1 in the upair array, and returns
    !! the index jpair of that pair in that aray. The variable jpair
    !! is returned as 0 if no match is found.
    !! This suboutine is called by:
    !!   EQPT/tprca.f
    !!   EQPT/tprn2.f
    !!   EQPT/tprnn.f
    !!   EQPT/tprnc.f
    !!   EQPT/tprna.f
    !! Principal input:
    !!   npx2t  = number of species pairs in the upair array
    !!   unam1  = the name of the first aqueous species in a pair
    !!   unam2  = the name of the second aqueous species in a pair
    !!   upair  = array containing the names of species in a pair
    !! Principal output:
    !!   jpair  = the index of the pair unam1, unam2 or unam2, unam1 in
    !!              the upair array, if that pair is present; else 0
    implicit none

    ! Calling sequence variable declarations.
    integer :: npx2mx

    integer :: jpair
    integer :: npx2t

    character(len=24) :: unam1
    character(len=24) :: unam2
    character(len=24) :: upair(2,npx2mx)

    ! Local variable declarations.
    integer :: j

    ! Search the upair array for an element matching (unam1, unam2).
    do j = 1,npx2t
        if (unam1(1:24) .eq. upair(1,j)(1:24)) then
            ! Look for unam2 in the matching position in the upair array.
            if (unam2(1:24) .eq. upair(2,j)(1:24)) then
                jpair = j
                go to 999
            end if
        end if
    end do

    ! Search the upair array for an element matching (unam2, unam1).
    do j = 1,npx2t
        if (unam2(1:24) .eq. upair(1,j)(1:24)) then
            ! Look for unam1 in the matching position in the upair array.
            if (unam1(1:24) .eq. upair(2,j)(1:24)) then
                jpair = j
                go to 999
            end if
        end if
    end do

    ! Did not find the pair in the upair array.
    jpair = 0

999 continue
end subroutine srch22