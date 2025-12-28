subroutine srch33(jtripl,unam1,unam2,unam3,utripl,npx3mx,npx3t)
    !! This subroutine searches for the species triplet corresponding
    !! to unam1, unam2, unam3 in the utripl array, and returns the
    !! index jtripl of that triplet in that array. The variable jtripl
    !! is returned as 0 if no match is found. The species names in
    !! the triplet array must be in the same order to generate a match.
    !! This subroutine is called by:
    !!   EQPT/wrpz3.f
    !! Principal input:
    !!   npx3t  = the number of aqueous species triplets in the utripl
    !!              array
    !!   unam1  = the name of the first aqueous species in a triplet
    !!   unam2  = the name of the second aqueous species in a triplet
    !!   unam3  = the name of the third aqueous species in a triplet
    !!   utripl = array of aqueous species triplets
    !! Principal output:
    !!   jtripl = the index of the triplet unam1, unam2, unam3 in the
    !!              utripl array, if that triplet is present; else 0
    implicit none

    ! Calling sequence variable declarations.
    integer :: npx3mx

    integer :: jtripl
    integer :: npx3t

    character(len=24) :: unam1
    character(len=24) :: unam2
    character(len=24) :: unam3
    character(len=24) :: utripl(3,npx3mx)

    ! Local variable declarations.
    integer :: j

    do j = 1,npx3t
        if (unam1(1:24) .eq. utripl(1,j)(1:24)) then
            if (unam2(1:24) .eq. utripl(2,j)(1:24)) then
                if (unam3(1:24) .eq. utripl(3,j)(1:24)) then
                    ! Found in the triplet in the utripl array.
                    jtripl = j
                    go to 999
                end if
            end if
        end if
    end do

    ! Did not find the triplet in the utripl array.
    jtripl = 0

999 continue
end subroutine srch33