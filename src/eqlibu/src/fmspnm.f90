subroutine fmspnm(jlen,unam48,uspn56)
    !! This subroutine formats a 48-character species name (unam48) into
    !! a 56-character string (uspn56) so that the phase part of the
    !! name (in the field composed of the second 24 characters) appears
    !! in parentheses following the actual species part (in the field
    !! composed of the first 24 characters). For example, ignoring
    !! trailing blanks, "Albite                  Plagioclase     " is
    !! reformatted as "Albite (Plagioclase)". However, the phase
    !! part (including the surrounding parentheses) is omitted from
    !! the formatted string if the phase name is not present in the
    !! 48-character string or if the phase name is identical to the
    !! species name. Thus, the returned string for pure Albite is
    !! 'Albite', not 'Albite (Albite)'.
    !! The formatted string is used mostly in error, warning, and note
    !! messages written by EQ3NR and EQ6.
    !! This subroutine is similar to fmspnx.f. However, that subroutine
    !! omits the phase part if the phase name is 'Aqueous solution'.
    !! This subroutine is called by:
    !!   Any
    !! Principal input:
    !!   unam48 = the species name (unformatted)
    !! Principal output:
    !!   jlen   = the length of the formatted species name
    !!   uspn56 = the formatted species name
    implicit none

    ! Calling sequence variable declarations.
    integer :: jlen

    character(len=48) :: unam48
    character(len=56) :: uspn56

    ! Local variable declarations.
    integer :: j2

    integer :: ilnobl

    character(len=24) :: uphnam

    uspn56 = unam48(1:24)
    jlen = ilnobl(uspn56)

    uphnam = unam48(25:48)
    j2 = ilnobl(uphnam)

    if (uphnam(1:24) .ne. unam48(1:24)) then
        if (j2 .gt. 0) then
            uspn56(jlen + 1:jlen + 2) = ' ('
            uspn56(jlen + 3:jlen + 2 + j2) = uphnam(1:j2)
            jlen = jlen + 3 + j2
            uspn56(jlen:jlen) = ')'
        end if
    end if
end subroutine fmspnm