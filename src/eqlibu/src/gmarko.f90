subroutine gmarko(nfldmx,nfldt,nmark,nttyo,ufield,uline1)
    !! This subroutine determines the choice for an option. This choice
    !! is marked with an asterisk on a line containing a field for
    !! each choice.
    !! This subroutine is called by:
    !!   Any
    !! Input:
    !!   nfldmx = the maximum number of fields on line uline1
    !!   nfldt  = the number of fields on line uline1
    !!   nttyo  = the unit number of the screen file
    !!   ufield = the array of field contents in uline1
    !!   uline1 = the line containing the option being examined
    !! Output:
    !!   nmark  = index of the first marked field on line uline1
    implicit none

    ! Calling sequence variable declarations.
    integer :: nfldmx

    integer :: nfldt
    integer :: nmark
    integer :: nttyo

    character(len=*) :: ufield(nfldmx)
    character(len=*) :: uline1

    ! Local variable declarations.
    integer :: j
    integer :: j2
    integer :: n
    integer :: ncount

    integer :: ilnobl

    ! Find the first option marked with an asterisk.
    ncount = 0
    nmark = 0

    do n = 2,nfldt
        j = index(ufield(n),'*')

        if (j .gt. 0) then
            ncount = ncount + 1

            if (ncount .eq. 1) then
                nmark = n
            end if
        end if
    end do

    j2 = ilnobl(uline1)
    j2 = min(j2,70)

    if (ncount .eq. 0) then
        write (nttyo,1030) uline1(1:j2)
1030 format(/' * Warning - (EQLIBU/gmarko) None of the options'  /7x,'is marked with an asterisk on the line beginning with',/7x,'"',a,'".')
    end if

    if (ncount .gt. 1) then
        write (nttyo,1040) uline1(1:j2)
1040 format(/' * Warning - (EQLIBU/gmarko) More than one of the'  /7x,'options is marked with an asterisk on the line',' beginning with',/7x,'"',a,'".')
    end if
end subroutine gmarko