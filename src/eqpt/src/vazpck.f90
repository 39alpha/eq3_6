subroutine vazpck(nat,natmax,nazt,naztmx,noutpt,nttyo,uaqsp,uazp)
    !! Validate the aqueous species names used to specify hard core
    !! diameters. Write a note if such a name does not appear on the
    !! main list of aqueous species.
    !! This subroutine is called by:
    !!   EQPT/eqpt.f
    !! Principal input:
    !!   nat    = the number of aqueous species
    !!   nazt   = the number of specified hard core diameters
    !!   uaqsp  = array of names of aqueous species
    !!   uazp   = array of aqueous species names used to specify
    !!              hard core diamters on the data file
    !! Principal output:
    !!   None
    implicit none

    ! Calling sequence variable declarations.
    integer :: natmax
    integer :: naztmx

    integer :: noutpt
    integer :: nttyo

    integer :: nat
    integer :: nazt

    character(len=24) :: uaqsp(natmax)
    character(len=24) :: uazp(naztmx)

    ! Local variable declarations.
    character(len=24) :: unam

    integer :: j2
    integer :: na
    integer :: naz

    integer :: ilnobl

    ! Loop over all aqueous species names used to specify hard core
    ! diameters.
    do naz = 1,nazt
        unam = uazp(naz)

        if (unam(1:7) .ne. '<blank>') then
            ! Check that each such name appears on the main list of
            ! aqueous species.
            do na = 1,nat
                if (unam(1:24) .eq. uaqsp(na)(1:24)) then
                    go to 100
                end if
            end do

            ! Did not find this name on the main list of aqueous species.
            j2 = ilnobl(unam)
            write (noutpt,1000) unam(1:j2)
            write (nttyo,1000) unam(1:j2)
1000 format(/' * Note - (EQPT/vazpck) A hard core diameter is',' specified on the',/7x,'data file for ',a,', but there is',' no species block for this',/7x,'in the aqueous species',' superblock. Therefore, this species will',/7x,'not appear',' in any computed model and the hard core diameter value',/7x,'will not be used. It is suggested to remove the',' hard core diameter',/7x,'in question or to add the',' corresponding species block.')

100 continue
        end if
    end do
end subroutine vazpck