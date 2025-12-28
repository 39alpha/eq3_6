subroutine dfbasp(nbaspa,nbta,nbta_asv,ndrsa,ndrsa_asv,ndrsra,nerr,noutpt,nsta,nsta_asv,nttyo,ubasp,uspeca)
    !! This subroutine decodes the ubasp array (list of basis species
    !! names built while reading the data file). It puts their
    !! species indices in the nbaspa array. The basis indices that
    !! were put in the ndrsa array are replaced by the corresponding
    !! species indices. The working basis set may be subsequently
    !! expanded when the input file is interpreted, as any species
    !! on the input file that is referenced as a basis species is
    !! added to that set if not already in it.
    !! This subroutine is called by:
    !!   EQLIB/indata.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: nbta_asv
    integer :: ndrsa_asv
    integer :: nsta_asv

    integer :: nbaspa(nbta_asv)
    integer :: ndrsa(ndrsa_asv)
    integer :: ndrsra(2,nsta_asv)
    integer :: nbta
    integer :: nerr
    integer :: noutpt
    integer :: nsta
    integer :: nttyo

    character(len=48) :: ubasp(nbta_asv)
    character(len=48) :: uspeca(nsta_asv)

    ! Local variable declarations.
    integer :: jlen
    integer :: n
    integer :: nb
    integer :: ns
    integer :: nr1
    integer :: nr1p1
    integer :: nr2
    integer :: nt

    character(len=56) :: uspn56

    ! Set up the nbaspa array and expand the elements of the ubasp array
    ! to the full 48 characters.
    do nb = 1,nbta
        do ns = 1,nsta
            if (uspeca(ns)(1:24) .eq. ubasp(nb)(1:24)) then
                nbaspa(nb) = ns
                ubasp(nb)(1:48) = uspeca(ns)(1:48)
                go to 105
            end if
        end do

        nbaspa(nb) = 0

        ! Calling sequence substitutions:
        !   ubasp(nb) for unam48
        call fmspnm(jlen,ubasp(nb),uspn56)
        write (noutpt,1000) uspn56(1:jlen)
        write (nttyo,1000) uspn56(1:jlen)
1000 format(/" * Error - (EQLIB/dfbasp) Couldn't find a basis",' species named',/7x,a,' among the species read from the',' data file.',/7x,'It must be referenced erroneously in',' the associated reaction',/7x,'for another species.')

        nerr = nerr + 1
105 continue
    end do

    ! Convert basis species indices in the ndrsa array to species
    ! indices.
    do ns = 1,nsta
        nt = ndrsra(2,ns) - ndrsra(1,ns) + 1

        if (nt .ge. 2) then
            nr1 = ndrsra(1,ns)
            nr1p1 = nr1 + 1
            nr2 = ndrsra(2,ns)

            do n = nr1p1,nr2
                nb = ndrsa(n)
                ndrsa(n) = nbaspa(nb)
            end do
        end if
    end do
end subroutine dfbasp