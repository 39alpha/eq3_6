subroutine bspchk(jsflag,nbaspd,nbtd,nbtmax,ndrsd,ndrsmx,ndrsrd,noutpt,nrdxsp,nstmax,nttyo,uspeca)
    !! This subroutine looks at each active auxiliary basis species.
    !! It prints a warning if any other species in the corresponding
    !! dissociation reaction is not present.
    !! This subroutine is called by:
    !!   EQ3NR/eq3nr.f
    !!   EQ6/eq6.f
    !! Principal input:
    !!   jsflag = array of status flags for species
    !!   nbaspd = array of indices of 'data file' basis species
    !! Principal output:
    !!   None
    implicit none

    ! Calling sequence variable declarations.
    integer :: nbtmax
    integer :: ndrsmx
    integer :: nstmax

    integer :: noutpt
    integer :: nttyo

    integer :: jsflag(nstmax)
    integer :: nbaspd(nbtmax)
    integer :: ndrsd(ndrsmx)
    integer :: ndrsrd(2,nstmax)
    integer :: nbtd
    integer :: nrdxsp

    character(len=48) :: uspeca(nstmax)

    ! Local variable declarations.
    integer :: jlen
    integer :: jlene
    integer :: n
    integer :: nb
    integer :: nr1
    integer :: nr2
    integer :: ns
    integer :: nse
    integer :: nt1

    character(len=56) :: uspe56
    character(len=56) :: uspn56

    do nb = 1,nbtd
        ns = nbaspd(nb)

        if (jsflag(ns) .le. 0) then
            nr1 = ndrsrd(1,ns)
            nr2 = ndrsrd(2,ns)
            nt1 = nr2 - nr1 + 1

            if (nt1 .ge. 2) then
                do n = nr1,nr2
                    nse = ndrsd(n)

                    if (jsflag(nse) .gt. 0) then
                        if (nse .ne. nrdxsp) then
                            ! Calling sequence substitutions:
                            !   jlene for jlen
                            !   uspeca(nse) for unam48
                            !   uspe56 for uspn56
                            call fmspnx(jlene,uspeca(nse),uspe56)

                            ! Calling sequence substitutions:
                            !   uspeca(ns) for unam48
                            call fmspnx(jlen,uspeca(ns),uspn56)

                            write (noutpt,1000) uspn56(1:jlen),uspe56(1:jlene)
                            write (nttyo,1000) uspn56(1:jlen),uspe56(1:jlene)
1000 format(/' The auxiliary basis species ',a,' is active even though',/3x,'it is a dependent',' species of ',a,', which is not',/3x,'present',' in the current model. This detached auxililary',' basis species',/3x,'will behave much like a',' strict basis species.')
                        end if
                    end if
                end do
            end if
        end if
    end do
end subroutine bspchk
