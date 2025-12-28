subroutine coasst(jassan,jassca,jassne,nat,natmax,uaqsp,zaqsp)
    !! Get the numbers of aqueous solute cations, anions, and neutral
    !! species. These will be used to generate lists of species pairs
    !! and triplets for use with Pitzer's equations.
    !! The fictive redox species aqueous e- and aqueous O2(g), if
    !! present, are excluded in computing these totals. So is solvent
    !! water.
    !! This subroutine is called by:
    !!   EQPT/eqpt.f
    !! Principal input:
    !!   nat    = the number of aqueous species
    !!   uaqsp  = array of names of aqueous species
    !!   zaqsp  = array of electrical charages of aqueous species
    !! Principal output:
    !!   jassca = number of aqueous cation species
    !!   jassan = number of aqueous anion species (excluding aqueous e-)
    !!   jassne = number of aqueous neutral species (excluding solvent
    !!              water and aqeuous O2(g))
    implicit none

    ! Calling sequence variable declarations.
    integer :: natmax

    integer :: jassan
    integer :: jassca
    integer :: jassne
    integer :: nat

    character(len=24) :: uaqsp(natmax)

    real(kind=8) :: zaqsp(natmax)

    ! Local variable declarations.
    integer :: n

    ! Count the aqueous solute cations, anions, and neutrals. Include
    ! in the counts only real solute species.
    jassne = 0
    jassca = 0
    jassan = 0

    do n = 2,nat
        if (zaqsp(n) .eq. 0.) then
            if (uaqsp(n)(1:6) .ne. 'O2(g) ') then
                jassne = jassne + 1
            end if
        else if (zaqsp(n) .gt. 0.) then
            jassca = jassca + 1
        else
            if (uaqsp(n)(1:3) .ne. 'e- ') then
                jassan = jassan + 1
            end if
        end if
    end do
end subroutine coasst