subroutine cpcomb(jassan,jassca,jassne,naapr,nat,ncapr,nccpr,nnapr,nncpr,nnnpr,nn2pr,naactr,na2ctr,nncatr,nccatr,nc2atr,nn2ntr,nn3tr)
    !! Compute all the pair and triplet combinations of aqueous
    !! solute species for use with Pitzer's equations. parameters.
    !! All pairs and triplets are distinct in the sense that the
    !! ordering of the members is not relevant; e.g., ca = ac,
    !! cc' = c'c, nca = anc = can = cna = acn = nac.
    !! This subroutine is called by:
    !!   EQPT/eqpt.f
    !! Principal input:
    !!   jassca = number of aqueous cation species
    !!   jassan = number of aqueous anion species (excluding aqueous e-)
    !!   jassne = number of aqueous neutral species (excluding solvent
    !!              water and aqeuous O2(g))
    !!   nat    = the number of aqueous species
    !! Principal output:
    !!   naapr  = number of aa' pairs
    !!   ncapr  = number of ca pairs
    !!   nccpr  = number of cc' pairs
    !!   nnapr  = number of na pairs
    !!   nncpr  = number of nc pairs
    !!   nnnpr  = number of nn' pairs
    !!   nn2pr  = number of nn pairs
    !!   naactr = number of aa'c triplets
    !!   na2ctr = number of aac triplets
    !!   nccatr = number of cc'a triplets
    !!   nc2atr = number of cca triplets
    !!   nncatr = number of nca triplets
    !!   nn3tr  = number of nnn triplets
    !!   nn2ntr = number of nnn' triplets
    implicit none

    ! Calling sequence variable declarations.
    integer :: jassan
    integer :: jassca
    integer :: jassne
    integer :: nat

    integer :: naapr
    integer :: ncapr
    integer :: nccpr
    integer :: nnapr
    integer :: nncpr
    integer :: nnnpr
    integer :: nn2pr
    integer :: naactr
    integer :: na2ctr
    integer :: nccatr
    integer :: nc2atr
    integer :: nncatr
    integer :: nn2ntr
    integer :: nn3tr

    ! Local variable declarations.
    integer :: jn
    integer :: jc
    integer :: ja

    jn = jassne
    jc = jassca
    ja = jassan

    ! Pair combinations.
    ! Get the number of ca pairs.
    ncapr = jc*ja

    ! Get the number of cc' pairs.
    nccpr = jc*(jc -1)/2

    ! Get the number of aa' pairs.
    naapr = ja*(ja -1)/2

    ! Get the number of nn pairs.
    nn2pr = jn

    ! Get the number of nn' pairs.
    nnnpr = jn*(jn -1)/2

    ! Get the number of nc pairs.
    nncpr = jn*jc

    ! Get the number of na pairs.
    nnapr = jn*ja

    ! Triplet combinations.
    ! Get the number of cca triplets.
    nc2atr = ncapr

    ! Get the number of aac triplets.
    na2ctr = ncapr

    ! Get the number of nnn triplets.
    nn3tr = nn2pr

    ! Get the number of nnn' triplets.
    nn2ntr = 2*nn2pr

    ! Get the number of cc'a triplets.
    nccatr = nccpr*ja

    ! Get the number of aa'c triplets.
    naactr = naapr*jc

    ! Get the number of nca triplets.
    nncatr = jn*jc*ja
end subroutine cpcomb