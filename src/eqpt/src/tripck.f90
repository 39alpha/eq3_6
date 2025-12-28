subroutine tripck(na,nc,nerr,nn,noutpt,nttyo,n1,n2,n3,qdup12,qdup13,qdup23,unam1,unam2,unam3)
    !! This suboutine checks the species triplets that were read
    !! from the DATA0 file for illegal combinations.
    !! The only legal combinations here correspond to ternary
    !! systems and maximal third order in the expansion of the
    !! Gibbs energy. These are:
    !!   nnn'      mu(nnn') maps to itself
    !!   nca       zeta(nca) maps to mu(nca)
    !!   cc'a      psi(cc'a) maps to mu(cc'a)
    !!   aa'c      psi(aa'c) maps to mu(aa'c)
    !! Here n = neutral, n' a different neutral, c = cation, c' = a
    !! different cation, a = anion, and a' = a different anion.
    !! Solvent water (w) may not appear in any combination in the
    !! normal Pitzer treatment of electrolyte solutions.
    !! Other possible combinations are not allowed because of one or
    !! more of the following reasons:
    !!   1. The combination corresponds to systems of higher order
    !!      (e.g., cc'c'')
    !!   2. The combination corresponds to systems of lesser order
    !!      (e.g., nnn and cca)
    !!   3. The combination corresponds to parameters that are defined
    !!      by convention to be zero (i.e., to unused parameters).
    !!      Examples include ccc and ncc.
    !! This suboutine is called by:
    !!   EQPT/rdpz3.f
    !! Principal input:
    !!   na     = number of anions in a species triplet
    !!   nc     = number of cations in a species triplet
    !!   nn     = number of neutral species in a species triplet
    !!   n1     = index of the first species in a triplet
    !!   n2     = index of the second species in a triplet
    !!   n3     = index of the third species in a triplet
    !!   unam1  = first name in a species triplet
    !!   unam2  = second name in a species triplet
    !!   unam3  = third name in a species triplet
    !! Principal output:
    !!   nerr   = error counter
    implicit none

    ! Calling sequence variable declarations.
    integer :: noutpt
    integer :: nttyo

    integer :: na
    integer :: nc
    integer :: nn
    integer :: nerr
    integer :: n1
    integer :: n2
    integer :: n3

    logical :: qdup12
    logical :: qdup13
    logical :: qdup23

    character(len=24) :: unam1
    character(len=24) :: unam2
    character(len=24) :: unam3

    ! Local variable declarations.
    integer :: j2
    integer :: j3
    integer :: j4

    integer :: ilnobl

    j2 = ilnobl(unam1)
    j3 = ilnobl(unam2)
    j4 = ilnobl(unam3)

    ! Check for illegal appearance of solvent water.
    if (unam1(1:4) .eq.'H2O ' .or. unam2(1:4) .eq.'H2O ' .or.  unam3(1:4).eq.'H2O ' .or. n1.eq.1 .or. n2.eq.1 .or.  n3.eq.1) then
        write (noutpt,1430) unam1(1:j2),unam2(1:j3),unam3(1:j4)
        write (nttyo,1430) unam1(1:j2),unam2(1:j3),unam3(1:j4)
1430 format(/' * Error - (EQPT/tripck) Have found an illegal data',' block for the',/7x,'species triplet ',a,', ',a,', ',a,' among the blocks for',/7x,'"mixture parameters".',' Solvent water may not appear in such a',/7x,'triplet in',' the normal Pitzer treatment of electrolyte',/7x,'solutions.')

        nerr = nerr + 1
    end if

    if (qdup12 .and. qdup13) then
        ! The same species appears three times.
        if (nc.gt.0 .or. na.gt.0) then
            ! The same cation or anion appears three times.(ccc or aaa).
            write (noutpt,1110) unam1(1:j2),unam2(1:j3),unam3(1:j4)
            write (nttyo,1110) unam1(1:j2),unam2(1:j3),unam3(1:j4)
1110 format(/' * Error - (EQPT/tripck) Have found an illegal data',' block for the',/7x,'species triplet ',a,', ',a,', ',a,' among the blocks',/7x,'for "mixture parameters". The same',' ion appears three times.',/7x,'This is not a valid',' combination, as the corrresponding lambda and mu',/7x,'data are defined to be zero by convention in the normal',' Pitzer',/7x,'treatment of electrolyte solutions.')
        end if

        if (nn .gt. 0) then
            ! The same neutral appears three times (nnn).
            write (noutpt,1112) unam1(1:j2),unam2(1:j3),unam3(1:j4)
            write (nttyo,1112) unam1(1:j2),unam2(1:j3),unam3(1:j4)
1112 format(/' * Error - (EQPT/tripck) Have found an illegal data',' block for the',/7x,'species triplet ',a,', ',a,', ',a,' among the blocks',/7x,'for "mixture parameters". The',' same neutral appears three times.',/7x,'This is not',' a valid combination for the present kind of block.',/7x,'Enter these data in a block for "single-salt',' parameters".')
        end if

        nerr = nerr + 1
        go to 999
    end if

    ! Check for other illegal combinations involving two or more
    ! neutrals.
    if (nn .eq. 3) then
        if (.not.(qdup12 .or. qdup13 .or. qdup23)) then
            ! Have three distinct neutrals (nn'n'').
            write (noutpt,1130) unam1(1:j2),unam2(1:j3),unam3(1:j4)
            write (nttyo,1130) unam1(1:j2),unam2(1:j3),unam3(1:j4)
1130 format(/' * Error - (EQPT/tripck) Have found an illegal data',' block for the',/7x,'species triplet ',a,', ',a,', ',a,' among the blocks',/7x,'for "mixture parameters". Three',' distinct neutrals appear. This is not',/7x,'a valid',' combination in the context of the normal Pitzer',/7x,'treatment of electrolyte solutions as it corresponds',' to a',/7x,'quaternary system.')

            nerr = nerr + 1
            go to 999
        end if
    end if

    if (nn .eq. 2) then
        if (nc.eq.1 .or. na.eq.1) then
            if (qdup12 .or. qdup23 .or. qdup13) then
                ! Have one of the following combinations: nnc or nna.
                write (noutpt,1132) unam1(1:j2),unam2(1:j3),unam3(1:j4)
                write (nttyo,1132) unam1(1:j2),unam2(1:j3),unam3(1:j4)
1132 format(/' * Error - (EQPT/tripck) Have found an illegal',' data block for the',/7x,'species triplet ',a,', ',a,', ',a,' among the blocks',/7x,'for "mixture parameters". One',' neutral appears twice with an ion.',/7x,'This is not a',' valid combination in the context of the normal Pitzer',/7x,'treatment of electrolyte solutions. The corresponding',' parameters',/7x,'are not used.')
            else
                ! Have one of the following combinations: nn'c or nn'a.
                write (noutpt,1136) unam1(1:j2),unam2(1:j3),unam3(1:j4)
                write (nttyo,1136) unam1(1:j2),unam2(1:j3),unam3(1:j4)
1136 format(/' * Error - (EQPT/tripck) Have found an illegal',' data block for the',/7x,'species triplet ',a,', ',a,', ',a,' among the blocks',/7x,'for "mixture parameters". Two',' neutrals appear with an ion.',/7x,'This is not a',' valid combination in the context of the normal Pitzer',/7x,'treatment of electrolyte solutions as it corresponds',' to a',/7x,'quaternary system.')
            end if

            nerr = nerr + 1
            go to 999
        end if
    end if

    ! Check for other illegal combinations involving two or more
    ! cations.
    if (nc .eq. 3) then
        if (.not.(qdup12 .or. qdup13 .or. qdup23)) then
            ! Have three distinct cations (cc'c'').
            write (noutpt,1140) unam1(1:j2),unam2(1:j3),unam3(1:j4)
            write (nttyo,1140) unam1(1:j2),unam2(1:j3),unam3(1:j4)
1140 format(/' * Error - (EQPT/tripck) Have found an illegal data',' block for the',/7x,'species triplet ',a,', ',a,', ',a,' among the blocks',/7x,'for "mixture parameters". Three',' distinct cations appear. This is',/7x,'not a valid',' combination in the context of the normal Pitzer',/7x,'treatment of electrolyte solutions as it corresponds',' to a',/7x,'quaternary system.')

            nerr = nerr + 1
            go to 999
        end if

        if (qdup12 .and. qdup13) then
            ! Have one cation appearing three times (ccc).
            write (noutpt,1142) unam1(1:j2),unam2(1:j3),unam3(1:j4)
            write (nttyo,1142) unam1(1:j2),unam2(1:j3),unam3(1:j4)
1142 format(/' * Error - (EQPT/tripck) Have found an illegal data',' block for the',/7x,'species triplet ',a,', ',a,', ',a,' among the blocks',/7x,'for "mixture parameters". A single',' cation appears three times. This',/7x,'is not a valid',' combination in the context of the normal Pitzer',/7x,'treatment of electrolyte solutions. The corresponding',' parameters are',/7x,'zero by convention.')

            nerr = nerr + 1
            go to 999
        end if

        if (qdup12 .or. qdup13 .or. qdup23) then
            ! Have one cation appearing two times (ccc').
            write (noutpt,1144) unam1(1:j2),unam2(1:j3),unam3(1:j4)
            write (nttyo,1144) unam1(1:j2),unam2(1:j3),unam3(1:j4)
1144 format(/' * Error - (EQPT/tripck) Have found an illegal data',' block for the',/7x,'species triplet ',a,', ',a,', ',a,' among the blocks',/7x,'for "mixture parameters". A single',' cation appears twice with another',/7x,'cation. This is',' not a valid combination in the context of the normal',/7x,'Pitzer treatment of electrolyte solutions. The',' corresponding',/7x,'parameters are not used.')

            nerr = nerr + 1
            go to 999
        end if
    end if

    if (nc .eq. 2) then
        if (qdup12 .or. qdup23 .or. qdup13) then
            if (nn .eq. 1) then
                ! Have a repeated cation and a neutral (ccn).
                write (noutpt,1150) unam1(1:j2),unam2(1:j3),unam3(1:j4)
                write (nttyo,1150) unam1(1:j2),unam2(1:j3),unam3(1:j4)
1150 format(/' * Error - (EQPT/tripck) Have found an illegal',' data block for the',/7x,'species triplet ',a,', ',a,', ',a,' among the blocks',/7x,'for "mixture parameters". One',' cation appears twice with a neutral.',/7x,'This is not a',' valid combination in the context of the normal Pitzer',/7x,'treatment of electrolyte solutions. The',' corresponding parameters',/7x,'are not used.')
            else if (na .eq. 1) then
                ! Have a repeated cation and an anion (cca).
                write (noutpt,1152) unam1(1:j2),unam2(1:j3),unam3(1:j4)
                write (nttyo,1152) unam1(1:j2),unam2(1:j3),unam3(1:j4)
1152 format(/' * Error - (EQPT/tripck) Have found an illegal',' data block for the',/7x,'species triplet ',a,', ',a,', ',a,' among the blocks',/7x,'for "mixture parameters". One',' cation appears twice with an anion.',/7x,'This is not a',' valid combination for the present type of block. The',/7x,'corresponding parameters properly appear in a block',' under',/7x,'"single-salt parameters".')
            end if

            nerr = nerr + 1
            go to 999
        else
            if (nn .eq. 1) then
                ! Have two cations and a neutral (c'cn).
                write (noutpt,1154) unam1(1:j2),unam2(1:j3),unam3(1:j4)
                write (nttyo,1154) unam1(1:j2),unam2(1:j3),unam3(1:j4)
1154 format(/' * Error - (EQPT/tripck) Have found an illegal',' data block for the',/7x,'species triplet ',a,', ',a,', ',a,' among the blocks',/7x,'for "mixture parameters". Two',' cations appear with a neutral.',/7x,'This is not a',' valid combination in the context of the normal Pitzer',/7x,'treatment of electrolyte solutions as it corresponds',' to a',/7x,'quaternary system.')

                nerr = nerr + 1
                go to 999
            end if

            ! Have two cations and an anion (cc'a). That is a legal
            ! combination.
        end if
    end if

    ! Check for other illegal combinations involving two or more
    ! anions.
    if (na .eq. 3) then
        if (.not.(qdup12 .or. qdup13 .or. qdup23)) then
            ! Have three distinct anions (aa'a'').
            write (noutpt,1160) unam1(1:j2),unam2(1:j3),unam3(1:j4)
            write (nttyo,1160) unam1(1:j2),unam2(1:j3),unam3(1:j4)
1160 format(/' * Error - (EQPT/tripck) Have found an illegal data',' block for the',/7x,'species triplet ',a,', ',a,', ',a,' among the blocks',/7x,'for "mixture parameters". Three',' distinct anions appear. This is',/7x,'not a valid',' combination in the context of the normal Pitzer',/7x,'treatment of electrolyte solutions as it corresponds',' to a',/7x,'quaternary system.')

            nerr = nerr + 1
            go to 999
        end if

        if (qdup12 .and. qdup13) then
            ! Have one anion appearing three times (aaa).
            write (noutpt,1162) unam1(1:j2),unam2(1:j3),unam3(1:j4)
            write (nttyo,1162) unam1(1:j2),unam2(1:j3),unam3(1:j4)
1162 format(/' * Error - (EQPT/tripck) Have found an illegal data',' block for the',/7x,'species triplet ',a,', ',a,', ',a,' among the blocks',/7x,'for "mixture parameters". A single',' anion appears three times. This',/7x,'is not a valid',' combination in the context of the normal Pitzer',/7x,'treatment of electrolyte solutions. The corresponding',' parameters are',/7x,'zero by convention.')

            nerr = nerr + 1
            go to 999
        end if

        if (qdup12 .or. qdup13 .or. qdup23) then
            ! Have one anion appearing two times (aaa').
            write (noutpt,1164) unam1(1:j2),unam2(1:j3),unam3(1:j4)
            write (nttyo,1164) unam1(1:j2),unam2(1:j3),unam3(1:j4)
1164 format(/' * Error - (EQPT/tripck) Have found an illegal data',' block for the',/7x,'species triplet ',a,', ',a,', ',a,' among the blocks',/7x,'for "mixture parameters". A single',' anion appears twice with another',/7x,'anion. This is',' not a valid combination in the context of the normal',/7x,'Pitzer treatment of electrolyte solutions. The',' corresponding',/7x,'parameters are not used.')

            nerr = nerr + 1
            go to 999
        end if
    end if

    if (na .eq. 2) then
        if (qdup12 .or. qdup23 .or. qdup13) then
            if (nn .eq. 1) then
                ! Have a repeated anion and a neutral (aan).
                write (noutpt,1170) unam1(1:j2),unam2(1:j3),unam3(1:j4)
                write (nttyo,1170) unam1(1:j2),unam2(1:j3),unam3(1:j4)
1170 format(/' * Error - (EQPT/tripck) Have found an illegal',' data block for the',/7x,'species triplet ',a,', ',a,', ',a,' among the blocks',/7x,'for "mixture parameters". One',' anion appears twice with a neutral.',/7x,'This is not a',' valid combination in the context of the normal Pitzer',/7x,'treatment of electrolyte solutions. The',' corresponding parameters',/7x,'are not used.')
            else if (nc .eq. 1) then
                ! Have a repeated anion and a cation (aac).
                write (noutpt,1173) unam1(1:j2),unam2(1:j3),unam3(1:j4)
                write (nttyo,1173) unam1(1:j2),unam2(1:j3),unam3(1:j4)
1173 format(/' * Error - (EQPT/tripck) Have found an illegal',' data block for the',/7x,'species triplet ',a,', ',a,', ',a,' among the blocks',/7x,'for "mixture parameters". One',' anion appears twice with a cation.',/7x,'This is not a',' valid combination for the present type of block. The',/7x,'corresponding parameters properly appear in a block',' under',/7x,'"single-salt parameters".')
            end if

            nerr = nerr + 1
            go to 999
        else
            if (nn .eq. 1) then
                ! Have two anions and a neutral (a'an).
                write (noutpt,1174) unam1(1:j2),unam2(1:j3),unam3(1:j4)
                write (nttyo,1174) unam1(1:j2),unam2(1:j3),unam3(1:j4)
1174 format(/' * Error - (EQPT/tripck) Have found an illegal',' data block for the',/7x,'species triplet ',a,', ',a,', ',a,' among the blocks',/7x,'for "mixture parameters". Two',' anions appear with a neutral.',/7x,'This is not a',' valid combination in the context of the normal Pitzer',/7x,'treatment of electrolyte solutions as it corresponds',' to a',/7x,'quaternary system.')

                nerr = nerr + 1
                go to 999
            end if

            ! Have two anions and a cation (aa'c). That is a legal
            ! combination.
        end if
    end if

999 continue
end subroutine tripck