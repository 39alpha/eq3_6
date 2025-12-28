subroutine flgchk(jpflag,jsflag,ncmpra,npta,nptmax,nstmax,qclnsa)
    !! This subroutine recalculates the jpflag array to insure
    !! consistency with the jsflag array. The jpflag value of a phase
    !! can't be lower than the lowest value of any of the species which
    !! comprise that phase. For example, if all of the species of
    !! a phase are not present (jsflag = 2), the phase can not be
    !! present (jpflag = 2). Also, if a solution phase (one with more
    !! than one component species) has only one species that is active
    !! (jsflag .lt. 2) and that species was created by cloning, the
    !! solution phase is purged (jpflag is set to 2). This subroutine
    !! should be called after calling other subroutines which may purge
    !! or suppress species. EQLIB/flgset.f is an example of such.
    !! This subroutine is called by:
    !!   EQ3NR/eq3nr.f
    !!   EQ6/eq6.f
    !! Principal input:
    !!   jpflag = array of status flags for phases
    !!   jsflag = array of status flags for species
    !!   ncmpra = array giving the range of indices of species
    !!              belonging to a given phase ('a' data set)
    !!   npta   = the number of phaes in the 'a' data set.
    !!   qclnsa = array of logical flags indicating if a species was
    !!              created by cloning.
    !! Principal output:
    !!   None
    implicit none

    ! Calling sequence variable declarations.
    integer :: nptmax
    integer :: nstmax

    integer :: jpflag(nptmax)
    integer :: jsflag(nstmax)
    integer :: ncmpra(2,nptmax)
    integer :: npta

    logical :: qclnsa(nstmax)

    ! Local variable declarations.
    integer :: ijs0
    integer :: ijs1
    integer :: jsflmn
    integer :: njs0
    integer :: npa
    integer :: nr1
    integer :: nr2
    integer :: nsa
    integer :: nt

    do npa = 1,npta
        if (jpflag(npa) .lt. 2) then
            nr1 = ncmpra(1,npa)
            nr2 = ncmpra(2,npa)
            jsflmn = 2

            do nsa = nr1,nr2
                jsflmn = min(jsflag(nsa),jsflmn)
            end do

            jpflag(npa) = max(jpflag(npa),jsflmn)
        end if
    end do

    do npa = 1,npta
        if (jpflag(npa) .lt. 2) then
            nr1 = ncmpra(1,npa)
            nr2 = ncmpra(2,npa)
            nt = nr2 - nr1 + 1

            if (nt .ge. 2) then
                ! Note the following local variables:
                !   njs0 = counter for non-cloned species with jsflag = 0
                !   ijs0 = counter for cloned and non-cloned species with
                !          jsflag = 0
                !   ijs1 = counter for cloned and non-cloned species with
                !          jsflag = 1
                njs0 = 0
                ijs1 = 0
                ijs0 = 0

                do nsa = nr1,nr2
                    if (jsflag(nsa) .eq. 0) then
                        ijs0 = ijs0 + 1

                        if (.not.qclnsa(nsa)) then
                            njs0 = njs0 + 1
                        end if
                    else if (jsflag(nsa) .eq. 1) then
                        ijs1 = ijs1 + 1
                    end if
                end do

                if (njs0.le.0 .and. ijs0.lt.2) then
                    jpflag(npa) = 1

                    if (ijs1 .le. 0) then
                        jpflag(npa) = 2
                    end if
                end if
            end if
        end if
    end do
end subroutine flgchk