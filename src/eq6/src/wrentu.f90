subroutine wrentu(actw,eh,fo2lg,iopg,iopt,kstep,nopgmx,noptmx,nttyo,ph,qredox,time1,xi1)
    !! This subroutine writes entertainment for the user while the
    !! run is underway. This output gives an idea of how the
    !! calculation is progressing. In that respect, it is somewhat
    !! like a progress bar.
    !! This subroutine is called by:
    !!   EQ6/path.f
    !! Principal input:
    !! Principal output:
    implicit none

    integer :: nopgmx
    integer :: noptmx

    integer :: nttyo

    integer :: iopg(nopgmx)
    integer :: iopt(noptmx)

    integer :: kstep

    logical :: qredox

    real(kind=8) :: actw
    real(kind=8) :: eh
    real(kind=8) :: fo2lg
    real(kind=8) :: ph
    real(kind=8) :: time1
    real(kind=8) :: xi1

    ! Local variable declarations.
    ! None
    if (iopt(2) .le. 0) then
        if (qredox) then
            if (iopg(1) .gt. 0) then
                write (nttyo,1000) kstep,xi1,ph,actw,fo2lg
1000 format(3x,'step=',i5,', Xi=',1pe10.3,', pH=',0pf7.3,', aw=',f6.3,', log fO2= ',f7.3)
            else
                write (nttyo,1010) kstep,xi1,ph,fo2lg
1010 format(3x,'step=',i5,', Xi=',1pe10.3,', pH=',0pf7.3,', log fO2= ',f7.3)
            end if
        else
            if (iopg(1) .gt. 0) then
                write (nttyo,1020) kstep,xi1,ph,actw
1020 format(3x,'step=',i5,', Xi=',1pe10.3,', pH=',0pf7.3,', aw=',f6.3)
            else
                write (nttyo,1030) kstep,xi1,ph
1030 format(3x,'step=',i5,', Xi=',1pe10.3,', pH=',0pf7.3)
            end if
        end if
    else
        if (qredox) then
            if (iopg(1) .gt. 0) then
                write (nttyo,1100) kstep,time1,ph,actw,fo2lg
1100 format(3x,'step=',i5,', time=',1pe10.3,' s, pH=',0pf7.3,', aw=',f6.3,', log fO2= ',f7.3)
            else
                write (nttyo,1110) kstep,time1,ph,fo2lg
1110 format(3x,'step=',i5,', time=',1pe10.3,' s, pH=',0pf7.3,', log fO2= ',f7.3)
            end if
        else
            if (iopg(1) .gt. 0) then
                write (nttyo,1120) kstep,time1,ph,actw
1120 format(3x,'step=',i5,', time=',1pe10.3,' s, pH=',0pf7.3,', aw=',f6.3)
            else
                write (nttyo,1130) kstep,time1,ph
1130 format(3x,'step=',i5,', time=',1pe10.3,' s, pH=',0pf7.3)
            end if
        end if
    end if
end subroutine wrentu