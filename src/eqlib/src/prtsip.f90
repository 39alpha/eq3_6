subroutine prtsip(affpd,iopr,ir1,ir2,jpflag,kpsat,kpsst,noprmx,noutpt,nptmax,sidrph,tolspf,ugroup,uphase)
    !! This subroutine prints a table of saturation indices and
    !! affinities of phases whose indices lie in the range ir1 to ir2.
    !! The group corresponding to this range is described by ugroup.
    !! The level of printing is controlled by the print control flag
    !! iopr(7):
    !!   -1 = Don't print
    !!    0 = Print for those phases not undersaturated by
    !!          more than 10 kcal
    !!    1 = Print for all phases
    !! This subroutine is called by:
    !!   EQLIB/prtsat.f
    !! Principal input:
    !!   affpd  = array of affinities of phases to form (precipitate),
    !!              computed from the reactions read from the data file
    !!   iopr   = array of print control options
    !!   ir1    = start of the range of pure liquid phases
    !!   ir2    = end of the range of pure liquid phases
    !!   jpflag = array of phase status flags
    !!   kpsat  = number of saturated phases
    !!   kpsst  = number of supersaturated phases
    !!   sidrph = array of saturation indices of the various phases,
    !!              computed from the reactions read from the data file
    !!   tolspf = saturation print flag tolerance, used to flag those
    !!              phases which are close to saturation
    !!   ugroup = string describing the group of phases for which
    !!              a table is to be printed
    !!   uphase = array of phase names
    !! Principal output:
    !!   None
    implicit none

    ! Calling sequence variable declarations.
    integer :: noprmx
    integer :: nptmax

    integer :: noutpt

    integer :: iopr(noprmx)
    integer :: jpflag(nptmax)
    integer :: ir1
    integer :: ir2
    integer :: kpsat
    integer :: kpsst

    character(len=24) :: uphase(nptmax)
    character(len=24) :: ugroup

    real(kind=8) :: affpd(nptmax)
    real(kind=8) :: sidrph(nptmax)
    real(kind=8) :: tolspf

    ! Local variable declarations.
    integer :: j2
    integer :: kount
    integer :: np

    integer :: ilnobl

    real(kind=8) :: afx

    j2 = ilnobl(ugroup)
    write (noutpt,1000) ugroup(1:j2)
1000 format(//11x,'--- Saturation States of ',a,' ---',//7x,'Phase   ',19x,'Log Q/K',4x,'Affinity, kcal',/)

    kount = 0

    do np = ir1,ir2
        if (jpflag(np) .lt. 2) then
            if (iopr(7).gt.0 .or. affpd(np).ge.-10.) then
                kount = kount + 1
                afx = abs(affpd(np))

                if (afx .le. tolspf) then
                    kpsat = kpsat + 1
                    write (noutpt,1010) uphase(np),sidrph(np),affpd(np)
1010 format(5x,a24,3x,f10.5,3x,f10.5,5x,'SATD')
                else if (affpd(np) .gt. tolspf) then
                    kpsst = kpsst + 1
                    write (noutpt,1020) uphase(np),sidrph(np),affpd(np)
1020 format(5x,a24,3x,f10.5,3x,f10.5,5x,'SSATD')
                else
                    write (noutpt,1030) uphase(np),sidrph(np),affpd(np)
1030 format(5x,a24,3x,f10.5,3x,f10.5)
                end if
            end if
        end if
    end do

    if (kount .le. 0) then
        write (noutpt,1040)
    end if

1040 format(5x,'None')

    if (iopr(7) .eq. 0) then
        write (noutpt,1050)
    end if

1050 format(/8x,'Phases with affinities less than -10 kcal are not',' listed.')

    write (noutpt,1060)
1060 format(1x)
end subroutine prtsip