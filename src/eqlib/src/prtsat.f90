subroutine prtsat(affpd,iern1,iern2,ilrn1,ilrn2,imrn1,imrn2,iopr,iopt,ixrn1,ixrn2,jpflag,noutpt,noprmx,noptmx,nptmax,sidrph,tolspf,uphase)
    !! This subroutine prints tables of saturation indices and affinities
    !! for the various non-aqueous phases. The level of printing is
    !! controlled by the print control flag iopr(7):
    !!   -1 = Don't print
    !!    0 = Print for those phases not undersaturated by
    !!          more than 10 kcal
    !!    1 = Print for all phases
    !! This subroutine is called by:
    !!   EQ3NR/scripx.f
    !!   EQ6/scripz.f
    !! Principal input:
    !!   affpd  = array of affinities of phases to form (precipitate),
    !!              computed from the reactions read from the data file
    !!   iern1  = start of the range of generic ion exchangers
    !!   iern2  = end of the range of generic ion exchangers
    !!   ilrn1  = start of the range of pure liquid phases
    !!   ilrn2  = end of the range of pure liquid phases
    !!   imrn1  = start of the range of pure mineral phases
    !!   imrn2  = end of the range of pure mineral phases
    !!   iopt   = array of model option switches
    !!   iopr   = array of print control options
    !!   ixrn1  = start of the range of solid solution phases
    !!   ixrn2  = end of the range of solid solution phases
    !!   jpflag = array of phase status flags
    !!   sidrph = array of saturation indices of the various phases,
    !!              computed from the reactions read from the data file
    !!   tolspf = saturation print flag tolerance, used to flag those
    !!              phases which are close to saturation
    !!   uphase = array of phase names
    !! Principal output:
    !!   None
    implicit none

    ! Calling sequence variable declarations.
    integer :: noprmx
    integer :: noptmx
    integer :: nptmax

    integer :: noutpt

    integer :: iopr(noprmx)
    integer :: iopt(noptmx)
    integer :: jpflag(nptmax)
    integer :: iern1
    integer :: iern2
    integer :: ilrn1
    integer :: ilrn2
    integer :: imrn1
    integer :: imrn2
    integer :: ixrn1
    integer :: ixrn2

    character(len=24) :: uphase(nptmax)

    real(kind=8) :: affpd(nptmax)
    real(kind=8) :: sidrph(nptmax)
    real(kind=8) :: tolspf

    ! Local variable declarations.
    integer :: ir1
    integer :: ir2
    integer :: kpsat
    integer :: kpsst

    character(len=24) :: ugroup
    character(len=24) :: upusol
    character(len=24) :: upuliq
    character(len=24) :: usosol
    character(len=24) :: usoliq
    character(len=24) :: ugexch
    character(len=24) :: ux24
    character(len=24) :: uy24

    data upusol /'Pure Solids             '/,upuliq /'Pure Liquids            '/,usosol /'Solid Solutions         '/,usoliq /'Liquid Solutions        '/,ugexch /'Generic Ion Exchangers  '/

    ! Note: the following statements don't really do anything except
    ! cause the compiler not to complain that usoliq is not used.
    ux24 = usoliq
    uy24 = ux24
    usoliq = uy24

    if (iopr(7) .le. -1) then
        go to 999
    end if

    ! Compute and print phase saturation states. these refer to
    ! the dissolution reactions as they are written in the
    ! 'd' set.
    !   kpsat = the number of saturated phases
    !   kpsst = the number of supersaturated phases
    kpsat = 0
    kpsst = 0

    ! Print a table for the pure solids.
    ugroup = upusol
    ir1 = imrn1
    ir2 = imrn2
    call prtsip(affpd,iopr,ir1,ir2,jpflag,kpsat,kpsst,noprmx,noutpt,nptmax,sidrph,tolspf,ugroup,uphase)

    ! Print a table for the pure liquids.
    ugroup = upuliq
    ir1 = ilrn1
    ir2 = ilrn2
    call prtsip(affpd,iopr,ir1,ir2,jpflag,kpsat,kpsst,noprmx,noutpt,nptmax,sidrph,tolspf,ugroup,uphase)

    if (iopt(4) .ge. 1) then
        ! Print a table for the solid solutions.
        ugroup = usosol
        ir1 = ixrn1
        ir2 = ixrn2
        call prtsip(affpd,iopr,ir1,ir2,jpflag,kpsat,kpsst,noprmx,noutpt,nptmax,sidrph,tolspf,ugroup,uphase)

        ! No table is presently printed for non-aqueous liquid
        ! solutions.
    end if

    if (iern2 .ge. iern1) then
        ! Print a table for the generic ion exchangers.
        ugroup = ugexch
        ir1 = iern1
        ir2 = iern2
        call prtsip(affpd,iopr,ir1,ir2,jpflag,kpsat,kpsst,noprmx,noutpt,nptmax,sidrph,tolspf,ugroup,uphase)
    end if

    write (noutpt,1150)
1150 format(/11x,'--- Summary of Saturated and Supersaturated',' Phases ---',/)

    if (kpsat .le. 0) then
        write (noutpt,1160)
1160 format(16x,'There are no saturated phases.')
    else if (kpsat .eq. 1) then
        write (noutpt,1170)
1170 format(16x,'There is 1 saturated phase.')
    else
        write (noutpt,1180) kpsat
1180 format(16x,'There are ',i4,' saturated phases.')
    end if

    if (kpsst .le. 0) then
        write (noutpt,1200)
1200 format(16x,'There are no supersaturated phases.')
    else if (kpsst .eq. 1) then
        write (noutpt,1210)
1210 format(16x,'There is 1 supersaturated phase.')
    else
        write (noutpt,1220) kpsst
1220 format(16x,'There are ',i4,' supersaturated phases.')
    end if

    write (noutpt,1230)
1230 format(1x)

999 continue
end subroutine prtsat