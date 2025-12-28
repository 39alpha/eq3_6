subroutine setffg(csts,iindx1,iffg,ipndx1,jpflag,jsflag,kbt,kdim,kmax,km1,kmt,kx1,kxt,losp,moffg,mtb,mtbaq,nbaspd,nbt,nbtmax,ncmpr,nffg,nffgmx,noutpt,nphasx,npt,nptmax,nstmax,nsts,nstsmx,nstsr,nttyo,qloffg,uffg,uspec,uzvec1,zvclg1,zvec1)
    !! This subroutine puts new fictive fugacity-fixing phases into the
    !! matrix if the corresponding masses to add to the equilibrium
    !! system (ES) are positive. It keep old such phases in if the
    !! corresponding net masses are positive. It deletes any fictive
    !! fugacity-fixing phases left over from a previous run for which
    !! the corresponding fixed fugacity options no longer apply.
    !! To carry out this last function, this subroutine must be called
    !! even if nffg is zero. If necessary this subroutine recalculates
    !! the mass balance totals.
    !! This subroutine is called by:
    !!   EQ6/eq6.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: kmax
    integer :: nbtmax
    integer :: nffgmx
    integer :: nptmax
    integer :: nstmax
    integer :: nstsmx

    integer :: noutpt
    integer :: nttyo

    integer :: iindx1(kmax)
    integer :: iffg(nffgmx)
    integer :: ipndx1(kmax)
    integer :: jpflag(nptmax)
    integer :: jsflag(nstmax)
    integer :: nbaspd(nbtmax)
    integer :: ncmpr(2,nptmax)
    integer :: nphasx(nstmax)
    integer :: nsts(nstsmx)
    integer :: nstsr(2,nstmax)

    integer :: kbt
    integer :: kdim
    integer :: km1
    integer :: kmt
    integer :: kx1
    integer :: kxt
    integer :: nbt
    integer :: nffg
    integer :: npt

    logical :: qloffg

    character(len=48) :: uspec(nstmax)
    character(len=48) :: uzvec1(kmax)
    character(len=24) :: uffg(nffgmx)

    real(kind=8) :: csts(nstsmx)
    real(kind=8) :: losp(nstmax)
    real(kind=8) :: moffg(nffgmx)
    real(kind=8) :: mtb(nbtmax)
    real(kind=8) :: mtbaq(nbtmax)
    real(kind=8) :: zvec1(kmax)
    real(kind=8) :: zvclg1(kmax)

    ! Local variable declarations.
    integer :: ier
    integer :: j2
    integer :: kcol
    integer :: n
    integer :: nb
    integer :: nerr
    integer :: np
    integer :: ns
    integer :: nse

    integer :: ilnobl

    logical :: qreclc

    real(kind=8) :: cx
    real(kind=8) :: mtbo
    real(kind=8) :: mtbx
    real(kind=8) :: mx
    real(kind=8) :: mz

    real(kind=8) :: coefst
    real(kind=8) :: texp
    real(kind=8) :: tlg

    qreclc = qloffg
    nerr = 0

    do kcol = km1,kxt
        ns = iindx1(kcol)
        losp(ns) = zvclg1(kcol)
    end do

    do n = 1,nffg
        mx = moffg(n)
        ns = iffg(n)
        np = nphasx(ns)

        if (jpflag(np) .ne. -1) then
            ! The fixed fugacity phase was not already in the matrix.
            losp(ns) = -99999.

            if (mx .gt. 0.) then
                qreclc = .true.
                jpflag(np) = -1
                jsflag(ns) = -1
                losp(ns) = tlg(mx)
            else if (mx .lt. 0.) then
                j2 = ilnobl(uffg(n))
                write (noutpt,1000) uffg(n)(1:j2)
                write (nttyo,1000) uffg(n)(1:j2)
1000 format(/' * Error - (EQ6/setffg) The new fixed fugacity',' phase',/7x,a," can't have negative moles added.")

                nerr = nerr + 1
            end if
        else if (mx .ne. 0.) then
            ! The fixed fugacity phase was already in the matrix.
            qreclc = .true.
            mz = texp(losp(ns))
            mz = mz + mx

            if (mz .gt. 0.) then
                losp(ns) = tlg(mz)
            else
                j2 = ilnobl(uffg(n))
                write (noutpt,1010) uffg(n)(1:j2),mz
                write (nttyo,1010) uffg(n)(1:j2),mz
1010 format(/' * Error - (EQ6/setffg) The adjusted number of',/7x,'moles of ',a,' is ',g12.5,'. A negative value is not',/7x,'permitted.')

                nerr = nerr + 1
            end if
        end if
    end do

    if (nerr .gt. 0) then
        stop
    end if

    if (.not.qreclc) then
        go to 200
    end if

    ! The following may may combine any of three functions:
    !   1. Add new fixed fugacity phases to the matrix.
    !   2. Delete left over fixed fugacity phases from the matrix.
    !   3. Reset the zvec1 and zvclg1 array elements for fixed
    !      fugacity phases.
    call miidxz(ier,iindx1,ipndx1,jpflag,jsflag,kbt,kdim,kmax,km1,kmt,kx1,kxt,losp,ncmpr,noutpt,npt,nptmax,nstmax,nttyo,uspec,uzvec1,zvclg1,zvec1)

    if (ier .gt. 0) then
        write (noutpt,1020)
        write (nttyo,1020)
1020 format(/" * Error - (EQ6/setffg) Can't recover from having",' exceeded the maximum ',i4,/7x,'elements of the iindx1',' array. Increase the dimensioning parameter kpar.')

        stop
    end if

    ! Recalculate the mass balance totals. This is done by adding to
    ! get the total, not by subtracting from the existing total.
    ! This purges the masses of any left over fixed fugacity phases
    ! from the mass balance totals and adds any new masses to these
    ! totals.
    if (nffg .gt. 0) then
        write (noutpt,1030)
    end if

1030 format(//6x,' --- Adding Fictive Fugacity-Fixing Phases to the',' ES ---',//13x,'Name',22x,'Moles Added',/)

    do n = 1,nffg
        ns = iffg(n)
        write (noutpt,1040) uspec(ns),moffg(n)
1040 format(11x,a24,3x,1pe12.5)
    end do

    write (noutpt,1050)
1050 format(//11x,' --- Component Total Numbers of Moles ---',//35x,'Initial',7x,'Adjusted',/)

    do nb = 1,nbt
        nse = nbaspd(nb)
        mtbo = mtb(nb)
        mtbx = mtbaq(nb)

        do kcol = km1,kxt
            ns = iindx1(kcol)
            cx = coefst(csts,nsts,nstsmx,nstsr,nb,ns,nstmax)
            mtbx = mtbx + cx*zvec1(kcol)
        end do

        mtb(nb) = mtbx
        write (noutpt,1060) uspec(nse),mtbo,mtbx
1060 format(5x,a24,3x,1pe12.5,3x,e12.5)
    end do

    write (noutpt,1070)
1070 format(1x)

200 continue

    do n = 1,nffg
        moffg(n) = 0.
    end do
end subroutine setffg