subroutine tstari(afrc1,alphar,betar,delxi,iodb,modr,nodbmx,noutpt,nrct,nrctmx,nrct1,nsscmx,qodeok,rirec1,rirecp,rrelr1,rrelrp,sscrew,time1,tistrt,ureac)
    !! This subroutine tests the accuracy of the rate law integration
    !! by comparing computed rates with predicted values.
    !! This subroutine is called by:
    !!   EQ6/path.f
    !! Principal input:
    !! Principal output:
    !!   qodeok = logical flag, .true. if the rate law integration
    !!              tolerances have been satisfied.
    implicit none

    ! Calling sequence variable declarations.
    integer :: nodbmx
    integer :: nrctmx
    integer :: nsscmx

    integer :: noutpt

    integer :: iodb(nodbmx)

    integer :: nrct
    integer :: nrct1

    logical :: qodeok

    character(len=24) :: ureac(nrctmx)

    real(kind=8) :: afrc1(nrctmx)
    real(kind=8) :: alphar(nrct1)
    real(kind=8) :: betar(nrct1)
    real(kind=8) :: modr(nrctmx)
    real(kind=8) :: rrelr1(nrctmx)
    real(kind=8) :: rrelrp(nrctmx)
    real(kind=8) :: sscrew(nsscmx)

    real(kind=8) :: delxi
    real(kind=8) :: rirec1
    real(kind=8) :: rirecp
    real(kind=8) :: time1
    real(kind=8) :: tistrt

    ! Local variable declarations.
    integer :: j2
    integer :: nrc

    integer :: ilnobl

    logical :: qbad

    real(kind=8) :: atx
    real(kind=8) :: etx
    real(kind=8) :: rdelmd
    real(kind=8) :: rdelti
    real(kind=8) :: ss4

    qodeok = .true.
    qbad = .false.
    ss4 = sscrew(4)**2

    do nrc = 1,nrct
        ! The error in a relative rate is acceptable is any of three
        ! tests is satisfied.
        qbad = .false.

        ! Test the fractional error in the relative rate.
        if (abs(betar(nrc)) .gt. sscrew(4)) then
            ! Test the absolute error in the relative rate.
            if (abs(alphar(nrc)) .gt. ss4) then
                if (abs(modr(nrc)) .le. 0.) then
                    qbad = .true.
                else
                    ! Test the error in the moles destroyed for this step
                    ! against the total moles destroyed.
                    rdelmd = (alphar(nrc)*delxi)/abs(modr(nrc))

                    if (abs(rdelmd) .gt. sscrew(3)) then
                        qbad = .true.
                    end if
                end if
            end if
        end if

        if (qbad) then
            qodeok = .false.

            if (iodb(2) .ge. 1) then
                j2 = ilnobl(ureac(nrc))
                write (noutpt,1010) ureac(nrc)(1:j2),rrelr1(nrc),rrelrp(nrc),betar(nrc),afrc1(nrc)
1010 format(/3x,'Relative rate for ',a,':',/5x,'Calculated value= ',1pe12.5,/5x,' Predicted value= ',e12.5,/5x,'Fractional error= ',e12.5,/5x,'Affinity= ',e12.5,/)
            end if
        end if
    end do

    ! The error in the inverse rate is acceptable is either of two
    ! tests is satisfied.
    qbad = .false.

    ! Test the fractional error in the inverse rate.
    if (abs(betar(nrct1)) .gt. sscrew(4)) then
        ! Test the absolute error in the time.
        atx = alphar(nrct1)*delxi

        ! XX
        if (abs(atx) .gt. 0.) then
            etx = time1 - tistrt

            if (etx.eq. 0.) then
                qbad = .true.
            else
                rdelti = atx/etx

                ! Test the error in the time against the elapsed time.
                if (abs(rdelti) .gt. sscrew(3)) then
                    qbad = .true.
                end if
            end if
        end if
    end if

    if (qbad) then
        qodeok = .false.

        if (iodb(2) .ge. 1) then
            write (noutpt,1000) rirec1,rirecp,betar(nrct1),rdelti
1000 format(/3x,'Inverse rate:',/5x,'Calculated value= ',1pe12.5,/5x,' Predicted value= ',e12.5,/5x,'Fractional error= ',e12.5,/5x,'Fractional error in time= ',e12.5,/)
        end if
    end if
end subroutine tstari