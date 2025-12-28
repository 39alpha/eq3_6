subroutine echgex(axlks,cdrs,cgexj,iern1,iern2,jern1,jern2,jetmax,jgext,jpflag,jsflag,narxmx,narxt,ndrs,ndrsmx,ndrsr,netmax,noutpt,nptmax,ntprmx,ntprt,nstmax,press,tempc,ugexj,ugexmo,uphase,uspec,xlks)
    !! This subroutine echoes a table for the generic ion exchangers,
    !! describing the setup of species, reactions, and corresponding
    !! thermodynamic data.
    !! This subroutine is called by:
    !!   EQ3NR/eq3nr.f
    !!   EQ6/eq6.f
    !! Principal input:
    !!   noutpt = unit number of the output file
    !! Principal output:
    !!   None
    implicit none

    ! Calling sequence variable declarations.
    integer :: jetmax
    integer :: narxmx
    integer :: ndrsmx
    integer :: netmax
    integer :: nptmax
    integer :: ntprmx
    integer :: nstmax

    integer :: noutpt

    integer :: jern1(jetmax,netmax)
    integer :: jern2(jetmax,netmax)
    integer :: jgext(netmax)
    integer :: jpflag(nptmax)
    integer :: jsflag(nstmax)
    integer :: narxt(ntprmx)
    integer :: ndrs(ndrsmx)
    integer :: ndrsr(2,nstmax)

    integer :: iern1
    integer :: iern2
    integer :: ntprt

    integer :: ilnobl

    character(len=48) :: uspec(nstmax)
    character(len=24) :: ugexmo(netmax)
    character(len=24) :: uphase(nptmax)
    character(len=8) :: ugexj(jetmax,netmax)

    real(kind=8) :: axlks(narxmx,ntprmx,nstmax)
    real(kind=8) :: cdrs(ndrsmx)
    real(kind=8) :: cgexj(jetmax,netmax)
    real(kind=8) :: xlks(nstmax)

    real(kind=8) :: press
    real(kind=8) :: tempc

    ! Local variable declarations.
    integer :: i
    integer :: j
    integer :: je
    integer :: j3
    integer :: j4
    integer :: j5
    integer :: ne
    integer :: np
    integer :: nr1
    integer :: nr2
    integer :: ns
    integer :: nt

    if (iern2 .lt. iern1) then
        go to 999
    end if

    write (noutpt,1000)
1000 format(/11x,' --- Generic Ion Exchange Phase Setup ---',/)

    write (noutpt,1010) tempc,press
1010 format(6x,'Log K data, etc., are for ',f10.3,'C and ',g12.5,' bars',/)

    do np = iern1,iern2
        ne = np - iern1 + 1
        j5 = ilnobl(ugexmo(ne))
        j4 = ilnobl(uphase(np))

        if (jpflag(np) .le. 0) then
            write (noutpt,1020) uphase(np)(1:j4)
1020 format(/16x,'--- ',a,' ---')
        else
            write (noutpt,1030) uphase(np)(1:j4)
1030 format(/16x,'--- ',a,' (Suppressed) ---')
        end if

        write (noutpt,1040) ugexmo(ne)(1:j5)
1040 format(/11x,'Model type= ',a)

        do je = 1,jgext(ne)
            j3 = ilnobl(ugexj(je,ne))
            write (noutpt,1050) ugexj(je,ne)(1:j3)
1050 format(/11x,'--- ',a,' ---',/)

            write (noutpt,1060) cgexj(je,ne)
1060 format(6x,'Site stoichiometric factor N= ',1pg12.5,/)

            do ns = jern1(je,ne),jern2(je,ne)
                ! Calling sequence substitutions:
                !   noutpt for nf
                call prreac(cdrs,ndrs,ndrsmx,ndrsr,noutpt,ns,nstmax,uspec)

                nr1 = ndrsr(1,ns)
                nr2 = ndrsr(2,ns)
                nt = nr2 - nr1 + 1

                if (nt .gt. 1) then
                    if (jsflag(ns) .le. 0) then
                        write (noutpt,1070) xlks(ns)
1070 format(/5x,'Log K= ',f12.4)

                        do j = 1,ntprt
                            write (noutpt,1100) j
1100 format(/7x,'Coefficients for temperature range ',i2,':')

                            write (noutpt,1110) (axlks(i,j,ns), i = 1,narxt(j))
1110 format( (4x,5(2x,g13.6)) )
                        end do

                        write (noutpt,1150)
1150 format(1x)
                    else
                        write (noutpt,1170)
1170 format(2x,'Suppressed',/)
                    end if
                end if
            end do
        end do
    end do

999 continue
end subroutine echgex