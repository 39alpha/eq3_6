subroutine evptzc(amu,aslm,ipbtmx,jpfcmx,jptffl,nmut,nmutmx,noutpt,nttyo,nslt,nsltmx,pmu,pslamn,tempc)
    !! This subroutine computes the S-lambda(n) (pslamn) and mu (pmu)
    !! coefficients of Pitzer's equations at the temperature tempc.
    !! This subroutine is called by:
    !!   EQLIB/evdata.f
    !! Principal input:
    !!   amu    = array of coefficients for calculating Pitzer mu
    !!            coefficients (pmu) as a function of temperature
    !!   aslm   = array of coefficients for calculating Pitzer
    !!              S-lambda(n) coefficients (pslamn) as a function
    !!              of temperature
    !!   jptffl = Pitzer parameter temperature coefficient flag:
    !!              -1 = 25C centric Taylor's series truncated at
    !!                     second order
    !!               0 = 25C centric LLNL 5TERM function
    !!                     (maximal 5th order)
    !!               1 = non-25C centric eight-term Greenberg and
    !!                     Moller (1989) function
    !!   nmut   = the number of species triples for which there are
    !!              mu coefficients
    !!   nslt   = the number of species pairs for which there are
    !!              S-lambda(n) coefficients
    !!   tempc  = the temperature (C)
    !! Principal output:
    !!   pmu    = array of Pitzer mu coefficients
    !!   pslamn = array of Pitzer S-lambda(n) coefficients
    implicit none

    ! Calling sequence variable declarations.
    integer :: ipbtmx
    integer :: jpfcmx
    integer :: nmutmx
    integer :: nsltmx

    integer :: noutpt
    integer :: nttyo

    integer :: jptffl
    integer :: nmut
    integer :: nslt

    real(kind=8) :: amu(jpfcmx,nmutmx)
    real(kind=8) :: aslm(jpfcmx,0:ipbtmx,nsltmx)
    real(kind=8) :: pmu(nmutmx)
    real(kind=8) :: pslamn(0:ipbtmx,nsltmx)

    real(kind=8) :: tempc

    ! Local variable declarations.
    integer :: j
    integer :: jpfclm
    integer :: j2
    integer :: k
    integer :: nmu
    integer :: nsl

    integer :: ilnobl

    character(len=8) :: ux8

    real(kind=8) :: px
    real(kind=8) :: t
    real(kind=8) :: tr

    ! Note: the dimension of dtf must match the maximum number of terms
    ! in any programmed temperature function.
    real(kind=8) :: dtf(2:8)

    t = tempc + 273.15
    tr = 298.15

    if (jptffl .eq. -1) then
        ! Classical 25C-centric Taylor's series truncated at second order.
        dtf(2) = t - tr
        dtf(3) = 0.5*dtf(2)*dtf(2)
        jpfclm = min(jpfcmx,3)
    else if (jptffl .eq. 0) then
        ! LLNL maximal 5-term equation.
        dtf(2) = (1./t) - (1./tr)
        dtf(3) = log(t/tr)
        dtf(4) = t - tr
        dtf(5) = t**2 - tr**2
        jpfclm = min(jpfcmx,5)
    else if (jptffl .eq. 1) then
        ! Greenberg and Moller (1989) eight-term equation.
        if ((abs(t - 227.) .le. 1.e-3)    .or. (abs(t - 263.) .le. 1.e-3)    .or. (abs(t - 680.) .le. 1.e-3)) then
            write (noutpt,1000) t,tempc
            write (nttyo,1000) t,tempc
1000 format(/' * Error - (EQLIBG/evptzc) The temperature of ',f7.4,'K (',f7.4,'C) is',/7x,'too close to a singularity',' in the TEQUIL eight-parameter Pitzer',/7x,'parameter',' temperature function.')

            stop
        end if

        dtf(2) = t
        dtf(3) = 1./t
        dtf(4) = log(t)
        dtf(5) = 1./(t - 263.)
        dtf(6) = t*t
        dtf(7) = 1./(680. - t)
        dtf(8) = 1./(t - 227.)
        jpfclm = min(jpfcmx,8)
    else
        write (ux8,'(i5)') jptffl
        call lejust(ux8)
        j2 = ilnobl(ux8)
        write (noutpt,1010) ux8(1:j2)
        write (nttyo,1010) ux8(1:j2)
1010 format(/' * Error - (EQLIBG/evptzc) Programming error trap:',' Have an',/7x,'illegal value of ',a,' for the variable',' jptffl, which determines',/7x,'the Pitzer parameter',' temperature function. Do not have such a',/7x,'function',' programmed that matches this value.')

        stop
    end if

    do nsl = 1,nslt
        do k = 0,2
            px = aslm(1,k,nsl)

            do j = 2,jpfclm
                px = px + aslm(j,k,nsl)*dtf(j)
            end do

            pslamn(k,nsl) = px
        end do
    end do

    do nmu = 1,nmut
        px = amu(1,nmu)

        do j = 2,jpfclm
            px = px + amu(j,nmu)*dtf(j)
        end do

        pmu(nmu) = px
    end do
end subroutine evptzc