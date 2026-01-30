subroutine setrcp(aftarg,dlxmax,dlxmin,dlxmx0,npslmx,nsscmx,nsslmx,prcinf,sscrew,tolaft,tolsat,tolsst,zkfac,zklgmn,zklogl,zklogu)
    use iso_fortran_env, only: dp => real64
    !! This subroutine sets the values of various run control parameters
    !! which are not read from the input file.
    !! This subroutine is called by:
    !!   EQ6/eq6.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: nsscmx

    integer :: npslmx
    integer :: nsslmx

    real(kind=8) :: sscrew(nsscmx)

    real(kind=8) :: aftarg
    real(kind=8) :: dlxmax
    real(kind=8) :: dlxmin
    real(kind=8) :: dlxmx0
    real(kind=8) :: prcinf
    real(kind=8) :: tolaft
    real(kind=8) :: tolsat
    real(kind=8) :: tolsst
    real(kind=8) :: zkfac
    real(kind=8) :: zklgmn
    real(kind=8) :: zklogl
    real(kind=8) :: zklogu

    ! Local variable declarations.
    integer :: n

    real(kind=8) :: tlg

    ! Supersaturation tolerance (if phase boundary searches are being
    ! made and this is exceeded, the step size is cut). Note that tolsst
    ! must be greater than tolsat (see aftarg and tolaft below).
    tolsst = 2.*tolsat

    ! Affinity target (the target for a search for a phase appearance
    ! boundary).
    aftarg = 0.5*(tolsst + tolsat)

    ! Affinity tolerance (the tolerance about the affinity target).
    ! This requires the affinity at the located phase appearnce
    ! boundary to lie between tolsat and tolsst.
    tolaft = tolsst - aftarg

    ! Setscrew parameters.
    !   sscrew(1) = Step size control parameter, bounds the absolute
    !                 value of the estimated error in any Taylor's
    !                 series for a basis species variable. This is used
    !                 to choose the order and the step size in normal
    !                 reaction path computational mode. It is used to
    !                 choose only the order in economy mode.
    !   sscrew(2) = Not used.
    !   sscrew(3) = Step size control parameter, bounds the absolute
    !                 value of the estimated error in Taylor's
    !                 series for rate functions (kinetic mode only).
    !                 It also serves a function similar to that of
    !                 sscrew(4) in testing the estimated error in the
    !                 absolute time or the reaction progress of any
    !                 individual irreversible reaction. If either of
    !                 the sscrew(3) or sscrew(4) tests for corrector
    !                 iteration or step size cuts is satisfied, no
    !                 corrective action is taken.
    !   sscrew(4) = Tolerance for corrector iteration or step size
    !                 cut for differences between actual and predicted
    !                 rate function values (kinetic mode only).
    !   sscrew(5) = Maximum magnitude of Newton-Raphson correction
    !                 term per iteration.
    !   sscrew(6) = Step size parameter for economy mode; this bounds
    !                 the change in a basis variable.
    do n = 1,nsscmx
        sscrew(n) = 0.
    end do

    sscrew(1) = 1.e-4

    ! if (sscrew(1) .lt. 1.e-6) sscrew(1) = 1.e-6
    ! if (sscrew(1) .gt. 1.e-4) sscrew(1) = 1.e-4
    sscrew(3) = 1.e-4

    ! if (sscrew(3) .lt. 1.e-6) sscrew(3) = 1.e-6
    ! if (sscrew(3) .gt. 1.e-4) sscrew(3) = 1.e-4
    sscrew(4) = 1.e-6

    ! if (sscrew(4) .lt. 1.e-8) sscrew(4) = 1.e-8
    ! if (sscrew(4) .gt. 1.e-4) sscrew(4) = 1.e-4
    sscrew(5) = 4.0

    ! if (sscrew(5) .lt. 1.0) sscrew(5) = 1.0
    ! if (sscrew(5) .gt. 6.0) sscrew(5) = 6.0
    sscrew(6) = 4.0

    ! if (sscrew(6) .lt. 0.5) sscrew(6) = 0.5
    ! if (sscrew(6) .gt. 1.e+38) sscrew(6) = 1.e+38
    ! Z parameters.
    !   zklogu = Threshhold and/or target value for the log number
    !              of moles of a phase; it is used in the following
    !              ways:
    !              (1) When the log number of moles of a phase is less
    !              than this value, EQ6 no longer tries to limit the
    !              step size to force the corresponding Taylor's
    !              series to meet the sscrew(1) accuracy criterion.
    !              try to keep the corresponding taylor's series
    !              accurate according to the sscrew(1) criterion.
    !              (2) This is the target value employed for the
    !              log number of moles when using Taylor's series
    !              to predict a phase disappearance boundary.
    !              (3) In the fluid-centered, flow-through open
    !              system model (iopt(1) = 2), this also determines
    !              (in conjunction with the zkfac parameter, see
    !              below) the maximum number of moles of a phase which
    !              is permitted to re-dissolve. It limits the number
    !              of moles which can be moved from the ES to the PRS
    !              (see zkfac below).
    !              Recommended values lie in the range -4.0 to -10.0.
    !   zklogl = The amount by which the log number of moles of a phase
    !              is decremented by a transfer from the ES to the PRS
    !              in the fluid-centered, flow-through open system
    !              model (iopt(1) = 2). A value of 2.0 means that 99%
    !              of the number of moles is transferred, a value of
    !              3.0, 99.9%. Recommended that values lie in the range
    !              of 2.0 to 4.0.
    !   zkfac  = Parameter which determines the minimum number of moles
    !              of a phase left in the ES after a shift to the PRS
    !              in the fluid-centered, flow-through open system
    !              model (iopt(1) = 2). This minimum value is calculated
    !              as zkfac*10**zklogu (equivalent to 10**zklgmn, see
    !              below). Recommended values lie in the range 0.9
    !              to 0.98.
    !   zklgmn = Secondary parameter equivalent to the minimum log
    !              number of moles of a phase left in the ES after a
    !              shift to the PRS in the fluid-centered, flow-through
    !              open system.
    zklogu = -7.0
    zklogl = 2.0
    zkfac = 0.98
    zklgmn = tlg(zkfac) + zklogu

    ! The minimum step size. No mechanism for cutting the step size will
    ! reduce it to less than this value. Typically, the minimum step
    ! size is less than the zero-order step size. It is undesirable to
    ! let the minimum step size get too big, as this determines the
    ! resolution of certain events of interest (such as hitting a
    ! desired pH value at the end of the run). The resolution of events
    ! in reaction progress space and time space is not limited by this
    ! parameter.
    dlxmin = 0.01*dlxmx0
    dlxmin = min(dlxmin,1.e-14_dp)

    ! Upper limit on the step size.
    dlxmax = prcinf

    ! Maximum number of attempts to slide over a phase boundary.
    npslmx = 8

    ! Maximum number of attempts to slide over a redox jump.
    nsslmx = 8
end subroutine setrcp
