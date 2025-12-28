module mod6pt
    !! This module contains data required to evaluate Pitzer's
    !! equations. It excludes the variables and arrays that belong
    !! to the 'a' set, which includes all the Pitzer data as read
    !! from the data file. Included here are the "compressed"
    !! equivalents.
    !! Integer scalars.
    integer :: napt
    integer :: nmut
    integer :: nslt

    ! xx    integer ifcphi1,ifcphi2,ifnnn,ifn2n,ifpsi1,ifpsi2,ifzeta,
    ! xx  $ ilcphi1,ilcphi2,ilnnn,iln2n,ilpsi1,ilpsi2,ilzeta
    !      Integer arrays. This are basically index arrays that allow
    !      for efficient storage of the data.
    integer, dimension(:), allocatable :: nalpha
    integer, dimension(:,:), allocatable :: nmux
    integer, dimension(:,:), allocatable :: nmxi
    integer, dimension(:,:), allocatable :: nmxx
    integer, dimension(:,:), allocatable :: nslx
    integer, dimension(:,:), allocatable :: nsxi
    integer, dimension(:,:), allocatable :: nsxx

    ! Coefficients for calculating empirical interaction coefficients
    ! as functions of temperature and pressure.
    real(kind=8), dimension(:,:), allocatable :: amu
    real(kind=8), dimension(:,:,:), allocatable :: aslm

    ! Empirical interaction coefficients and related functions at
    ! temperature and pressure.
    real(kind=8), dimension(:), allocatable :: pmu
    real(kind=8), dimension(:), allocatable :: pslm
    real(kind=8), dimension(:,:), allocatable :: dpslm
    real(kind=8), dimension(:,:), allocatable :: gpit
    real(kind=8), dimension(:,:), allocatable :: palpha
    real(kind=8), dimension(:,:), allocatable :: pslamn
    real(kind=8), dimension(:,:,:), allocatable :: dgpit

    ! Theoretical higher-order electrostatic interaction coefficients
    ! and related functions at temperature and pressure.
    real(kind=8), dimension(:,:), allocatable :: elam
    real(kind=8), dimension(:,:), allocatable :: pelm
    real(kind=8), dimension(:,:,:), allocatable :: delam
    real(kind=8), dimension(:,:,:), allocatable :: dpelm

    real(kind=8), dimension(:), allocatable :: selm
    real(kind=8), dimension(:,:), allocatable :: dselm
end module mod6pt