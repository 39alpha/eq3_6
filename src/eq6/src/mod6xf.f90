module mod6xf
    !! This module contains standard state thermodynamic data and
    !! standard grid pressure data. This excludes interpolating
    !! polynomial coefficient arrays belonging to the 'a' set (e.g.,
    !! axhfsa, axlksa, axvfsa) which includes all such data read
    !! from the data file. Included here are the "compressed"
    !! equivalents.,
    !! Pressure parameters.
    real(kind=8) :: prehw
    real(kind=8) :: presg

    real(kind=8), dimension(:,:), allocatable :: aprehw
    real(kind=8), dimension(:,:), allocatable :: apresg

    ! "Eh reaction" parameters.
    real(kind=8) :: xhfe
    real(kind=8) :: xlke
    real(kind=8) :: xvfe

    real(kind=8), dimension(:), allocatable :: dhfe
    real(kind=8), dimension(:), allocatable :: dvfe
    real(kind=8), dimension(:,:), allocatable :: axlke
    real(kind=8), dimension(:,:), allocatable :: axhfe
    real(kind=8), dimension(:,:), allocatable :: axvfe
    real(kind=8), dimension(:,:,:), allocatable :: adhfe
    real(kind=8), dimension(:,:,:), allocatable :: advfe

    ! Regular reaction parameters.
    real(kind=8), dimension(:), allocatable :: xhfs
    real(kind=8), dimension(:), allocatable :: xhfsd
    real(kind=8), dimension(:), allocatable :: xlks
    real(kind=8), dimension(:), allocatable :: xlksd
    real(kind=8), dimension(:), allocatable :: xvfs
    real(kind=8), dimension(:), allocatable :: xvfsd

    real(kind=8), dimension(:,:), allocatable :: dhfs
    real(kind=8), dimension(:,:), allocatable :: dhfsd
    real(kind=8), dimension(:,:), allocatable :: dvfs
    real(kind=8), dimension(:,:), allocatable :: dvfsd

    real(kind=8), dimension(:,:,:), allocatable :: axhfs
    real(kind=8), dimension(:,:,:), allocatable :: axhfsd
    real(kind=8), dimension(:,:,:), allocatable :: axhfsx
    real(kind=8), dimension(:,:,:), allocatable :: axlks
    real(kind=8), dimension(:,:,:), allocatable :: axlksd
    real(kind=8), dimension(:,:,:), allocatable :: axlksx
    real(kind=8), dimension(:,:,:), allocatable :: axvfs
    real(kind=8), dimension(:,:,:), allocatable :: axvfsd
    real(kind=8), dimension(:,:,:), allocatable :: axvfsx

    real(kind=8), dimension(:,:,:,:), allocatable :: adhfs
    real(kind=8), dimension(:,:,:,:), allocatable :: adhfsd
    real(kind=8), dimension(:,:,:,:), allocatable :: adhfsx
    real(kind=8), dimension(:,:,:,:), allocatable :: advfs
    real(kind=8), dimension(:,:,:,:), allocatable :: advfsd
    real(kind=8), dimension(:,:,:,:), allocatable :: advfsx
end module mod6xf