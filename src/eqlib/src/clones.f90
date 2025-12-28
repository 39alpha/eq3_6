subroutine clones(axlksa,cdrsa,cessa,mwtspa,narx_asv,ndrsa,ndrsa_asv,ndrsn,ndrsra,nessa,nessa_asv,nessn,nessra,np,npta_asv,ns,nsc,nsta_asv,ntpr_asv,uphasa,uspeca,zchara)
    !! This subroutine clones the nsc-th species into the ns-th. The
    !! ns-th species belongs to the np-th phase. Typically, the nsc-th
    !! species is a pure mineral or liquid and the ns-th species is the
    !! corresponding component species of a solid or liquid solution.
    !! This subroutine is called by:
    !!   EQLIB/indata.f
    !!   EQLIB/indatp.f
    !! Principal input:
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: narx_asv
    integer :: ndrsa_asv
    integer :: nessa_asv
    integer :: npta_asv
    integer :: nsta_asv
    integer :: ntpr_asv

    integer :: ndrsa(ndrsa_asv)
    integer :: ndrsra(2,nsta_asv)
    integer :: nessa(nessa_asv)
    integer :: nessra(2,nsta_asv)

    integer :: ndrsn
    integer :: nessn
    integer :: np
    integer :: ns
    integer :: nsc

    character(len=24) :: uphasa(npta_asv)
    character(len=48) :: uspeca(nsta_asv)

    real(kind=8) :: axlksa(narx_asv,ntpr_asv,nsta_asv)
    real(kind=8) :: cessa(nessa_asv)
    real(kind=8) :: cdrsa(ndrsa_asv)
    real(kind=8) :: mwtspa(nsta_asv)
    real(kind=8) :: zchara(nsta_asv)

    ! Local variable declarations.
    integer :: i
    integer :: j
    integer :: ndrsn1
    integer :: nessn1
    integer :: nr1
    integer :: nr2
    integer :: nse
    integer :: nt

    ! Copy the name, molecular weight, and electrical charge.
    ! Note that the phase part of the name is different.
    uspeca(ns)(1:24) = uspeca(nsc)(1:24)
    uspeca(ns)(25:48) = uphasa(np)(1:24)
    mwtspa(ns) = mwtspa(nsc)
    zchara(ns) = zchara(nsc)

    ! Copy the elemental composition.
    nessra(1,ns) = nessn + 1
    nr1 = nessra(1,nsc)
    nr2 = nessra(2,nsc)

    do nessn1 = nr1,nr2
        nessn = nessn + 1
        cessa(nessn) = cessa(nessn1)
        nessa(nessn) = nessa(nessn1)
    end do

    nessra(2,ns) = nessn

    ! Copy the associated reaction.
    ndrsra(1,ns) = ndrsn + 1
    nr1 = ndrsra(1,nsc)
    nr2 = ndrsra(2,nsc)
    nt = nr2 - nr1 + 1

    if (nt .lt. 2) then
        ! The nsc-th species is a data file basis species that
        ! has only an identity reaction. Write a linking reaction.
        ndrsn = ndrsn + 1
        cdrsa(ndrsn) = -1.
        ndrsa(ndrsn) = ns
        ndrsn = ndrsn + 1
        cdrsa(ndrsn) = 1.
        ndrsa(ndrsn) = nsc
    else
        ! The nsc-th species has a reaction which can be copied.
        do ndrsn1 = nr1,nr2
            ndrsn = ndrsn + 1
            cdrsa(ndrsn) = cdrsa(ndrsn1)
            nse = ndrsa(ndrsn1)

            if (nse .eq. nsc) then
                nse = ns
            end if

            ndrsa(ndrsn) = nse
        end do
    end if

    ndrsra(2,ns) = ndrsn

    ! Copy the coefficients of the polynomials which represent
    ! the thermodynamic functions of the reaction. Note that
    ! this is valid whether the reaction is an identity reaction
    ! or a real one.
    do j = 1,ntpr_asv
        do i = 1,narx_asv
            axlksa(i,j,ns) = axlksa(i,j,nsc)
        end do
    end do
end subroutine clones