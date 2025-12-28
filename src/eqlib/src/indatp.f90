subroutine indatp(apxa,axlksa,bpxa,cdrsa,cessa,iapxa_asv,iapxta,ibpxa_asv,ibpxta,ikta_asv,jsola,mwtspa,nad1,narx_asv,ncmpra,ndrsa,ndrsa_asv,ndrsn,ndrsra,nerr,nessa,nessa_asv,nessn,nessra,nmrn1a,nmrn2a,noutpt,np,npta_asv,ns,nsta_asv,ntpr_asv,nttyo,nxta,nxta_asv,qclnsa,uendit,uspeca,uphasa,uptsld,uptypa,zchara)
    !! This subroutine reads the solid solution superblock from the
    !! supporting data file "data1".
    !! This subroutine is called by:
    !!   EQLIB/indata.f
    !! Principal input:
    !!   nad1   = unit number of the data file
    !! Principal output:
    implicit none

    ! Calling sequence variable declarations.
    integer :: iapxa_asv
    integer :: ibpxa_asv
    integer :: ikta_asv
    integer :: narx_asv
    integer :: ndrsa_asv
    integer :: nessa_asv
    integer :: npta_asv
    integer :: nsta_asv
    integer :: ntpr_asv
    integer :: nxta_asv

    integer :: noutpt
    integer :: nttyo

    integer :: iapxta(nxta_asv)
    integer :: ibpxta(nxta_asv)
    integer :: jsola(nxta_asv)
    integer :: ncmpra(2,npta_asv)
    integer :: ndrsa(ndrsa_asv)
    integer :: ndrsra(2,nsta_asv)
    integer :: nessa(nessa_asv)
    integer :: nessra(2,nsta_asv)

    integer :: nad1
    integer :: ndrsn
    integer :: nerr
    integer :: nessn
    integer :: nmrn1a
    integer :: nmrn2a
    integer :: np
    integer :: ns
    integer :: nxta

    logical :: qclnsa(nsta_asv)

    character(len=24) :: uphasa(npta_asv)
    character(len=24) :: uptypa(npta_asv)
    character(len=24) :: uptsld
    character(len=48) :: uspeca(nsta_asv)
    character(len=8) :: uendit

    real(kind=8) :: apxa(iapxa_asv,nxta_asv)
    real(kind=8) :: axlksa(narx_asv,ntpr_asv,nsta_asv)
    real(kind=8) :: bpxa(ibpxa_asv,nxta_asv)
    real(kind=8) :: cdrsa(ndrsa_asv)
    real(kind=8) :: cessa(nessa_asv)
    real(kind=8) :: mwtspa(nsta_asv)
    real(kind=8) :: zchara(nsta_asv)

    ! Local variable declarations.
    integer :: j2
    integer :: j3
    integer :: jsolv
    integer :: n
    integer :: ncmptv
    integer :: nsc
    integer :: nx

    integer :: ilnobl

    character(len=24) :: uphasv

    ! Local variable declarations with global dimensioning.
    ! The following need not be SAVEd.
    character(len=24), dimension(:), allocatable :: ucompv

    ! Allocate  a scratch array for the names of the end-member
    ! components.
    ALLOCATE(ucompv(ikta_asv))

    ! The label below defines a loop for reading the blocks in the
    ! superblock.
100 continue

    ! Read a species block.
    ! uphasv  = name of solid solution
    ! ncmptv  = number of component species
    ! jsolv   = code for solid solution activity coefficients
    read (nad1) uphasv,ncmptv,jsolv
    j3 = ilnobl(uphasv)

    ! Check for end of this current superblock.
    if (uphasv(1:8) .eq. uendit(1:8)) then
        go to 200
    end if

    np = np + 1
    uphasa(np) = uphasv
    uptypa(np) = uptsld

    nxta = nxta + 1
    nx = nxta
    jsola(nx) = jsolv

    if (ncmptv .gt. 0) then
        ! Read the names of component species.
        read (nad1) (ucompv(n), n = 1,ncmptv)
    else
        write (noutpt,1020) uphasv(1:j3)
        write (noutpt,1020) uphasv(1:j3)
1020 format(/' * Error - (EQLIB/indatp) The solid solution ',a,/7x,'has no components specified on the data file.')

        nerr = nerr + 1
    end if

    ! Read the activity coefficient parameters, first the ordinary
    ! parameters (apx), then the site-mixing parameters (bpx).
    read (nad1) iapxta(nx)

    if (iapxta(nx) .gt. 0) then
        read (nad1) (apxa(n,nx), n = 1,iapxta(nx))
    end if

    read (nad1) ibpxta(nx)

    if (ibpxta(nx) .gt. 0) then
        read (nad1) (bpxa(n,nx), n = 1,ibpxta(nx))
    end if

    ! Clone the component species.
    ncmpra(1,np) = ns + 1
    ncmpra(2,np) = ns + ncmptv

    do n = 1,ncmptv
        ns = ns + 1

        do nsc = nmrn1a,nmrn2a
            if (ucompv(n)(1:24) .eq. uspeca(nsc)(1:24)) then
                call clones(axlksa,cdrsa,cessa,mwtspa,narx_asv,ndrsa,ndrsa_asv,ndrsn,ndrsra,nessa,nessa_asv,nessn,nessra,np,npta_asv,ns,nsc,nsta_asv,ntpr_asv,uphasa,uspeca,zchara)
                qclnsa(ns) = .true.
                go to 130
            end if
        end do

        ! The listed component was not found as a pure mineral species.
        j2 = ilnobl(ucompv(n))
        write (noutpt,1040) ucompv(n)(1:j2),uphasv(1:j3)
        write (nttyo,1040) ucompv(n)(1:j2),uphasv(1:j3)
1040 format(/" * Error - (EQLIB/indatp) Couldn't find the pure",' phase equivalent',/7x,'of the species',/7x,a,' (',a,'),',' therefore',/7x,'could not create it by cloning.')

        nerr = nerr + 1
130 continue
    end do

    ! Loop back to read another species block.
    go to 100

    ! Deallocate  the scratch array.
200 continue
    DEALLOCATE(ucompv)

999 continue
end subroutine indatp