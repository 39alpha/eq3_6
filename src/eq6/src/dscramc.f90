subroutine dscramc(ntabx,ntabs,nllnmx,ulinex)
    !! This routine descrambles a file of tables whose lines are
    !! interspersed, but which are marked by an id string that
    !! apears first on any such line. The id string may contain up
    !! to 8 characters and is followed by a comma.
    !! The present routine is intended to support the creation of
    !! .csv file to be opened in Excel or other spreadsheet or plotting
    !! software, rather than a simple text file (such as the original
    !! EQ6 TAB file).
    !! The id string conceptually consists of eight characters.
    !! Trailing blanks may be omitted on a line read by this
    !! subroutine. The id string starts with a six-character tag
    !! that identifies the corresponding table. The remaining two
    !! characters may be used for other purposes. The string "vh'
    !! identifies a variable header. The id string for a "Table A"
    !! might be "A       ".
    !! The contents of the scrambled file are copied to the descrambled
    !! file as the descrambling takes place. The descrambled file must
    !! already be open.
    !! This subroutine is called by:
    !!   EQ6/eq6.f
    !! Input:
    !!   ntabx  = unit number of the scrambled file
    !!   ntabs  = unit number of the descrambled file (here the
    !!            TAB scratch file)
    !!   nllnmx = the maximum character length of a line on the
    !!            TABX, TABS, and TAB files
    !!   ulinex = variable holding a line, including the id string
    !! Output:
    !!   None
    implicit none

    ! Calling sequence variable declarations.
    integer :: ntabx
    integer :: ntabs

    integer :: nllnmx

    character(len=nllnmx) :: ulinex

    ! Local variable declarations.
    integer :: j1
    integer :: j2
    integer :: j
    integer :: jc
    integer :: n
    integer :: ntable
    integer :: ntablx

    integer, dimension(:), allocatable :: iftabl
    integer, dimension(:), allocatable :: iltabl

    logical :: qtablf
    logical :: qvhfix

    character(len=nllnmx), dimension(:), allocatable :: ulvhtb

    character(len=8), dimension(:), allocatable :: utable
    character(len=8), dimension(:), allocatable :: utablx

    character(len=8) :: uidstr
    character(len=8) :: uidstx
    character(len=8) :: utx
    character(len=8) :: ux8

    ! Determine the six-character table identifiers (utable)
    ! that are present on the scrambled file.
    ntablx = 10
    ALLOCATE(utable(ntablx))

    do n = 1,ntablx
        utable(n) = ' '
    end do

    rewind ntabx
    ntable = 0

    uidstr = ' '
    read (ntabx,'(a8)',end=110) ux8
    jc = index(ux8,',')

    if (jc .gt. 0) then
        uidstr(1:jc - 1) = ux8(1:jc - 1)
    else
        uidstr(1:8) = ux8(1:8)
    end if

    ntable = 1
    utable(1)(1:6) = uidstr(1:6)

100 continue
    uidstr = ' '
    read (ntabx,'(a8)',end=110) ux8
    jc = index(ux8,',')

    if (jc .gt. 0) then
        uidstr(1:jc - 1) = ux8(1:jc - 1)
    else
        uidstr(1:8) = ux8(1:8)
    end if

    ! Is the current line's id string in the list?
    do n = 1,ntable
        if (uidstr(1:6) .eq. utable(n)(1:6)) then
            go to 100
        end if
    end do

    ! Does the list of table identifiers need to be redimensioned?
    if (ntable .eq. ntablx) then
        ALLOCATE(utablx(ntablx))

        do n = 1,ntable
            utablx(n) = utable(n)
        end do

        do n = ntable + 1,ntablx
            utablx(n) = ' '
        end do

        DEALLOCATE(utable)
        ntablx = ntablx + 10
        ALLOCATE(utable(ntablx))

        do n = 1,ntable
            utable(n) = utablx(n)
        end do

        do n = ntable + 1,ntablx
            utable(n) = ' '
        end do

        DEALLOCATE(utablx)
    end if

    ! Add the new table identifier to the list.
    ntable = ntable + 1
    utable(ntable)(1:6) = uidstr(1:6)
    go to 100

110 continue

    ! Organize the scrambled file so that all the lines belonging
    ! to a table are contiguous. Descramble to a scratch file (TABS).
    ! In the process, find the last instance of a variable header
    ! in each table, if there is one. A variable header instance is
    ! marked by the string 'vh' in the last two characters of the
    ! id string.
    ALLOCATE(ulvhtb(ntable))

    do n = 1,ntable
        ulvhtb(n) = 'None'
    end do

    rewind ntabs

    do n = 1,ntable
        ! Loop over table designators.
        utx = utable(n)
        rewind ntabx

140 continue
        uidstr = ' '
        read (ntabx,'(a)',end=150) ulinex

        jc = index(ulinex,',')

        if (jc .gt. 0) then
            uidstr(1:jc - 1) = ulinex(1:jc - 1)
        else
            uidstr(1:8) = ulinex(1:8)
        end if

        if (uidstr(1:6) .eq. utx(1:6)) then
            j2 = len_trim(ulinex)
            write (ntabs,'(a)') ulinex(1:j2)

            if (uidstr(7:8) .eq. 'vh') then
                ux8 = uidstr
                ux8(7:8) = '  '
                j1 = len_trim(ux8)
                jc = index(ulinex,",")
                ulvhtb(n)(1:j1) = ux8(1:j1)
                ulvhtb(n)(j1 + 1:nllnmx) = ulinex(jc:nllnmx)
            end if
        end if

        go to 140
150 continue
    end do

    ! Copy the partially descrambled file back over the
    ! original descrambled file. After this has been done,
    ! the scrambled file will still be suitable for extension
    ! in a possible subsequent EQ6 run.
    rewind ntabx
    rewind ntabs

170 continue
    read (ntabs,'(a)',end=180) ulinex
    j2 = len_trim(ulinex)
    write (ntabx,'(a)') ulinex(1:j2)
    go to 170

180 continue

    ! The last instance of a possible variable header line
    ! has already been found. It is the only instance that
    ! is desired, and it will be placed in the relative position
    ! of the first instance on the descrambled file.
    ! A variable header line defines quantities associated with
    ! some of the table columns, but whose contents are built up
    ! during the run and cannot be known in advance. For example,
    ! a table might give the number of moles of precipitated
    ! minerals. The identities of these minerals cannot be known
    ! until the run is complete. A first instance is written
    ! in the desired relative position for the full variable header.
    ! This will usually be incomplete. As more minerals are
    ! precipitated, updated instances of the variable header will
    ! be written. Only the last one is complete (as far as the
    ! current run is concerned).
    ! Loop over the tables.
    rewind ntabs

    do n = 1,ntable
        ! Loop over table designators.
        qtablf = .false.
        qvhfix = .false.
        utx = utable(n)
        rewind ntabx
240 continue
        uidstr = ' '
        read (ntabx,'(a)',end=250) ulinex

        jc = index(ulinex,',')

        if (jc .gt. 0) then
            uidstr(1:jc - 1) = ulinex(1:jc - 1)
        else
            uidstr(1:8) = ulinex(1:8)
        end if

        if (uidstr(1:6) .eq. utx(1:6)) then
            qtablf = .true.

            if (uidstr(7:8) .eq. 'vh') then
                if (.not.qvhfix) then
                    qvhfix = .true.
                    ulinex = ulvhtb(n)
                else
                    go to 240
                end if
            end if

            j2 = len_trim(ulinex)
            write (ntabs,'(a)') ulinex(1:j2)
        else if (qtablf) then
            go to 250
        end if

        go to 240
250 continue
    end do

    DEALLOCATE(utable)
    DEALLOCATE(ulvhtb)
end subroutine dscramc