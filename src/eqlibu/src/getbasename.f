      subroutine getbasename(filename, indices)
            implicit none
            character(*), intent(in) :: filename
            integer, intent(out) :: indices(2)
            integer :: n, i, extindex, slashindex

            n = len(filename)
            slashindex = 0
            extindex = n + 1
            do i = 1, n
                if (filename(i:i).eq.'/') then
                    slashindex = i
                else if (filename(i:i).eq.'.') then
                    extindex = i
                end if
            end do

            indices(1) = slashindex + 1
            indices(2) = extindex - 1
      end subroutine
