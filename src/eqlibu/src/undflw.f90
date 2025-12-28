subroutine undflw
    !! This subroutine disables underflow trapping. This is not
    !! necessary on most modern platform/OS/compiler combinations, as
    !! they do not engage in underflow trapping by default at the
    !! system level. More commonly, underflow trapping is now a compiler
    !! option (which should be turned off in the case of EQ3/6).
    !! On some (mainly older) systems, underflow trapping is built into
    !! the system in a way that can only be overridden for a given
    !! program by putting a platform-dependent system subroutine call
    !! in that program. If such a call is necessary, it should be made
    !! in the present subroutine. If it appears that such a call is
    !! necessary, consult your local system documentation, compiler
    !! documentation, compiler vendor, or local computer support person
    !! to find out what this call should be. Do not call LLNL to find
    !! out.
    !! This subroutine is called by:
    !!   EQPT/eqpt.f
    !!   EQ3NR/eq3nr.f
    !!   EQ6/eq6.f
    !! Input:
    !!   None
    !! Output:
    !!   None
    implicit none

    ! om   BEGIN_MACHINE_DEPENDENT_CODE
    ! om
    ! om     If a system subroutine call is necessary to disable underflow
    ! om     trapping, add it here.
    ! om
    ! om
    ! om   END_MACHINE_DEPENDENT_CODE
end subroutine undflw