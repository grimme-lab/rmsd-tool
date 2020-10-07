! This file is part of mctc-rmsd.
!
! mctc-rmsd is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! mctc-rmsd is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with mctc-rmsd.  If not, see <https://www.gnu.org/licenses/>.

module test_rmsd_filter
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check
   use rmsd
   use rmsd_filter
   implicit none
   private

   public :: collect_rmsd_filter


contains


!> Collect all exported unit tests
subroutine collect_rmsd_filter(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      & new_unittest("valid1-filter", test_filter1), &
      & new_unittest("valid2-filter", test_filter2), &
      & new_unittest("valid3-filter", test_filter3) &
      & ]

end subroutine collect_rmsd_filter


subroutine test_filter1(error)
   use rmsd_toml

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(rmsd_filter_type) :: filter
   type(toml_table) :: table
   type(toml_array), pointer :: array
   integer :: i, stat

   table = toml_table()
   table%key = "test"
   call add_array(table, "include", array, stat)
   do i = 1, 8
      call set_value(array, i, i+2, stat)
   end do

   call new_rmsd_filter(filter, table, error)
   if (allocated(error)) return

   call check(error, size(filter%num), 8)
   if (allocated(error)) return
   call check(error, size(filter%sym), 0)
   if (allocated(error)) return

end subroutine test_filter1


subroutine test_filter2(error)
   use mctc_io_symbols, only : to_symbol
   use rmsd_toml

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(rmsd_filter_type) :: filter
   type(toml_table) :: table
   type(toml_array), pointer :: array
   integer :: i, stat

   table = toml_table()
   table%key = "test"
   call add_array(table, "exclude", array, stat)
   do i = 1, 16
      call set_value(array, i, trim(to_symbol(i+2)), stat)
   end do

   call new_rmsd_filter(filter, table, error)
   if (allocated(error)) return

   call check(error, size(filter%num), 0)
   if (allocated(error)) return
   call check(error, size(filter%sym), 16)
   if (allocated(error)) return

end subroutine test_filter2


subroutine test_filter3(error)
   use mctc_io_symbols, only : to_symbol
   use rmsd_toml

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(rmsd_filter_type) :: filter
   type(toml_table) :: table
   type(toml_array), pointer :: array
   integer :: i, stat

   table = toml_table()
   table%key = "test"
   call add_array(table, "include", array, stat)
   do i = 1, 20
      if (mod(i, 2) == 0) then
         call set_value(array, i, i, stat)
      else
         call set_value(array, i, trim(to_symbol(i)), stat)
      end if
   end do

   call new_rmsd_filter(filter, table, error)
   if (allocated(error)) return

   call check(error, size(filter%num), 10)
   if (allocated(error)) return
   call check(error, size(filter%sym), 10)
   if (allocated(error)) return

end subroutine test_filter3


end module test_rmsd_filter
