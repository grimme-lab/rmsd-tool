! This file is part of mctc-rmsd.
! SPDX-Identifier: LGPL-3.0-or-later
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

!> Configuration data for calculations of RMSDs between structures
module rmsd_config
   use mctc_env_accuracy, only : wp
   use mctc_env_error, only : error_type, fatal_error
   use mctc_io_convert, only : autoaa
   use rmsd_filter, only : rmsd_filter_type, new_rmsd_filter
   use rmsd_toml, only : toml_table, toml_array, toml_key, add_table, add_array, &
      & get_value, set_value, toml_stat
   implicit none
   private

   public :: rmsd_config_type, new_rmsd_config


   !> Configuration data for RMSD calculations
   type :: rmsd_config_type

      !> Error on mismatching structures
      logical :: strict

      !> Conversion factor for output
      real(wp) :: conv

      !> Length unit for RMSD output
      character(len=:), allocatable :: length_unit

      !> Available RMSD filters
      type(rmsd_filter_type), allocatable :: filter(:)

   contains

      !> Select a filter from the configuration data
      procedure :: get_filter

   end type rmsd_config_type


contains


!> Create new configuration data from TOML data structure
subroutine new_rmsd_config(self, table, error)

   !> Instance of the configuration data
   type(rmsd_config_type), intent(out) :: self

   !> TOML data structure
   type(toml_table), intent(inout) :: table

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_table), pointer :: child, child2
   type(toml_array), pointer :: array
   integer :: stat

   call get_value(table, "strict", self%strict, .true., stat=stat)
   if (stat /= toml_stat%success) then
      call fatal_error(error, "Could not error policy value from strict entry")
      return
   end if

   call get_value(table, "unit", self%length_unit, "AA")
   if (.not.allocated(self%length_unit)) then
      call fatal_error(error, "Could not retrieve length unit")
      return
   end if

   call get_unit_conversion(self%conv, self%length_unit, error)
   if (allocated(error)) return

   call get_value(table, "filter", child, requested=.false.)
   if (associated(child)) then
      call new_rmsd_filter(self%filter, child, error)
      if (allocated(error)) return
   else
      call add_table(table, "filter", child)
      call add_table(child, "heavy", child2)
      call add_array(child2, "exclude", array)
      call set_value(array, 1, "H")
      call set_value(array, 2, "h")
      call new_rmsd_filter(self%filter, child, error)
      if (allocated(error)) return
   end if

end subroutine new_rmsd_config


!> Get unit conversion factor from string
subroutine get_unit_conversion(conv, length_unit, error)

   !> Conversion factor
   real(wp), intent(out) :: conv

   !> Length unit for RMSD output
   character(len=*), intent(in) :: length_unit

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   select case(length_unit)
   case default
      call fatal_error(error, "Unknown unit '"//length_unit//"' encountered")
   case("AA", "Angstrom")
      conv = autoaa
   case("a0", "Bohr")
      conv = 1.0_wp
   end select

end subroutine get_unit_conversion


!> Retrieve filter from configuration data
subroutine get_filter(self, name, filter)

   !> Instance of the configuration data
   class(rmsd_config_type), intent(in), target :: self

   !> Name of the filter
   character(len=*), intent(in) :: name

   !> Filter to use
   type(rmsd_filter_type), pointer :: filter

   integer :: i

   nullify(filter)
   if (.not.allocated(self%filter)) return

   do i = 1, size(self%filter)
      if (name == self%filter(i)%name) then
         filter => self%filter(i)
         exit
      end if
   end do

end subroutine get_filter


end module rmsd_config
