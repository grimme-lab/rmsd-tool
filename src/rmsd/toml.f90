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

!> Proxy interface to TOML-Fortran library
module rmsd_toml
   use mctc_env, only : error_type, fatal_error
   use tomlf, only : toml_table, toml_array, toml_key, toml_stat, toml_error, &
      & toml_parse, toml_serializer, add_table, add_array, get_value, set_value, len
   implicit none
   private

   public :: read_config_file
   public :: toml_table, toml_array, toml_key, toml_stat, toml_serializer
   public :: add_table, add_array, get_value, set_value, len


contains


!> Process the configuration file to a TOML data structure
subroutine read_config_file(table, file, error)

   !> TOML data structure
   type(toml_table), allocatable, intent(out) :: table

   !> Name of the configuration file
   character(len=*), intent(in) :: file

   !> Error status of the operation
   type(error_type), allocatable, intent(out) :: error

   type(toml_error), allocatable :: parse_error
   integer :: unit
   logical :: exist

   inquire(file=file, exist=exist)

   if (.not.exist) then
      call fatal_error(error, "'"//file//"' not found")
      return
   end if

   open(file=file, newunit=unit)
   call toml_parse(table, unit, parse_error)
   close(unit)

   if (allocated(parse_error)) then
      allocate(error)
      call move_alloc(parse_error%message, error%message)
      return
   end if

end subroutine read_config_file


end module rmsd_toml
