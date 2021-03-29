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

program rmsd_main
   use, intrinsic :: iso_fortran_env, only : output_unit, error_unit
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io, only : structure_type, read_structure
   use targ
   use rmsd
   use rmsd_toml
   use rmsd_config
   use rmsd_filter
   implicit none
   character(len=*), parameter :: prog_name = "mctc-rmsd"

   interface read_structure
      procedure :: read_structure_from_arg
   end interface read_structure

   integer :: iarg, stat
   type(arg_type), allocatable :: args(:)
   type(error_type), allocatable :: error
   type(structure_type) :: ref, mol
   type(rmsd_config_type) :: config
   type(rmsd_filter_type), pointer :: filter
   type(toml_table), allocatable :: table, opts
   type(toml_table), pointer :: child
   type(toml_serializer) :: ser
   character(len=:), allocatable :: rcfile, filter_name
   real(wp) :: rmsd_val
   type(targ_type) :: argparse
   logical :: show_version, show_help, show_rc
   logical, allocatable :: mask(:)

   argparse = new_argument_parser()
   call argparse%add_option("help")
   call argparse%add_option("version")
   call argparse%add_option("rc")
   call argparse%add_option("filter", require=1)
   call get_arguments(args)
   call get_options(argparse, args, opts)

   call get_value(opts, "filter", filter_name)
   call get_value(opts, "version", show_version, .false.)
   call get_value(opts, "rc", show_rc, .false.)
   call get_value(opts, "help", show_help, size(args) <= 1 .and. .not.show_rc)
   call opts%destroy

   if (show_version) then
      call version(output_unit)
      stop
   end if

   if (show_help) then
      call help(output_unit)
      if (size(args) <= 1) error stop
      stop
   end if

   call get_config_file(rcfile)
   if (allocated(rcfile)) then
      call read_config_file(table, rcfile, error)
      if (allocated(error)) then
         write(error_unit, '(a)') error%message
         error stop
      end if
   else
      table = toml_table()
   end if

   call get_value(table, "rmsd", child, requested=.true.)
   if (associated(child)) then
      call new_rmsd_config(config, child, error)
   else
      call fatal_error(error, "Type mismatch, rmsd must be a table")
   end if
   if (allocated(error)) then
      write(error_unit, '(a)') error%message
      error stop
   end if

   if (show_rc) then
      if (allocated(rcfile)) then
         write(output_unit, '(*(a:, 1x))') "#", "configuration file:", rcfile
      else
         write(output_unit, '(*(a:, 1x))') "#", "internal defaults"
      end if
      ser = toml_serializer(output_unit)
      call table%accept(ser)
      stop
   end if

   nullify(filter)
   if (allocated(filter_name)) then
      call config%get_filter(filter_name, filter)
      if (.not.associated(filter)) then
         write(error_unit, '(a)') &
            & "No filter with name '"//filter_name//"' present"
         error stop
      end if
   end if

   call read_structure(ref, args(1), error)
   if (allocated(error)) then
      write(error_unit, '(a)') error%message
      error stop
   end if

   if (associated(filter)) then
      call filter%get_mask(ref, mask)
   end if

   stat = 0
   do iarg = 2, size(args)
      call read_structure(mol, args(iarg), error)
      if (allocated(error)) then
         write(error_unit, '(a)') error%message
         stat = stat + 1
         cycle
      end if

      if (config%strict) then
         if (ref%nat /= mol%nat) then
            write(error_unit, '(a)') &
               "Number of atoms for '"//args(iarg)%arg//"' mismatch with reference"
            stat = stat + 1
            cycle
         end if
         if (any(ref%num(ref%id) /= mol%num(mol%id))) then
            write(error_unit, '(a)') &
               "Atomic number of '"//args(iarg)%arg//"' mismatch with reference"
            stat = stat + 1
            cycle
         end if
      end if

      if (allocated(error)) then
         stat = stat + 1
         cycle
      end if

      call get_rmsd(ref, mol, rmsd_val, mask=mask)
      write(output_unit, '(a, t40, es20.10, 1x, a)') &
         & args(iarg)%arg, rmsd_val*config%conv, config%length_unit
   end do

   if (stat > 0) then
      error stop
   end if


contains


subroutine help(unit)
   integer, intent(in) :: unit

   write(unit, '(a, *(1x, a))') &
      "Usage: "//prog_name//" [options] <file> <file>..."

   write(unit, '(a)') &
      "", &
      "Requires at least two input files. Supported formats are:", &
      "- xyz, mol, sdf, coord (0D), gen (C), pdb", &
      "", &
      "Configuration data is read from [rmsd] table in mctc.toml.", &
      "Just place the configuration file mctc.toml (or .mctc.toml) in your home directory.", &
      "Example:", &
      "", &
      "    [rmsd]", &
      "    unit = ""AA""", &
      "    [rmsd.filter]", &
      "    heavy.exclude = [ ""H"", ""h"" ]", &
      "", &
      "Options", &
      "-------", &
      ""

   write(unit, '(3x, a, t25, a)') &
      "--filter <name>", "Use <name> filter from configuration data to apply mask", &
      "--rc", "check configuration data and print it to standard out", &
      "--version", "Print program version and exit", &
      "--help", "Show this help message"

   write(unit, '(a)') &
      "", &
      "Filter", &
      "------", &
      "", &
      "Filters can be defined in the [rmsd.filter] section, they take a list of", &
      "atomic numbers and/or element symbols to define the allow-/deny-list.", &
      "For example, to only check all carbon, nitrogen and oxygen atoms create", &
      "a filter named organic with:", &
      "", &
      "    organic.include = [6, 7, 8]", &
      "", &
      "Similarly, to create a filter for all heavy elements, effectively just", &
      "excluding hydrogen with standard symbols, use:", &
      "", &
      "    heavy.exclude = [""H"", ""h""]", &
      "", &
      "Note that this approach will still consider deuterium labeled as D,", &
      "which would be excluded as well when using the atomic number instead.", &
      "", &
      "To create a PDB specific filter use the four character PDB identifier", &
      "of the atoms and enable the PDB functionality.", &
      "To match only the proteine backbone use", &
      "", &
      "    c-alpha.include = ["" CA "", "" N  "", "" C  "", "" O  ""]", &
      "    c-alpha.pdb = true", &
      "", &
      "Atomic numbers and element symbols can be included here as well.", &
      ""

end subroutine help


subroutine read_structure_from_arg(self, arg, error)
   type(structure_type), intent(out) :: self
   type(arg_type), intent(in) :: arg
   type(error_type), allocatable, intent(out) :: error
   call read_structure(self, arg%arg, error)
end subroutine read_structure_from_arg


subroutine get_variable(var, val)

   !> Name of variable
   character(len=*), intent(in) :: var

   !> Value of variable
   character(len=:), allocatable, intent(out) :: val

   integer :: length, stat

   call get_environment_variable(var, length=length, status=stat)
   if (stat /= 0) then
      return
   endif

   allocate(character(len=length) :: val, stat=stat)
   if (stat /= 0) then
      return
   endif

   if (length > 0) then
      call get_environment_variable(var, val, status=stat)
      if (stat /= 0) then
         deallocate(val)
         return
      end if
   end if

end subroutine get_variable


subroutine get_config_file(config)

   !> Name of the configuration file
   character(len=:), allocatable, intent(out) :: config

   character(len=*), parameter :: rc = 'mctc.toml'
   character(len=:), allocatable :: tmp, prefix
   character :: sep
   logical :: exist

   if (is_windows()) then
      sep = '\'
      call get_variable('HOMEDRIVE', prefix)
      call get_variable('HOMEDIR', tmp)
      prefix = prefix // tmp
   else
      sep = '/'
      call get_variable('HOME', prefix)
   end if
   if (allocated(prefix)) then
      tmp = prefix // sep // rc
      inquire(file=tmp, exist=exist)
      if (exist) then
         config = tmp
         return
      end if

      tmp = prefix // sep // '.' // rc
      inquire(file=tmp, exist=exist)
      if (exist) then
         config = tmp
         return
      end if
   end if

   inquire(file=rc, exist=exist)
   if (exist) then
      config = rc
      return
   end if

end subroutine get_config_file


!> Try to determine if we run on Windows and don't have POSIX compliance around
function is_windows() result(windows)

   !> Operating system seems to be Windows
   logical :: windows

   character(len=:), allocatable :: tmp

   windows = .false.
   call get_variable('OS', tmp)
   if (allocated(tmp)) then
      windows = index(tmp, 'Windows_NT') > 0
   end if
   if (.not.windows) then
      call get_variable('OSTYPE', tmp)
      if (allocated(tmp)) then
         windows = index(tmp, 'win') > 0 .or. index(tmp, 'msys') > 0
      end if
   end if

end function is_windows


subroutine version(unit)
   integer, intent(in) :: unit
   character(len=:), allocatable :: version_string
   call get_rmsd_version(string=version_string)
   write(unit, '(a, *(1x, a))') &
      & prog_name, "version", version_string
end subroutine version


end program rmsd_main
