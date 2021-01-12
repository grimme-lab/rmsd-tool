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

!> Implementation of a filter type to store selection rules for molecular
!> structure components based on atomic numbers and element symbols.
module rmsd_filter
   use mctc_env_error, only : error_type, fatal_error
   use mctc_io_structure, only : structure_type
   use mctc_io_symbols, only : symbol_length, to_number
   use rmsd_toml, only : toml_table, toml_array, toml_key, get_value, toml_stat, &
      & len
   implicit none
   private

   public :: rmsd_filter_type, new_rmsd_filter


   !> Filter to allow selecting of atoms from a molecular structure type
   type :: rmsd_filter_type

      !> Name of the filter
      character(len=:), allocatable :: name

      !> Atomic number
      integer, allocatable :: num(:)

      !> Element symbols
      character(symbol_length), allocatable :: sym(:)

      !> Whether the filter contains an allow or a deny list
      logical :: allow

      !> Use PDB identifiers for filtering
      logical :: pdb

   contains

      !> Create mask from filter and structure
      procedure :: get_mask

   end type rmsd_filter_type


   !> Overloaded constructor for RMSD filters
   interface new_rmsd_filter
      module procedure :: new_rmsd_filter_all
      module procedure :: new_rmsd_filter_one
      module procedure :: new_rmsd_filter_tbl
   end interface new_rmsd_filter


contains


!> Create a list of RMSD filters from a TOML data structure
subroutine new_rmsd_filter_all(self, table, error)

   !> List of new species filters
   type(rmsd_filter_type), allocatable, intent(out) :: self(:)

   !> TOML data structure
   type(toml_table), intent(inout) :: table

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_key), allocatable :: list(:)
   type(toml_table), pointer :: child
   integer :: i

   call table%get_keys(list)
   allocate(self(size(list)))

   do i = 1, size(list)
      call get_value(table, list(i)%key, child)
      if (.not.associated(child)) then
         call fatal_error(error, "Entry '"//list(i)%key//"' is not a table")
         exit
      end if

      call new_rmsd_filter(self(i), child, error)
      if (allocated(error)) exit
   end do
   if (allocated(error)) return

end subroutine new_rmsd_filter_all


!> Create a new RMSD filter from a TOML data structure
subroutine new_rmsd_filter_tbl(self, table, error)

   !> Instance of species filter
   type(rmsd_filter_type), intent(out) :: self

   !> TOML data structure
   type(toml_table), intent(inout) :: table

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(toml_array), pointer :: array
   integer, allocatable :: num(:)
   character(symbol_length), allocatable :: sym(:)
   character(len=:), allocatable :: name, cval
   integer :: ndim, inum, isym, ival, i, stat
   logical :: allow, pdb

   call table%get_key(name)

   call get_value(table, "pdb", pdb, .false.)
   call get_value(table, "include", array, requested=.false.)
   allow = associated(array)
   if (.not.allow) then
      call get_value(table, "exclude", array, requested=.false.)
      if (.not.associated(array)) then
         call fatal_error(error, "Filter '"//name//"' requires either an include or an exclude list")
         return
      end if
   end if
   ndim = len(array)
   allocate(num(ndim))
   allocate(sym(ndim))
   isym = 0
   inum = 0
   do i = 1, ndim
      call get_value(array, i, ival, stat=stat)
      if (stat == toml_stat%success) then
         if (ival > 0) then
            inum = inum + 1
            num(inum) = ival
            cycle
         else
            call fatal_error(error, "Atomic numbers must be larger than zero")
            exit
         end if
      endif
      call get_value(array, i, cval, stat=stat)
      if (stat == toml_stat%success) then
         if (to_number(cval) > 0) then
            isym = isym + 1
            sym(isym) = cval
            cycle
         else
            call fatal_error(error, "Unknown element symbol '"//cval//"'")
            exit
         end if
      end if
      call fatal_error(error, "Type mismatch for '"//name//"' filter")
      exit
   end do
   if (allocated(error)) return

   call new_rmsd_filter(self, name=name, num=num(:inum), sym=sym(:isym), &
      & allow=allow, pdb=pdb)

end subroutine new_rmsd_filter_tbl


!> Create a new RMSD filter from parts
subroutine new_rmsd_filter_one(self, name, num, sym, allow, pdb)

   !> Instance of species filter
   type(rmsd_filter_type), intent(out) :: self

   !> Name of the filter
   character(len=*), intent(in) :: name

   !> Whether the filter contains an allow or a deny list
   logical, intent(in) :: allow

   !> Filter is specific for PDB identifiers
   logical, intent(in) :: pdb

   !> Atomic number
   integer, intent(in) :: num(:)

   !> Element symbols
   character(symbol_length), intent(in) :: sym(:)

   self%name = name
   self%sym = sym
   self%num = num
   self%allow = allow
   self%pdb = pdb

end subroutine new_rmsd_filter_one


!> Create a logical mask for a given molecular structure type
subroutine get_mask(self, mol, mask)

   !> Instance of species filter
   class(rmsd_filter_type), intent(in) :: self

   !> Instance of the molecular structure data
   class(structure_type), intent(in) :: mol

   !> Filter make
   logical, allocatable, intent(out) :: mask(:)

   logical, allocatable :: tmp(:)
   integer :: iat, izp

   allocate(tmp(mol%nid))

   if (self%allow) then
      do izp = 1, mol%nid
         tmp(izp) = any(mol%num(izp) == self%num) &
            &  .or. any(mol%sym(izp) == self%sym)
      end do
   else
      do izp = 1, mol%nid
         tmp(izp) = all(mol%num(izp) /= self%num) &
            & .and. all(mol%sym(izp) /= self%sym)
      end do
   end if

   allocate(mask(mol%nat))
   do iat = 1, mol%nat
      izp = mol%id(iat)
      mask(iat) = tmp(izp)
   end do

   if (self%pdb .and. allocated(mol%pdb)) then
      if (self%allow) then
         do iat = 1, mol%nat
            mask(iat) = mask(iat) .or. any(mol%pdb(iat)%name == self%sym)
         end do
      else
         do iat = 1, mol%nat
            mask(iat) = mask(iat) .and. all(mol%pdb(iat)%name /= self%sym)
         end do
      end if
   end if
      
end subroutine get_mask


end module rmsd_filter
