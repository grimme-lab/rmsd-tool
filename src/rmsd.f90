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

module rmsd
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use rmsd_ls, only : get_rmsd
   use rmsd_version, only : get_rmsd_version
   implicit none
   private

   public :: get_rmsd, get_rmsd_version

   interface get_rmsd
      module procedure :: get_rmsd_for_structure
   end interface get_rmsd


contains


pure subroutine get_rmsd_for_structure(struc1, struc2, rmsd, gradient, trafo, mask)

   !> Reference structure
   type(structure_type), intent(in) :: struc1

   !> Molecular structure to compare against
   type(structure_type), intent(in) :: struc2

   !> Root mean square deviation between the two structures
   real(wp), intent(out) :: rmsd

   !> Gradient of the RMSD w.r.t. the coordinates
   real(wp), intent(out), optional :: gradient(:, :)

   !> Rotation matrix between the two structures
   real(wp), intent(out), optional :: trafo(:, :)

   !> Atoms to include in the RMSD calculation
   logical, intent(in), optional :: mask(:)

   call get_rmsd(struc1%xyz, struc2%xyz, rmsd, gradient, trafo, mask)

end subroutine get_rmsd_for_structure


end module rmsd
