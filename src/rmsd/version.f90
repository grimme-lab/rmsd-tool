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

!> Versioning information on this library.
module rmsd_version
   implicit none
   private

   public :: rmsd_version_string, rmsd_version_compact
   public :: get_rmsd_version


   !> String representation of the mctc-rmsd version
   character(len=*), parameter :: rmsd_version_string = "0.1.2"

   !> Numeric representation of the mctc-rmsd version
   integer, parameter :: rmsd_version_compact(3) = [0, 1, 2]


contains


!> Getter function to retrieve mctc-rmsd version
subroutine get_rmsd_version(major, minor, patch, string)

   !> Major version number of the mctc-rmsd version
   integer, intent(out), optional :: major

   !> Minor version number of the mctc-rmsd version
   integer, intent(out), optional :: minor

   !> Patch version number of the mctc-rmsd version
   integer, intent(out), optional :: patch

   !> String representation of the mctc-rmsd version
   character(len=:), allocatable, intent(out), optional :: string

   if (present(major)) then
      major = rmsd_version_compact(1)
   end if
   if (present(minor)) then
      minor = rmsd_version_compact(2)
   end if
   if (present(patch)) then
      patch = rmsd_version_compact(3)
   end if
   if (present(string)) then
      string = rmsd_version_string
   end if

end subroutine get_rmsd_version


end module rmsd_version
