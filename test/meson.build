# This file is part of mctc-rmsd.
# SPDX-Identifier: LGPL-3.0-or-later
#
# mctc-rmsd is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# mctc-rmsd is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with mctc-rmsd.  If not, see <https://www.gnu.org/licenses/>.

tests = [
  'rmsd',
  'rmsd-filter',
]

test_srcs = files(
  'main.f90',
)
foreach t : tests
  test_srcs += files('test_@0@.f90'.format(t.underscorify()))
endforeach

tester = executable(
  'tester',
  sources: test_srcs,
  dependencies: rmsd_dep,
)

foreach t : tests
  test(t, tester, args: t)
endforeach
