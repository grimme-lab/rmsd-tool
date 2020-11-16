# This file is part of mctc-rmsd.
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

project(
  'mctc-rmsd',
  'fortran',
  version: '0.1.0',
  license: 'LGPL-3.0-or-later',
  meson_version: '>=0.53',
  default_options: [
    'buildtype=debugoptimized',
    'default_library=both',
  ],
)
install = not (meson.is_subproject() and get_option('default_library') == 'static')

# General configuration information
inc_dirs = []
lib_deps = []
subdir('config')

# Collect source of the project
srcs = []
subdir('src')

# Library target
rmsd_lib = library(
  meson.project_name(),
  sources: srcs,
  version: meson.project_version(),
  dependencies: lib_deps,
  include_directories: inc_dirs,
  install: install,
)

# Export dependency for other projects and test suite
rmsd_inc = rmsd_lib.private_dir_include()
rmsd_dep = declare_dependency(
  link_with: rmsd_lib,
  include_directories: rmsd_inc,
  dependencies: lib_deps,
)

# Add applications
subdir('app')

# Package the license files
rmsd_lic = files(
  'COPYING',
  'COPYING.LESSER',
)

if install
  # Distribute the license files in share/licenses/<name>
  install_data(
    rmsd_lic,
    install_dir: get_option('datadir')/'licenses'/meson.project_name()
  )

  pkg = import('pkgconfig')
  pkg.generate(
    rmsd_lib,
    description: 'Structure comparison by root mean square deviation',
  )

  asciidoc = find_program('asciidoctor', required: false)
  if asciidoc.found()
    install_man(
      configure_file(
        command: [asciidoc, '-b', 'manpage', '@INPUT@', '-o', '@OUTPUT@'],
        input: files('man/mctc-rmsd.1.adoc'),
        output: '@BASENAME@',
      )
    )
  endif
endif

# add the testsuite
subdir('test')
