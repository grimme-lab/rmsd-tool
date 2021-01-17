# Root mean square deviation for structural comparison

[![LGPL-3.0](https://img.shields.io/github/license/grimme-lab/rmsd-tool)](COPYING.LESSER)
[![Release](https://img.shields.io/github/v/release/grimme-lab/rmsd-tool)](https://github.com/grimme-lab/rmsd-tool/releases/latest)
[![CI](https://github.com/grimme-lab/rmsd-tool/workflows/CI/badge.svg)](https://github.com/grimme-lab/rmsd-tool/actions)
[![docs](https://github.com/grimme-lab/rmsd-tool/workflows/docs/badge.svg)](https://grimme-lab.github.io/rmsd-tool)

## Installation

To build this project from the source code in this repository you need to have
- a Fortran compiler supporting Fortran 2008
- [meson](https://mesonbuild.com) version 0.53 or newer
- a build-system backend, *i.e.* [ninja](https://ninja-build.org) version 1.7 or newer

Setup a build with

```
meson setup _build
```

You can select the Fortran compiler by the `FC` environment variable, currently this project supports GCC and Intel compilers.
To compile the project run

```
meson compile -C _build
```

You can run the projects testsuite with

```
meson test -C _build --print-errorlogs
```

If the testsuite passes you can install with

```
meson configure _build --prefix=/path/to/install
meson install -C _build
```

This might require administrator access depending on the chosen install prefix.
Now you are ready to use ``mctc-rmsd``.


## Usage

To calculate the RMSD between two structures use the [``mctc-rmsd(1)``](man/mctc-rmsd.1.adoc) program with

```
mctc-rmsd struc-a.xyz struc-b.xyz
```

The RMSD is calculated relative to the first provided structure.
Multiple structures can be provided to perform additional RMSD calculations,

```
mctc-rmsd coord conf-1.xyz conf-2.xyz conf-3.xyz
```

To select a subset of atoms use a filter

```
mctc-rmsd xrd.pdb opt.pdb --filter heavy
```

The heavy atom RMSD (excluding all H atoms) is built-in by default.
New filter can be defined in the ``mctc.toml`` configuration file in the current directory or in your home directory with:

```toml
[rmsd.filter]
# Definition of the built-in heavy atom filter
heavy.exclude = ["H", "h"]
# Include all carbon, nitrogen and oxygen atoms
organic.include = [6, 7, 8]
```

For more details look up [the manual page](man/mctc-rmsd.1.adoc).


## License

This project is free software: you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This project is distributed in the hope that it will be useful,
but without any warranty; without even the implied warranty of
merchantability or fitness for a particular purpose.  See the
GNU Lesser General Public License for more details.

Unless you explicitly state otherwise, any contribution intentionally
submitted for inclusion in this project by you, as defined in the
GNU Lesser General Public license, shall be licensed as above, without any
additional terms or conditions.
