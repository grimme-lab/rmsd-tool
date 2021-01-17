---
project: mctc-rmsd
summary: Root mean square deviation for structural comparison
project_github: https://github.com/grimme-lab/rmsd-tool
project_download: https://github.com/grimme-lab/rmsd-tool/releases
author: Grimme group, Bonn
github: https://github.com/grimme-lab
src_dir: ./src
         ./app
output_dir: ./_docs
exclude_dir: ./test
docmark: <
predocmark: >
source: true
graph: false
sort: alpha
print_creation_date: true
extra_mods: iso_fortran_env:https://gcc.gnu.org/onlinedocs/gfortran/ISO_005fFORTRAN_005fENV.html
creation_date: %Y-%m-%d %H:%M %z
md_extensions: markdown.extensions.toc
               markdown.extensions.smarty
---

Tool to calculate root mean square deviations between molecular structures.


## Getting Started

### Meson

Create a new meson project and include `rmsd` either as git-submodule in your subprojects directory or create a wrap file to fetch it from upstream:

```ini
[wrap-git]
directory = mctc-rmsd
url = https://github.com/grimme-lab/rmsd-tool
revision = head
```

The `rmsd` library depends on the [MCTC-library](https://github.com/grimme-lab/mctc-lib) and [TOML-Fortran](https://github.com/toml-f/toml-f).
You might have to promote those dependencies to your subprojects by copying the wrap files from the `rmsd` subprojects.

To load the project the necessary boilerplate code for subprojects is just

<!--pygments doesn't know about meson, python highlighting looks okayish-->
```python
rmsd_prj = subproject(
  'mctc-rmsd',
  version: '>=0.1',
  default_options: [
    'default_library=static',
  ],
)
rmsd_dep = rmsd_prj.get_variable('rmsd_dep')
```

Now you can add `rmsd_dep` to your dependencies and access the public API by the `rmsd` module.

We recommend to set the default library type of `rmsd` to static when linking your applications or library against it.
Note for library type both and shared `rmsd` will install itself along with your project.

For more fine-tuned control you can access:

- the library target with `rmsd_lib`
- the private include dir of this target, containing the Fortran module files, with `rmsd_inc`
- the license files of `rmsd` with `rmsd_lic`

If you are linking your application statically against `rmsd` and still want to distribute the license files of `rmsd` (thank you), just use

```python
install_data(
  rmsd_prj.get_variable('rmsd_lic'),
  install_dir: get_option('datadir')/'licenses'/meson.project_name()/'mctc-rmsd',
)
```


### Fortran Package Manager (fpm)

This project supports [fpm](https://github.com/fortran-lang/fpm) as build system as well.
Just add it to the dependencies in your `fpm.toml` file:

```toml
[dependencies]
[dependencies.mctc-rmsd]
git = "https://github.com/grimme-lab/rmsd-tool"
```
