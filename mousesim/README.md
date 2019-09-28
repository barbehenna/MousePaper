# mousesim package

This package is basically a wrapper of all the simulation functionality I use in this project. 

## Notes

- To update the R wrappers use the command `Rcpp::compileAttributes()`
- To build documentaion, you may need to first compile a `.so` library with `pkgbuild::compile_dll()`
- `devtools::document()`, `devtools::build()`, `devtools::check()`, `devtools::install()` all help the build and install process along

## To-do

- add unit tests based on the R file in `cpp` that was used for initial debugging
- maybe hide some unneeded functions
 