*qcsim* is an R package to simulate the power of clinical laboratory QC schemes to detect changes in assay bias or imprecision.  It simulates a variety of traditional and nontraditional QC rules, and it effectively scales to highly multiplexed assays.

To install using the *remotes* library (available from CRAN), use the following command:
remotes::install_github("srmaster90/qcsim")

Documentation and examples for the *qcsim()* function are found in the man files (available via *help()*).

Since this package contains C++, a compiler will be required for installation as for similar R packages.  **To build on Mac**: install Xcode if compilers are not already present.  **To build on Windows:** Install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) before installing this package. Choose the version of Rtools that matches your R version (e.g., Rtools44 for R 4.4.x). Rtools provides the C++ compiler needed to build packages with compiled code from source. Once Rtools is installed, the `remotes::install_github()` command above will work as expected. (thanks to R. Julian for Windows instructions)

