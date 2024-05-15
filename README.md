# APECSS 

<p align="left">
  <a href="https://github.com/polycfd/apecss/actions/workflows/test_build.yml">
    <img src="https://github.com/polycfd/apecss/actions/workflows/test_build.yml/badge.svg" alt="Issues">
  </a>
  <a href="https://github.com/polycfd/apecss/actions/workflows/test_run.yml">
    <img src="https://github.com/polycfd/apecss/actions/workflows/test_run.yml/badge.svg" alt="License">
  </a> 
  <a href="https://doi.org/10.5281/zenodo.7249297">
    <img src="https://img.shields.io/badge/DOI-10.5281/zenodo.7249297-blue" alt="Latest version">
  </a>
  <a href="https://joss.theoj.org/papers/27166cd5496c33d62e74132712efec8a">
    <img src="https://joss.theoj.org/papers/27166cd5496c33d62e74132712efec8a/status.svg">
  </a>
</p>

APECSS is a software toolbox to compute pressure-driven bubble dynamics and the resulting acoustic emissions. It is written in C and has been developed with simplicity, versatility and performance in mind. The acronym APECSS stands for "Acoustic Pulse Emitted by Cavitation in Spherical Symmetry".

Key features of APECSS are:
- Bubble dynamics using widely-used models (Rayleigh-Plesset, Keller-Miksis, Gilmore), solved using an in-built 5th-order Runge-Kutta scheme with adaptive time stepping.
- Acoustic emissions of the bubble under different assumptions (incompressible, quasi-acoustic, fully compressible).
- *(Work in progress)* Interbubble interactions within a cavitation cluster
- Prediction of the formation and attenuation of shock fronts emitted by the bubble.
- Viscoelastic media (Kelvin-Voigt, Zener, Oldroyd-B).
- Lipid monolayer coating of the bubble as used for ultrasound contrast agents.
- APECSS has, aside from the C standard library, no external dependencies.

For more information, please consult the [JOSS paper](https://doi.org/10.21105/joss.05435) discussing the features of APECSS as well as the [PoF paper](https://doi.org/10.1063/5.0131930) introducing the theoretical foundation.

## Developers
- [Fabian Denner](mailto:fabian.denner@polymtl.ca) (principal developer, maintainer)
- [Sören Schenke](mailto:soeren.schenke@ovgu.de)

## License and Copyright
APECSS is under the copyright of its developers and made available as open-source software under the terms of the [Mozilla Public License Version 2.0](LICENSE.txt).

## Requirements

APECSS has been developed and tested on Unix systems, using Linux and MacOS operating systems. APECSS has also been installed and executed successfully on Windows 11. The only mandatory requirement to compile and run APECSS is a standard C compiler, such as [````gcc````](https://gcc.gnu.org).

Optionally, [````cmake````](https://cmake.org) (version 3.12 or higher) is required if you would like to make use of the provided compilation and test scripts. To use the provided Python scripts for visualizing the output of APECSS, [````Python````](https://www.python.org) (version 3), [````numpy````](https://numpy.org) and [````matplotlib````](https://matplotlib.org) are required. 

## Quick Start Guide

### **Installation on Unix systems**
Getting started with APECSS using a Unix operating system, e.g. Linux or MacOS, is easy. After downloading APECSS in the directory ````<path to APECSS>````, define the following environment variables:
- ````APECSS_DIR```` to the directory in which APECSS is located. Using bash, for instance, simply execute the command ````export APECSS_DIR=<path to APECSS>```` or, even better, add this command to your bash profile.
- ````USRLIB_DIR```` to the directory in which libm.a or libm.dylib (the standard _math_ library) is located. This may, for instance, be ````/usr/lib64/```` on Linux systems or ````/usr/lib/```` on MacOS systems.

Now, navigate into the folder ````$APECSS_DIR/lib```` and execute ````./compile_lib.sh````. This shell script will compile the APECSS library using _cmake_ with the ````CMakeLists.txt```` file provided in this folder. By default, APECSS is compiled with double precision and in _Release_ mode, meaning all optimization flags are enabled. That's it, you've successfully compiled APECSS!

### **Installation on Windows 11**
Installing APECSS on Windows is possible, but a little more involved than the installation on Unix. A step-by-step installation guide on Windows 11, kindly provided by Maxime Montois, can be found in the [pdf-documentation](/documentation/APECSS-Documentation.pdf) (see Section 2.2.2). The installation uses `MinGW`, `MSYS2` and `Git Bash` and the instructions also detail small changes to the compile and cmake scripts, as well as the source code, that are required to make APECSS work smoothly on Windows.

### **Examples**
There are several ways in which you can use the APECSS library. You can either incorporate selected features of APECSS into your own software or you can program an interface to use APECSS as a standalone program. Some representative examples are given in the [````$APECSS_DIR/examples````](/examples/) directory. Each directory contains the following:
- A ````README.md```` file explaining the purpose and specificities of this/these example(s).
- A ````src```` folder with a file called ````*_apecss.c```` that acts as the standalone interface to the APECSS library. This file contains the ````main()```` function and any additional functionality required to simulate a specific scenario.
- A ````build```` folder containing the ````CMakeLists.txt```` file and a shell script ````compile.sh```` with which this example can be compiled using the command ````./compile.sh````.
- One or several ````*.apecss```` files in which the options for a specific case are defined.
- One or several ````plot_*.py```` Python scripts to plot the results of the example. These Python scripts require ````numpy```` and ````matplotlib````, and plot the results into a pdf-file. 
- Some examples have a ````reference-data```` folder, which contains results of previous studies to quickly validate the overall correct working of APECSS and which are plotted alongside the results produced by APECSS using the Python script(s).

You can run **all** examples by executing ````./run_all.sh```` in the [````examples````](/examples/) folder. This shell script compiles, runs and plots the results of all test cases. This is a great way of testing the installation of APECSS and see what sort of results can be output and how these results may be visualized. The script also compares the results computed for a respresentative example (the [Rayleigh collapse with emissions](/examples/rayleighcollapse/)) with reference results of APECSS, checking whether the installation produces correct results. Note that ````run_all.sh```` assumes that *Python* can be executed using the command ````python3````; please adjust this if necessary.

Alternatively, if you're, for instance, specifically interested in an ultrasound-driven bubble, navigate to [````$APECSS_DIR/examples/ultrasound````](/examples/ultrasound/). There you can find different test cases, described in more detail in the accompanying [````README.md````](/examples/ultrasound/README.md). Execute the ````compile.sh```` script in the [````build````](/examples/ultrasound/build/) directory. You can now choose which test case to run, by using the corresponding execution command suggested in the [````README.md````](/examples/ultrasound/README.md). Also, the settings and options given by the execution command and the corresponding ````*.apecss```` file typically reproduce a case previously presented in the literature, as indicated in the ````README.md```` file.

More information about each example, how to run it and what the results might be compared to can be found in the accompanying ````README.md```` file.

## Repository Structure
The APECSS repository is structured as follows:
- The [documentation](/documentation/) folder contains a short [pdf](/documentation/APECSS-Documentation.pdf) documentation of APECSS. The documentation discusses the theory behind APECSS, explains the code structure and how to use APECSS. The documentation will be amended and expanded over time.
- The [examples](/examples/) folder contains representative examples of how to use APECSS and to demonstrate the most important features of APECSS. A short explanation on how to run the examples is given in the [Quick Start Guide](#quick-start-guide) above.
- The [include](/include/) folder contains the [apecss.h](/include/apecss.h) header file, in which all variables, macros and functions of APECSS are defined.
- The [lib](/lib/) folder in which the APECSS library is compiled (at least if you follow the [Quick Start Guide](#quick-start-guide) above).
- The [src](/src/) folder contains all source files (*.c) of APECSS.
- The [.clang-format](.clang-format) file, which defines the formatting rules for the source code.
- The [.gitignore](.gitignore) file tells _git_ which folders and files to ignore.
- The [CITATION.cff](CITATION.cff) file contains information on how to cite APECSS.
- The [CODE_OF_CONDUCT.md](CODE_OF_CONDUCT.md) file outlines the code of conduct for this repository.
- The [CONTRIBUTING.md](CONTRIBUTING.md) file contains brief guidelines on how to contribute to APECSS.
- The [JOSS-Paper.pdf](JOSS-Paper.pdf) file is the paper about APECSS in the Journal of Open Source Software.
- The [LICENSE.txt](LICENSE.txt) file contains the Mozilla Public License Version 2.0.
- The [README.md](README.md) file is the file you are currently reading.

## How to cite us
If you use APECSS for your scientific work, please consider citing the [paper](https://doi.org/10.21105/joss.05435) introducing the features and capabilities of APECSS

    F. Denner and S. Schenke, APECSS: A software library for cavitation bubble dynamics and acoustic emissions. Journal of Open Source Software 8 (2023), 5435. https://doi.org/10.21105/joss.05435 

and, in the interest of reproducibility, the version of APECSS you've used for your work, e.g.

    F. Denner and S. Schenke, APECSS (v1.2), Zenodo (2022). https://doi.org/10.5281/zenodo.7465050

All releases of APECSS and the corresponding DOIs can be found on the [Zenodo page](https://doi.org/10.5281/zenodo.7249297) of APECSS.

## Acknowledgements
The development of APECSS has directly benefitted from research funding provided by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation), grant number 441063377.
