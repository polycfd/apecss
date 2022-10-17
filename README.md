# APECSS
APECSS is a software toolbox to compute pressure-driven bubble dynamics and the resulting acoustic emissions. It is written in C and has been developed with simplicity, versatility and performance in mind. The acronym APECSS stands for "Acoustic Pulse Emitted by Cavitation in Spherical Symmetry".

Key features of APECSS are:
- Includes widely-used models for the bubble dynamics (Rayleigh-Plesset, Keller-Miksis, Gilmore).
- Acoustic emissions of the bubble under different assumptions (incompressible, quasi-acoustic, fully compressible).
- Prediction of the formation and attenuation of shock fronts emitted by the bubble.
- Viscoelastic media (Kelvin-Voigt, Zener, Oldroyd-B).
- Lipid monolayer coating of the bubble as used for ultrasound contrast agents.
- All ODEs are solved with an in-built 5th-order Runge-Kutta scheme with 4th-order error estimate, based on Dormand and Prince (1980).
- APECSS has no external dependencies, aside from some common C headers (math.h, stdio.h, stdlib.h, string.h).

## Developers
- [Fabian Denner](mailto:fabian.denner@ovgu.de) (principal developer, maintainer)
- [Sören Schenke](mailto:soeren.schenke@ovgu.de)

## License and Copyright
APECSS is under the copyright of its developers and made available as open-source software under the terms of the [Mozilla Public License Version 2.0](LICENSE.txt).

## Repository Structure
The APECSS repository is structured as follows
- The [documentation](/documentation/) folder contains a short documentation of APECSS, written in Latex.
- The [examples](/examples/) folder contains representative examples of how to use APECSS and to demonstrate the most important features of APECSS. These examples also serve to validate APECSS against results reported in the literature.
- The [include](/include/) folder contains [apecss.h](/include/apecss.h) header file, in which all variable, macro and function of APECSS are defined.
- The [lib](/lib/) folder in which the APECSS library is compiled (at least if you follow the [Quick Start Guide](#quick-start-guide) below).
- The [src](/src/) folder contains all source files (*.c) of APECSS.
- The [.clang-format](.clang-format) file, which defines the formatting rules for the source code.
- The [.gitignore](.gitignore) file telling _git_ which folders and files to ignore.
- The [LICENSE.txt](LICENSE.txt) file containing the Mozilla Public License Version 2.0.
- The [README.txt](README.txt) file is the file you are currently reading.

## Quick Start Guide
Getting started with APECSS is easy. After downloading APECSS in the directory ````<path to APECSS>````, define the following environment variables:
- ````APECSS_DIR```` to the directory in which APECSS is located. Using bash, for instance, simply execute the command ````export APECSS_DIR=<path to APECSS>```` or, even better, add this command to your bash profile.
- ````USRLIB_DIR```` to the directory in which libm.a or libm.dylib (the standard _math_ library) is located. This may, for instance, be ````/usr/lib64/```` on Linux systems or ````/usr/lib/```` on MacOS systems.

Now navigate into the folder ````$APECSS_DIR/lib```` and execute ````./compile_lib.sh````. This command will compile the APECSS library using the cmake with the ````CMakeLists.txt```` file provided in this folder. By default, APECSS is compiled with double precision and in _Release_ mode, meaning all optimization flags are enabled. That's it, you've succesfully compiled the APECSS library!

There are several ways in which you can use the APECSS library. You can either incorporate selected features of APECSS into your own software or you can program an interface to use APECSS as a standalone program. Some representative examples are given in the ````$APECSS_DIR/examples```` directory. Each directory contains the following:
- A ````README.md```` file explaining the purpose and specificities of this/these example(s).
- A ````src```` folder with a file called ````*_apecss.c```` that acts as the standalone interface to the APECSS library. This file contains the ````main()```` function and any additional functionality required to simulate a specific scenario.
- A ````build```` folder containing the ````CMakeLists.txt```` file and a shell script ````compile.sh```` with which this example can be compiled using the command ````./compile.sh````.
- One or several ````*.apecss```` files in which the options for a specific case are defined.

## Acknowledgements
The development of APECSS has directly benefitted from research funding provided by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation), grant number 441063377.