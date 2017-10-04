*******************************************************************************
*******************************************************************************
# PowerI4
A power spectrum estimator with interlacing and fourth-order interpolation
*******************************************************************************
v1 Oct 2017
*******************************************************************************

Latest version on https://github.com/sefusatti/PowerI4
DOI https://

*******************************************************************************

Contents of this file:
--- Files and dependencies
--- Installation
--- Usage and testing

*******************************************************************************
Files and dependencies:
*******************************************************************************

The code is organised in the following files:
- PowerI4.f90 : Main code, contains also the subroutines for density grid
normalization and file outputs
- grid.f90 : specifies the FFT interpolated grid
- input_catalog.f90 : contains all the subroutines for reading different input
catalog formats. You can add here a subroutine reading your own format.
- interpolation.f90 : contains all the interpolation subroutines (all different
mass assignment schemes for assignment on single and double grids)
- parbox.f90 : contains global variables

PowerI4 requires the use of FFTW (www.fftw.org). Current version uses FFTW3.
Make sure you download and compile FFTW and that you put the correct path to
FFTW libraries in the Makefile (in the variables FFTW_LIBR and FFTW_INCL). The
fftw3.f is included in the first line of PowerI4.f90 .

*******************************************************************************
Installation:
*******************************************************************************

A Makefile to compile PowerI4 is provided. Detail your system configuration
in the variable SYSTYPE. Running make should produce the executable power.e

*******************************************************************************
Usage and testing:
*******************************************************************************

We provide a parameters.ini file for a test run, and a set of particles to
measure the spectrum from. They are both located in the /test/ folder. To run

./power.e < test/parameters.ini

This run should produce a power spectrum file ps_n0256_4i.out in the folder
/test/ that should match the file ps_n0256_4i_ref.out in the same directory.
For further runs, set-up your own parameters.ini.

A brief description of each entry in the parameters.ini file is provided in the
file itself, while a more detailed description can be found in the user
manual PowerI4.pdf

*******************************************************************************
*******************************************************************************
