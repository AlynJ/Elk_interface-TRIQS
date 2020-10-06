                       +------------------------------+
                       |    The Elk interface Code    |
                       +------------------------------+

This is an altered elk-6.2.8 version which outputs the required files for the 
Elk-TRIQS interface. This is distributed under the terms of the GNU General 
Public License. Details of the original Elk README file is at the bottom of 
this document.

This is an initial version written by Alyn D. N. James from the University of 
Bristol (UoB), United Kingdom. Note that this version is subject to changes. 
The intent is to provide an open source FP-LAPW code interfaced with the TRIQS 
library. Therefore, DFT+DMFT will be more accessible to the community. 
We would like to thank the Elk and TRIQS developers for their help with this 
interface.

This interface is open source, but we ask you to cite "TBA" if you use it. 

Full documentaion of using the interface can be found at 
https://triqs.github.io/dft_tools/.

There are examples of the Elk part of the interface in the 
./examples/Wannier-Projector directory.

If you come across any issue with this interface, please create an issue on 
this GitHub with a concise description about the problem and how we can 
reproduce it. 

We hope you enjoy using the interface!



--------------------------------------------------------------------------------


                       +------------------------------+
                       |     The Elk FP-LAPW Code     |
                       +------------------------------+

This code is distributed under the terms of the GNU General Public License.
See the file COPYING for license details.

    Elk is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Elk is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
    License for more details.

    You should have received a copy of the GNU General Public License
    along with Elk. If not, see http://www.gnu.org/licenses/. 

Elk can be compiled by first running "./setup" in this directory followed by
"make all". This will compile the main code as well as several auxiliary
programs. For optimal performance we strongly recommend that you tune the
Fortran compiler options in the file "make.inc" and use machine-optimised
BLAS/LAPACK libraries. Setting the OpenMP and MPI options of your compiler will
enable Elk to run in parallel mode on multiprocessor systems and clusters.

A test suite is available: entering "make test" will check the output of your
executable against a standard set. This may take some time to complete.

Auxiliary programs include "spacegroup" for producing crystal geometries from
spacegroup data, and "eos" for fitting equations of state to energy-volume data.

Elk is updated regularly with new features and bug fixes. Features not listed as
"experimental" may be used for production but, as with any code, please check
the consistency of your results carefully.

--------------------------------------------------------------------------------
J. K. Dewhurst, S. Sharma
L. Nordstrom,  F. Cricchio, O. Granas
E. K. U. Gross

