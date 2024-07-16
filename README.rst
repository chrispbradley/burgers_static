==============
Burgers Static
==============

Burgers equation is the simplified version of the Navier-Stokes equations (by dropping the pressure term). This equation can be investigated in one spatial dimension.
The static form of the Burgers equation is implemented in this example.


Building the example
====================

The fortran version of the example can be configured and built with CMake::

  git clone https://github.com/OpenCMISS-Examples/burgers_static.git
  cd burgers_static
  mkdir build
  cd build
  cmake -DOpenCMISS_INSTALL_ROOT=/path/to/opencmiss/install/ ../.
  make


Running the example
===================

Fortran version::

  cd src/fortran/
  ./burgers_static_fortran

The results can be visualised by running `visualiseBurgersStatic.cmgui <./src/fortran/visualiseBurgersStatic.cmgui>`_ with the `Cmgui visualiser <http://physiomeproject.org/software/opencmiss/cmgui/download>`_.


References
==========

Heinrich, J.C. and Pepper, D.W., 1999. Intermediate finite element method: fluid flow and heat transfer applications. (pg. 216)


Prerequisites
=============

There are no additional input files required for this example as it is self-contained.


License
=======

License applicable to this example is described in `LICENSE <./LICENSE>`_.
