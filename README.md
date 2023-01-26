c-urchin
=======

c-urchin is a small, C++-based cyclotron radiation power calculation for 90° betas in a waveguide using analytical solutions for circular and rectangular waveguides (publication upcoming).

Dependencies
------------

**External**
- CMake (3.6 or better)
- Boost (filesystem, math, program options. 1.46 or better)

Operating System Support
------------------------

* Mac OS X (tested on OS X 13.1)


Tips on Installing the Dependencies
-----------------------------------

MacOS: [Homebrew](https://brew.sh/) is a convenient package manager. Before using it for the installation, make sure that the version available on brew is compatible with what is listed in the dependency list above. [CMake](http://brewformulas.org/Cmake), [Boost](http://brewformulas.org/Boost).

Installing
----------


1. Download the repository and create a build directory:
  ```
  $ git clone "https://github.com/Helium6CRES/c-urchin"
  ```

2. Unzip the Data Files:
  ```
  $ bunzip2 Data/*bz2
  ```

3. Make your Build Directory and Configure the installation: 

  ```
  $ mkdir build
  $ cd build
  $ cmake ..
  ```

  The install prefix is specified by the CMake variable `CMAKE_INSTALL_PREFIX`.
  The library, binaries, and header files will be installed in the
  lib, bin, and include subdirectories. The default install prefix is the build directory.

You should set the CMake variable `CMAKE_BUILD_TYPE` to either `RELEASE`, `STANDARD`, or `DEBUG` (default), in order of how much text output you would like (from least to most) and how much compiler optimization should be performed (from most to least).
`RELEASE` is recommended when running in production, as it is >2x faster than `DEBUG`.

4. Build and install.
  ```
  $ make install
  ```

If you made a change to the dependencies, you may have to wipe the build directory and start again from step 1; simply writing `make install` again will not always work. 

To run from any directory, you may have to set your path variable in your .bashrc/.zshrc/ etc.
  ```
    export PATH=/path/to/c-urchin/build/bin:$PATH
  ```

Instructions for Use
--------------------
To produce one beta with a particular set of parameters:

```
  > c-urchin --start-frequency 19.1e9 --start-rho 1e-3 --start-field 2.25
```

To see all command options, use:

```
  > c-urchin --help
```

The start, end, and number parameters function like the numpy [linspace](https://numpy.org/doc/stable/reference/generated/numpy.linspace.html)

**Please be careful**, by default the simulation output writes to the directory you call c-urchin from, appending to file. If files aren't managed, they will become quite large.

Development
-----------

Issues should be posted via [GitHub](https://github.com/Helium6CRES/c-urchin/issues)
