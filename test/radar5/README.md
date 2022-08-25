The file `waltman.csv` contains a reference solution of the stiff DDE problem described in Section 5.1, eq. (7), in the paper at time points 0, 0.1, ..., 300.
The solution was computed with the DDE solver RADAR5 by Hairer and Guglielmi.

If you want to re-compute it, you should:

- Download `waltman.patch` that you find in this directory.

- Download [RADAR5](http://www.unige.ch/~hairer/radar5-v2.1.tar) and extract it on your computer.

- Navigate to the folder `RADAR5-V2.1` that you just extracted.

- Patch the code for the DDE problem and its build script by applying the patch `waltman.patch`:
  ```shell
  $ patch -ruN -p1 < path/to/waltman.patch
  ```

- Navigate to the subfolder `WALTMAN`:
  ```shell
  cd WALTMAN
  ```

- Compile the solver (this requires [`make`](https://www.gnu.org/software/make/) and [`gfortran`](https://gcc.gnu.org/fortran/)):
  ```shell
  make prog
  ```

- Compute the solution:
  ```shell
  ./waltman
  ```

The last step will create the file `waltman.csv` with the desired solution in the current directory.