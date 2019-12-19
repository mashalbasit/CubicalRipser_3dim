# CubicalRipser : Persistent homology for 2D image and 3D voxel data

copyright by Takeki Sudo and Kazushi Ahara, Meiji University, 2018

modified by Shizuo Kaji, Kyushu University, 2019

## Description
CubicalRipser is an adaptation of [Ripser](http://ripser.org) by Ulrich Bauer for a filtered cubical complex.
For 2 and 3 dimensional cubical complexes, we believe CubicalRipser is currently the fastest program for computing persistent homology.

## License
CubicalRipser is free software: you can redistribute it and/or modify it under
the terms of the GNU Lesser General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your option)
any later version.

## How to use
To build the software:

    % make

NOTE: We recommend gcc. For some reasons, Apple Clang deteriorates performance drastically.

To see the command-line options:

    % ./cubicalripser

Example:

    % ./cubicalripser --print --location --maxdim 2 --output out.csv demo/3dimsample.txt

## Input file format
CubicalRipser takes three types of input files: NUMPY, TEXT, DIPHA.

- 2D or 3D Numpy array (recommended). The filename should end with ".npy". DType must be float64. see the Jupyter Notebook example found under the demo directory.
- Text file. The filename should end with ".txt"
```
3
max_x
max_y
max_z
val[1,1,1]
val[2,1,1]
...
val[max_x,max_y,max_z]
```
- [DIPHA binary format](https://github.com/DIPHA/dipha#file-formats) 


## Difference with the original version
I added the following functionality:
- input/output numpy array file (.npy): see the Jupyter Notebook example found under the demo directory.
- no limitation on the input size
- cleaned up/optimised codes
- output birth location
- bug fixes (including memory leaks)
