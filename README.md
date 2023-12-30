
Code from "Combinatorial and Hodge Laplacians: Similarity and Difference." Computes the multidimensional combinatorial and Hodge Laplacians of the grid representation of an elementary volume with normal boundary conditions. Output is the smallest eigenvalues of both Laplacians.

[arXiv:2204.12218](https://doi.org/10.48550/arXiv.2204.12218)

Dependencies: Matlab

Build and run the project using CMake.

1. clone the repo

2. create build folder within repo

`mkdir build`
`cd build`

3. generate build files using cmake (example creates an Xcode project)

`cmake -G Xcode ..`

4. open Xcode project NormalEulerianMesh3D.xcodeproj

5. change build target from ALL_BUILD to NormalEulerianMesh3D in drop down menu next to run/stop buttons

6. click Edit Scheme... in drop down menu

7. add args in Run/Arguments

- arg1, the length of the model (the distance of the outermost boundary from the origin)
- arg2, grid length
- arg3, model type ('c' for cube, 's' for sphere/ball, 't' for torus, 'l' for spherical shell).
- arg4, dimension of Laplacian (0, 1, or 3)