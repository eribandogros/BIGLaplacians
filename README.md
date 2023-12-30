
Code from "Combinatorial and Hodge Laplacians: Similarity and Difference." Computes the multidimensional combinatorial and Hodge Laplacians of the grid representation of an elementary volume with normal boundary conditions. Output is the smallest eigenvalues of both Laplacians.

Paper: [arXiv:2204.12218](https://doi.org/10.48550/arXiv.2204.12218)

Laplacian eigenfields displayed in Mathematica with ListStreamPlot3D.

<img height="200" alt="nonzeroU" src="https://github.com/eribandogros/BIGLaplacians/assets/14114157/b103e6c4-078e-4dfb-a079-a4cf708d18d1">
<img height="200" alt="field" src="https://github.com/eribandogros/BIGLaplacians/assets/14114157/f96e0443-4c59-40d2-94b9-dfeab8e9a44f">

---
**Build and run the project using CMake.**

Dependencies: Matlab

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

<img height="200" alt="cube" src="https://github.com/eribandogros/BIGLaplacians/assets/14114157/db8e7561-e377-4581-8d2f-644b63295f16">
<img height="200" alt="ball" src="https://github.com/eribandogros/BIGLaplacians/assets/14114157/979e4324-3f54-4592-9512-8206ccc72dd6">
<img height="200" alt="torus" src="https://github.com/eribandogros/BIGLaplacians/assets/14114157/a0cd0522-1e73-4df9-9172-340fd3320772">
<img height="200" alt="shell" src="https://github.com/eribandogros/BIGLaplacians/assets/14114157/3a6d0873-7f0c-4ceb-8757-4434ed2465e7">

