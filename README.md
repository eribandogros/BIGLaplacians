Dependencies: Matlab

1. clone the repo

2. create build folder

cd EulerianMesh3D

mkdir build

3. generate build files using cmake (creates an Xcode project)

cmake -G Xcode EulerianMesh3D ..

4. open Xcode project EulerianMesh3D.xcodeproj

5. change build target from ALL_BUILD to EulerianMesh3D in drop down menu next to run/stop buttons

6. click Edit Scheme... in drop down menu

7. add args in Run/Arguments

The first arg is the length of the model, the second is the grid length, and the third is which model ('c' for cube, 's' for sphere and 't' for torus).
