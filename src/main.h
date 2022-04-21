//
//  main.h

#ifndef main_hpp
#define main_hpp

#include <stdio.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include "Vector3d.h"

CVector3d return_grid_pos(int i,int j,int k,double grid_size,CVector3d start_pos);

double return_cube_distance(CVector3d pos,double cube_length);
double return_sphere_distance(CVector3d pos,double sphere_radius);
double return_ellipsoid_distance(CVector3d pos,double sphere_radius);
double return_inside_cuboid(CVector3d pos,double cube_legnth);
double return_torus_distance(CVector3d pos,double radius);

void compute_level_set_cube(double cube_length,double grid_size);
void compute_level_set_cube_inside_cube(double cube_length,double grid_size);
void compute_level_set_sphere(double radius,double grid_size);
void compute_level_set_ellipsoid(double radius,double grid_size);
void compute_level_set_torus(double radius,double grid_size);

#endif /* main_h */
