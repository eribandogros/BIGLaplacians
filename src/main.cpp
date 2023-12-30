//
//  main.cpp

#include "main.h"
#include "grid_mesh.h"
#include "Vector3d.h"
#include <algorithm>
#include <string>

CVector3d return_grid_pos(int i, int j, int k, double grid_length, CVector3d start_pos){
    CVector3d diff = CVector3d(i*grid_length,j*grid_length,k*grid_length);
    return (diff+start_pos);
}

double return_cube_distance(CVector3d pos, double cube_length){
    double tmp = std::max(fabs(pos[0])-cube_length/2.0,fabs(pos[1])-cube_length/2.0);
    double max_value = std::max(tmp,fabs(pos[2])-cube_length/2.0);

    return (max_value);
}

double return_sphere_distance(CVector3d pos, double sphere_radius){
    return (pos.LengthSquared() - sphere_radius*sphere_radius);
}

double return_ellipsoid_distance(CVector3d pos, double sphere_radius){
    return (pos[0]*pos[0]*2+pos[1]*pos[1]*3+pos[2]*pos[2] - sphere_radius*sphere_radius);
}

double return_torus_distance(CVector3d pos, double outer_radius){
    double inner_radius = outer_radius-0.5;
    CVector3d pos_xd(pos[0], 0, pos[2]);
    double x = pos_xd.Length() - outer_radius;
    double y = pos[1];
    CVector3d new_pos(x, y, 0);
    return (new_pos.Length()-inner_radius);
}

double return_spherical_shell_distance(CVector3d pos, double outer_radius){
    double inner_radius = 0.5;
    double r = pos.Length(); //distance from origin
    if (r < inner_radius){ //inner cavity
        return inner_radius - r;
    }
    else if (r >= inner_radius && r <= outer_radius){ //shell interior
        double d = -std::min(r - inner_radius, outer_radius - r);
        return d;
    }
    else if (r > outer_radius){ //outside shell
        return r - outer_radius;
    }
    return 0.0;
}

void compute_level_set_cube(double cube_length, double grid_length){
    int x_num_times = ceil((cube_length/2.)/grid_length)*2+1;
    int z_num_times = ceil((cube_length/2.)/grid_length)*2+1;
    int y_num_times = ceil((cube_length/2.)/grid_length)*2+1;
        
    double start = -(x_num_times-1)/2.*grid_length;
    CVector3d start_pos = CVector3d(start, start, start);

    int x_num_dimension = x_num_times;
    int y_num_dimension = y_num_times;
    int z_num_dimension = z_num_times;

    std::ofstream file("levelset_cube");
    file << x_num_dimension << std::endl;
    file << y_num_dimension << std::endl;
    file << z_num_dimension << std::endl;
    file << grid_length << std::endl;
    
    for(int k=0;k<z_num_dimension;k++)
        for(int j=0;j<y_num_dimension;j++)
            for(int i=0;i<x_num_dimension;i++){
                CVector3d grid_pos = return_grid_pos(i,j,k,grid_length,start_pos);
                double distance_val = return_cube_distance(grid_pos,cube_length);
                file << distance_val << std::endl;
            }
    file.close();
}

void compute_level_set_sphere(double radius, double grid_length){
    int num_times = ceil(radius/grid_length)+1;
    CVector3d start_pos = CVector3d(-num_times*grid_length,-num_times*grid_length,
        -num_times*grid_length);

    int num_dimension = 2*num_times+1;

    std::ofstream file("levelset_sphere");
    file << num_dimension << std::endl;
    file << num_dimension << std::endl;
    file << num_dimension << std::endl;
    file << grid_length << std::endl;

    for(int k=0;k<num_dimension;k++)
        for(int j=0;j<num_dimension;j++)
            for(int i=0;i<num_dimension;i++){
                CVector3d grid_pos = return_grid_pos(i,j,k,grid_length,start_pos);
                double distance_val = return_sphere_distance(grid_pos,radius);
                file << distance_val << std::endl;

            }
    file.close();
}

void compute_level_set_torus(double radius, double grid_length){
    double inner_radius = radius-0.5;
    int num_times = ceil((radius+inner_radius)/grid_length)+1;
    CVector3d start_pos = CVector3d(-num_times*grid_length,-num_times*grid_length,
        -num_times*grid_length);

    int num_dimension = 2*num_times;

    std::ofstream file("levelset_torus");
    file << num_dimension << std::endl;
    file << num_dimension << std::endl;
    file << num_dimension << std::endl;
    file << grid_length << std::endl;

    for(int k=0;k<num_dimension;k++)
        for(int j=0;j<num_dimension;j++)
            for(int i=0;i<num_dimension;i++){
                CVector3d grid_pos = return_grid_pos(i,j,k,grid_length,start_pos);
                double distance_val = return_torus_distance(grid_pos,radius);
                file << distance_val << std::endl;
            }
    file.close();
}

void compute_level_set_spherical_shell(double radius, double grid_length){

    int num_times = ceil(radius/grid_length)+1; //radius plus one more grid point
    CVector3d start_pos = CVector3d(-num_times*grid_length,-num_times*grid_length,
        -num_times*grid_length);

    int num_dimension = 2*num_times+1; //each half of grid plus axes

    std::ofstream file("levelset_spherical_shell");
    file << num_dimension << std::endl;
    file << num_dimension << std::endl;
    file << num_dimension << std::endl;
    file << grid_length << std::endl;

    for(int k=0;k<num_dimension;k++)
        for(int j=0;j<num_dimension;j++)
            for(int i=0;i<num_dimension;i++){
                CVector3d grid_pos = return_grid_pos(i,j,k,grid_length,start_pos);
                double distance_val = return_spherical_shell_distance(grid_pos,radius);
                file << distance_val << std::endl;
            }
    file.close();
}

int main(int argc, char* argv[]){
    grid_mesh g;
    g.m_grid_size = std::atof(argv[1]);
    g.m_grid_length = std::atof(argv[2]);
    std::cout << "model size: " << g.m_grid_size << std::endl;
    std::cout << "grid length: " << g.m_grid_length << std::endl;
    
    std::string fname;
    std::cout << "num of args (should be 5): " << argc << std::endl;
    
    if(argc > 3){
        g.model = *argv[3];
        if(*argv[3] == 'c'){
            std::cout << "cube model" << std::endl;
            compute_level_set_cube(g.m_grid_size, g.m_grid_length);
            fname = "levelset_cube";
        }
        else if(*argv[3] == 's'){
            std::cout << "sphere model" << std::endl;
            compute_level_set_sphere(g.m_grid_size, g.m_grid_length); //radius 1st param
            fname = "levelset_sphere";
        }
        else if (*argv[3] == 't'){
            std::cout << "torus model" << std::endl;
            compute_level_set_torus(g.m_grid_size, g.m_grid_length);
            fname = "levelset_torus";
        }
        else if (*argv[3] == 'l'){
            std::cout << "spherical shell model" << std::endl;
            compute_level_set_spherical_shell(g.m_grid_size, g.m_grid_length);
            fname = "levelset_spherical_shell";
        }
    }
    else{
        g.model = 'c';
        compute_level_set_cube(g.m_grid_size, g.m_grid_length);
        fname = "levelset_cube";
    }
    
    if(argc > 4){
        g.option = *argv[4];
    }
    
    g.read_level_set(fname);
    g.init();
    
    return 0;
}
