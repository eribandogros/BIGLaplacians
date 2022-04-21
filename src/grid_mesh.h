//
//  grid_mesh.h

#ifndef grid_mesh_h
#define grid_mesh_h

#include <stdio.h>
#include "SparseMatrix.h"
#include "Vector3d.h"
#include <string>
#include <vector>

#include "Eigen/Dense"
#include "Eigen/Sparse"

#include "MatlabEngine.hpp"
#include "MatlabDataArray.hpp"

#define PI 3.1415926535897932384

class grid_mesh{
public:
    char model;
    char option = '0';
    double m_grid_length;
    double m_grid_size;
    int m_x_num_dimension,m_y_num_dimension,m_z_num_dimension; //number of dimensions in x,y,z direction
    
    int m_estimated_num_faces;//note that it is larger than the real number of faces
    int m_num_interior_faces,m_num_interior_edges;
    int m_num_grid_points;
    int m_num_vertices_new_t;
    int m_num_vertices_new_n;
    int m_estimated_num_cells;
    int m_estimated_num_edges;
    int m_num_edges_new_t;
    int m_num_edges_new_n;
    int m_num_faces_new_n;
    bool *m_is_counted_face; //indicate whether this face is counted in the structure
    bool *m_is_interior_face; //indicate whether this face is an interior face or not
    bool *m_is_cell_included; //each cell is considered included if at least one of the grid point is non-positive
    int *m_face_real_index,*m_edge_real_index; //-1 for faces not included or boundary
    double *m_cell_volume; //the true volume for each cell
    
    double *m_star2_array, *m_star1_array;
    double *m_star1_array_rm;
    double *m_star0_array, *m_star0_array_rm_n;
    double *m_star3_array;
    double *m_star2_inv;

    void initial_mesh_edge_vertex_L0n();
    void initial_mesh_edge_vertex_L1n();
    
    std::vector<int> exterior_vertices;
    std::vector<int> special_vertices;
    std::vector<int> faces;
    std::vector<int> new_vertex_indices;
    std::vector<int> edge_removal;
    std::vector<int> face_removal;
    std::vector<int> edge_indices;
    std::vector<double> eig_vectors;
    std::vector<double> mdivergence;
    
    double *m_faces_size;
    bool *m_is_edge_counted,*m_is_boundary_edge, *m_is_boundary_vertex;
    
    int m_num_basis;
    double *m_vector_bases,*m_vorticity_bases; //note that m_vector_bases is normalized while vorticity_bases are not
    double *m_specific_mode_basis;
    CVector3d *m_cell_velcoity;
    double *m_cell_vf_bases; //store the vector bases as a constant vf per cell
    double *m_flux; //store the current flux
    double *m_coeff;
    double *m_cell_density;
    
    CSparseMatrix *m_d0_normal_rm;
    CSparseMatrix *m_d2_normal;
    CSparseMatrix *m_d0;
    CSparseMatrix *m_d1;
    
    CSparseMatrix *projection0;
    CSparseMatrix *projection1;
    CSparseMatrix *projection2;
    CSparseMatrix *projection3;
    
    int m_num_real_faces;
    bool *m_is_face_duplicate;

    struct volHeader {
        char     ID[3]; // VOL
        char     version; // 5
        int        encoding; // 1
        int        dimX, dimY, dimZ;
        int        channels; // 1
    };
    
    struct LVL {
        char txt[20];
        int nx,ny,nz;
        int num_estimated_faces;
        int version;
    };
    
    struct grid_point{
        int index[3]; //index for x,y,z direction
        float val; //value of the level set
    };
    float *m_grid_points;
    
    void read_level_set(std::string filename);
    void init();
    
    void initial_mesh();
    void initial_mesh_cell();
    void initial_mesh_face();
    void initial_mesh_edge_vertex();
    
    int return_face_index(int index_x,int index_y,int index_z);
    int return_edge_index(int index_x,int index_y,int index_z);
    int return_face_index_with_type(int index_x,int index_y,int index_z,int type);
    int return_point_index(int index_x,int index_y,int index_z);
    int return_cell_index(int index_x,int index_y,int index_z);
    bool is_face_included(int i,int j,int k,int type,bool &is_interior, double &face_size, int &num_np);
    void decode_face_index(int face_index,int &x,int &y,int &z,int &type);
    bool is_edge_counted(int i,int j,int k, int type,int *coeff,int *adj_face_index,double &inv_s1_value,int &p1_index,int &p2_index);
    void tag_faces_edges(int c_i,int c_j,int c_k);
    void check_cell_info();
    
    CVector3d return_grid_pos(int i,int j,int k,double grid_size,CVector3d start_pos);
    
    double return_length_ratio_included(double negative_val,double pos_val);
    double compute_area_marching_square(double *val,int &num_np);

    void matlabEIGS();
    
};

#endif /* grid_mesh_h */
