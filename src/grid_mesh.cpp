//
//  grid_mesh.cpp

#include "grid_mesh.h"
#include <iostream>
#include <fstream>
#include "MarchingCubes.h"

//read level set info
void grid_mesh::read_level_set(std::string filename){

    std::ifstream file(filename);
    file >> m_x_num_dimension;
    file >> m_y_num_dimension;
    file >> m_z_num_dimension;
    file >> m_grid_length;

    std::cout<<"grid dimensions: "<<m_x_num_dimension<<' '<<m_y_num_dimension<<' '<<m_z_num_dimension<<std::endl;

    m_grid_points = new float[m_x_num_dimension*m_y_num_dimension*m_z_num_dimension];
    
    int track_index = 0;
    //read level set and construct grid points
    for(int k=0;k<m_z_num_dimension;k++){
        for(int j=0;j<m_y_num_dimension;j++){
            for(int i=0;i<m_x_num_dimension;i++){
                double current_value;
                file >> current_value;
                CVector3d pos = CVector3d(i*m_grid_length,j*m_grid_length,k*m_grid_length);
                m_grid_points[track_index] = current_value;
                
                if(current_value > 0){
                    exterior_vertices.push_back(track_index);
                }
                track_index++;
            }
        }
    }
//    std::cout << "num vertices " << track_index << std::endl;
    file.close();
}

CVector3d grid_mesh::return_grid_pos(int i,int j,int k,double grid_length,CVector3d start_pos){
    CVector3d diff = CVector3d(i*grid_length,j*grid_length,k*grid_length);
    return (diff+start_pos);
}

//helper functions for mesh initialization

int grid_mesh::return_point_index(int index_x,int index_y,int index_z){
    //out of range
    if(index_x<0 || index_y<0 || index_z<0)
        return -1;
    if(index_x>=m_x_num_dimension || index_y>=m_y_num_dimension || index_z>=m_z_num_dimension)
        return -1;

    return (index_z*m_x_num_dimension*m_y_num_dimension + index_y*m_x_num_dimension + index_x);
}

//return the cell index, where index_x, index_y, index_z is its minimal point index
int grid_mesh::return_cell_index(int index_x,int index_y,int index_z){
    if(index_x<0 || index_x>=m_x_num_dimension-1)
        return -1;
    if(index_y<0 || index_y>=m_y_num_dimension-1)
        return -1;
    if(index_z<0 || index_z>=m_z_num_dimension-1)
        return -1;

    return (index_z*(m_y_num_dimension-1)*(m_x_num_dimension-1) +
        index_y*(m_x_num_dimension-1) + index_x);
}

//return one index of the face, where index_x, index_y, index_z is its minimal point index
int grid_mesh::return_face_index(int index_x,int index_y,int index_z){
    return 3*(index_z*m_x_num_dimension*m_y_num_dimension +
        index_y*m_x_num_dimension + index_x);
}

int grid_mesh::return_edge_index(int index_x,int index_y,int index_z){
    return 3*(index_z*m_x_num_dimension*m_y_num_dimension +
        index_y*m_x_num_dimension + index_x);
}

int grid_mesh::return_face_index_with_type(int index_x,int index_y,int index_z,int type){
    if(index_x<0 || index_y<0 || index_z<0)
        return -1;
    if(index_x>=m_x_num_dimension || index_y>=m_y_num_dimension || index_z>=m_z_num_dimension)
        return -1;
    return (return_face_index(index_x,index_y,index_z) + type);
}

void grid_mesh::decode_face_index(int face_index,int &x,int &y,int &z,int &type){
    type = face_index%3;

    int tmp1 = (face_index-type)/3;
    z = tmp1/(m_x_num_dimension*m_y_num_dimension);

    int tmp2 = tmp1 - z*m_x_num_dimension*m_y_num_dimension;
    y = tmp2/m_x_num_dimension;

    x = tmp2%m_x_num_dimension;
}

void grid_mesh::tag_faces_edges(int c_i,int c_j,int c_k){

    int adj_face_index[6];
    adj_face_index[0] = return_face_index_with_type(c_i,  c_j,  c_k,  0);
    adj_face_index[1] = return_face_index_with_type(c_i,  c_j,  c_k+1,0);
    adj_face_index[2] = return_face_index_with_type(c_i,  c_j,  c_k,  1);
    adj_face_index[3] = return_face_index_with_type(c_i,  c_j+1,c_k,  1);
    adj_face_index[4] = return_face_index_with_type(c_i,  c_j,  c_k,  2);
    adj_face_index[5] = return_face_index_with_type(c_i+1,c_j,  c_k,  2);

    int adj_edge_index[12];
    adj_edge_index[0] = return_edge_index(c_i,  c_j,  c_k  );
    adj_edge_index[1] = return_edge_index(c_i+1,c_j,  c_k  )+2;
    adj_edge_index[2] = return_edge_index(c_i,  c_j,  c_k+1);
    adj_edge_index[3] = return_edge_index(c_i,  c_j,  c_k  )+2;

    adj_edge_index[4] = return_edge_index(c_i,  c_j+1,c_k  );
    adj_edge_index[5] = return_edge_index(c_i+1,c_j+1,c_k  )+2;
    adj_edge_index[6] = return_edge_index(c_i,  c_j+1,c_k+1);
    adj_edge_index[7] = return_edge_index(c_i,  c_j+1,c_k  )+2;

    adj_edge_index[8]  = return_edge_index(c_i,  c_j,  c_k  )+1;
    adj_edge_index[9]  = return_edge_index(c_i+1,c_j,  c_k  )+1;
    adj_edge_index[10] = return_edge_index(c_i,  c_j,  c_k+1)+1;
    adj_edge_index[11] = return_edge_index(c_i+1,c_j,  c_k+1)+1;

    for(int i=0;i<6;i++)
        m_is_counted_face[adj_face_index[i]] = true;
    for(int i=0;i<12;i++)
        m_is_edge_counted[adj_edge_index[i]] = true;

}

void grid_mesh::check_cell_info(){
    //HS3 read in from file to save time, be careful that the correct file is used!
    bool read_from_file = false;
    if (read_from_file && option == '3'){
        std::ifstream file("hs3.txt"); //update
        for(int i = 0; i < m_estimated_num_cells; i++){
            double current_value;
            file >> current_value;
            m_cell_volume[i] = current_value;
        }
        file.close();
    }
    
    for(int k=0;k<m_z_num_dimension-1;k++)
        for(int j=0;j<m_y_num_dimension-1;j++)
            for(int i=0;i<m_x_num_dimension-1;i++){
                int cell_index = return_cell_index(i,j,k);

                double cell_vol;

                //data[i][j][k]
                std::vector<std::vector<std::vector<double>>> data(2, std::vector<std::vector<double> > (2, std::vector<double> (2, 0)));
                data[0][0][0] = m_grid_points[return_point_index(i,  j,  k  )];
                data[1][0][0] = m_grid_points[return_point_index(i+1,j,  k  )];
                data[0][1][0] = m_grid_points[return_point_index(i,  j+1,k  )];
                data[0][0][1] = m_grid_points[return_point_index(i,  j,  k+1)];
                data[1][1][0] = m_grid_points[return_point_index(i+1,j+1,k  )];
                data[1][0][1] = m_grid_points[return_point_index(i+1,j,  k+1)];
                data[0][1][1] = m_grid_points[return_point_index(i,  j+1,k+1)];
                data[1][1][1] = m_grid_points[return_point_index(i+1,j+1,k+1)];

                //use marching cube to compute primal volume
                if (!read_from_file && option == '3') {
                    //cell is inside model
                    if (data[0][0][0]<=0 && data[1][0][0]<=0 && data[0][1][0]<=0 && data[0][0][1]<=0 &&
                        data[1][1][0]<=0 && data[1][0][1]<=0 && data[0][1][1]<=0 && data[1][1][1]<=0){
                        m_cell_volume[cell_index] = pow(m_grid_length,3);;
                    }
                    //cell is outside model
                    else if (data[0][0][0]>0 && data[1][0][0]>0 && data[0][1][0]>0 && data[0][0][1]>0 &&
                             data[1][1][0]>0 && data[1][0][1]>0 && data[0][1][1]>0 && data[1][1][1]>0){
                        m_cell_volume[cell_index] = 0.0;
                    }
                    //cell crosses the boundary, marching cubes
                    else{
                        MarchingCubes mc;
                        mc.set_resolution(2,2,2);
                        mc.init_all(data,0.);
                        mc.run(m_grid_length,cell_vol);
                        m_cell_volume[cell_index] = cell_vol; //
                        mc.clean_temps();
                    }
                }
                //std::cout << "index " << cell_index << " vol " << cell_vol << std::endl;
                
                if(data[0][0][0]<=0){m_is_cell_included[cell_index] = true; tag_faces_edges(i,j,k); continue;}
                if(data[1][0][0]<=0){m_is_cell_included[cell_index] = true; tag_faces_edges(i,j,k); continue;}
                if(data[0][1][0]<=0){m_is_cell_included[cell_index] = true; tag_faces_edges(i,j,k); continue;}
                if(data[0][0][1]<=0){m_is_cell_included[cell_index] = true; tag_faces_edges(i,j,k); continue;}
                if(data[1][1][0]<=0){m_is_cell_included[cell_index] = true; tag_faces_edges(i,j,k); continue;}
                if(data[1][0][1]<=0){m_is_cell_included[cell_index] = true; tag_faces_edges(i,j,k); continue;}
                if(data[0][1][1]<=0){m_is_cell_included[cell_index] = true; tag_faces_edges(i,j,k); continue;}
                if(data[1][1][1]<=0){m_is_cell_included[cell_index] = true; tag_faces_edges(i,j,k); continue;}

                m_is_cell_included[cell_index] = false;

            }
}

//return the edge length included
double grid_mesh::return_length_ratio_included(double negative_val,double pos_val){
    return (-negative_val/(pos_val-negative_val));
}

double grid_mesh::compute_area_marching_square(double *val,int &num_np){
    num_np = 0;
    for(int i=0;i<4;i++)
        if(val[i]<=0)
            num_np++;

    //only one non-positive point
    if(num_np==1){
        int current_np_index = 0;
        for(int i=0;i<4;i++)
            if(val[i]<=0){
                current_np_index = i;
                break;
            }
        double l_1 = return_length_ratio_included(val[current_np_index],val[(current_np_index+1)%4]);
        double l_2 = return_length_ratio_included(val[current_np_index],val[(current_np_index+3)%4]);
        double tmp = l_1*l_2/2.;
        return tmp*(m_grid_length*m_grid_length);
    }

    //two non-positive points
    if(num_np==2){
        int current_np_index = 0;
        for(int i=0;i<4;i++){
            if(val[i]<=0){
                current_np_index = i;
                break;
            }
        }

        //case 1:
        if(val[(current_np_index+1)%4]<=0){
            double l_1 = return_length_ratio_included(val[(current_np_index+1)%4],val[(current_np_index+2)%4]);
            double l_2 = return_length_ratio_included(val[current_np_index],val[(current_np_index+3)%4]);
            return ((l_1+l_2)/2.*m_grid_length*m_grid_length);
        }
        if(val[(current_np_index+3)%4]<=0){
            double l_1 = return_length_ratio_included(val[(current_np_index+3)%4],val[(current_np_index+2)%4]);
            double l_2 = return_length_ratio_included(val[current_np_index],val[(current_np_index+1)%4]);
            return ((l_1+l_2)/2.*m_grid_length*m_grid_length);
        }

        //case 2:
        if(val[(current_np_index+2)%4]<=0){
            double l_01 = return_length_ratio_included(val[current_np_index],val[(current_np_index+1)%4]);
            double l_03 = return_length_ratio_included(val[current_np_index],val[(current_np_index+3)%4]);

            double l_23 = return_length_ratio_included(val[(current_np_index+2)%4],val[(current_np_index+3)%4]);
            double l_21 = return_length_ratio_included(val[(current_np_index+2)%4],val[(current_np_index+1)%4]);

            return (l_01*l_03+l_21*l_23)/2.*(m_grid_length*m_grid_length);
        }
    }

    //three non-positive points
    if(num_np==3){
        int current_positive_index = 0;
        for(int i=0;i<4;i++)
            if(val[i]>0){
                current_positive_index = i;
                break;
            }
        double l_1 = return_length_ratio_included(val[(current_positive_index+1)%4],val[current_positive_index]);
        double l_2 = return_length_ratio_included(val[(current_positive_index+3)%4],val[current_positive_index]);
        double tmp = l_1*l_2/2;
        return (1 - tmp)*(m_grid_length*m_grid_length);
    }
    
    return m_grid_length*m_grid_length; //
}

//check whether the face existed
bool grid_mesh::is_face_included(int i,int j,int k,int type,bool &is_interior, double &face_size, int &num_np){
    is_interior = false;
    CVector3d origin_pos = CVector3d(0.0,0.0,0.0);
    double val[4];
    CVector3d pos[4];

    if(type==0){
        //centered point index (i+0.5, j+0.5, k)
        //the face is out of range
        if(i==m_x_num_dimension-1 || j==m_y_num_dimension-1)
            return false;

        //note: in an counterclockwise order
        val[0] = m_grid_points[return_point_index(i,  j,  k)];
        val[1] = m_grid_points[return_point_index(i+1,j,  k)];
        val[2] = m_grid_points[return_point_index(i+1,j+1,k)];
        val[3] = m_grid_points[return_point_index(i,  j+1,k)];

        pos[0] = return_grid_pos(i,  j,  k,  m_grid_length,origin_pos);
        pos[1] = return_grid_pos(i+1,j,  k,  m_grid_length,origin_pos);
        pos[2] = return_grid_pos(i+1,j+1,k,  m_grid_length,origin_pos);
        pos[3] = return_grid_pos(i,  j+1,k,  m_grid_length,origin_pos);

        if(k==0 || k==m_z_num_dimension-1)
            is_interior = false;
        else{
            if(m_is_cell_included[return_cell_index(i,j,k)] && m_is_cell_included[return_cell_index(i,j,k-1)])
                is_interior = true;
            else
                is_interior = false;
        }
    
    }

    if(type==1){
        //centered point index (i+0.5, j, k+0.5)
        //the face is out of range
        if(i==m_x_num_dimension-1 || k==m_z_num_dimension-1)
            return false;

        //the values of the four points are all positive
        val[0] = m_grid_points[return_point_index(i,  j,  k  )];
        val[1] = m_grid_points[return_point_index(i,  j,  k+1)];
        val[2] = m_grid_points[return_point_index(i+1,j,  k+1)];
        val[3] = m_grid_points[return_point_index(i+1,j,  k  )];

        pos[0] = return_grid_pos(i,  j,  k,  m_grid_length,origin_pos);
        pos[1] = return_grid_pos(i,  j,  k+1,m_grid_length,origin_pos);
        pos[2] = return_grid_pos(i+1,j,  k+1,m_grid_length,origin_pos);
        pos[3] = return_grid_pos(i+1,j,  k,  m_grid_length,origin_pos);

        if(j==0 || j==m_y_num_dimension-1)
            is_interior = false;
        else{
            if(m_is_cell_included[return_cell_index(i,j,k)] && m_is_cell_included[return_cell_index(i,j-1,k)])
                is_interior = true;
            else
                is_interior = false;
        }
    }

    if(type==2){
        //centered point index (i, j+0.5, k+0.5)

        //the face is out of range
        if(j==m_y_num_dimension-1 || k==m_z_num_dimension-1)
            return false;

        //the values of the four points are all positive
        val[0] = m_grid_points[return_point_index(i,  j,  k  )];
        val[1] = m_grid_points[return_point_index(i,  j+1,k  )];
        val[2] = m_grid_points[return_point_index(i,  j+1,k+1)];
        val[3] = m_grid_points[return_point_index(i,  j,  k+1)];

        pos[0] = return_grid_pos(i,  j,  k,  m_grid_length,origin_pos);
        pos[1] = return_grid_pos(i,  j+1,k,  m_grid_length,origin_pos);
        pos[2] = return_grid_pos(i,  j+1,k+1,m_grid_length,origin_pos);
        pos[3] = return_grid_pos(i,  j,  k+1,m_grid_length,origin_pos);

        if(i==0 || i==m_x_num_dimension-1)
            is_interior = false;
        else{
            if(m_is_cell_included[return_cell_index(i,j,k)] && m_is_cell_included[return_cell_index(i-1,j,k)])
                is_interior = true;
            else
                is_interior = false;
        }
    
    }

    //if all vertices of the face are on the bdy or outside, return false (tmp in caller fn)
    if(val[0]>-1e-12 && val[1]>-1e-12 && val[2]>-1e-12 && val[3]>-1e-12){
        face_size = 0;
        return false;
    }
    //all inside
    if(val[0]<=1e-12 && val[1]<=1e-12 && val[2]<=1e-12 && val[3]<=1e-12){
        num_np = 4;
        face_size = m_grid_length*m_grid_length;
        return true;
    }

    //compute the primal face size
    face_size = compute_area_marching_square(val,num_np);
    if (face_size < 1e-12){
        face_size = 1e-12;
    }

    return true;

}

void grid_mesh::initial_mesh_cell(){
    m_is_cell_included = new bool[m_estimated_num_cells];
    m_cell_volume = new double[m_estimated_num_cells];
    m_is_counted_face = new bool[m_estimated_num_faces]; //indicate whether this face is counted in the structure
    m_is_edge_counted = new bool[m_estimated_num_edges];
    m_d2_normal = new CSparseMatrix(m_estimated_num_cells,m_estimated_num_faces);
    m_star3_array = new double[m_estimated_num_cells];
    projection3 = new CSparseMatrix(m_estimated_num_cells,m_estimated_num_cells);

    memset(m_is_cell_included,false,sizeof(bool)*m_estimated_num_cells);
    memset(m_cell_volume,0.0,sizeof(double)*m_estimated_num_cells);
    memset(m_is_counted_face,false,sizeof(bool)*m_estimated_num_faces);
    memset(m_is_edge_counted,false,sizeof(bool)*m_estimated_num_edges);
    memset(m_star3_array,0.0,sizeof(double)*m_estimated_num_cells);

    check_cell_info();
    
    int count = 0;
    double val[4];
    for(int k=0;k<m_z_num_dimension-1;k++)
        for(int j=0;j<m_y_num_dimension-1;j++)
            for(int i=0;i<m_x_num_dimension-1;i++){
                int cell_index = return_cell_index(i,j,k);
                if(!m_is_cell_included[cell_index])
                    continue;

                int adj_face_index[6];
                int coeff[6];
                adj_face_index[0] = return_face_index_with_type(i,  j,  k,  0); coeff[0] = -1;
                adj_face_index[1] = return_face_index_with_type(i,  j,  k+1,0); coeff[1] = 1;
                adj_face_index[2] = return_face_index_with_type(i,  j,  k,  1); coeff[2] = -1;
                adj_face_index[3] = return_face_index_with_type(i,  j+1,k,  1); coeff[3] = 1;
                adj_face_index[4] = return_face_index_with_type(i,  j,  k,  2); coeff[4] = -1;
                adj_face_index[5] = return_face_index_with_type(i+1,j,  k,  2); coeff[5] = 1;

                for(int l=0;l<6;l++){
                    if(adj_face_index[l]==-1)
                        continue;
                    //check face
                    if(l == 0 || l == 1){
                        //grab vertices of face
                        val[0] = m_grid_points[return_point_index(i,  j,  k)];
                        val[1] = m_grid_points[return_point_index(i+1,j,  k)];
                        val[2] = m_grid_points[return_point_index(i+1,j+1,k)];
                        val[3] = m_grid_points[return_point_index(i,  j+1,k)];
                    }
                    else if (l == 2 || l == 3){
                        val[0] = m_grid_points[return_point_index(i,  j,  k  )];
                        val[1] = m_grid_points[return_point_index(i,  j,  k+1)];
                        val[2] = m_grid_points[return_point_index(i+1,j,  k+1)];
                        val[3] = m_grid_points[return_point_index(i+1,j,  k  )];
                    }
                    else {
                        val[0] = m_grid_points[return_point_index(i,  j,  k  )];
                        val[1] = m_grid_points[return_point_index(i,  j+1,k  )];
                        val[2] = m_grid_points[return_point_index(i,  j+1,k+1)];
                        val[3] = m_grid_points[return_point_index(i,  j,  k+1)];
                    }
                    
                    m_d2_normal->add1Value(cell_index,adj_face_index[l],coeff[l]);
                }
                
                if (option == '3' && m_cell_volume[cell_index] > std::pow(1e-12,3)){ //threshold
                    m_star3_array[cell_index] = 1./m_cell_volume[cell_index];
                    projection3->add1Value(count++,cell_index,1.0);
                }
            }
}

void grid_mesh::initial_mesh_face(){
    m_is_interior_face = new bool[m_estimated_num_faces]; //indicate whether this face is an interior face or not
    m_faces_size = new double[m_estimated_num_faces]; //the true/exact face size
    m_star2_array = new double[m_estimated_num_faces];
    m_star2_inv = new double[m_estimated_num_faces];
    projection2 = new CSparseMatrix(m_estimated_num_faces,m_estimated_num_faces);

    memset(m_faces_size,0.0,sizeof(double)*m_estimated_num_faces);
    memset(m_star2_array,0.0,sizeof(double)*m_estimated_num_faces);
    memset(m_star2_inv,0.0,sizeof(double)*m_estimated_num_faces);

    m_num_interior_faces = 0;

    //initialization
    for(int i=0;i<m_estimated_num_faces;i++){
        m_is_interior_face[i] = false;
    }

    //check for each face
    int count = 0;
    for(int k=0;k<m_z_num_dimension;k++)
        for(int j=0;j<m_y_num_dimension;j++)
            for(int i=0;i<m_x_num_dimension;i++){
                int start_face_index = return_face_index(i,j,k);
                for(int type=0;type<3;type++){
                    bool is_interior;
                    double face_size;
                    int num_np = 4;
                    int current_face_index = start_face_index+type;

                    if(m_is_counted_face[current_face_index]){  //taken from cell inclusion
                        
                        bool tmp = is_face_included(i,j,k,type,is_interior,face_size,num_np);
                        m_is_interior_face[current_face_index] = is_interior;
                        m_faces_size[current_face_index] = face_size; //face_size = m_grid_length*m_grid_length if inside

                        if(is_interior) {
                            m_num_interior_faces++;
                        }
                            
                        if(tmp && face_size>1e-12){
                            m_star2_inv[current_face_index] = face_size/m_grid_length;
                            m_star2_array[current_face_index] = m_grid_length/face_size; //pick one, rm other
                            projection2->add1Value(count++,current_face_index,1.0);
                        }
                    }
                }
            }
}

//return false: inv_s1 is zero either it is because it is a boundary edge or both of the points are positive
bool grid_mesh::is_edge_counted(int i,int j,int k, int type,int *coeff,int *adj_face_index,double &inv_s1_value,int &p1_index,int &p2_index){
    p1_index = return_point_index(i,j,k);
    double val_1 = m_grid_points[p1_index];
    double val_2;
    double primal_edge_length;
    
    if(type==0){
        if(i==m_x_num_dimension-1)
            return false;

        p2_index = return_point_index(i+1,j,k);
        val_2 = m_grid_points[p2_index];
    }
    if(type==1){
        if(j==m_y_num_dimension-1)
            return false;

        p2_index = return_point_index(i,j+1,k);
        val_2 = m_grid_points[p2_index];
        
    }
    if(type==2){
        if(k==m_z_num_dimension-1)
            return false;
        

        p2_index = return_point_index(i,j,k+1);
        val_2 = m_grid_points[p2_index];
    }

    //even if this edge does exist and both of the points are non-positive,
    //the edges related with boundary face is still not included
    if(type==0){
        adj_face_index[0] = return_face_index_with_type(i,j,  k,  0); coeff[0] = 1;
        adj_face_index[1] = return_face_index_with_type(i,j-1,k,  0); coeff[1] = -1;
        adj_face_index[2] = return_face_index_with_type(i,j,  k-1,1); coeff[2] = 1;
        adj_face_index[3] = return_face_index_with_type(i,j,  k,  1); coeff[3] = -1;
    }

    if(type==1){
        adj_face_index[0] = return_face_index_with_type(i-1,j,k,  0); coeff[0] = 1;
        adj_face_index[1] = return_face_index_with_type(i,  j,k,  0); coeff[1] = -1;
        adj_face_index[2] = return_face_index_with_type(i,  j,k,  2); coeff[2] = 1;
        adj_face_index[3] = return_face_index_with_type(i,  j,k-1,2); coeff[3] = -1;
    }

    if(type==2){
        adj_face_index[0] = return_face_index_with_type(i-1,j,k,1); coeff[0] = -1;
        adj_face_index[1] = return_face_index_with_type(i,  j,k,1); coeff[1] = 1;
        adj_face_index[2] = return_face_index_with_type(i,  j,k,2); coeff[2] = -1;
        adj_face_index[3] = return_face_index_with_type(i,j-1,k,2); coeff[3] = 1;
    }


    if(val_1>=-1e-12 && val_2>=-1e-12){
        inv_s1_value = 0.0;
        return false;
    }

    if(val_1<=0 && val_2<=0)
        primal_edge_length = m_grid_length;
    else{
        if(val_1<=0 && val_2>0)
            primal_edge_length = (1-val_2/(val_2-val_1))*m_grid_length;
        else
            primal_edge_length = (1-val_1/(val_1-val_2))*m_grid_length;
        
    }
    
    inv_s1_value = primal_edge_length/(m_grid_length*m_grid_length);

    return true;
}

void grid_mesh::initial_mesh_edge_vertex(){
    m_is_boundary_edge = new bool[m_estimated_num_edges];
    m_star1_array = new double[m_estimated_num_edges];
    m_is_boundary_vertex = new bool[m_num_grid_points];
    m_star0_array = new double[m_num_grid_points];
    m_d0 = new CSparseMatrix(m_estimated_num_edges,m_num_grid_points);
    m_d1 = new CSparseMatrix(m_estimated_num_faces,m_estimated_num_edges);
    projection0 = new CSparseMatrix(m_num_grid_points,m_num_grid_points);
    projection1 = new CSparseMatrix(m_estimated_num_edges,m_estimated_num_edges);
    
    m_num_interior_edges = 0;

    memset(m_star1_array,0.0,sizeof(double)*m_estimated_num_edges);
    
    //initialization
    for(int i=0;i<m_estimated_num_edges;i++)
        m_is_boundary_edge[i] = false;
    for(int i=0;i<m_num_grid_points;i++)
        m_is_boundary_vertex[i] = false;

    //check for all the interior faces
    for(int i=0;i<m_estimated_num_faces;i++)
       if(m_is_counted_face[i] && !m_is_interior_face[i]){
            int x,y,z,type;
            decode_face_index(i,x,y,z,type);
            if(type==0){
                m_is_boundary_edge[return_edge_index(x,  y,  z  )  ] = true;
                m_is_boundary_edge[return_edge_index(x+1,y,  z  )+1] = true;
                m_is_boundary_edge[return_edge_index(x,  y+1,z  )  ] = true;
                m_is_boundary_edge[return_edge_index(x,  y,  z  )+1] = true;

                m_is_boundary_vertex[return_point_index(x,  y,  z)] = true;
                m_is_boundary_vertex[return_point_index(x+1,y,  z)] = true;
                m_is_boundary_vertex[return_point_index(x+1,y+1,z)] = true;
                m_is_boundary_vertex[return_point_index(x,  y+1,z)] = true;
            }

            if(type==1){
                m_is_boundary_edge[return_edge_index(x,  y,  z  )  ] = true;
                m_is_boundary_edge[return_edge_index(x+1,y,  z  )+2] = true;
                m_is_boundary_edge[return_edge_index(x,  y,  z+1)  ] = true;
                m_is_boundary_edge[return_edge_index(x,  y,  z  )+2] = true;


                m_is_boundary_vertex[return_point_index(x,  y,  z  )] = true;
                m_is_boundary_vertex[return_point_index(x+1,y,  z  )] = true;
                m_is_boundary_vertex[return_point_index(x+1,y,  z+1)] = true;
                m_is_boundary_vertex[return_point_index(x,  y,  z+1)] = true;
            }

            if(type==2){
                m_is_boundary_edge[return_edge_index(x,  y,  z  )+1] = true;
                m_is_boundary_edge[return_edge_index(x,  y+1,z  )+2] = true;
                m_is_boundary_edge[return_edge_index(x,  y,  z+1)+1] = true;
                m_is_boundary_edge[return_edge_index(x,  y,  z  )+2] = true;

                m_is_boundary_vertex[return_point_index(x,  y,  z  )] = true;
                m_is_boundary_vertex[return_point_index(x,  y,  z+1)] = true;
                m_is_boundary_vertex[return_point_index(x,  y+1,z+1)] = true;
                m_is_boundary_vertex[return_point_index(x,  y+1,z  )] = true;
            }
       }
    
    if(option == '0'){
        initial_mesh_edge_vertex_L0n();
    }
    else if (option == '1'){
        initial_mesh_edge_vertex_L1n();
    }
    
}

void grid_mesh::initial_mesh_edge_vertex_L0n(){
    std::vector<int> edge_removal_n; //edges to keep
    std::vector<int> vertex_removal_n; //vertices to keep
    std::vector<int> i_0n;
    std::vector<int> j_0n;
    std::vector<double> k_0n;
    
    for(int k=0;k<m_z_num_dimension;k++)
        for(int j=0;j<m_y_num_dimension;j++)
            for(int i=0;i<m_x_num_dimension;i++){
                int general_edge_id = return_edge_index(i,j,k);
                for(int type=0;type<3;type++){
                    int coeff[4],adj_face_index[4];
                    double inv_s1_value;
                    int p1_index,p2_index;
                    if(m_is_edge_counted[general_edge_id+type]){
                        bool tmp = is_edge_counted(i,j,k,type,coeff,adj_face_index,inv_s1_value,p1_index,p2_index);
                        
                        //HS1
                        if(tmp){
                            bool add = false;
                            if (m_grid_points[p1_index] < 1e-12){
                                edge_removal_n.push_back(general_edge_id+type); //keep
                                bool add = true;
                                new_vertex_indices.push_back(p1_index);
                            }
                            else if (m_grid_points[p2_index] < 1e-12){
                                if (!add)
                                    edge_removal_n.push_back(general_edge_id+type);
                                new_vertex_indices.push_back(p2_index);
                            }
                            m_star1_array[general_edge_id+type] = 1.0/inv_s1_value;
                        }

                        //L0n
                        if (tmp && !m_is_boundary_vertex[p1_index] && m_grid_points[p1_index] < 1e-12){
                            vertex_removal_n.push_back(p1_index);
                            i_0n.push_back(edge_removal_n.size()-1);  //new edge index
                            j_0n.push_back(p1_index);
                            k_0n.push_back(1.0);
                            m_star0_array[p1_index] = m_grid_length*m_grid_length*m_grid_length;
                        }
                        if (tmp && !m_is_boundary_vertex[p2_index] && m_grid_points[p2_index] < 1e-12){
                            vertex_removal_n.push_back(p2_index);
                            i_0n.push_back(edge_removal_n.size()-1);  //new edge index
                            j_0n.push_back(p2_index);
                            k_0n.push_back(-1.0);
                            m_star0_array[p2_index] = m_grid_length*m_grid_length*m_grid_length;
                        }
                    }
                }
            }
    

    //clean up
    std::sort(vertex_removal_n.begin(), vertex_removal_n.end());
    std::sort(new_vertex_indices.begin(), new_vertex_indices.end());
    vertex_removal_n.erase(unique(vertex_removal_n.begin(), vertex_removal_n.end() ), vertex_removal_n.end());
    new_vertex_indices.erase(unique(new_vertex_indices.begin(), new_vertex_indices.end() ), new_vertex_indices.end());
    
    int num_entries_n = i_0n.size();
    int i0n[num_entries_n-1];
    int j0n[num_entries_n-1];
    double k0n[num_entries_n-1];
    for(int i = 0; i < num_entries_n; i++){
        i0n[i] = i_0n[i]; //edge indices
        k0n[i] = k_0n[i];

        //find new index in vertex_removal
        std::vector<int>::iterator it = std::find(vertex_removal_n.begin(), vertex_removal_n.end(), j_0n[i]);
        if(it != vertex_removal_n.end()){
            //new vertex index for j[i] is its index in vertex_removal
            j0n[i] = std::distance(vertex_removal_n.begin(), it);
        }
    }
    
    //hs0 reindex
    m_num_vertices_new_n = vertex_removal_n.size();
    //std::cout << m_num_vertices_new_n << std::endl;
    m_star0_array_rm_n = new double[m_num_vertices_new_n];
    for(int i = 0; i < m_num_vertices_new_n; i++){
        m_star0_array_rm_n[i] = m_star0_array[vertex_removal_n[i]];
    }

    m_num_edges_new_n = edge_removal_n.size();
    
    //hs1 reindex
    m_star1_array_rm = new double[m_num_edges_new_n];
    for(int i = 0; i < m_num_edges_new_n; i++){
        m_star1_array_rm[i] = m_star1_array[edge_removal_n[i]];
    }
    
    m_d0_normal_rm = new CSparseMatrix(edge_removal_n.size(),m_num_vertices_new_n);
    m_d0_normal_rm->setValues(num_entries_n,i0n,j0n,k0n);
}

void grid_mesh::initial_mesh_edge_vertex_L1n(){
    int count_edge = 0;
    
    for(int k=0;k<m_z_num_dimension;k++)
        for(int j=0;j<m_y_num_dimension;j++)
            for(int i=0;i<m_x_num_dimension;i++){
                int general_edge_id = return_edge_index(i,j,k);
                for(int type=0;type<3;type++){
                    int coeff[4],adj_face_index[4];
                    double inv_s1_value;
                    int p1_index,p2_index;
                    if(m_is_edge_counted[general_edge_id+type]){
                        
                        bool tmp = is_edge_counted(i,j,k,type,coeff,adj_face_index,inv_s1_value,p1_index,p2_index);
                        
                        if(tmp){
                            m_star1_array[general_edge_id+type] = 1.0/inv_s1_value;
                            projection1->add1Value(count_edge++,general_edge_id+type,1.0);
                        }
                        
                        for(int l=0;l<4;l++)
                            if(adj_face_index[l]!=-1)
                                m_d1->add1Value(adj_face_index[l],general_edge_id+type,coeff[l]);
                        
                        m_d0->add1Value(general_edge_id+type,p1_index,1);
                        m_d0->add1Value(general_edge_id+type,p2_index,-1);

                    }
                }
            }
    
    //need to build the vertex projection separately
    int count_vertex = 0;
    for(int k=0;k<m_z_num_dimension;k++)
        for(int j=0;j<m_y_num_dimension;j++)
            for(int i=0;i<m_x_num_dimension;i++){
                int p_index = return_point_index(i,j,k);
                if(m_grid_points[p_index] < -1e-12) //
                    projection0->add1Value(count_vertex++,p_index,1.0);
            }
}


void grid_mesh::initial_mesh(){
    initial_mesh_cell();
    initial_mesh_face();
    initial_mesh_edge_vertex();
}

void grid_mesh::matlabEIGS(){
    using namespace matlab::engine;

    std::cout << "Start MATLAB" << std::endl;
    std::unique_ptr<MATLABEngine> mat = startMATLAB();
    matlab::data::ArrayFactory factory;
    matlab::data::ArrayDimensions ad;
    
    std::vector<double> gl = {m_grid_length}; //switch to grid length
    matlab::data::TypedArray<double> mgl = factory.createArray<std::vector<double>::iterator>({ 1, 1 }, gl.begin(), gl.end());
    
    mat->setVariable(u"grid_length", std::move(mgl)); //switch to grid length

    if (option == '0') {
        //HS0
        ad.push_back(m_num_vertices_new_n);
        ad.push_back(1);
        matlab::data::TypedArray<double> mstar0_array_rm
                = factory.createArray(ad, m_star0_array_rm_n, m_star0_array_rm_n+m_num_vertices_new_n);
        ad.clear();
        
        //HS1
        ad.push_back(m_num_edges_new_n);
        ad.push_back(1);
        matlab::data::TypedArray<double> mstar1_array_rm
                = factory.createArray(ad, m_star1_array_rm, m_star1_array_rm+m_num_edges_new_n);
        ad.clear();
        
        //D0
        int numEld0 = 0;
        CMatrixElement* theElem;
        for(int i = 0; i < m_d0_normal_rm->numRows; i++){
            for(theElem = m_d0_normal_rm->rowList[i]; theElem != NULL; theElem =
                theElem->rowNext){
                    numEld0++;
            }
        }
        double rowd0[numEld0];
        double cold0[numEld0];
        double vald0[numEld0];
        int count = 0;
        for(int i = 0; i < m_d0_normal_rm->numRows; i++){
            for(theElem = m_d0_normal_rm->rowList[i]; theElem != NULL; theElem =
                theElem->rowNext){
                rowd0[count] = theElem->i+1;
                cold0[count] = theElem->j+1;
                vald0[count] = theElem->value;
                count++;
            }
        }
        
        ad.push_back(numEld0);
        ad.push_back(1);
        matlab::data::TypedArray<double> mrowd0 = factory.createArray(ad, rowd0, rowd0+numEld0);
        matlab::data::TypedArray<double> mcold0 = factory.createArray(ad, cold0, cold0+numEld0);
        matlab::data::TypedArray<double> mvald0 = factory.createArray(ad, vald0, vald0+numEld0);
        ad.clear();
        std::vector<double> vsized0 = {static_cast<double>(m_d0_normal_rm->numRows), static_cast<double>(m_d0_normal_rm->numCols)};
        matlab::data::TypedArray<double> msize0 = factory.createArray<std::vector<double>::iterator>({1, 2}, vsized0.begin(), vsized0.end());
        
        mat->setVariable(u"star0_array_rm", std::move(mstar0_array_rm));
        mat->setVariable(u"star1_array_rm", std::move(mstar1_array_rm));
        
        mat->setVariable(u"row0", std::move(mrowd0));
        mat->setVariable(u"col0", std::move(mcold0));
        mat->setVariable(u"val0", std::move(mvald0));
        mat->setVariable(u"size0", std::move(msize0));
    }
    
    if (option == '1'){
        //is_interior_face
        ad.push_back(m_estimated_num_faces);
        ad.push_back(1);
        matlab::data::TypedArray<bool> mis_interior_face
                = factory.createArray(ad, m_is_interior_face, m_is_interior_face+m_estimated_num_faces);
        ad.clear();
        
        //is_boundary_edge
        ad.push_back(m_estimated_num_edges);
        ad.push_back(1);
        matlab::data::TypedArray<bool> mis_boundary_edge
                = factory.createArray(ad, m_is_boundary_edge, m_is_boundary_edge+m_estimated_num_edges);
        ad.clear();
        
        //is_edge_counted
        ad.push_back(m_estimated_num_edges);
        ad.push_back(1);
        matlab::data::TypedArray<bool> mis_edge_counted
                = factory.createArray(ad, m_is_edge_counted, m_is_edge_counted+m_estimated_num_edges);
        ad.clear();
        
        //HS1
        ad.push_back(m_estimated_num_edges);
        ad.push_back(1);
        matlab::data::TypedArray<double> mstar1_array
                = factory.createArray(ad, m_star1_array, m_star1_array+m_estimated_num_edges);
        ad.clear();
        
        //HS2
        ad.push_back(m_estimated_num_faces);
        ad.push_back(1);
        matlab::data::TypedArray<double> mstar2_array
                = factory.createArray(ad, m_star2_array, m_star2_array+m_estimated_num_faces);
        ad.clear();
        
        //D0
        int numEld0 = 0;
        CMatrixElement* theElem;
        for(int i = 0; i < m_d0->numRows; i++){
            for(theElem = m_d0->rowList[i]; theElem != NULL; theElem =
                theElem->rowNext){
                    numEld0++;
            }
        }
        double rowd0[numEld0];
        double cold0[numEld0];
        double vald0[numEld0];
        int count = 0;
        for(int i = 0; i < m_d0->numRows; i++){
            for(theElem = m_d0->rowList[i]; theElem != NULL; theElem =
                theElem->rowNext){
                rowd0[count] = theElem->i+1;
                cold0[count] = theElem->j+1;
                vald0[count] = theElem->value;
                count++;
            }
        }
        
        //D1
        int numEld1 = 0;
        for(int i = 0; i < m_d1->numRows; i++){
            for(theElem = m_d1->rowList[i]; theElem != NULL; theElem =
                theElem->rowNext){
                    numEld1++;
            }
        }
        double rowd1[numEld1];
        double cold1[numEld1];
        double vald1[numEld1];
        count = 0;
        for(int i = 0; i < m_d1->numRows; i++){
            for(theElem = m_d1->rowList[i]; theElem != NULL; theElem = theElem->rowNext){
                rowd1[count] = theElem->i+1;
                cold1[count] = theElem->j+1;
                vald1[count] = theElem->value;
                count++;
            }
        }
        
        //P0
        int numElp0 = 0;
        for(int i = 0; i < projection0->numRows; i++){
            for(theElem = projection0->rowList[i]; theElem != NULL; theElem =
                theElem->rowNext){
                numElp0++;
            }
        }
        double rowp0[numElp0];
        double colp0[numElp0];
        double valp0[numElp0];
        count = 0;
        for(int i = 0; i < projection0->numRows; i++){
            for(theElem = projection0->rowList[i]; theElem != NULL; theElem =
                theElem->rowNext){
                rowp0[count] = theElem->i+1;
                colp0[count] = theElem->j+1;
                valp0[count] = theElem->value;
                count++;
            }
        }
        
        //P1
        int numElp1 = 0;
        for(int i = 0; i < projection1->numRows; i++){
            for(theElem = projection1->rowList[i]; theElem != NULL; theElem =
                theElem->rowNext){
                numElp1++;
            }
        }
        double rowp1[numElp1];
        double colp1[numElp1];
        double valp1[numElp1];
        count = 0;
        for(int i = 0; i < projection1->numRows; i++){
            for(theElem = projection1->rowList[i]; theElem != NULL; theElem = theElem->rowNext){
                rowp1[count] = theElem->i+1;
                colp1[count] = theElem->j+1;
                valp1[count] = theElem->value;
                count++;
            }
        }
        
        int numElp2 = 0;
        for(int i = 0; i < projection2->numRows; i++){
            for(theElem = projection2->rowList[i]; theElem != NULL; theElem =
                theElem->rowNext){
                    numElp2++;
            }
        }
        double rowp2[numElp2];
        double colp2[numElp2];
        double valp2[numElp2];
        count = 0;
        for(int i = 0; i < projection2->numRows; i++){
            for(theElem = projection2->rowList[i]; theElem != NULL; theElem =
                theElem->rowNext){
                rowp2[count] = theElem->i+1;
                colp2[count] = theElem->j+1;
                valp2[count] = theElem->value;
                count++;
            }
        }

        ad.push_back(numEld0);
        ad.push_back(1);
        matlab::data::TypedArray<double> mrowd0 = factory.createArray(ad, rowd0, rowd0+numEld0);
        matlab::data::TypedArray<double> mcold0 = factory.createArray(ad, cold0, cold0+numEld0);
        matlab::data::TypedArray<double> mvald0 = factory.createArray(ad, vald0, vald0+numEld0);
        ad.clear();
        std::vector<double> vsized0 = {static_cast<double>(m_d0->numRows), static_cast<double>(m_d0->numCols)};
        matlab::data::TypedArray<double> msized0 = factory.createArray<std::vector<double>::iterator>({1, 2}, vsized0.begin(), vsized0.end());

        ad.push_back(numEld1);
        ad.push_back(1);
        matlab::data::TypedArray<double> mrowd1 = factory.createArray(ad, rowd1, rowd1+numEld1);
        matlab::data::TypedArray<double> mcold1 = factory.createArray(ad, cold1, cold1+numEld1);
        matlab::data::TypedArray<double> mvald1 = factory.createArray(ad, vald1, vald1+numEld1);
        ad.clear();
        std::vector<double> vsized1 = {static_cast<double>(m_d1->numRows), static_cast<double>(m_d1->numCols)};
        matlab::data::TypedArray<double> msized1 = factory.createArray<std::vector<double>::iterator>({ 1, 2 }, vsized1.begin(), vsized1.end());
        
        ad.push_back(numElp0);
        ad.push_back(1);
        matlab::data::TypedArray<double> mrowp0 = factory.createArray(ad, rowp0, rowp0+numElp0);
        matlab::data::TypedArray<double> mcolp0 = factory.createArray(ad, colp0, colp0+numElp0);
        matlab::data::TypedArray<double> mvalp0 = factory.createArray(ad, valp0, valp0+numElp0);
        ad.clear();
        std::vector<double> vsizep0 = {static_cast<double>(projection0->numRows), static_cast<double>(projection0->numCols)};
        matlab::data::TypedArray<double> msizep0 = factory.createArray<std::vector<double>::iterator>({ 1, 2 }, vsizep0.begin(), vsizep0.end());
        
        ad.push_back(numElp1);
        ad.push_back(1);
        matlab::data::TypedArray<double> mrowp1 = factory.createArray(ad, rowp1, rowp1+numElp1);
        matlab::data::TypedArray<double> mcolp1 = factory.createArray(ad, colp1, colp1+numElp1);
        matlab::data::TypedArray<double> mvalp1 = factory.createArray(ad, valp1, valp1+numElp1);
        ad.clear();
        std::vector<double> vsizep1 = {static_cast<double>(projection1->numRows), static_cast<double>(projection1->numCols)};
        matlab::data::TypedArray<double> msizep1 = factory.createArray<std::vector<double>::iterator>({ 1, 2 }, vsizep1.begin(), vsizep1.end());
        
        ad.push_back(numElp2);
        ad.push_back(1);
        matlab::data::TypedArray<double> mrowp2 = factory.createArray(ad, rowp2, rowp2+numElp2);
        matlab::data::TypedArray<double> mcolp2 = factory.createArray(ad, colp2, colp2+numElp2);
        matlab::data::TypedArray<double> mvalp2 = factory.createArray(ad, valp2, valp2+numElp2);
        ad.clear();
        std::vector<double> vsizep2 = {static_cast<double>(projection2->numRows), static_cast<double>(projection2->numCols)};
        matlab::data::TypedArray<double> msizep2 = factory.createArray<std::vector<double>::iterator>({1, 2}, vsizep2.begin(), vsizep2.end());
        
        mat->setVariable(u"star1_array", std::move(mstar1_array));
        mat->setVariable(u"star2_array", std::move(mstar2_array));
        
        mat->setVariable(u"rowd0", std::move(mrowd0));
        mat->setVariable(u"cold0", std::move(mcold0));
        mat->setVariable(u"vald0", std::move(mvald0));
        mat->setVariable(u"sized0", std::move(msized0));
        
        mat->setVariable(u"rowd1", std::move(mrowd1));
        mat->setVariable(u"cold1", std::move(mcold1));
        mat->setVariable(u"vald1", std::move(mvald1));
        mat->setVariable(u"sized1", std::move(msized1));
        
        mat->setVariable(u"rowp0", std::move(mrowp0));
        mat->setVariable(u"colp0", std::move(mcolp0));
        mat->setVariable(u"valp0", std::move(mvalp0));
        mat->setVariable(u"sizep0", std::move(msizep0));
        
        mat->setVariable(u"rowp1", std::move(mrowp1));
        mat->setVariable(u"colp1", std::move(mcolp1));
        mat->setVariable(u"valp1", std::move(mvalp1));
        mat->setVariable(u"sizep1", std::move(msizep1));
        
        mat->setVariable(u"rowp2", std::move(mrowp2));
        mat->setVariable(u"colp2", std::move(mcolp2));
        mat->setVariable(u"valp2", std::move(mvalp2));
        mat->setVariable(u"sizep2", std::move(msizep2));
    }
    
    if (option == '3'){
        //HS2
        ad.push_back(m_estimated_num_faces);
        ad.push_back(1);
        matlab::data::TypedArray<double> mstar2_inv
            = factory.createArray(ad, m_star2_inv, m_star2_inv+m_estimated_num_faces);
        ad.clear();
        
        //HS3
        ad.push_back(m_estimated_num_cells);
        ad.push_back(1);
        matlab::data::TypedArray<double> mstar3_array
            = factory.createArray(ad, m_star3_array, m_star3_array+m_estimated_num_cells);
        ad.clear();
        
        //D2
        int numEld2 = 0;
        CMatrixElement* theElem;
        for(int i = 0; i < m_d2_normal->numRows; i++){
            for(theElem = m_d2_normal->rowList[i]; theElem != NULL; theElem =
                theElem->rowNext){
                    numEld2++;
            }
        }
        double rowd2[numEld2];
        double cold2[numEld2];
        double vald2[numEld2];
        int count = 0;
        for(int i = 0; i < m_d2_normal->numRows; i++){
            for(theElem = m_d2_normal->rowList[i]; theElem != NULL; theElem =
                theElem->rowNext){
                rowd2[count] = theElem->i+1;
                cold2[count] = theElem->j+1;
                vald2[count] = theElem->value;
                count++;
            }
        }
        
        //Projection matrix
        int numElp = 0;
        for(int i = 0; i < projection3->numRows; i++){
            for(theElem = projection3->rowList[i]; theElem != NULL; theElem =
                theElem->rowNext){
                    numElp++;
            }
        }
        double rowp[numElp];
        double colp[numElp];
        double valp[numElp];
        count = 0;
        for(int i = 0; i < projection3->numRows; i++){
            for(theElem = projection3->rowList[i]; theElem != NULL; theElem =
                theElem->rowNext){
                rowp[count] = theElem->i+1;
                colp[count] = theElem->j+1;
                valp[count] = theElem->value;
                count++;
            }
        }
        
        int numElp2 = 0;
        for(int i = 0; i < projection2->numRows; i++){
            for(theElem = projection2->rowList[i]; theElem != NULL; theElem =
                theElem->rowNext){
                    numElp2++;
            }
        }
        double rowp2[numElp2];
        double colp2[numElp2];
        double valp2[numElp2];
        count = 0;
        for(int i = 0; i < projection2->numRows; i++){
            for(theElem = projection2->rowList[i]; theElem != NULL; theElem =
                theElem->rowNext){
                rowp2[count] = theElem->i+1;
                colp2[count] = theElem->j+1;
                valp2[count] = theElem->value;
                count++;
            }
        }
        
        ad.push_back(numEld2);
        ad.push_back(1);
        matlab::data::TypedArray<double> mrowd2 = factory.createArray(ad, rowd2, rowd2+numEld2);
        matlab::data::TypedArray<double> mcold2 = factory.createArray(ad, cold2, cold2+numEld2);
        matlab::data::TypedArray<double> mvald2 = factory.createArray(ad, vald2, vald2+numEld2);
        ad.clear();
        std::vector<double> vsized2 = {static_cast<double>(m_d2_normal->numRows), static_cast<double>(m_d2_normal->numCols)};
        matlab::data::TypedArray<double> msize2 = factory.createArray<std::vector<double>::iterator>({1, 2}, vsized2.begin(), vsized2.end());
        
        ad.push_back(numElp2);
        ad.push_back(1);
        matlab::data::TypedArray<double> mrowp2 = factory.createArray(ad, rowp2, rowp2+numElp2);
        matlab::data::TypedArray<double> mcolp2 = factory.createArray(ad, colp2, colp2+numElp2);
        matlab::data::TypedArray<double> mvalp2 = factory.createArray(ad, valp2, valp2+numElp2);
        ad.clear();
        std::vector<double> vsizep2 = {static_cast<double>(projection2->numRows), static_cast<double>(projection2->numCols)};
        matlab::data::TypedArray<double> msizep2 = factory.createArray<std::vector<double>::iterator>({1, 2}, vsizep2.begin(), vsizep2.end());
        
        ad.push_back(numElp);
        ad.push_back(1);
        matlab::data::TypedArray<double> mrowp3 = factory.createArray(ad, rowp, rowp+numElp);
        matlab::data::TypedArray<double> mcolp3 = factory.createArray(ad, colp, colp+numElp);
        matlab::data::TypedArray<double> mvalp3 = factory.createArray(ad, valp, valp+numElp);
        ad.clear();
        std::vector<double> vsizep = {static_cast<double>(projection3->numRows), static_cast<double>(projection3->numCols)};
        matlab::data::TypedArray<double> msizep3 = factory.createArray<std::vector<double>::iterator>({1, 2}, vsizep.begin(), vsizep.end());

        mat->setVariable(u"star2_inv", std::move(mstar2_inv));
        mat->setVariable(u"star3_array", std::move(mstar3_array));
        
        mat->setVariable(u"row2", std::move(mrowd2));
        mat->setVariable(u"col2", std::move(mcold2));
        mat->setVariable(u"val2", std::move(mvald2));
        mat->setVariable(u"size2", std::move(msize2));
        
        mat->setVariable(u"rowp2", std::move(mrowp2));
        mat->setVariable(u"colp2", std::move(mcolp2));
        mat->setVariable(u"valp2", std::move(mvalp2));
        mat->setVariable(u"sizep2", std::move(msizep2));
        
        mat->setVariable(u"rowp3", std::move(mrowp3));
        mat->setVariable(u"colp3", std::move(mcolp3));
        mat->setVariable(u"valp3", std::move(mvalp3));
        mat->setVariable(u"sizep3", std::move(msizep3));
    }
    
    std::cout << "Calculate eigenvalues" << std::endl;
    
    //options: '0' L0n, '1' L1n, '3' L3n
    
    // L0n
    if (option == '0'){
        mat->eval(u"my_s1_rm=spdiags(star1_array_rm(:),0,length(star1_array_rm),length(star1_array_rm));");
        mat->eval(u"s0_matrix=spdiags(star0_array_rm(:),0,length(star0_array_rm),length(star0_array_rm));");
        mat->eval(u"D0=sparse(row0,col0,val0,size0(1),size0(2));");
        mat->eval(u"L0_H=D0'*my_s1_rm*D0;");
        mat->eval(u"L0_H=L0_H+speye(size0(2))*1e-12;");
        mat->eval(u"[V0_Hn,D0_Hn]=eigs(L0_H, s0_matrix, 40, 'smallestabs');");
        
        //graph
        mat->eval(u"L0_G=D0'*D0;");
        mat->eval(u"L0_G=L0_G+speye(size0(2))*1e-12;");
        mat->eval(u"[V0_Gn,D0_Gn]=eigs(L0_G, speye(size0(2))*grid_length*grid_length, 40, 'smallestabs');");
    }
    // L1n
    if (option == '1'){
        mat->eval(u"s0_inv_value=1.0/(grid_length*grid_length*grid_length);");
        mat->eval(u"s0_inv=s0_inv_value*speye(sized0(2));");
        mat->eval(u"s1=spdiags(star1_array(:),0,length(star1_array),length(star1_array));");
        mat->eval(u"s2=spdiags(star2_array(:),0,length(star2_array),length(star2_array));");
        mat->eval(u"D0=sparse(rowd0,cold0,vald0,sized0(1),sized0(2));");
        mat->eval(u"D1=sparse(rowd1,cold1,vald1,sized1(1),sized1(2));");
        mat->eval(u"P0=sparse(rowp0,colp0,valp0,sizep0(1),sizep0(2));");
        mat->eval(u"P1=sparse(rowp1,colp1,valp1,sizep1(1),sizep1(2));");
        mat->eval(u"P2=sparse(rowp2,colp2,valp2,sizep2(1),sizep2(2));");
        
        mat->eval(u"s0_inv_p=P0*s0_inv*P0';");
        mat->eval(u"s1_p=P1*s1*P1';");
        mat->eval(u"s2_p=P2*s2*P2';");
        mat->eval(u"D0_p = P1*D0*P0';");
        mat->eval(u"D1_p = P2*D1*P1';");
        
        mat->eval(u"L1_right=s1_p*D0_p*s0_inv_p*D0_p'*s1_p;"); // s1 D0 s0^-1 D0^T s1
        mat->eval(u"L1_left=D1_p'*s2_p*D1_p;"); // D1^T S2 D1
        mat->eval(u"L1_H=L1_left+L1_right;"); //Hodge L1n
        mat->eval(u"L1_H=(L1_H+L1_H')/2.;");
        mat->eval(u"L1_H=L1_H+speye(size(L1_H,1))*1e-10;");
        mat->eval(u"[V1_Hn,D1_Hn]=eigs(L1_H, s1_p, 40, 'sm');");
        mat->eval(u"D1_Hn=real(D1_Hn);");
        mat->eval(u"V1_Hn=real(V1_Hn);");
        
        mat->eval(u"L1_G=D1_p'*D1_p+D0_p*D0_p';"); //Graph L1n
        mat->eval(u"L1_G=L1_G+speye(size(L1_G,1))*1e-12;");
        mat->eval(u"[V1_Gn,D1_Gn]=eigs(L1_G, P1*speye(size(L1_G,1))*grid_length*grid_length*P1', 40, 'sm');");
    }
    // L3n
    if (option == '3'){
        mat->eval(u"my_s2_inv=spdiags(star2_inv(:),0,length(star2_inv),length(star2_inv));");
        mat->eval(u"my_s3=spdiags(star3_array(:),0,length(star3_array),length(star3_array));");
        mat->eval(u"D2=sparse(row2,col2,val2,size2(1),size2(2));");
        mat->eval(u"P2=sparse(rowp2,colp2,valp2,sizep2(1),sizep2(2));");
        mat->eval(u"P3=sparse(rowp3,colp3,valp3,sizep3(1),sizep3(2));");
        
        mat->eval(u"s2_inv_p=P2*my_s2_inv*P2';");
        mat->eval(u"s3_p=P3*my_s3*P3';");
        mat->eval(u"D2_p=P3*D2*P2';");
    
        mat->eval(u"L3_H=s3_p*D2_p*s2_inv_p*D2_p'*s3_p;"); //S3 D2 S2^-1 D2^T S3
        mat->eval(u"L3_H=L3_H+speye(size(L3_H,1))*1e-12;");
        mat->eval(u"[V3_Hn,D3_Hn]=eigs(L3_H, P3*(my_s3+speye(size(my_s3,1))*1e-12)*P3', 100, 'sm');");

        //graph
        mat->eval(u"L3_G=P3*P3'*D2_p*D2_p'*P3*P3';");
        mat->eval(u"L3_G=L3_G+speye(size(L3_G,1))*1e-12;");
        mat->eval(u"[V3_Gn,D3_Gn]=eigs(L3_G, P3*(speye(size(my_s3,1))*(grid_length*grid_length)+speye(size(my_s3,1))*1e-12)*P3', 100, 'sm');");
    }

    std::cout << "Write eigenvalues to file" << std::endl;

    //L0n output
    if (option == '0'){
        matlab::data::TypedArray<double> d0hn = mat->getVariable(u"D0_Hn");
        matlab::data::TypedArray<double> d0gn = mat->getVariable(u"D0_Gn");
        std::vector<double> eval0hn;
        std::vector<double> eval0gn;
        eval0hn.reserve(d0hn.getDimensions()[0]);
        eval0gn.reserve(d0gn.getDimensions()[0]);
        for(int i=0; i<d0hn.getDimensions()[0]; ++i){
            eval0hn.push_back(d0hn[i][i]);
            eval0gn.push_back(d0gn[i][i]);
        }
        std::string fname0_hn = "ev0_Hn_" + std::to_string(m_grid_length) + "_" + model + ".txt";
        std::string fname0_gn = "ev0_Gn_" + std::to_string(m_grid_length) + "_" + model + ".txt";
        std::ofstream eigval0_Hn(fname0_hn);
        std::ofstream eigval0_Gn(fname0_gn);
        for(int i=0; i<d0hn.getDimensions()[0]; i++){
            eigval0_Hn << eval0hn[i] << std::endl;
            eigval0_Gn << eval0gn[i] << std::endl;
        }
        eigval0_Hn.close();
        eigval0_Gn.close();
    }
    //L1n output
    if (option == '1'){
        matlab::data::TypedArray<double> d1hn = mat->getVariable(u"D1_Hn"); //
        matlab::data::TypedArray<double> d1gn = mat->getVariable(u"D1_Gn");
        
        std::vector<double> eval1hn, eval1gn;
        eval1hn.reserve(d1hn.getDimensions()[0]);
        eval1gn.reserve(d1gn.getDimensions()[0]);
        for(int i=0; i<d1gn.getDimensions()[0]; ++i){ //assumes have same num of eigenvalues
            eval1hn.push_back(d1hn[i][i]);
            eval1gn.push_back(d1gn[i][i]);
        }
        std::sort(eval1hn.begin(), eval1hn.end());
        std::sort(eval1gn.begin(), eval1gn.end());
        std::string fname1_hn = "ev1_Hn_" + std::to_string(m_grid_length) + "_" + model + ".txt";
        std::string fname1_gn = "ev1_Gn_" + std::to_string(m_grid_length) + "_" + model + ".txt";
        std::ofstream eigval1_Hn(fname1_hn);
        std::ofstream eigval1_Gn(fname1_gn);
        for(int i=0; i<d1gn.getDimensions()[0]; i++){
            eigval1_Hn << eval1hn[i] << std::endl;
            eigval1_Gn << eval1gn[i] << std::endl;
        }
        eigval1_Hn.close();
        eigval1_Gn.close();
    }
    //L3n output
    if (option == '3'){
        matlab::data::TypedArray<double> d3hn = mat->getVariable(u"D3_Hn");
        matlab::data::TypedArray<double> d3gn = mat->getVariable(u"D3_Gn");
        std::vector<double> eval3hn, eval3gn;
        eval3hn.reserve(d3hn.getDimensions()[0]);
        eval3gn.reserve(d3gn.getDimensions()[0]);
        for(int i=0; i<d3gn.getDimensions()[0]; ++i){ //assumes have same num of eigenvalues
            eval3hn.push_back(d3hn[i][i]);
            eval3gn.push_back(d3gn[i][i]);
        }
        std::string fname3_hn = "ev3_Hn_" + std::to_string(m_grid_length) + "_" + model + ".txt";
        std::string fname3_gn = "ev3_Gn_" + std::to_string(m_grid_length) + "_" + model + ".txt";
        std::ofstream eigval3_Hn(fname3_hn);
        std::ofstream eigval3_Gn(fname3_gn);
        for(int i=0; i<d3gn.getDimensions()[0]; i++){
                eigval3_Hn << eval3hn[i] << std::endl;
                eigval3_Gn << eval3gn[i] << std::endl;
        }
        eigval3_Hn.close();
        eigval3_Gn.close();
    }
}

void grid_mesh::init(){
    
    m_estimated_num_faces = 3*m_x_num_dimension*m_y_num_dimension*m_z_num_dimension;
    m_estimated_num_cells = (m_x_num_dimension-1)*(m_y_num_dimension-1)*(m_z_num_dimension-1);
    m_estimated_num_edges = m_estimated_num_faces;
    m_num_grid_points = m_x_num_dimension*m_y_num_dimension*m_z_num_dimension;

    std::cout << "num vertices: " << m_num_grid_points << std::endl;
    std::cout << "num edges/faces: " << m_estimated_num_edges << std::endl;
    std::cout << "num cells: " << m_estimated_num_cells << std::endl;
    
    initial_mesh();
    
    matlabEIGS();
    
}
