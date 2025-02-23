/* 
This program generates a cubed sphere mesh.
*/

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cassert>
#include <vector>
#include <algorithm>
#include <cmath>
#include "cmdline.h"

using namespace std;

typedef enum 
{
    SHELL,
    SPHEROID
} DomainType;

typedef struct 
{
    double x_min, x_max;
    double y_min, y_max;
    double z_min, z_max;
} Cube;

typedef struct 
{
    double Coords[3];
} Point;

typedef struct 
{
    int Nodes[8];
} Element;

// a map that maps (i,j,k)
class MapCube
{
    protected:
        int N_cube;
        
    public:
    
        MapCube(int N_cube)
        {
            this->N_cube = N_cube;
        }
        
        int get_index(int i, int j, int k)
        {
            return (k * (N_cube + 1) + j) * (N_cube + 1) + i;
        }
};

class MapSurface
{
    protected:
        int N_face;
        int N_radial;
        int N_cube;
    public:
        MapSurface(int N_radial, int N_cube)
        {
            this->N_face = 6;
            this->N_radial = N_radial;
            this->N_cube = N_cube;
        }
        
        int get_index(int face, int layer, int i, int j)
        {
            return layer*(
                6 * (N_cube + 1) * (N_cube + 1)
            ) + face*(N_cube + 1)*(N_cube + 1) + i*(N_cube + 1) + j;
        }
};


void mesh_in_cube(Cube cube, int N_cube, vector<Point> &nodes, vector<Element> &elements)
{
    double dx = (cube.x_max - cube.x_min) / N_cube;
    double dy = (cube.y_max - cube.y_min) / N_cube;
    double dz = (cube.z_max - cube.z_min) / N_cube;

    int N_nodes = (N_cube + 1) * (N_cube + 1) * (N_cube + 1);
    int N_elements = N_cube * N_cube * N_cube;
    
    MapCube map(N_cube);

    int node_index = 0;
    for (int k = 0; k <= N_cube; k++)
    {
        for (int j = 0; j <= N_cube; j++)
        {
            for (int i = 0; i <= N_cube; i++)
            {
                node_index = map.get_index(i, j, k);
                nodes[node_index].Coords[0] = cube.x_min + i * dx;
                nodes[node_index].Coords[1] = cube.y_min + j * dy;
                nodes[node_index].Coords[2] = cube.z_min + k * dz;
                // printf("node %d: %f %f %f\n", node_index, nodes[node_index].Coords[0], nodes[node_index].Coords[1], nodes[node_index].Coords[2]);
            }
        }
    }
    
    int elem_idx = 0;
    for (int k = 0; k < N_cube; k++)
    {
        for (int j = 0; j < N_cube; j++)
        {
            for (int i = 0; i < N_cube; i++)
            {
                elements[elem_idx].Nodes[0] = map.get_index(i, j, k);
                elements[elem_idx].Nodes[1] = map.get_index(i + 1, j, k);
                elements[elem_idx].Nodes[2] = map.get_index(i + 1, j + 1, k);
                elements[elem_idx].Nodes[3] = map.get_index(i, j + 1, k);
                elements[elem_idx].Nodes[4] = map.get_index(i, j, k + 1);
                elements[elem_idx].Nodes[5] = map.get_index(i + 1, j, k + 1);
                elements[elem_idx].Nodes[6] = map.get_index(i + 1, j + 1, k + 1);
                elements[elem_idx].Nodes[7] = map.get_index(i, j + 1, k + 1);
                // printf("element %d: %d %d %d %d %d %d %d %d\n", elem_idx, elements[elem_idx].Nodes[0], elements[elem_idx].Nodes[1], elements[elem_idx].Nodes[2], elements[elem_idx].Nodes[3], elements[elem_idx].Nodes[4], elements[elem_idx].Nodes[5], elements[elem_idx].Nodes[6], elements[elem_idx].Nodes[7]);
                elem_idx++;
            }
        }   
    }
}
                

void WriteToVTK(vector<Point> &nodes, vector<Element> &elements, string filename)
{
    ofstream out(filename);
    // set precision
    out.precision(16);
    out << "# vtk DataFile Version 3.0\n";
    out << "Cubed Sphere\n";
    out << "ASCII\n";
    out << "DATASET UNSTRUCTURED_GRID\n";
    out << "POINTS " << nodes.size() << " double\n";
    for (int i = 0; i < nodes.size(); i++)
    {
        out << nodes[i].Coords[0] << " " << nodes[i].Coords[1] << " " << nodes[i].Coords[2] << "\n";
    }
    out << "CELLS " << elements.size() << " " << elements.size() * 9 << "\n";
    for (int i = 0; i < elements.size(); i++)
    {
        out << "8 " << elements[i].Nodes[0] << " " << elements[i].Nodes[1] << " " << elements[i].Nodes[2] << " " << elements[i].Nodes[3] << " " << elements[i].Nodes[4] << " " << elements[i].Nodes[5] << " " << elements[i].Nodes[6] << " " << elements[i].Nodes[7] << "\n";
    }
    out << "CELL_TYPES " << elements.size() << "\n";
    for (int i = 0; i < elements.size(); i++)
    {
        out << "12\n";
    }
    
    out.close();
    
    printf("Wrote %s\n", filename.c_str());
}

int main(int argc, char **argv)
{
    cmdline::parser a;
    
    a.add<int>("N_cube", 'n', "number of elements on each face of the cube is N_cube x N_cube", false, 30, cmdline::range(1, 1000));
    a.add<int>("N_radial", 'm', "number of radial layers", false, 30, cmdline::range(1, 1000));
    a.add<double>("R_outer", 'R', "outer radius of the sphere", false, 1.0, cmdline::range(0.0, 1000.0));
    a.add<double>("R_inner", 'r', "inner radius of the sphere", false, 0.5, cmdline::range(0.0, 1000.0));
    a.add<double>("frac", 'f', "the length of the cube is frac * R_outer", false, 0.3, cmdline::range(0.0, 1.0));
    a.add<string>("output", 'o', "output file name", false, "cubed-sphere.vtk");
    a.add<int>("type", 't', "type of the domain", false, SPHEROID, cmdline::oneof<int>(SHELL, SPHEROID));
    a.add<double>("Eccentricity", 'E', "eccentricity of the ellipsoid", false, 0.0, cmdline::range(0.0, 1.0));
    a.add<double>("Boundary thickness", 'b', "thickness of elements near boundary", false, 0.0, cmdline::range(0.0, 1.0));
    a.add<double>("Interior thickness", 'i', "thickness of elements near interior", false, 0.0, cmdline::range(0.0, 1.0));
    a.add<bool>("Evenly spaced", 'e', "whether the elements are evenly spaced", false, false);
    
    a.parse_check(argc, argv);
    
    int N_cube = a.get<int>("N_cube");
    int N_radial = a.get<int>("N_radial");
    double R_outer = a.get<double>("R_outer");
    double R_inner = a.get<double>("R_inner");
    double frac = a.get<double>("frac");
    DomainType type = (DomainType)a.get<int>("type");
    string output = a.get<string>("output");
    double Ec = a.get<double>("Eccentricity");
    double thick = a.get<double>("Boundary thickness");
    double thick2 = a.get<double>("Interior thickness");
    bool evenly_spaced = a.get<bool>("Evenly spaced");
    
    double S_cube = type==SPHEROID? frac * R_outer * 2.0 : frac * R_inner * 2.0;
    
    int N_elem = 0;
    if(type == SPHEROID)
    {
        N_elem += N_cube * N_cube * N_cube; // Elements in the cube
        N_elem += N_radial * (N_cube * N_cube * 6); // Elements in the radial layers
    }
    else if(type == SHELL)
    {
        N_elem = N_radial * (N_cube * N_cube * 6); // Elements in the radial layers
    }
    
    int N_surface = 0;
    N_surface += 6 * (N_cube + 1) * (N_cube + 1); 
    N_surface -= 12 * (N_cube + 1); // Nodes on the 6 faces of the cube
    N_surface += 8; // 8 corners of the cube
   
    int N_node = 0;
    if(type == SPHEROID)
    {
        N_node += (N_cube + 1) * (N_cube + 1) * (N_cube + 1); // Nodes in the cube
        N_node += N_surface * N_radial; // Nodes in the radial layers
    }
    else if (type == SHELL)
    {
        N_node = N_surface * (N_radial + 1); // Nodes in the radial layers
    }
    
    Cube cube;
    cube.x_min = -S_cube / 2.0;
    cube.x_max = S_cube / 2.0;
    cube.y_min = -S_cube / 2.0;
    cube.y_max = S_cube / 2.0;
    cube.z_min = -S_cube / 2.0/ (sqrt(1-Ec*Ec));
    cube.z_max = S_cube / 2.0/ (sqrt(1-Ec*Ec));
    
    vector<Point> nodes(N_node);
    vector<Element> elements(N_elem);
    
    vector<Point> cube_nodes((N_cube + 1) * (N_cube + 1) * (N_cube + 1));
    vector<Element> cube_elements(N_cube * N_cube * N_cube);
    
    // fill the elements and nodes in the cube
    mesh_in_cube(cube, N_cube, cube_nodes, cube_elements);
    
    if(type == SPHEROID)
    {
        for(int i = 0; i < cube_nodes.size(); i++)
        {
            nodes[i] = cube_nodes[i];
        }
        for(int i = 0; i < cube_elements.size(); i++)
        {
            elements[i] = cube_elements[i];
        }
    }

    vector<double> t_radial(N_radial);
    if(evenly_spaced)
    {
        for (int i = 0; i < N_radial; i++)
        {
            t_radial[i] = (i + 1.0) / (N_radial);
        }
    }
    else
    {
        if(type==SPHEROID)
        {
            double dist = (1.0-frac)*R_outer;
            double a = (N_radial*thick)/dist -1.0;
            double b = 1-a;
            for (int i = 0; i < N_radial; i++)
            {
                double t = (i + 1.0) / (N_radial);
                t_radial[i] = a*t*t + b*t;
            }
            int N_recommend = (int)(b * dist / thick2);
            printf("N_recommend = %d for interior thickness\n", N_recommend);
        }
        else if(type==SHELL)
        {
            double c = (N_radial*thick)/(R_outer-R_inner);
            double b = 3*(1-c);
            double a = 1 - b - c;
            for (int i = 0; i < N_radial; i++)
            {
                double t = (i + 1.0) / (N_radial);
                t_radial[i] = a*t*t*t + b*t*t + c*t;
            }
            int N_recommend = (int)((3.0/4.0*a + b + c) * (R_outer-R_inner) / thick2);
            printf("N_recommend = %d for interior thickness\n", N_recommend);
        }
    }
    
    vector<int> face_node_idx;
    face_node_idx.resize((N_radial + 1) * 6 * (N_cube + 1) * (N_cube + 1));
    MapSurface map_surf(N_radial, N_cube);
    MapCube map_cube(N_cube);
    
    for(int i = 0; i < N_cube + 1; i++)
    {
        for(int j = 0; j < N_cube + 1; j++)
        {
            face_node_idx[map_surf.get_index(0, 0, i, j)] = map_cube.get_index(i, j, 0);
            face_node_idx[map_surf.get_index(1, 0, i, j)] = map_cube.get_index(i, j, N_cube);            
            face_node_idx[map_surf.get_index(2, 0, i, j)] = map_cube.get_index(i, 0, j);
            face_node_idx[map_surf.get_index(3, 0, i, j)] = map_cube.get_index(i, N_cube, j);
            face_node_idx[map_surf.get_index(4, 0, i, j)] = map_cube.get_index(0, i, j);
            face_node_idx[map_surf.get_index(5, 0, i, j)] = map_cube.get_index(N_cube, i, j);
        }
    }
    
    // get the nodes on the surface of the cube
    vector<int> nodes_surface(6*(N_cube + 1)*(N_cube + 1));
    for (int i=0; i<nodes_surface.size(); i++)
    {
        nodes_surface[i] = face_node_idx[i];
    }
    sort(nodes_surface.begin(), nodes_surface.end());
    auto ip = unique(nodes_surface.begin(), nodes_surface.end());
    nodes_surface.resize(distance(nodes_surface.begin(), ip));
    
    vector<int> cube_surface_idx((N_cube + 1) * (N_cube + 1) * (N_cube + 1), -1);
    for(int i = 0; i < nodes_surface.size(); i++)
    {
        cube_surface_idx[nodes_surface[i]] = i;
    }
    
    // inner layer of the shell
    if(type == SHELL)
    {
        for(int i = 0; i < nodes_surface.size(); i++)
        {
            Point node_on_cube = cube_nodes[nodes_surface[i]];
            double r = sqrt(node_on_cube.Coords[0] * node_on_cube.Coords[0] + node_on_cube.Coords[1] * node_on_cube.Coords[1] + node_on_cube.Coords[2] * node_on_cube.Coords[2]/(1.0-Ec*Ec));
            nodes[i].Coords[0] = node_on_cube.Coords[0] * R_inner / r;
            nodes[i].Coords[1] = node_on_cube.Coords[1] * R_inner / r;
            nodes[i].Coords[2] = node_on_cube.Coords[2] * R_inner / r;
        }
        for(int i = 0; i < N_cube + 1; i++)
        {
            for(int j = 0; j < N_cube + 1; j++)
            {
                for(int k = 0; k < 6; k++)
                {
                    face_node_idx[map_surf.get_index(k, 0, i, j)] = cube_surface_idx[face_node_idx[map_surf.get_index(k, 0, i, j)]];
                }
            }
        }
    }
    
    // fill the nodes in the radial layers
    for (int i=0; i<N_radial; i++)
    {
        int start;
        if (type == SPHEROID)
        {
            start = (N_cube + 1) * (N_cube + 1) * (N_cube + 1);
            
            for (int j=0; j<N_surface; j++)
            {
                Point node_on_cube = cube_nodes[nodes_surface[j]];
                double r = sqrt(node_on_cube.Coords[0] * node_on_cube.Coords[0] + node_on_cube.Coords[1] * node_on_cube.Coords[1] + node_on_cube.Coords[2] * node_on_cube.Coords[2]/(1.0-Ec*Ec));
                
                for(int k = 0; k < 3; k++)
                {
                    nodes[start + i * N_surface + j].Coords[k] = node_on_cube.Coords[k] / r * (R_outer - r) * t_radial[i] + node_on_cube.Coords[k];
                }
            }
        }
        else if (type == SHELL)
        {
            start = N_surface;
            for (int j=0; j<N_surface; j++)
            {
                Point node_on_inner = nodes[j];
                double r = sqrt(node_on_inner.Coords[0] * node_on_inner.Coords[0] + node_on_inner.Coords[1] * node_on_inner.Coords[1] + node_on_inner.Coords[2] * node_on_inner.Coords[2]/(1.0-Ec*Ec));

                for(int k = 0; k < 3; k++)
                {
                    nodes[start + i * N_surface + j].Coords[k] = node_on_inner.Coords[k] / r * (R_outer - R_inner) * t_radial[i] + node_on_inner.Coords[k] * R_inner / r;
                }
            }
        }
    }
    
    for(int layer = 1; layer <= N_radial; layer++)
    {
        int start;
        if (type == SPHEROID)
        {
            start = (N_cube + 1) * (N_cube + 1) * (N_cube + 1) + (layer-1) * N_surface;

            for(int i = 0; i < N_cube + 1; i++)
            {
                for(int j = 0; j < N_cube + 1; j++)
                {
                    face_node_idx[map_surf.get_index(0, layer, i, j)] = start + cube_surface_idx[map_cube.get_index(i, j, 0)];
                    face_node_idx[map_surf.get_index(1, layer, i, j)] = start + cube_surface_idx[map_cube.get_index(i, j, N_cube)];
                    face_node_idx[map_surf.get_index(2, layer, i, j)] = start + cube_surface_idx[map_cube.get_index(i, 0, j)];
                    face_node_idx[map_surf.get_index(3, layer, i, j)] = start + cube_surface_idx[map_cube.get_index(i, N_cube, j)];
                    face_node_idx[map_surf.get_index(4, layer, i, j)] = start + cube_surface_idx[map_cube.get_index(0, i, j)];
                    face_node_idx[map_surf.get_index(5, layer, i, j)] = start + cube_surface_idx[map_cube.get_index(N_cube, i, j)];
                }
            }
        }
        else if (type == SHELL)
        {
            start = layer * N_surface;
            for(int i = 0; i < N_cube + 1; i++)
            {
                for(int j = 0; j < N_cube + 1; j++)
                {
                    for(int k = 0; k < 6; k++)
                    {
                        face_node_idx[map_surf.get_index(k, layer, i, j)] = start + face_node_idx[map_surf.get_index(k, 0, i, j)];
                    }
                }
            }
        }
    }
    
    // fill the elements in the radial layers
    bool reverse = false;
    for(int i = 0; i<N_radial; i++)
    {
        for(int j = 0; j<6; j++)
        {
            for(int k = 0; k<N_cube; k++)
            {
                for(int l = 0; l<N_cube; l++)
                {
                    int elem_idx;
                    if (type == SPHEROID)
                    {
                        elem_idx = N_cube * N_cube * N_cube + i * 6 * N_cube * N_cube + j * N_cube * N_cube + k * N_cube + l;
                    }
                    else if (type == SHELL)
                    {
                        elem_idx = i * 6 * N_cube * N_cube + j * N_cube * N_cube + k * N_cube + l;
                    }
                    reverse = (j==0 || j==3 || j==4) ? true : false;
                    if(!reverse)
                    {
                        elements[elem_idx].Nodes[0] = face_node_idx[map_surf.get_index(j, i, k, l)];
                        elements[elem_idx].Nodes[1] = face_node_idx[map_surf.get_index(j, i, k + 1, l)];
                        elements[elem_idx].Nodes[2] = face_node_idx[map_surf.get_index(j, i, k + 1, l + 1)];
                        elements[elem_idx].Nodes[3] = face_node_idx[map_surf.get_index(j, i, k, l + 1)];
                        elements[elem_idx].Nodes[4] = face_node_idx[map_surf.get_index(j, i + 1, k, l)];
                        elements[elem_idx].Nodes[5] = face_node_idx[map_surf.get_index(j, i + 1, k + 1, l)];
                        elements[elem_idx].Nodes[6] = face_node_idx[map_surf.get_index(j, i + 1, k + 1, l + 1)];
                        elements[elem_idx].Nodes[7] = face_node_idx[map_surf.get_index(j, i + 1, k, l + 1)];
                    }
                    else
                    {
                        elements[elem_idx].Nodes[0] = face_node_idx[map_surf.get_index(j, i, k, l + 1)];
                        elements[elem_idx].Nodes[1] = face_node_idx[map_surf.get_index(j, i, k + 1, l + 1)];
                        elements[elem_idx].Nodes[2] = face_node_idx[map_surf.get_index(j, i, k + 1, l)];
                        elements[elem_idx].Nodes[3] = face_node_idx[map_surf.get_index(j, i, k, l)];
                        elements[elem_idx].Nodes[4] = face_node_idx[map_surf.get_index(j, i + 1, k, l + 1)];
                        elements[elem_idx].Nodes[5] = face_node_idx[map_surf.get_index(j, i + 1, k + 1, l + 1)];
                        elements[elem_idx].Nodes[6] = face_node_idx[map_surf.get_index(j, i + 1, k + 1, l)];
                        elements[elem_idx].Nodes[7] = face_node_idx[map_surf.get_index(j, i + 1, k, l)];
                    }
                }
            }
        }
    }
    
    // output
    WriteToVTK(nodes, elements, output);
    
    printf("N_node = %d\n", N_node);
    printf("N_elem = %d\n", N_elem);
    printf("N_surface_pts = %d\n", N_surface);
    printf("N_surface_elems = %d\n", 6 * N_cube * N_cube);
    
    return 0;
}