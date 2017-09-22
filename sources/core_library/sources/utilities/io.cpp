/*
    This file is part of the LibTri library.

    Author(s): Riccardo Fellegara (riccardo.fellegara@gmail.com)

    This project has been supported by the Italian Ministry of Education and
    Research under the PRIN 2009 program, and by the National Science Foundation
    under grant number IIS-1116747.

    The LibTri library is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    The LibTri library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with the LibTri library.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "io.h"
#include "utilities/string_management.h"

bool IO::read_mesh(Spatial_Mesh &mesh, string path)
{
    string extension = string_management::get_file_extension(path);
    if(extension == "tri")
        return IO::read_mesh_tri(mesh,path);
    else if(extension == "off")
        return IO::read_mesh_off(mesh,path);
    else
    {
        cerr << "[ERROR] unsopported file format. " << endl;
        return false;
    }
}

bool IO::read_mesh_off(Spatial_Mesh &mesh, string path)
{
    ifstream input(path.c_str());

    if (input.is_open() == false) {
        cerr << "Error in file " << path << "\nThe file could not exist, be unreadable or incorrect." << endl;
        return false;
    }

    string l;
    getline(input,l); // trow away the first line
    getline(input,l);
    vector<string> lt;
    string_management::tokenize(l,lt," ");
    itype num_vertices = atol(lt.at(0).c_str());
    itype num_triangles = atol(lt.at(1).c_str());

    if (num_vertices == 0 || num_triangles == 0)
    {
        cerr << "This is not a valid .off file: " << path << endl;
        return false;
    }

    mesh.reserve(num_vertices, num_triangles);

    IO::read_vertices_list(mesh,input,num_vertices);

    itype v[3], index;
    itype num_v;
    for (itype i = 0; i < num_triangles; i++)
    {
        input >> num_v;

        if(num_v != 3)
        {
            cerr << "[ERROR] the input mesh must be a pure triangle mesh. read a simplex with "<< num_v << "vertices." << endl;
            return false;
        }

        for (itype j = 0; j < num_v; j++)
        {
            input >> index;
            v[j] = index;
        }
        Triangle t = Triangle(v[0], v[1], v[2]);
        mesh.add_triangle(t);
    }

    return true;
}

bool IO::read_mesh_tri(Spatial_Mesh& mesh, string path)
{
    ifstream input(path.c_str());

    if (input.is_open() == false) {
        cerr << "Error in file " << path << "\nThe file could not exist, be unreadable or incorrect." << endl;
        return false;
    }

    string l;
    getline(input,l);
    vector<string> lt;
    string_management::tokenize(l,lt," ");
    itype num_vertices = atoi(lt.at(0).c_str());

    if (num_vertices == 0)
    {
        cerr << "This is not a valid .tri file: " << path << endl;
        return false;
    }

    mesh.reserve_vertices_space(num_vertices);

    IO::read_vertices_list(mesh,input,num_vertices);

    itype num_triangles;
    input >> num_triangles;

    if (num_triangles == 0)
    {
        cerr << "This is not a valid .tri file: " << path << endl;
        return false;
    }

    mesh.reserve_triangles_space(num_triangles);

    itype v[3], index;
    for (itype i = 0; i < num_triangles; i++)
    {
        for (int j = 0; j < 3; j++) {
            input >> index;
            v[j] = index;
        }
        Triangle t = Triangle(v[0], v[1], v[2]);
        mesh.add_triangle(t);
    }

    return true;
}

void IO::read_vertices_list(Spatial_Mesh &mesh, ifstream &input, itype num_vertices)
{
    bool is2D = false;
    string line;
    //legge i vertici aggiustando il dominio..
    for (itype i = 0; i < num_vertices; i++)
    {
        getline(input,line);
        vector<string> line_tokens;
        string_management::tokenize(line,line_tokens," ");

        if (input.eof())
            break;
        Vertex v;
        if(line_tokens.size() == 3) // 3D
            v = Vertex(atof(line_tokens[0].c_str()),atof(line_tokens[1].c_str()),atof(line_tokens[2].c_str()));
        else if(line_tokens.size() == 2) // 2D
        {
            v = Vertex(atof(line_tokens[0].c_str()),atof(line_tokens[1].c_str()));
            is2D = true;
        }
        mesh.add_vertex(v);
    }

    if(is2D)
        cerr<<"[NOTA] The points are embedded in a 2D space."<<endl;
    else
        cerr<<"[NOTA] The points are embedded in a 3D space."<<endl;
}

bool IO::write_mesh_connectivity(Spatial_Mesh& mesh, string path)
{
    stringstream ss; ss<<path<<".ia_conn";
    ofstream output(ss.str().c_str());

    output<<"IA_CONNECTIVITY"<<endl;
    output<<mesh.get_vertices_num()<<" 1 0,2,"<<mesh.get_triangles_num()<<endl;
    for(itype v=0; v<mesh.get_vertices_num(); v++)
        output<<"2,"<<mesh.get_vertex(v).get_VTstar()<<endl;
    for(itype t=0; t<mesh.get_triangles_num(); t++)
    {
        Triangle &tri = mesh.get_triangle(t);
        output<<tri.TT(0)<<" "<<tri.TT(1)<<" "<<tri.TT(2)<<endl;
    }
    output<<"0"<<endl; /// for non-manifold adjacencies
    output.close();
    return true;
}
