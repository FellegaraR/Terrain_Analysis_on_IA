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

#include "quad_mesh.h"

#include <sstream>
#include <fstream>

void Quad_Mesh::loadQuad(Spatial_Mesh &mesh)
{
    map<dvect, itype> vertici;

    dvect vert;
    itype num_vert=0;
    for(itype i=0; i<mesh.get_triangles_num(); i++)
    {
        Triangle &tri = mesh.get_triangle(i);

        //crea vertici
        dvect barycenter = compute_barycenter(i,mesh); //baricentro del tetraedro


        pair<map<dvect, itype>::iterator,bool> ret = vertici.insert(make_pair(barycenter,num_vert));
        if(ret.second)
            num_vert++;

        dvect original;
        for(int pos = 0; pos<tri.vertices_num(); pos++)
        {
            itype v = tri.TV(pos);
            original = mesh.get_vertex(v).get_coordinates();

            pair<map<dvect, itype>::iterator,bool> ret = vertici.insert(make_pair(original,num_vert));
            if(ret.second)
                num_vert++;

            ivect vertici_quad; vertici_quad.assign(4,-1);
            vertici_quad[0] = ret.first->second;
            vertici_quad[2] = vertici.find(barycenter)->second;

            int count = 1;

            for(int pos2 = 0; pos2<tri.vertices_num(); pos2++)
            {
                if(pos == pos2) continue;

                vert = compute_barycenter(v, tri.TV(pos2), mesh);
                ret = vertici.insert(make_pair(vert,num_vert));
                if(ret.second)
                    num_vert++;
                vertici_quad[count] = ret.first->second;
                count += 2;
            }

            QuadTri q = QuadTri(vertici_quad, v, i);
            this->quads.push_back(q);
        }
    }

    vertices = vector<QuadVertex>(vertici.size());
    for(map<dvect, itype>::iterator it = vertici.begin(); it != vertici.end(); it++){
        vertices[it->second].setCoords(it->first);
    }
}

dvect Quad_Mesh::compute_barycenter(int t, Spatial_Mesh &mesh)
{
    Triangle &tri = mesh.get_triangle(t);

    dvect coords = mesh.get_vertex(tri.TV(0)).get_coordinates();
    for(int pos = 1; pos<tri.vertices_num(); pos++)
    {
        Vertex &v=mesh.get_vertex(tri.TV(pos));
        for(int i=0; i<v.get_dimension(); i++)
            coords[i] += v.get_c(i);
    }

    coords[0] = coords[0]/3.0;
    coords[1] = coords[1]/3.0;
    coords[2] = coords[2]/3.0;

    return coords;
}

dvect Quad_Mesh::compute_barycenter(int v1, int v2, Spatial_Mesh &mesh)
{
    dvect coords;
    Vertex &v_1=mesh.get_vertex(v1);
    Vertex &v_2=mesh.get_vertex(v2);
    for(int i=0; i<v_1.get_dimension(); i++)
        coords.push_back((v_1.get_c(i)+v_2.get_c(i))/2.0);
    return coords;
}

void Quad_Mesh::save_quad_mesh(string file_name)
{
    cout<<"save the quad mesh with V "<<this->vertices.size()<<" and H "<<this->quads.size()<<endl;
    std::stringstream ss;
    ss << file_name << "_quad.off";
    ofstream output(ss.str().c_str());
    output << "OFF" << endl;
    output<<this->vertices.size()<<" "<<this->quads.size()<<" 0"<<endl;
    for(uint i=0; i<this->vertices.size(); i++)
    {
        for(int v=0; v<3; v++)
            output << this->vertices[i].get_i_coord(v) << " ";
        output << endl;
    }

    for(uint i=0; i<this->quads.size(); i++)
    {
        output << "4 ";
        for(int v=0; v<this->quads[i].get_num_vertices(); v++)
            output << this->quads[i].getIndex(v) << " ";
        output << endl;
    }
    output.close();
}

//void Quad_Mesh::write_vtk(char* filename)
//{
//    std::stringstream ss;
//    ss << filename << ".vtk";

//    FILE* file;
//    file = fopen(ss.str().c_str(), "w");

//    fprintf(file, "# vtk DataFile Version 2.0\n\n");
//    fprintf(file, "ASCII \n");
//    fprintf(file, "DATASET UNSTRUCTURED_GRID\n\n");
//    fprintf(file, "POINTS %d double\n", (int)vertices.size());

//    for(uint i=0; i<vertices.size(); i++){
//        //        assert(vertices[i].getCoords()->size() == 4);
//        for(int j=0; j<3; j++)
//            fprintf(file, "%f ", vertices[i].get_i_coord(j));
//        fprintf(file, "\n");
//    }
//    fprintf(file, "\n\n");

//    fprintf(file, "CELLS %d %d\n", (int)this->quads.size(),(int)quads.size()*9);
//    for(uint i=0; i<quads.size(); i++){
//        fprintf(file, "8 ");
//        fprintf(file, "%d %d %d %d ", quads[i].getIndex(0), quads[i].getIndex(1), quads[i].getIndex(6), quads[i].getIndex(3));
//        fprintf(file, "%d %d %d %d ", quads[i].getIndex(5), quads[i].getIndex(4), quads[i].getIndex(7), quads[i].getIndex(2));
//        fprintf(file, "\n");
//    }
//    fprintf(file, "\n");

//    fprintf(file, "CELL_TYPES %d\n", (int)quads.size());

//    for(uint i=0; i<quads.size(); i++)
//        fprintf(file, "%d ", 12);
//    fprintf(file, "\n\n");

//    fclose(file);
//}
