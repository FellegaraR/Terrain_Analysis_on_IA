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

#ifndef QUAD_MESH_H
#define QUAD_MESH_H

#include <vector>
#include <map>
#include "ia/mesh.h"
#include "ia/vertex.h"
#include "ia/triangle.h"

using namespace std;

class QuadVertex
{
    dvect coords;

public:
    inline QuadVertex(){ coords= dvect(3,0);}
    inline QuadVertex(coord_type x, coord_type y, coord_type z){
        coords.push_back(x);
        coords.push_back(y);
        coords.push_back(z);
    }
    inline dvect* getCoords(){return &coords;}
    inline void setCoords(dvect co){coords = co;}
    inline coord_type get_i_coord(itype i){return coords[i];}

    inline bool operator==(QuadVertex v){
        for(itype i=0; i<3; i++){
            if(this->get_i_coord(i) == v.get_i_coord(i)) continue;
            else return false;
        }
        return true;
    }

    inline bool operator<(QuadVertex v){
        for(itype i=0; i<3; i++){
            if(this->get_i_coord(i) < v.get_i_coord(i)) return true;
            else if(this->get_i_coord(i) == v.get_i_coord(i)) continue;
            else return true;
        }
        return false;
    }

    inline bool operator>(QuadVertex v){
        for(itype i=0; i<3; i++){
            if(this->get_i_coord(i) > v.get_i_coord(i)) return true;
            else if(this->get_i_coord(i) == v.get_i_coord(i)) continue;
            else return true;
        }
        return false;
    }
};

class QuadTri
{
    ivect vertices;
    itype original_vertex;
    itype original_tetra;

public:
    inline QuadTri(ivect vertices, itype vert, itype tetra){
        this->vertices = vertices;
        original_vertex = vert;
        original_tetra = tetra;
    }

    inline itype getIndex(itype i){return vertices[i];}
    inline itype getOriginalVertex(){return original_vertex;}
    inline itype get_num_vertices() { return vertices.size(); }

    inline friend std::ostream& operator<<(std::ostream& out, const QuadTri& p)
    {
        out <<"Q[ ";
        for(utype i=0; i<p.vertices.size(); i++)
            out << p.vertices[i] << " ";
        out << "]";
        return out;
    }
};


class Quad_Mesh
{
private:
    vector<QuadVertex> vertices;
    vector<QuadTri> quads;

private:
    dvect compute_barycenter(itype v1, itype v2, Spatial_Mesh &mesh);
    dvect compute_barycenter(itype t, Spatial_Mesh &mesh);
public:
    inline Quad_Mesh() { }

    void loadQuad(Spatial_Mesh &mesh);
    void save_quad_mesh(string file_name);
//    void write_vtk(char* filename);
};

#endif // QUAD_MESH_H
