/*
    This file is part of the LibTri library.

    Author(s): Riccardo Fellegara (riccardo.fellegara@gmail.com)
               Paola Magillo (paola.magillo@unige.it)

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

#include "concentrated_curvature.h"
#include "geometry_curvature.h"

//function for computing mean angle deficit curvature at vertex v
coord_type ConcentratedCurvature::compute(itype v, Spatial_Mesh &mesh)
{
    itype v1,v2, pos=-1;
    coord_type totang = 0.0;
    coord_type curva;

    bool is_border = false;
    ivect vt = mesh.VT(v,is_border);

    //t==-1 cannot happen if the mesh has no isolated vertices
    if (vt.size() == 0) return 0.0;

    for(auto t : vt)
    {
        Triangle &current = mesh.get_triangle(t);
        pos = current.vertex_index(v);
        v1 = current.TV((pos+1)%current.vertices_num());
        v2 = current.TV((pos+2)%current.vertices_num());
        coord_type res = Geometry_Curvature::cos_angle(mesh.get_vertex(v1),mesh.get_vertex(v),mesh.get_vertex(v2));
        totang += acos(res);
    }

    if(is_border)
        curva = PI-totang; // boundary vertex
    else
        curva = 2.0*PI-totang;   // interior vertex

    if (this->divide_by_area) // gaussian angle deficit
        curva /= Geometry_Curvature::voronoi_barycentric_area(v,mesh);
    // else concentrated curvature

    return curva;
}
