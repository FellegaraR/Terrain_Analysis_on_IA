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

#ifndef ABSTRACT_CURVATURE_H
#define ABSTRACT_CURVATURE_H

#include "ia/mesh.h"
#include "ia/vertex.h"
#include "ia/triangle.h"

class AbstractCurvature
{
public:
    //compute for all vertices
    void compute_values(Spatial_Mesh &mesh);
    ///
    inline coord_type get_curvature(itype pos) { return curvature[pos]; }
    inline void set_curvature(int pos, coord_type curv) { this->curvature[pos] = curv; }

    void print_curvature_stats(Spatial_Mesh &mesh);

protected:
    dvect curvature;
    //flag if we want to divide by area of the neighborhood or not
    bool divide_by_area;

    //constructor
    AbstractCurvature(Spatial_Mesh &mesh, bool div_area=false)
    {
        divide_by_area = div_area;
        curvature.assign(mesh.get_vertices_num(),0);
    }

    virtual coord_type compute(itype v, Spatial_Mesh &mesh)=0;

};

#endif // ABSTRACT_CURVATURE_H
