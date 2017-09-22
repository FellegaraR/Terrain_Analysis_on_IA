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

#include "abstract_curvature.h"

//compute for all vertices
void AbstractCurvature::compute_values(Spatial_Mesh &mesh)
{
    for(itype i=0; i<mesh.get_vertices_num(); i++)
    {
        curvature[i]=compute(i,mesh); //mean angle deficit
    }
}

void AbstractCurvature::print_curvature_stats(Spatial_Mesh &mesh)
{
    coord_type min=INFINITY, max=-INFINITY, avg=0;
    for(itype v=0; v<mesh.get_vertices_num(); v++)
    {
        if(curvature[v] < min)
            min = curvature[v];
        if(curvature[v] > max)
            max = curvature[v];
        avg += curvature[v];
    }
    cerr<<"[STATS] curvature min: "<<min<<" avg: "<<avg/(coord_type)mesh.get_vertices_num()<<" max: "<<max<<endl;
}
