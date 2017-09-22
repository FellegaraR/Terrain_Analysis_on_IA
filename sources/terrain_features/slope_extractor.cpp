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

#include "slope_extractor.h"

void Slope_Extractor::compute_triangles_slopes(Spatial_Mesh& mesh)
{
    t_slopes.assign(mesh.get_triangles_num(),0);

    for(itype t = 0; t < mesh.get_triangles_num(); t++)
    {
        Triangle& tr = mesh.get_triangle(t);
        t_slopes[t] = Geometry_Slope::compute_triangle_slope(tr,mesh);

        if(t_slopes[t] < this->min)
            this->min = t_slopes[t];
        if(t_slopes[t] > this->max)
            this->max = t_slopes[t];
        this->avg += t_slopes[t];
    }

    print_slopes_stats(mesh.get_triangles_num());
    reset_stats();
    t_slopes.clear();
}

void Slope_Extractor::compute_edges_slopes(Spatial_Mesh& mesh)
{
    for(itype t = 0; t < mesh.get_triangles_num(); t++)
    {
        Triangle& tr = mesh.get_triangle(t);
        for(int v=0; v<tr.vertices_num(); v++)
        {
            Edge e = tr.TE(v);
            if(e_slopes.find(e)==e_slopes.end())
            {
                coord_type s = Geometry_Slope::compute_edge_slope(e,mesh);
                e_slopes.insert(make_pair(e,s));
                if(s < this->min)
                    this->min = s;
                if(s > this->max)
                    this->max = s;
                this->avg += s;
            }
        }
    }

    print_slopes_stats(e_slopes.size());
    reset_stats();
    e_slopes.clear();
}
