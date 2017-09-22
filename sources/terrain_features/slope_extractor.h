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

#ifndef SLOPE_EXTRACTOR_H
#define SLOPE_EXTRACTOR_H

#include <climits>
#include "terrain_features/geometry_slope.h"
#include "ia/mesh.h"

class Slope_Extractor
{
public:

    Slope_Extractor()
    {
        min = INT_MAX;
        max = INT_MIN;
        avg = 0;
    }

    //compute for all triangles the slope
    //storing the values in a global array
    void compute_triangles_slopes(Spatial_Mesh &mesh);

    void compute_edges_slopes(Spatial_Mesh &mesh);

    inline void print_slopes_stats(utype num) { cerr<<"   min: "<<min<<" avg: "<<avg/(coord_type)num<<" max: "<<max<<endl; }

    inline void reset_stats()
    {
        min = INT_MAX;
        max = INT_MIN;
        avg = 0;
    }


private:
    map<Edge,coord_type> e_slopes;
    dvect t_slopes;
    coord_type min, avg, max;
};



#endif // SLOPE_EXTRACTOR_H
