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

#ifndef GEOMETRY_SLOPE_H
#define GEOMETRY_SLOPE_H

#include "utilities/basic_wrappers.h"
#include "ia/mesh.h"

class Geometry_Slope
{
public:
    static coord_type compute_triangle_slope(Triangle &t, Spatial_Mesh &mesh);

    static coord_type compute_edge_slope(Edge &e, Spatial_Mesh &mesh);

private:
    Geometry_Slope() {}
};

#endif // GEOMETRY_SLOPE_H
