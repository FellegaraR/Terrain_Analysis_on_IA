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

#include "geometry_slope.h"

coord_type Geometry_Slope::compute_triangle_slope(Triangle &t, Spatial_Mesh &mesh)
{
    Vertex &v1 = mesh.get_vertex(t.TV(0));
    Vertex &v2 = mesh.get_vertex(t.TV(1));
    Vertex &v3 = mesh.get_vertex(t.TV(2));

    // get the two vectors in the triangle
    Vertex u = v2 - v1;
    Vertex v = v3 - v1;

    // get the cross product to get the normal
    dvect n = { u.get_c(1)*v.get_c(2) - u.get_c(2)*v.get_c(1) ,
                u.get_c(2)*v.get_c(0) - u.get_c(0)*v.get_c(2),
                u.get_c(0)*v.get_c(1) - u.get_c(1)*v.get_c(0) };
    // get the magnitude of the normal
    coord_type magnitude_n = sqrt( n[0]*n[0] + n[1]*n[1] + n[2]*n[2] );

    // the dot product between n and the ground normal is only formed by the normalized y component of n
    n[1] /= magnitude_n;
    // thus, we do not have to compute the real dot product but we use the y part of the norm
    coord_type angle = acos(n[1]);

    return angle;
}

coord_type Geometry_Slope::compute_edge_slope(Edge &e, Spatial_Mesh &mesh)
{
    Vertex &v1 = mesh.get_vertex(e.EV(0));
    Vertex &v2 = mesh.get_vertex(e.EV(1));
    Vertex v3 = Vertex(v1.get_c(0),v1.get_c(1),v2.get_c(2)); // temp vertex to get a ground

    // get the two vectors in the triangle
    Vertex u = v2 - v1;
    Vertex v = v3 - v1;

    coord_type u_mag = u.norm();
    u /= u_mag;

    coord_type v_mag = v.norm();
    v /= v_mag;

    coord_type prod = u.scalar_product(v);

    coord_type angle = acos(prod);

    return angle;
}
