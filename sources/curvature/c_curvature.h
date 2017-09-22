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

#ifndef C_CURVATURE_H
#define C_CURVATURE_H
#include <cmath>

#include "abstract_curvature.h"
#include "geometry_curvature.h"

enum CCurvatureType {GAUSS_CCURVATURE = 0, MEAN_CCURVATURE = 1};

class C_Curvature : public AbstractCurvature
{
protected:
    //type of Ccurvature: MEAN_CCURVATURE or GAUSS_CCURVATURE
    CCurvatureType type;
    //in case of mean Ccurvature, if 1 compute mean as sum of all values
    //divided by number of values, instead of min+max divided by two
    bool take_mean_of_all;
    //curvature estimation function: depending on type
    coord_type compute(itype v, Spatial_Mesh &mesh);
    //compute curvature of polyline segment p1 v2 p3
    //and set its sign according to the direction of
    //surface normal v2norm at vertex v2
    coord_type wedge_curvature(Vertex &p1, Vertex &v2, Vertex &p3, Vertex &v2norm);
    //Given vertex v, its normal vnorm, and a point w, compute the curvature
    //in v of the line obtained by intersecting the surface with the plane
    //through v, vnorm, and w.
    coord_type compute_curve_and_curvature(itype v_id, Vertex &v, ivect &vt, Vertex &vnorm, Vertex &w, Spatial_Mesh &mesh);

public:
    //constructor given the type (see above for values), the flag whether
    //we compute mean by summing all values
    C_Curvature(CCurvatureType tp, Spatial_Mesh &mesh, bool mn=false);
};

#endif //C_CURVATURE_H
