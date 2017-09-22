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

#ifndef MEANCURVATURE_H
#define MEANCURVATURE_H

#include "abstract_curvature.h"

class MeanCurvature : public AbstractCurvature
{
public:
    //constructor:
    MeanCurvature(bool div_by_area, Spatial_Mesh &mesh) : AbstractCurvature(mesh,div_by_area) {}

protected:
    //curvature estimation function: mean angle deficit
    coord_type compute(itype v, Spatial_Mesh &mesh);
    //compute contribution of a pair of consecutive triangles incident
    //in vertex v, and put the contribution in vector somma
    void triangle_pair_contrib(itype v, itype t, itype tnext, Vertex &summation, Spatial_Mesh &mesh);
};

#endif // MEANCURVATURE_H
