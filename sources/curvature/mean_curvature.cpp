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

#include "mean_curvature.h"
#include "geometry_curvature.h"

//function for computing mean angle deficit curvature at vertex v
coord_type MeanCurvature::compute(itype v, Spatial_Mesh &mesh)
{
    Vertex summation = Vertex(3);
    coord_type curva;

    itype pred = -1;
    itype current = mesh.get_vertex(v).get_VTstar();

    itype k = mesh.get_triangle(current).vertex_index(v);

    //choose one side to visit
    pred = current;
    current = mesh.get_triangle(current).TT((k+1)%3);

    bool is_border = false;
    while(current != mesh.get_vertex(v).get_VTstar())
    {
        if(current != -1)
        {
            triangle_pair_contrib(v,pred,current,summation,mesh);
            mesh.get_triangle(current).next_triangle_around_v(v,current,pred);
        }
        else if(current == -1) // border
        {
            if(is_border) // if it is the second time that I get to a border I exit
                break;
            else
            {
                // otherwise I hit a triangle on the border for the first time
                // then, I visit in the opposite direction starting again from VTstar
                is_border = true;
                pred = mesh.get_vertex(v).get_VTstar();
                k = mesh.get_triangle(pred).vertex_index(v);
                current = mesh.get_triangle(pred).TT((k+2)%3);
                if(current==-1) // no need visit the star in the opposite direction
                    break;
            }
        }
    }

    if (this->divide_by_area)
    {
        coord_type Atot = Geometry_Curvature::voronoi_barycentric_area(v,mesh);

        if(Atot>0.00001)//approximation error
        {
            summation /= (2*Atot);
        }
        else
        {
            summation = Vertex(3);
        }
    }
    else
    {
        summation /= 2.0;
    }

    curva = 0.5*(summation.norm());

    return(curva);
}


//compute contribution of a pair of consecutive triangles incident
//in vertex v, and put the contribution in vector somma
void MeanCurvature::triangle_pair_contrib(itype v, itype t, itype tnext, Vertex &summation, Spatial_Mesh &mesh)
{
    itype i, inext;
    itype w, v1, v2;
    coord_type cos_a, sin_a, cot_a;
    coord_type cos_b, sin_b, cot_b;

    Vertex &vert = mesh.get_vertex(v);

    Triangle &tr = mesh.get_triangle(t);
    Triangle &tr_next = mesh.get_triangle(tnext);

    i = tr.vertex_index(v);
    inext = tr_next.vertex_index(v);
    w = tr.TV((i+1)%3);
    v1 = tr.TV((i+2)%3);
    v2 = tr_next.TV((inext+1)%3);

    if(v2==v1) //se sono uguali inverto l'ordine dei due triangoli
    {
        w = tr_next.TV((inext+1)%3);
        v1 = tr_next.TV((inext+2)%3);
        v2 = tr.TV((i+1)%3);
    }

    Vertex &vert1 = mesh.get_vertex(v1);
    Vertex &vert2 = mesh.get_vertex(v2);
    Vertex &vertW = mesh.get_vertex(w);

    cos_a = Geometry_Curvature::cos_angle(vert,vert1,vertW);
    cos_b = Geometry_Curvature::cos_angle(vert,vert2,vertW);

    /* if angle is zero the triangles are degenerate */
    if ( (cos_a==1.0) || (cos_a==-1.0) )
    {
        fprintf(stderr,"Degenerate triangle %d %d %d\n", v,v1,w);
        return;
    }
    if ( (cos_b==1.0) || (cos_b==-1.0) )
    {
        fprintf(stderr,"Degenerate triangle %d %d %d\n", v,v2,w);
        return;
    }

    sin_a = sin(acos(cos_a));
    cot_a = cos_a / sin_a;

    sin_b = sin(acos(cos_b));
    cot_b = cos_b / sin_b;

    summation += (vert - vertW) * (cot_a+cot_b);
}
