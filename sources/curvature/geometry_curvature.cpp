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

#include "geometry_curvature.h"

//compute cosinus of angle formed by 3 vertices
coord_type Geometry_Curvature::cos_angle(Vertex &v1, Vertex &v2, Vertex &v3)
{
    coord_type norm21=v2.norm(v1);
    coord_type norm23=v2.norm(v3);
    coord_type product=v2.scalar_product(v1,v3);
    coord_type costeta=product/(norm21*norm23);
    if (costeta>1.0) costeta = 1.0;
    if (costeta<-1.0) costeta = -1.0;
    return costeta;
}

//compute total area of triangles incident in v
coord_type Geometry_Curvature::fan_area(ivect& vt, Spatial_Mesh &mesh)
{
    coord_type a = 0.0;

    for(utype i=0;i<vt.size();i++)
    {
        a += triangle_area(mesh.get_triangle(vt[i]),mesh);
    }

    return a;
}

//compute triangle area
coord_type Geometry_Curvature::triangle_area(Triangle& t, Spatial_Mesh &mesh)
{
    coord_type prodscaluv;
    coord_type normau;
    coord_type normav;
    coord_type cosalpha;
    coord_type senalpha;

    Vertex &v1 = mesh.get_vertex(t.TV(0));
    Vertex &v2 = mesh.get_vertex(t.TV(1));
    Vertex &v3 = mesh.get_vertex(t.TV(2));

    prodscaluv=v1.scalar_product(v2,v3);
    normau=v1.norm(v2);
    normav=v1.norm(v3);
    cosalpha = prodscaluv / (normau*normav);
    senalpha=sqrt(1-(cosalpha*cosalpha));
    if(std::isnan(senalpha)) senalpha=0.0001;

    return (normau*normav*senalpha)/2;
}

void Geometry_Curvature::triangle_normal(Triangle &t, Vertex &tnorm, Spatial_Mesh &mesh)
{
    Vertex &v1 = mesh.get_vertex(t.TV(0));
    Vertex &v2 = mesh.get_vertex(t.TV(1));
    Vertex &v3 = mesh.get_vertex(t.TV(2));

    coord_type a[3], b[3];
    coord_type norm=0;

    for(int i=0; i < v1.get_dimension(); i++)
    {
        a[i] = v1.get_c(i)-v2.get_c(i);
        b[i] = v1.get_c(i)-v3.get_c(i);
    }

    tnorm = Vertex(a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]);

    norm = tnorm.norm();
    tnorm /= norm;
}

//test if segment v1-v2 intersect plane of equation
// a(x-vx)+b(y-vy)+c(z-vz)=0 where v=(vx,vy,vz), three possible results:
// 0 : no intersection
// 1 : intersection is v1
// 2 : intersection is v2
// 3 : proper intersection
int Geometry_Curvature::intersect_plane(Vertex &v1, Vertex &v2, dvect &abc, Vertex &v)
{
    coord_type res1=0, res2=0;
    for(int i=0; i < v1.get_dimension(); i++)
    {
        res1 += abc[i]*(v1.get_c(i)-v.get_c(i));
        res2 += abc[i]*(v2.get_c(i)-v.get_c(i));
    }

    if ((-TOLER<=res1)&&(res1<=TOLER)) return 1;
    if ((-TOLER<=res2)&&(res2<=TOLER)) return 2;
    if ( (res1>0.0)&&(res2<0.0) ) return 3;
    if ( (res2>0.0)&&(res1<0.0) ) return 3;
    return 0;
}

int Geometry_Curvature::same_point(Vertex &p1, Vertex &p2)
{
    for(int i=0; i < p1.get_dimension(); i++)
        if (fabs(p1.get_c(i)-p2.get_c(i))>SMALL_TOLER)
            return 0;
    return 1;
}

//Find the triangle, around V v, that is intersected by the plane through
//v, the normal vnorm to v, and point w, the equation of such plane is
//a(x-vx)+b(y-vy)+c(z-vz)=0 dove v=(vx,vy,vz):
//compute the intersection point and return it.
Vertex Geometry_Curvature::find_plane_intersection(itype v_id, Vertex &v, ivect &vt, Vertex &vnorm, Vertex &w, dvect &abc, Spatial_Mesh &mesh)
{
    //    int t1;
    bool found = false;
    itype v1_id, v2_id;
    Vertex w1;

    for(ivect_iter it=vt.begin(); it!=vt.end(); ++it)
    {
        Triangle &t = mesh.get_triangle(*it);
        int v_pos = t.vertex_index(v_id);

        v1_id = t.TV((v_pos+1)%t.vertices_num());
        v2_id = t.TV((v_pos+2)%t.vertices_num());
        Vertex &v1 = mesh.get_vertex(v1_id);
        Vertex &v2 = mesh.get_vertex(v2_id);

        if ( !Geometry_Curvature::same_point(v1,w) && !Geometry_Curvature::same_point(v2,w) )
        {
            switch(Geometry_Curvature::intersect_plane(v1,v2,abc,v))
            {
            case 1: // intersection is v1
                w1 = v1;
                found = true;
                break;
            case 2: // intersection is v2
                w1 = v2;
                found = true;
                break;
            case 3: // proper intersection, compute intersection
                coord_type d=0,n=0;
                for(int i=0; i < v.get_dimension(); i++)
                {
                    d += -abc[i]*v1.get_c(i) + abc[i]*v.get_c(i);
                    n +=  abc[i]*v2.get_c(i) - abc[i]*v1.get_c(i);
                }
                coord_type s = d / n;
                w1 = Vertex(v.get_dimension());
                for(int i=0; i < v1.get_dimension(); i++)
                    w1.set_c(i,s*(v2.get_c(i)-v1.get_c(i)) + v1.get_c(i));
                found = true;
                break;
            }

            // se ho trovato w stesso devo cercare l'altro punto di intersezione
            if (found && same_point(w1,w))
                found = false;
        }

        if(found)
            break;
    }

    if(!found)
        // v is a boundary V, we have visited all triangles
        // and not found an edge which intersects...
    {
        // simulate a vertical wall
        w1 = Vertex(v.get_dimension());
        for(int i=0; i < v.get_dimension(); i++)
            w1.set_c(i,v.get_c(i)-vnorm.get_c(i));
    }

    return w1;
}

//compute barycentric area of triangle t centered at
//vertex v, where v is a vertex of t.
//compute Voronoi-barycentric area around vertex v
coord_type Geometry_Curvature::voronoi_barycentric_area(itype v, Spatial_Mesh &mesh)
{
    coord_type a = 0.0;
    itype pred = -1;
    itype current = mesh.get_vertex(v).get_VTstar();
    a += voronoi_barycentric_area(v, current, mesh);

    itype k = mesh.get_triangle(current).vertex_index(v);

    //choose one side to visit
    pred = current;
    current = mesh.get_triangle(current).TT((k+1)%3);

    bool is_border = false;
    while(current != mesh.get_vertex(v).get_VTstar())
    {
        if(current == -1) // border
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

        a += voronoi_barycentric_area(v, current, mesh);
        mesh.get_triangle(current).next_triangle_around_v(v,current,pred);
    }
    return a;
}

//compute the contribution of triangle t to the mixed area around
//vertex v, where v is a vertex of t.
//the contribution of t is the Voronoi area if
//all angles are acute, or the barycentric area otherwise
coord_type Geometry_Curvature::voronoi_barycentric_area(itype v, itype t, Spatial_Mesh &mesh)
{
    itype i, v1, v2;
    coord_type cos_a, cos_b, cos_c;

    Triangle &tr = mesh.get_triangle(t);

    // triangle t is v v1 v2
    // a,b,c are the angles in v, v1, v2
    i = tr.vertex_index(v);
    v1 = tr.TV((i+1)%3);
    v2 = tr.TV((i+2)%3);

    Vertex &vert = mesh.get_vertex(v);
    Vertex &vert1 = mesh.get_vertex(v1);
    Vertex &vert2 = mesh.get_vertex(v2);

    cos_a = Geometry_Curvature::cos_angle(vert2,vert,vert1);
    cos_b = Geometry_Curvature::cos_angle(vert,vert1,vert2);
    cos_c = Geometry_Curvature::cos_angle(vert1,vert2,vert);

    if ((cos_a<0.0)||(cos_b<0.0)||(cos_c<0.0)) // one is obtuse
    {
        if (cos_a<0.0) // the obtuse angle is in v
            return 0.5 * Geometry_Curvature::triangle_area(tr,mesh);
        else
            return 0.25 * Geometry_Curvature::triangle_area(tr,mesh);
    }
    else // Voronoi area
    {
        coord_type sin_b, cot_b, sin_c, cot_c, e1, e2;
        /* if angles are zero area is zero */
        if (cos_b==1.0) return 0.0;
        if (cos_c==1.0) return 0.0;

        sin_b = sin(acos(cos_b));
        sin_c = sin(acos(cos_c));

        cot_b = cos_b / sin_b;
        cot_c = cos_c / sin_c;

        e1 = vert1.norm(vert);
        e2 = vert2.norm(vert);
        return (cot_c*e1*e1 + cot_b*e2*e2) / 8.0;
    }
}
