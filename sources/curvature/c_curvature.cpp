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

#include <cstdio>

#include "c_curvature.h"

//compute curvature of polyline segment p1 v2 p3 and set its sign 
//according to the direction of surface normal v2norm at vertex v2
coord_type C_Curvature::wedge_curvature(Vertex &p1, Vertex &v2, Vertex &p3, Vertex &v2norm)
{
    coord_type cosinus, curva;

    cosinus = Geometry_Curvature::cos_angle(p1,v2,p3);
    if (cosinus<=-1.0) curva = 0.0;
    else curva = PI - acos(cosinus);

    if (curva==0.0)
        return curva;

    Vertex v2_plus_norm = v2 + v2norm;
    coord_type cos2, cos3, ang1, ang2, ang3;

    if (cosinus<=-1.0)
        ang1 = PI;
    else
        ang1 = acos(cosinus);
    cos2 = Geometry_Curvature::cos_angle(p1,v2,v2_plus_norm);
    if (cos2<=-1.0)
        ang2 = PI;
    else
        ang2 = acos(cos2);
    cos3 = Geometry_Curvature::cos_angle(p3,v2,v2_plus_norm);
    if (cos3<=-1.0)
        ang3 = PI;
    else
        ang3 = acos(cos3);
    // Note that ang1 is the angle <180 among the two angle formed by p1,v2,p3
    // so if surface normal lies inside this angle then surface is concave.
    // Therefore: change curvature sign iff ang1 == ang2+ang3,
    // but numerical errors may be present in the angles,
    // so do the following check:
    if ( ((ang2+ang3)<PI) && (ang2<ang1) && (ang3<ang1) )
        return -curva;
    else
        return curva;
}



//Given vertex v, its normal vnorm, and a point w, compute the curvature
//in v of the line obtained by intersecting the surface with the plane
//through v, vnorm, and w.
coord_type C_Curvature::compute_curve_and_curvature(itype v_id, Vertex &v, ivect &vt, Vertex &vnorm, Vertex &w, Spatial_Mesh &mesh)
{
//    cout<<"compute_curve_and_curvature"<<endl;
    dvect abc= { 0, 0, 0 };// coefficients of plane through v,w and parallel to normal

    Vertex sub = v - w;
    abc[0] =  +vnorm.get_c(2)*(sub.get_c(1)) -vnorm.get_c(1)*(sub.get_c(2));
    abc[1] =  -vnorm.get_c(2)*(sub.get_c(0)) +vnorm.get_c(0)*(sub.get_c(2));
    abc[2] =  -vnorm.get_c(0)*(sub.get_c(1)) +vnorm.get_c(1)*(sub.get_c(0));
    // il piano ha equazione:
    // a(x-vx)+b(y-vy)+c(z-vz)=0 dove v=(vx,vy,vz)

    Vertex w1 = Geometry_Curvature::find_plane_intersection(v_id, v, vt, vnorm, w, abc, mesh);

    //take curvature of wedge w v w1
    return wedge_curvature(w,v,w1, vnorm);
}

//function for computing mean and gaussian Ccurvature
coord_type C_Curvature::compute(itype v, Spatial_Mesh &mesh)
{
    Vertex &vert = mesh.get_vertex(v);

    itype v1,v2;
    Vertex vnorm(vert.get_dimension()), tnorm(vert.get_dimension());
    coord_type totang;
    coord_type cc;
    bool is_border = false; //forse non serve
    coord_type l;
    coord_type cc_sum = 0.0; // really used only if (take_mean_of_all)

    ivect vt = mesh.VT(v,is_border);
    ivect vv = mesh.VV(v);

    if (vt.size()==1) return 0.0; //one incident, we say that v is flat vertex

    //compute total angle at v and vertex normal as a weighted
    //sum of triangle normals, the weight of each triangle is
    //its angle in v.
    totang = 0.0;

    for(auto t : vt)
    {
        coord_type ang;
        Triangle &tr = mesh.get_triangle(t);
        int v_pos = tr.vertex_index(v);
        v1 = tr.TV((v_pos+1)%3);
        v2 = tr.TV((v_pos+2)%3);
        ang = acos(Geometry_Curvature::cos_angle(mesh.get_vertex(v1),vert,mesh.get_vertex(v2)));
        totang += ang;
        Geometry_Curvature::triangle_normal(tr,tnorm,mesh);
        vnorm += tnorm * ang;
    }
    vnorm /= (vt.size()*totang);
    l = vnorm.norm();
    vnorm /= l;

    coord_type max_cur=-INFINITY;
    coord_type min_cur=INFINITY;

    //outer loop on neighbor vertices w represented by incident triangles
    //(consider w = the first vertex following v in each triangle),
    //inner loop finds the triangle incident in v such that its edge opposite
    // to v intersects the plane through v and w and parallel to the surface
    // normal in v; compute the plane-edge intersection poits w1;
    //and compute curvature of curve w v w1, store min and max
    for(auto vid : vv)
    {
        cc = compute_curve_and_curvature(v,vert,vt,vnorm,mesh.get_vertex(vid),mesh);

        if (cc<min_cur)
            min_cur = cc;
        if (cc>max_cur)
            max_cur = cc;
        cc_sum += cc; // really used only if (take_mean_of_all)
    }

    //mean Ccurvature
    if (type==MEAN_CCURVATURE)
    {
        if (take_mean_of_all)
            return (cc_sum/vv.size());
        else
            return ((min_cur+max_cur)/2.0);
    }
    //gaussian Ccurvature
    else if (type==GAUSS_CCURVATURE)
        return (min_cur*max_cur);
    //this will never happen
    return 0.0;
}

//constructor: read mesh and compute curvature for each vertex
C_Curvature :: C_Curvature(CCurvatureType tp, /*, coord_type as*/ Spatial_Mesh &mesh, bool mn)
    : AbstractCurvature(mesh)
{    
    type = tp;
    take_mean_of_all = mn;
    if (type==MEAN_CCURVATURE)
    {
        printf("Mean Ccurvature with mean as ");
        if (mn==0) printf("(min+max)/2 \n");
        else printf("(sum of all N)/N \n");
    }
    else if (type==GAUSS_CCURVATURE)
    {
        printf("Gaussian Ccurvature \n");
    }
    else printf("Error: type of Ccurvature is not valid\n");
}
