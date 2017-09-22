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

#ifndef MESH_H
#define MESH_H

#include <vector>
#include <string>
#include <cstdlib>
#include <iostream>
#include <queue>
#include <cmath>

#include "utilities/sorting.h"
#include "ia/edge.h"
#include "ia/triangle.h"
#include "ia/vertex.h"

/*
Return true iff v is an extreme of the edge v1-v2
*/
#define is_endpoint(v,v1,v2) ( ( (v)==(v1) ) || ( (v)==(v2) ) )

using namespace std;
///A class representing a generic mesh parametrized by the type of top simplexes
template<class V> class Mesh
{
public:
    ///A constructor method
    Mesh();
    ///A constructor method
    Mesh(const Mesh& orig);
    ///A destructor method
    virtual ~Mesh();
    ///A public method that returns the vertex at the i-th position in the mesh list
    /*!
     * \param id an integer argument, representing the position in the list
     * \return a Vertex&, the vertex at the id-th position in the list
     */
    inline V& get_vertex(itype id) { return this->vertices[id]; }
    ///A public method that returns the number of mesh vertices
    /*!
     * \return an integer, representing the number of vertices
     */
    inline itype get_vertices_num() { return this->vertices.size(); }
    ///A public method that adds a vertex to the vertices list
    /*!
     * \param v a Vertex& argument, representing the vertex to add
     */
    inline void add_vertex(V& v) { this->vertices.push_back(v); }
    ///A public method that returns the tetrahedron at the i-th position in the mesh list
    /*!
     * \param id an integer argument, representing the position in the list
     * \return a T&, the tetrahedron at the id-th position in the list
     */
    inline Triangle& get_triangle(itype id) { return this->triangles[id]; }
    ///A public method that returns the number of triangles
    /*!
     * \return an integer, representing the number of triangles
     */
    inline itype get_triangles_num() { return this->triangles.size(); }
    ///A public method that adds a top simplex to the triangles list
    /*!
     * \param t a Triangle& argument, representing the triangle to add
     */
    inline void add_triangle(Triangle& t) { this->triangles.push_back(t); }

    bool build();

    ivect VT(itype center, bool &is_border);
    vector<Edge> VE(itype center);
    ivect VV(itype center);

    ivect ET(Edge &e);
    vector<Edge> EE(Edge &e);

    // check if the center is on the border
    bool is_boundary(int center);

    ///A public method that initializes the space needed by the vertices and triangles arrays
    /*!
     * \param numV an itype, represents the number of mesh vertices
     * \param numT an itype, represents the number of mesh triangles
     */
    inline void reserve(itype numV, itype numT)
    {
        this->vertices.reserve(numV);
        this->triangles.reserve(numT);
    }
    ///A public method that initializes the space needed by the vertices array
    /*!
     * \param numV an itype, represents the number of mesh vertices
     */
    inline void reserve_vertices_space(itype numV) { this->vertices.reserve(numV); }
    ///A public method that initializes the space needed by the triangles array
    /*!
     * \param numT an itype, represents the number of mesh triangles
     */
    inline void reserve_triangles_space(itype numT) { this->triangles.reserve(numT); }

protected:
    ///A private varible representing the vertices list of the mesh
    vector<V> vertices;
    ///A private varible representing the top simplexes list of the mesh
    vector<Triangle> triangles;

    void link_adj (itype t1, itype t2);
};

template<class V> Mesh<V>::Mesh()
{
    vertices = vector<V>();
    triangles = vector<Triangle>();
}

template<class V> Mesh<V>::Mesh(const Mesh& orig)
{
    this->vertices = orig.vertices;
    this->triangles = orig.triangles;
}

template<class V> Mesh<V>::~Mesh()
{
    vertices.clear();
    triangles.clear();
}

#include <climits>

template<class V> bool Mesh<V>::build()
{
    aux *tr_vec ;
    utype k;
    itype t1, t2;

    utype num = (this->get_triangles_num()*3);
    cout<<"edges: "<<num<<endl;

    tr_vec = (aux*) calloc( num , sizeof(aux) ) ;
    if (!tr_vec)
    {
        cerr << "failed malloc in find_adj" <<endl;
        return false;
    }

    k = 0;

    Edge e;

    for (itype j=0; j<this->get_triangles_num(); j++)
    {
        Triangle &t = this->get_triangle(j);
        for (itype i=0;i<t.vertices_num();i++)
        {
            // initialize (if unset) the partial VT of a vertex
            if(this->get_vertex(t.TV(i)).get_VTstar()==-1)
                this->get_vertex(t.TV(i)).set_VTstar(j);

            e = t.TE(i);

            tr_vec[k].t = j;
            tr_vec[k].v1 = e.EV(0);
            tr_vec[k].v2 = e.EV(1);
            k++;
        }
    }

    qsort(tr_vec,num,sizeof(aux),cmp_aux) ;

    for(utype k=0;k<num-1;k++)
    {
        if ( is_endpoint(tr_vec[k].v1,tr_vec[k+1].v1,tr_vec[k+1].v2) &&
             is_endpoint(tr_vec[k].v2,tr_vec[k+1].v1,tr_vec[k+1].v2) )
        {  /* i due triangoli hanno lo stesso lato */
            t1 = tr_vec[k].t;
            t2 = tr_vec[k+1].t;
            link_adj(t1,t2);
        }
    }
    free(tr_vec) ;

    return true;
}

template<class V> void Mesh<V>::link_adj(itype t1, itype t2)
        /* Lega t1 come adiacente di t2 e viceversa */
{
    Triangle &tr1 = this->get_triangle(t1);
    Triangle &tr2 = this->get_triangle(t2);

    itype i, j, k, pos1[2], pos2[2];
    {
        k = 0;
        for (i=0; ((i<3)&&(k<2)); i++)
        {
            for (j=0; ((j<3)&&(k<2)); j++)
            {
                if(tr1.TV(i) == tr2.TV(j))
                {
                    pos1[k] = i;
                    pos2[k] = j;
                    k++;
                }
            }
        }
        if (k<2)
        {
            cerr << "error in link_adj" <<endl;
        }
        else
        {
            tr1.set_TT(3-pos1[0]-pos1[1],t2);
            tr2.set_TT(3-pos2[0]-pos2[1],t1);
        }
    }
}

template<class V> ivect Mesh<V>::VT(itype center, bool &is_border)
{
    ivect triangles;

    itype pred = -1;
    itype current = this->get_vertex(center).get_VTstar();
    triangles.push_back(current);

    itype k = this->get_triangle(current).vertex_index(center);

    //choose one side to visit
    pred = current;
    current = this->get_triangle(current).TT((k+1)%3);

    while(current != this->get_vertex(center).get_VTstar())
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
                pred = this->get_vertex(center).get_VTstar();
                k = this->get_triangle(pred).vertex_index(center);
                current = this->get_triangle(pred).TT((k+2)%3);
                if(current==-1) // no need visit the star in the opposite direction
                    break;
            }
        }

        triangles.push_back(current);
        this->get_triangle(current).next_triangle_around_v(center,current,pred);
    }
    return triangles;
}

template<class V> vector<Edge> Mesh<V>::VE(itype center)
{
    vector<Edge> edges;

    itype pred = -1;
    itype current = this->get_vertex(center).get_VTstar();

    itype k = this->get_triangle(current).vertex_index(center);

    edges.push_back(this->get_triangle(current).TE((k+1)%3));

    //scelgo un giro a caso da prendere
    pred = current;
    current = this->get_triangle(current).TT((k+1)%3);

    bool is_border = false;

    while(current != this->get_vertex(center).get_VTstar())
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
                pred = this->get_vertex(center).get_VTstar();
                Triangle &tri = this->get_triangle(pred);
                k = tri.vertex_index(center);
                //
                edges.push_back(tri.TE((k+2)%3));
                current = tri.TT((k+2)%3);
                if(current==-1) // no need visit the star in the opposite direction
                    break;
            }
        }

        Triangle &tri = this->get_triangle(current);
        k = tri.vertex_index(center);

        if(tri.TT((k+1)%3) == pred)
        {
            edges.push_back(tri.TE((k+2)%3));
            pred = current;
            current = tri.TT((k+2)%3);
        }
        else if(tri.TT((k+2)%3) == pred)
        {
            edges.push_back(tri.TE((k+1)%3));
            pred = current;
            current = tri.TT((k+1)%3);
        }
    }

    return edges;
}

template<class V> ivect Mesh<V>::VV(int center)
{
    ivect vertices;
    itype pred = -1;
    itype current = this->get_vertex(center).get_VTstar();

    itype k = this->get_triangle(current).vertex_index(center);
    vertices.push_back(this->get_triangle(current).TV((k+2)%3));

    //scelgo un giro a caso da prendere
    pred = current;
    current = this->get_triangle(current).TT((k+1)%3);

    bool is_border = false;

    while(current != this->get_vertex(center).get_VTstar())
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
                pred = this->get_vertex(center).get_VTstar();
                Triangle &tri = this->get_triangle(pred);
                k = tri.vertex_index(center);
                //
                vertices.push_back(tri.TV((k+1)%3));
                current = tri.TT((k+2)%3);
                if(current==-1) // no need visit the star in the opposite direction
                    break;
            }
        }

        Triangle &tri = this->get_triangle(current);
        k = tri.vertex_index(center);

        if(tri.TT((k+1)%3) == pred)
        {
            vertices.push_back(tri.TV((k+1)%3));
            pred = current;
            current = tri.TT((k+2)%3);
        }
        else if(tri.TT((k+2)%3) == pred)
        {
            vertices.push_back(tri.TV((k+2)%3));
            pred = current;
            current = tri.TT((k+1)%3);
        }
    }

    return vertices;
}

template<class V> ivect Mesh<V>::ET(Edge &e)
{
    ivect triangles;
    bool is_border=false;
    ivect vt = this->VT(e.EV(0),is_border);
    for(auto t_id : vt)
    {
        Triangle &t = this->get_triangle(t_id);
        for(int j=0;j<t.vertices_num();j++)
        {
            if(t.has_vertex(e.EV(1)))
            {
                triangles.push_back(t_id);
                break;
            }
        }
    }

    return triangles;
}

template<class V> vector<Edge> Mesh<V>::EE(Edge &e)
{
    vector<Edge> ve0 = this->VE(e.EV(0));
    vector<Edge> ve1 = this->VE(e.EV(1));

    ve0.insert(std::end(ve0), std::begin(ve1), std::end(ve1));
    ve0.erase(std::remove(ve0.begin(), ve0.end(), e), ve0.end());
    return ve0;
}

template<class V> bool Mesh<V>::is_boundary(int center)
{
    itype pred = -1;
    itype current = this->get_vertex(center).VTstar();
    itype k = this->get_triangle(current).vertex_index(center);
    pred = current;
    current = this->get_triangle(current).TT((k+1)%3);

    while(current != this->get_vertex(center).VTstar())
    {
        if(current == -1)
            return true;
        this->get_triangle(current).next_triangle_around_v(center,current,pred);
    }

    return false;
}

/** @typedef Abs_Mesh
* @brief A Mesh<V> specialization, specific for abstract triangle meshes
*
*/
typedef Mesh<VT_star> Abs_Mesh;
/** @typedef Spatial_Mesh
* @brief A Mesh<V> specialization, specific for spatially embedded triangle meshes
*
*/
typedef Mesh<Vertex> Spatial_Mesh;

#endif // MESH_H
