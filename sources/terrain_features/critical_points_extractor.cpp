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

#include "critical_points_extractor.h"

void Critical_Points_Extractor::compute_critical_points(Spatial_Mesh &mesh)
{
    this->critical_points.assign(mesh.get_vertices_num(),Point_Type::REGULAR); //the default is a regular vertex
    this->extract_critical_points(mesh);
}

void Critical_Points_Extractor::extract_critical_points(Spatial_Mesh &mesh, flat_areas &fa)
{
    for(itype i=0; i<mesh.get_vertices_num(); i++)
    {
        Vertex &v = mesh.get_vertex(i);
        ivect vv = mesh.VV(i);
        bool isb=false;
        ivect vt = mesh.VT(i,isb);

        bool into_flat_area = false;
        map<itype,itype> upper, lower;

        for(auto vid : vv)
        {
            Vertex &v2 = mesh.get_vertex(vid);

            if(v.get_c(2) == v2.get_c(2))
            {
                flat_areas::iterator ret = fa.find(v.get_c(2));
                if(ret == fa.end())
                {
                    iset triangles;
                    triangles.insert(vt.begin(),vt.end());
                    fa.insert(make_pair(v.get_c(2),triangles));
                }
                else
                {
                    ret->second.insert(vt.begin(),vt.end());
                }
                into_flat_area = true;
                break;
            }
            else if(v.get_c(2) > v2.get_c(2))
            {
                lower.insert(make_pair(vid,vid)); //initially the component is formed only by the vertex itself
            }
            else
            {
                upper.insert(make_pair(vid,vid)); //initially the component is formed only by the vertex itself
            }
        }

        if(!into_flat_area)
        {
            if(upper.empty()) //the vertex is a maximum
                this->critical_points[i] = Point_Type::MAXIMUM;
            else if(lower.empty()) //the vertex is a minimum
                this->critical_points[i] = Point_Type::MINIMUM;
            else //we have to check the connected components of the edges in the link
            {
                map<itype,iset> adj_upper;
                map<itype,iset> adj_lower;

                // get the adjacent vertices that have greater or lower field values w.r.t. v
                init_adjacent_vertices(i,v,vt,mesh,adj_upper,adj_lower);

                // merge upper components
                utype uc_num = get_components_num(adj_upper,upper);
                // merge lower components
                utype lc_num = get_components_num(adj_lower,lower);

                if(uc_num == 2 && lc_num == 2) // simple saddle
                    this->critical_points[i] = Point_Type::SADDLE;
                else if(uc_num > 2 && lc_num > 2) // multiple saddle
                    this->critical_points[i] = Point_Type::MULTIPLE_SADDLE;
                // else the point is already flagged as regular point
            }
        }
    }
}

void Critical_Points_Extractor::extract_critical_points(Spatial_Mesh &mesh)
{

    for(itype i=0; i<mesh.get_vertices_num(); i++)
    {
        Vertex &v = mesh.get_vertex(i);
        ivect vv = mesh.VV(i);
        bool isb=false;
        ivect vt = mesh.VT(i,isb);

        bool into_flat_area = false;
        map<itype,itype> upper, lower;

        for(auto vid : vv)
        {
            Vertex &v2 = mesh.get_vertex(vid);

            if(v.get_c(2) == v2.get_c(2))
            {
                into_flat_area = true;
                break;
            }
            else if(v.get_c(2) > v2.get_c(2))
            {
                lower.insert(make_pair(vid,vid)); //initially the component is formed only by the vertex itself
            }
            else
            {
                upper.insert(make_pair(vid,vid)); //initially the component is formed only by the vertex itself
            }
        }

        if(!into_flat_area)
        {
            if(upper.empty()) //the vertex is a maximum
                this->critical_points[i] = Point_Type::MAXIMUM;
            else if(lower.empty()) //the vertex is a minimum
                this->critical_points[i] = Point_Type::MINIMUM;
            else //we have to check the connected components of the edges in the link
            {
                map<itype,iset> adj_upper;
                map<itype,iset> adj_lower;

                // get the adjacent vertices that have greater or lower field values w.r.t. v
                init_adjacent_vertices(i,v,vt,mesh,adj_upper,adj_lower);

                // merge upper components
                utype uc_num = get_components_num(adj_upper,upper);
                // merge lower components
                utype lc_num = get_components_num(adj_lower,lower);

                if(uc_num == 2 && lc_num == 2) // simple saddle
                    this->critical_points[i] = Point_Type::SADDLE;
                else if(uc_num > 2 && lc_num > 2) // multiple saddle
                    this->critical_points[i] = Point_Type::MULTIPLE_SADDLE;
                // else the point is already flagged as regular point
            }
        }
    }
}

void Critical_Points_Extractor::extract_critical_points_from_flat_areas(flat_areas &fa, Spatial_Mesh &mesh)
{
    map<itype,itype> upper, lower;
    map<itype,iset> adj_upper;
    map<itype,iset> adj_lower;

    for(flat_areas::iterator it=fa.begin(); it!=fa.end(); ++it)
    {
        itype key_vertex = -1;
        double field = it->first;
        iset &vt = it->second;

        for(auto tid : vt)
        {
            Triangle &t = mesh.get_triangle(tid);
            ivect link_v;
            for(int i=0; i<t.vertices_num(); i++)
            {
                Vertex &v = mesh.get_vertex(t.TV(i));
                if(v.get_c(2) != field)
                {
                    if(v.get_c(2) < field)
                        lower[t.TV(i)] = t.TV(i);
                    else
                        upper[t.TV(i)] = t.TV(i);
                    link_v.push_back(t.TV(i));
                }
                else if(key_vertex == -1)//same field and still unitialized critical vertex
                    key_vertex = t.TV(i);
            }
            if(link_v.size() == 2)
            {
                if(field < mesh.get_vertex(link_v[0]).get_c(2) &&
                        field < mesh.get_vertex(link_v[1]).get_c(2))
                {
                    adj_upper[link_v[0]].insert(link_v[1]);
                    adj_upper[link_v[1]].insert(link_v[0]);
                }
                else if(field > mesh.get_vertex(link_v[0]).get_c(2) &&
                        field > mesh.get_vertex(link_v[1]).get_c(2))
                {
                    adj_lower[link_v[0]].insert(link_v[1]);
                    adj_lower[link_v[1]].insert(link_v[0]);
                }
            }
        }

        if(upper.empty()) //the vertex is a maximum
            this->critical_points[key_vertex] = Point_Type::MAXIMUM;
        else if(lower.empty()) //the vertex is a minimum
            this->critical_points[key_vertex] = Point_Type::MINIMUM;
        else //we have to check the connected components of the edges in the link
        {
            // merge upper components
            utype uc_num = get_components_num(adj_upper,upper);
            // merge lower components
            utype lc_num = get_components_num(adj_lower,lower);

            if(uc_num == 2 && lc_num == 2) // simple saddle
                this->critical_points[key_vertex] = Point_Type::SADDLE;
            else if(uc_num > 2 && lc_num > 2) // multiple saddle
                this->critical_points[key_vertex] = Point_Type::MULTIPLE_SADDLE;
            // else the point is already flagged as regular point
        }

        // reset the local d.s.
        upper.clear();
        lower.clear();
        adj_lower.clear();
        adj_upper.clear();
    }
}

void Critical_Points_Extractor::init_adjacent_vertices(itype v_id, Vertex &v, ivect &vt, Spatial_Mesh &mesh,
                                                       map<itype,iset> &adj_upper, map<itype,iset> &adj_lower)
{
    for(ivect::iterator it=vt.begin(); it!=vt.end(); ++it)
    {
        Triangle &t = mesh.get_triangle(*it);
        int v_pos = t.vertex_index(v_id);
        Edge e = t.TE(v_pos);

        // if the other two vertices have the fields higher/lower than v
        // -> then we flag them as adjacent
        if(v.get_c(2) < mesh.get_vertex(e.EV(0)).get_c(2) &&
                v.get_c(2) < mesh.get_vertex(e.EV(1)).get_c(2))
        {
            adj_upper[e.EV(0)].insert(e.EV(1));
            adj_upper[e.EV(1)].insert(e.EV(0));
        }
        else if(v.get_c(2) > mesh.get_vertex(e.EV(0)).get_c(2) &&
                v.get_c(2) > mesh.get_vertex(e.EV(1)).get_c(2))
        {
            adj_lower[e.EV(0)].insert(e.EV(1));
            adj_lower[e.EV(1)].insert(e.EV(0));
        }
    }
}

utype Critical_Points_Extractor::get_components_num(map<itype,iset> &adj_map, map<itype,itype> &v_flag)
{
    for(auto am : adj_map)
    {
        if(am.first == v_flag[am.first]) // not set yet
        {
            itype flag = v_flag[am.first];
            iqueue q;
            q.push(am.first);

            while(!q.empty())
            {
                itype current = q.front();
                q.pop();
                if(v_flag[current] == current) // not set
                {
                    v_flag[current] = flag;
                    iset &adj = adj_map[current];
                    for(auto a : adj)
                    {
                        q.push(a);
                    }
                }
            }
        }
    }

    iset component;
    for(map<itype,itype>::iterator it=v_flag.begin(); it!=v_flag.end(); ++it)
        component.insert(it->second);
    return component.size();
}
