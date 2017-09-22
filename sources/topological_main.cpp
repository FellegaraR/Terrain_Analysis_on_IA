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

#ifndef TOPOLOGICAL_MAIN_H
#define TOPOLOGICAL_MAIN_H

#include <map>

#include "ia/mesh.h"
#include "utilities/usage.h"
#include "utilities/timer.h"

template<class M> void VT_ALL(M& mesh);
template<class M> void ALL(M& mesh);

template<class M> void ALL(M& mesh)
{
    Timer time;
    time.start();
    bool is_border;
    for(int vid=0; vid<mesh.get_vertices_num(); vid++)
    {
        ivect ret = mesh.VT(vid,is_border);
        ret.clear();
    }
    time.stop();
    time.print_elapsed_time("[TIME] batched VT extraction: ");

    cerr << "[MEMORY] peak for extracting the VT relations: " <<
            to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " Mbs" << std::endl;

    time.start();
    for(int vertexId=0; vertexId<mesh.get_vertices_num(); vertexId++)
    {
        vector<Edge> ret = mesh.VE(vertexId);
        ret.clear();
    }
    time.stop();
    time.print_elapsed_time("[TIME] batched VE extraction: ");

    cerr << "[MEMORY] peak for extracting the VE relations: " <<
            to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " Mbs" << std::endl;

    time.start();
    for(int vid=0; vid<mesh.get_vertices_num(); vid++)
    {
        ivect ret = mesh.VV(vid);
        ret.clear();
    }
    time.stop();
    time.print_elapsed_time("[TIME] batched VV extraction: ");

    cerr << "[MEMORY] peak for extracting the VV relations: " <<
            to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " Mbs" << std::endl;

    time.start();
    set<Edge> edges;
    for(int tId=0; tId<mesh.get_triangles_num(); tId++)
    {
        Triangle &triangle = mesh.get_triangle(tId);
        for(int edgePos=0; edgePos<triangle.vertices_num(); edgePos++)
        {
            Edge edge = triangle.TE(edgePos);
            pair<set<Edge>::iterator,bool> iter = edges.insert(edge);
            if(iter.second)
            {
                ivect ret = mesh.ET(edge);
                ret.clear();
            }
        }
    }
    time.stop();
    time.print_elapsed_time("[TIME] batched ET extraction: ");
    edges.clear();

    cerr << "[MEMORY] peak for extracting the ET relations: " <<
            to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " Mbs" << std::endl;

    time.start();
    for(int tId=0; tId<mesh.get_triangles_num(); tId++)
    {
        Triangle &triangle = mesh.get_triangle(tId);
        for(int edgePos=0; edgePos<triangle.vertices_num(); edgePos++)
        {
            Edge edge = triangle.TE(edgePos);
            pair<set<Edge>::iterator,bool> iter = edges.insert(edge);
            if(iter.second)//I don't have computed the EE for the current edge
            {
                vector<Edge> ret = mesh.EE(edge);
                ret.clear();
            }
        }
    }
    time.stop();
    time.print_elapsed_time("[TIME] batched EE extraction: ");
    edges.clear();

    cerr << "[MEMORY] peak for extracting the EE relations: " <<
            to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " Mbs" << std::endl;

    return;
}

template<class M> void VT_ALL(M& mesh)
{
    Timer time;
    time.start();
    for(itype vid=0; vid<mesh.get_vertices_num(); vid++)
    {
        bool is_border=false;
        ivect ret = mesh.VT(vid,is_border);
        ret.clear();
    }
    time.stop();
    time.print_elapsed_time("[TIME] batched VT extraction: ");

    cerr << "[MEMORY] peak for extracting the VT relations: " <<
            to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " Mbs" << std::endl;

    return;
}

#endif // TOPOLOGICAL_MAIN_H
