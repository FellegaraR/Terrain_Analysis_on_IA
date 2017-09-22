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

#ifndef CRITICAL_POINTS_EXTRACTOR_H
#define CRITICAL_POINTS_EXTRACTOR_H

#include <map>

#include "ia/mesh.h"
#include "utilities/basic_wrappers.h"

// the key corresponds to the field value of the flat area
// the value correspond to the vertices and triangles belonging to the flat area
typedef map<double,iset> flat_areas;

enum class Point_Type: short {REGULAR=0, MINIMUM=1, SADDLE=2, MULTIPLE_SADDLE=3, MAXIMUM=4};

class Critical_Points_Extractor
{
public:
    //
    Critical_Points_Extractor() { }
    //
    void compute_critical_points(Spatial_Mesh &mesh);
    //
    inline void print_stats()
    {
        int num_reg = 0, num_min = 0, num_max = 0, num_saddle = 0, num_multisaddle = 0;
        for(vector<Point_Type>::iterator it=critical_points.begin(); it!=critical_points.end(); ++it)
        {
            switch(*it)
            {
            case Point_Type::REGULAR:
                num_reg++;
                break;
            case Point_Type::MINIMUM:
                num_min++;
                break;
            case Point_Type::SADDLE:
                num_saddle++;
                break;
            case Point_Type::MULTIPLE_SADDLE:
                num_multisaddle++;
                break;
            case Point_Type::MAXIMUM:
                num_max++;
                break;
            }
        }
        cerr<<"[STAT] Extracted Critical points"<<endl;
        cerr<<"   mesh vertices: "<<critical_points.size()<<endl;
        cerr<<"   regular_points: "<<num_reg<<" -- minima: "<<num_min<<endl;
        cerr<<"   saddles: "<<num_saddle<<" -- multi-saddles: "<<num_multisaddle<<endl;
        cerr<<"   maxima: "<<num_max<<endl;
    }
    inline vector<Point_Type>& get_critical_points() { return this->critical_points; }

private:
    vector<Point_Type> critical_points;

    void extract_critical_points(Spatial_Mesh &mesh, flat_areas &fa); /// buggy
    void extract_critical_points(Spatial_Mesh &mesh);
    void extract_critical_points_from_flat_areas(flat_areas &fa, Spatial_Mesh &mesh); /// buggy
    void init_adjacent_vertices(itype v_id, Vertex &v, ivect &vt, Spatial_Mesh &mesh, map<itype, iset> &adj_upper,
                                map<itype, iset> &adj_lower);
    utype get_components_num(map<itype, iset> &adj_map, map<itype, itype> &v_flag);

};

#endif // CRITICAL_POINTS_EXTRACTOR_H
