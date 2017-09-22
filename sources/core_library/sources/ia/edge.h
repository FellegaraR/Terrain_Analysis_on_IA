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

#ifndef EDGE_H
#define EDGE_H

#include <algorithm>
#include "utilities/basic_wrappers.h"

class Edge
{
public:
    Edge()
    {
        this->vertices = { -1, -1 };
    }

    Edge(itype v1, itype v2)
    {
        this->vertices = { std::min(v1,v2), std::max(v1,v2) };
    }

    inline itype EV(itype pos) { return this->vertices[pos]; }

    inline friend bool operator== (const Edge &p, const Edge &q)
    {
        if(p.vertices[0]!=q.vertices[0]) return false;
        if(p.vertices[1]!=q.vertices[1]) return false;
        return true;
    }

    inline friend bool operator!= (const Edge &p, const Edge &q) { return !(p==q); }

    inline bool operator<(const Edge& s) const
    {
        return ((this->vertices[0] < s.vertices[0]) || (this->vertices[0] == s.vertices[0] && this->vertices[1] < s.vertices[1]));
    }
    inline bool operator>(const Edge& s) const
    {
        return ((this->vertices[0] > s.vertices[0]) || (this->vertices[0] == s.vertices[0] && this->vertices[1] > s.vertices[1]));
    }

private:
    ivect vertices;
};


#endif // EDGE_H
