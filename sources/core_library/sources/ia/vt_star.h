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

#ifndef VT_STAR_H
#define VT_STAR_H

#include "utilities/basic_wrappers.h"

class VT_star
{
public:
    ///A constructor method
    VT_star() { this->vtstar = -1; }
    ///A constructor method
    VT_star(const VT_star& orig) { this->vtstar = orig.vtstar; }
    ///A destructor method
    virtual ~VT_star() {}
    ///
    friend bool operator== (const VT_star& p, const VT_star &q) { return (p.vtstar == q.vtstar); }
    ///
    friend bool operator!= (const VT_star& p, const VT_star &q) { return !(p == q); }
    ///
    itype get_VTstar() { return this->vtstar; }
    ///
    void set_VTstar(itype vtstar) { this->vtstar = vtstar; }
protected:
    ///
    itype vtstar;
};

#endif // VT_STAR_H
