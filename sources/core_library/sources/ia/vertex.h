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

#ifndef _VERTEX_H
#define	_VERTEX_H

#include <math.h>
#include "ia/vt_star.h"

///An inner-class, extending Vertex2D, representing a vertex in a tetrahedral mesh
class Vertex : public VT_star
{
public:
    ///A constructor method
    Vertex() : VT_star() {  }
    ///A constructor method
    Vertex(const Vertex& orig) : VT_star(orig) { this->coords=orig.coords; }
    ///A constructor method for a 2D point
    /*!
     * \param x a coord_type argument, representing the x coordinate
     * \param y a coord_type argument, representing the y coordinate
     */
    Vertex(coord_type x, coord_type y) : VT_star() { coords = { x, y }; }
    ///A constructor method for a 3D point
    /*!
     * \param x a coord_type argument, representing the x coordinate
     * \param y a coord_type argument, representing the y coordinate
     * \param z a coord_type argument, representing the z coordinate
     */
    Vertex(coord_type x, coord_type y, coord_type z) : VT_star() { coords = { x, y, z }; }
    ///A constructor method for a nD point
    Vertex(dvect &coords) : VT_star() { this->coords = coords; }
    ///A constructor method for an unset nD point
    Vertex(int dim) : VT_star() { this->coords.assign(dim,0); }
    ///A destructor method
    virtual ~Vertex() {}
    ///
    friend bool operator== (const Vertex& p, const Vertex &q)
    {
        if(p.get_dimension() != q.get_dimension())
            return false;
        for(int i=0; i<p.get_dimension(); i++)
        {
            if(p.get_c(i) != q.get_c(i))
                return false;
        }
        return true;
    }
    ///
    friend bool operator!= (const Vertex& p, const Vertex &q) { return !(p == q); }

    inline Vertex operator+(const Vertex& other)
    {
        Vertex copy = *this;
        copy+=other;
        return copy;
    }
    inline Vertex& operator+=(const Vertex& other)
    {
        for(int i=0; i<this->get_dimension(); i++)
            this->coords[i] += other.coords[i];
        return *this;
    }
    inline Vertex operator-(const Vertex& other)
    {
        Vertex copy = *this;
        copy -= other;
        return copy;
    }
    inline Vertex& operator-=(const Vertex& other)
    {
        for(int i=0; i<this->get_dimension(); i++)
            this->coords[i] -= other.coords[i];
        return *this;
    }
    inline Vertex operator*(const Vertex& other)
    {
        Vertex copy = *this;
        copy*=other;
        return copy;
    }
    inline Vertex& operator*=(const Vertex& other)
    {
        for(int i=0; i<this->get_dimension(); i++)
            this->coords[i] *= other.coords[i];
        return *this;
    }
    inline Vertex operator*(const coord_type c)
    {
        Vertex copy = *this;
        copy *= c;
        return copy;
    }
    inline Vertex& operator*=(const coord_type c)
    {
        for(int i=0; i<this->get_dimension(); i++)
            this->coords[i] *= c;
        return *this;
    }
    inline Vertex operator/(const Vertex& other)
    {
        Vertex copy = *this;
        copy/=other;
        return copy;
    }
    inline Vertex& operator/=(const Vertex& other)
    {
        for(int i=0; i<this->get_dimension(); i++)
            this->coords[i] /= other.coords[i];
        return *this;
    }
    inline Vertex operator/(const coord_type c)
    {
        Vertex copy = *this;
        copy /= c;
        return copy;
    }
    inline Vertex& operator/=(const coord_type c)
    {
        for(int i=0; i<this->get_dimension(); i++)
            this->coords[i] /= c;
        return *this;
    }

    inline int get_dimension() const { return coords.size(); }
    // no range check!
    inline void set_c(int pos, coord_type c) { this->coords[pos] = c; }
    // no range check!
    inline coord_type get_c(int pos) const { return this->coords[pos]; }

    //this function returns the norm of vector vec-v
    coord_type norm(Vertex& v)
    {
        coord_type sum = 0;
        for(int i=0; i<get_dimension(); i++)
        {
            sum += (v.get_c(i)-this->get_c(i))*(v.get_c(i)-this->get_c(i));
        }
        return sqrt(sum);
    }

    //this function returns the scalar products between vectors v1-vec and v2-vec
    coord_type scalar_product(Vertex& v1,Vertex& v2)
    {
        coord_type sum = 0;
        for(int i=0; i<get_dimension(); i++)
        {
            sum += (v1.get_c(i)-this->get_c(i))*(v2.get_c(i)-this->get_c(i));
        }
        return sum;
    }

    //norm of this vector
    coord_type norm()
    {
        coord_type sum = 0;
        for(int i=0; i<get_dimension(); i++)
        {
            sum += this->get_c(i)*this->get_c(i);
        }
        return sqrt(sum);
    }

    //scalar products between vectors vec and v1
    coord_type scalar_product(Vertex& v1)
    {
        coord_type sum = 0;
        for(int i=0; i<get_dimension(); i++)
        {
            sum += (v1.get_c(i)*this->get_c(i));
        }
        return sum;
    }

    inline dvect& get_coordinates() { return coords; }

    /**
     * @brief operator <<
     * @param out
     * @param p
     * @return
     */
    inline friend std::ostream& operator<<(std::ostream& out, const Vertex& p)
    {
        for(auto c : p.coords)
            out << c << " ";
        return out;
    }

private:    
    ///A protected variable representing the z coordinate of the point
    dvect coords;

};

#endif	/* _VERTEX_H */

