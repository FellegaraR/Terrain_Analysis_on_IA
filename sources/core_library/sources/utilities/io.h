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

#ifndef _IO_H
#define	_IO_H

#include <vector>
#include <string>

#include <sstream>
#include <algorithm>
#include <iostream>
#include <fstream>

#include "ia/vt_star.h"
#include "ia/vertex.h"
#include "ia/mesh.h"
#include "ia/triangle.h"

using namespace std;
///A class that handles the input/output to initialize the main library structures
class IO {
public:
    ///A public method that reads a file containing a triangle mesh (TRI format)
    /*!
     * \param mesh a Mesh& argument, representing the mesh to initialize
     * \param path a string argument, representing the path to the mesh file
     * \return a boolean value, true if the file is correctly readed, false otherwise
     */
    static bool read_mesh(Spatial_Mesh& mesh, string path);

    static bool write_mesh_connectivity(Spatial_Mesh& mesh, string path);
private:
    ///A constructor method
    IO() {}
    ///A constructor method
    IO(const IO& ) {}
    ///A destructor method
    virtual ~IO() {}

    ///A public method that reads a file containing a triangle mesh (TRI format)
    /*!
     * \param mesh a Mesh& argument, representing the mesh to initialize
     * \param path a string argument, representing the path to the mesh file
     * \return a boolean value, true if the file is correctly readed, false otherwise
     */
    static bool read_mesh_off(Spatial_Mesh& mesh, string path);
    ///A public method that reads a file containing a triangle mesh (TRI format)
    /*!
     * \param mesh a Mesh& argument, representing the mesh to initialize
     * \param path a string argument, representing the path to the mesh file
     * \return a boolean value, true if the file is correctly readed, false otherwise
     */
    static bool read_mesh_tri(Spatial_Mesh& mesh, string path);
    /**
     * @brief A private method that reads a list of vertices with just the coordinates
     *
     * @param mesh represents the mesh to initialize
     * @param input represents the input stream
     * @param num_vertices represents the number of vertices of the mesh
     *
     * NOTA: the procedure also set the mesh domain
     */
    static void read_vertices_list(Spatial_Mesh& mesh, ifstream &input, itype num_vertices);
};

#endif	/* _IO_H */

