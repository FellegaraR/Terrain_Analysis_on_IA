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

#include <cstdlib>
#include <iostream>
#include <string.h>
#include <stdio.h>
#include <fstream>

#include "ia/mesh.h"
#include "curvature/mean_curvature.h"
#include "curvature/concentrated_curvature.h"
#include "curvature/c_curvature.h"

#include "terrain_features/critical_points_extractor.h"
#include "terrain_features/slope_extractor.h"

#include "topological_main.cpp"

#include "utilities/quad_mesh.h"
#include "utilities/string_management.h"
#include "utilities/usage.h"
#include "utilities/io.h"
#include "utilities/timer.h"

using namespace std;
using namespace string_management;

#define BOLD  "\033[1m\033[33m" //for dark background shell
//#define BOLD "\033[1m\033[31m"  //for white background shell
#define RESET   "\033[0m"

void print_help();
void print_paragraph(string stringa, int cols);

int main(int argc, char *argv[])
{
    if(argc == 1)
    {
        print_help();
        return 0;
    }

    if(argc < 2)
    {
        cout << "[ERROR] too few arguments" << endl;
        return -1;
    }
    if(argc > 4)
    {
        cout << "[ERROR] too many arguments" << endl;
        return -1;
    }

    Timer time;
    Spatial_Mesh mesh = Spatial_Mesh();
    IO::read_mesh(mesh,argv[2]);
    cerr << argv[2] << " vertices: " << mesh.get_vertices_num() << " triangles: " << mesh.get_triangles_num() << endl;
    cerr << "[MEMORY] peak for loading the terrain: " << to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;
    time.start();
    bool generated = mesh.build();
    if(!generated)
    {
        cerr<<"[ERROR] Failed to generate the IA data structure."<<endl;
        return -1;
    }
    time.stop();
    time.print_elapsed_time("[TIME] IA generation: ");
    cerr << "[MEMORY] peak for generating the IA: " << to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;

    if(strcmp(argv[1],"concurv")==0)
    {
        ConcentratedCurvature ccurv = ConcentratedCurvature(true,mesh);
        time.start();
        ccurv.compute_values(mesh);
        time.stop();
        time.print_elapsed_time("[TIME] Computing Concentrated Curvatures: ");
        cerr << "[MEMORY] peak for extracting the Concentrated Curvatures: " <<
                to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;
    }
    else if(strcmp(argv[1],"meancurv")==0)
    {
        MeanCurvature mcurv = MeanCurvature(true,mesh);
        time.start();
        mcurv.compute_values(mesh);
        time.stop();
        time.print_elapsed_time("[TIME] Computing Mean Curvatures: ");
        cerr << "[MEMORY] peak for extracting the Mean Curvatures: " <<
                to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;
    }
    else if(strcmp(argv[1],"gccurv")==0 || strcmp(argv[1],"mccurv")==0)
    {
        CCurvatureType type = (strcmp(argv[1],"gccurv")==0)? GAUSS_CCURVATURE : MEAN_CCURVATURE;
        C_Curvature ccurv = C_Curvature(type,mesh,(type==MEAN_CCURVATURE));
        time.start();
        ccurv.compute_values(mesh);
        time.stop();
        time.print_elapsed_time("[TIME] Computing CCurvatures: ");
        cerr << "[MEMORY] peak for extracting the C Curvatures: " <<
                to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;
    }
    else if(strcmp(argv[1],"vtall")==0)
        VT_ALL(mesh);
    else if(strcmp(argv[1],"all")==0)
        ALL(mesh);
    else if(strcmp(argv[1],"quad")==0)
    {
        Quad_Mesh quad_mesh;
        quad_mesh.loadQuad(mesh);
        cerr << "[MEMORY] peak for extracting the dual Quad mesh: " <<
                to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;
        cout<<mesh.get_vertices_num()<<" "<<mesh.get_triangles_num()<<endl;
        quad_mesh.save_quad_mesh(string_management::get_path_without_file_extension(argv[2]));
    }
    else if(strcmp(argv[1],"eslope")==0)
    {
        Slope_Extractor se;
        time.start();
        se.compute_edges_slopes(mesh);
        time.stop();
        time.print_elapsed_time("[TIME] Computing Edge Slopes: ");
        cerr << "[MEMORY] peak for extracting the edge slopes: " <<
                to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;
    }
    else if(strcmp(argv[1],"tslope")==0)
    {
        Slope_Extractor se;
        time.start();
        se.compute_triangles_slopes(mesh);
        time.stop();
        time.print_elapsed_time("[TIME] Computing Triangle Slopes: ");
        cerr << "[MEMORY] peak for extracting the triangles slopes: " <<
                to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;
    }
    else if(strcmp(argv[1],"crit")==0)
    {
        Critical_Points_Extractor cpe;
        time.start();
        cpe.compute_critical_points(mesh);
        time.stop();
        time.print_elapsed_time("[TIME] Computing Critical Points: ");
        cerr << "[MEMORY] peak for extracting the critical points: " <<
                to_string(MemoryUsage().get_Virtual_Memory_in_MB()) << " MBs" << std::endl;
        cpe.print_stats();
    }
    else if(strcmp(argv[1],"save")==0)
    {
        cout<<"[NOTA] Saving mesh connectivity."<<endl;
        IO::write_mesh_connectivity(mesh,string_management::get_path_without_file_extension(argv[2]));
    }
    else
        cout<<"error reading arguments"<<endl;

    return 0;
}

void print_help(){

    //annoying stuff to get the dimension of the output shell (!!! not sure it works on Mac,
    //everything based on the command tput. If it doesn't work the dimension si setted to 80 by default)
    FILE* fp;
    char path[1035];

    int cols;
    fp = popen("tput cols", "r");
    if(fp != NULL){
        fgets(path, sizeof(path)-1, fp);
        cols = atoi(path);
    }
    else{
        cols = 80;
    }

    //printf start
    printf(BOLD "\n  NAME:\n\n" RESET);
    printf("\tLibTri library - Adjacency-based data structure for representing and analyzing triangle meshes.\n\n" RESET);

    printf(BOLD "  USAGE: \n\n" RESET);
    printf(BOLD "    .\\libtri [operation] [mesh_name]\n\n" RESET);
    print_paragraph("NOTA: the arguments order is fixed.", cols);

    printf(BOLD "    [operation]\n\n" RESET);
    print_paragraph("the operation argument can be vtall, all, meancurv, concurv, gcurv, mccurv, eslope, tslope, crit.",cols);
    printf(BOLD "        vtall\n" RESET); print_paragraph(" extracts all the VT relations of the input mesh (prints timings - no output).",cols);
    printf(BOLD "        all\n" RESET); print_paragraph(" extracts all the topological relations of the input mesh (prints timings - no output).",cols);
    printf(BOLD "        meancurv\n" RESET); print_paragraph(" computes the Mean Curvature for all the mesh vertices.",cols);
    printf(BOLD "        concurv\n" RESET); print_paragraph(" computes the Concentrated Curvature for all the mesh vertices.",cols);
    printf(BOLD "        mccurv\n" RESET); print_paragraph(" computes the Mean CCurvature for all the mesh vertices.",cols);
    printf(BOLD "        gccurv\n" RESET); print_paragraph(" computes the Gauss CCurvature for all the mesh vertices.",cols);
    printf(BOLD "        quad\n" RESET); print_paragraph(" extracts the dual quad mesh from the input mesh and saves it in off format.",cols);
    printf(BOLD "        eslope\n" RESET); print_paragraph(" computes the the slope values for each edge of the mesh.",cols);
    printf(BOLD "        tslope\n" RESET); print_paragraph(" computes the the slope values for each triangle of the mesh.",cols);
    printf(BOLD "        crit\n" RESET); print_paragraph(" computes the critical points of the mesh.",cols);

    printf(BOLD "    [mesh_name]\n\n" RESET);
    print_paragraph("the mesh_name argument represents the triangular mesh (in .tri format).",cols);

    printf(BOLD "  EXAMPLE: \n\n" RESET);
    printf("          .\\libtri vtall mesh.tri\n\n");
    print_paragraph("read as input file the mesh [mesh.tri] and the vt relations for all vertices.", cols);


    printf(BOLD "  IMPLEMENTATION:\n\n" RESET);
    printf("          Author: Riccardo Fellegara\n");
    printf("          Group: UMD GeoVisualization Group\n");
    printf("          Last Update: September 2017\n\n");

    printf(BOLD "  DESCRIPTION: \n\n" RESET);
    print_paragraph("A 2D implementation of the IA data structure for the representation "
                    "and analysis triangle meshes. We have developed a tool based on the LibTri for "
                    "terrain analysis, which includes state-of-the-art estimators for slope "
                    "and curvature, and for the extraction of critical points. .", cols);
}

void print_paragraph(string stringa, int cols){
    if((int)stringa.size() < cols-20){
        printf("          %s\n\n",stringa.c_str());
    }
    else{
        float dim = (float)(stringa.size()/((float)cols-20));
        int dim_int = dim;
        if(dim > dim_int) dim_int++;
        for(int i=0; i<dim_int; i++){
            if(stringa.at((cols-20)*i) == ' ')
                printf("         %s\n", stringa.substr( (cols-20)*i, cols-20).c_str());
            else printf("          %s\n", stringa.substr( (cols-20)*i, cols-20).c_str());
        }
        printf("\n");
    }
}



