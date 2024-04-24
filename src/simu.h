#pragma once

#include "mesh.h"
#include "fem.h"
#include <math.h>
#include <cmath>
#include <iostream>

namespace FEM2A {
    namespace Simu {

        //#################################
        //  Useful functions
        //#################################

        double unit_fct( vertex v )
        {
            return 1.;
        }

        double zero_fct( vertex v )
        {
            return 0.;
        }

        double xy_fct( vertex v )
        {
            return v.x + v.y;
        }

        //#################################
        //  Simulations
        //#################################

        void pure_dirichlet_pb( const std::string& mesh_filename, bool verbose, bool source)
        {
// Sans terme source :

            if (not source){
		    std::cout << "Solving a pure Dirichlet problem without source" << std::endl;

		    if ( verbose ) {
		        std::cout << " with lots of printed details..." << std::endl;
		    }

		    Mesh mesh;
		    mesh.load(mesh_filename);

		    Quadrature quad = Quadrature::get_quadrature(2);

		    ShapeFunctions shp_func_triangle(2, 1);

		    SparseMatrix K(mesh.nb_vertices());// déclaration de K

		    std::vector<double> F(mesh.nb_vertices(), 0);// init de F
		    
		    std::vector<bool> attribute_dirichlet(2, false);
		    attribute_dirichlet[1] = true;
		    mesh.set_attribute(unit_fct, 1, true);
		    
		    std::vector< double > values(mesh.nb_vertices());
		    
		    for (int elem_ind = 0; elem_ind < mesh.nb_triangles(); elem_ind ++){//triangle

		    	DenseMatrix Ke;
			Ke.set_size(3,3);

		    	ElementMapping elt_map(mesh, false, elem_ind);
		    	assemble_elementary_matrix(elt_map, shp_func_triangle, quad, unit_fct, Ke);
		    	local_to_global_matrix( mesh, elem_ind, Ke, K);
		    }
		    for (int k = 0; k < mesh.nb_vertices(); k ++){
		   	values[k] = xy_fct(mesh.get_vertex(k));		
		    }
		    apply_dirichlet_boundary_conditions(mesh, attribute_dirichlet, values, K, F);
		    std::vector<double> u(mesh.nb_vertices());
		    solve(K,F, u);
		    std::string export_name = "Dirichlet_pur";
		    mesh.save(export_name+".mesh");
		    save_solution(u, export_name+".bb");
	    
	    }
	    
// Avec terme source :
       	    if (source){
		    std::cout << "Solving a pure Dirichlet problem with source" << std::endl;

		    if ( verbose ) {
		        std::cout << " with lots of printed details..." << std::endl;
		    }

		    Mesh mesh;
		    mesh.load(mesh_filename);

		    Quadrature quad = Quadrature::get_quadrature(2);

		    ShapeFunctions shp_func_triangle(2, 1);

		    SparseMatrix K(mesh.nb_vertices());// déclaration de K

		    std::vector<double> F(mesh.nb_vertices(), 0);// init de F
		    
		    std::vector<bool> attribute_dirichlet(2, false);
		    attribute_dirichlet[1] = true;
		    mesh.set_attribute(unit_fct, 1, true);
		    
		    std::vector< double > values(mesh.nb_vertices());
		    
		    for (int elem_ind = 0; elem_ind < mesh.nb_triangles(); elem_ind ++){//triangle		

		    	ElementMapping elt_map(mesh, false, elem_ind);
		    	// Assemblage K
		    	DenseMatrix Ke;
			Ke.set_size(3,3);
		    	assemble_elementary_matrix(elt_map, shp_func_triangle, quad, unit_fct, Ke);
		    	local_to_global_matrix( mesh, elem_ind, Ke, K);
		    	//Assemblage F
		    	std::vector< double > Fe;
		    	assemble_elementary_vector(elt_map, shp_func_triangle, quad, unit_fct, Fe);
		    	std::cout<<Fe[0]<<Fe[1]<<Fe[2]<<std::endl;
		    	local_to_global_vector(mesh, false, elem_ind, Fe, F);		    	
		    }
		    
		    for (int k = 0; k < mesh.nb_vertices(); k ++){
		   	values[k] = zero_fct(mesh.get_vertex(k));

		    }
		    
		    apply_dirichlet_boundary_conditions(mesh, attribute_dirichlet, values, K, F);
		    std::vector<double> u(mesh.nb_vertices());
		    solve(K,F, u);
		    std::string export_name = "Dirichlet_pur_source";
		    mesh.save(export_name+".mesh");
		    save_solution(u, export_name+".bb");
	     }

        }

    }

}
