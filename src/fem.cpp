#include "fem.h"
#include "mesh.h"

#include <iomanip>
#include <iostream>
#include <cmath>
#include <limits>
#include <stdlib.h>
#include <assert.h>

namespace FEM2A {

    void print( const std::vector<double>& x )
    {
        for ( int i = 0; i < x.size(); ++i ) {
            std::cout << x[i] << " ";
        }
        std::cout << std::endl;
    }

    /****************************************************************/
    /* Implementation of Quadrature */
    /****************************************************************/
    int Quadrature::nb_points() const
    {
        return wxy_.size() / 3 ;
    }

    vertex Quadrature::point( int i ) const
    {
        assert( i < nb_points() ) ;
        vertex v ;
        v.x = wxy_[3 * i + 1] ;
        v.y = wxy_[3 * i + 2] ;
        return v ;
    }

    double Quadrature::weight( int i ) const
    {
        assert( i < nb_points() ) ;
        return wxy_[3 * i + 0] ;
    }

    const double triangle_P0[3] = {
        0.5, 0.333333333333333, 0.333333333333333
    };

    const double triangle_P2[9] = {
        0.166666666666667, 0.166666666666667, 0.166666666666667,
        0.166666666666667, 0.166666666666667, 0.666666666666667,
        0.166666666666667, 0.666666666666667, 0.166666666666667
    };

    const double triangle_P4[18] = {
        0.0549758718276609, 0.0915762135097707, 0.0915762135097707,
        0.0549758718276609, 0.0915762135097707, 0.816847572980459,
        0.0549758718276609, 0.816847572980459, 0.0915762135097707,
        0.111690794839006, 0.445948490915965, 0.445948490915965,
        0.111690794839006, 0.445948490915965, 0.10810301816807,
        0.111690794839006, 0.10810301816807, 0.445948490915965
    };

    const double triangle_P6[36] = {
        0.0254224531851034, 0.0630890144915022, 0.0630890144915022,
        0.0254224531851034, 0.0630890144915022, 0.873821971016996,
        0.0254224531851034, 0.873821971016996, 0.0630890144915022,
        0.0583931378631897, 0.24928674517091, 0.24928674517091,
        0.0583931378631897, 0.24928674517091, 0.501426509658179,
        0.0583931378631897, 0.501426509658179, 0.24928674517091,
        0.0414255378091868, 0.0531450498448169, 0.310352451033784,
        0.0414255378091868, 0.310352451033784, 0.0531450498448169,
        0.0414255378091868, 0.0531450498448169, 0.636502499121399,
        0.0414255378091868, 0.636502499121399, 0.0531450498448169,
        0.0414255378091868, 0.310352451033784, 0.636502499121399,
        0.0414255378091868, 0.636502499121399, 0.310352451033784
    };

    const double segment_P0[2] = {
        1., 0.5
    };

    const double segment_P2[4] = {
        0.5, 0.21132486540518708,
        0.5, 0.7886751345948129
    };

    Quadrature Quadrature::get_quadrature( int order, bool border )
    {
        double* pts = NULL;
        int nb_pts = 0;
        Quadrature Q;
        if ( order == 0 && !border ) {
            pts = const_cast<double*>(triangle_P0);
            nb_pts = 1;
        } else if ( order == 2 && !border ) {
            pts = const_cast<double*>(triangle_P2);
            nb_pts = 3;
        } else if ( order == 4 && !border ) {
            pts = const_cast<double*>(triangle_P4);
            nb_pts = 6;
        } else if ( order == 6 && !border ) {
            pts = const_cast<double*>(triangle_P6);
            nb_pts = 12;
        } else if ( order == 0 && border ) {
            pts = const_cast<double*>(segment_P0);
            nb_pts = 1;
        } else if ( order == 2 && border ) {
            pts = const_cast<double*>(segment_P2);
            nb_pts = 2;
        } else {
            //std::cout << "Quadrature not implemented for order " << order << std::endl;
            assert( false );
        }
        Q.wxy_.resize(nb_pts * 3);
        for ( int i = 0; i < nb_pts; ++i ) {
            if ( !border ) {
                Q.wxy_[3*i+0] = pts[3*i+0];
                Q.wxy_[3*i+1] = pts[3*i+1];
                Q.wxy_[3*i+2] = pts[3*i+2];
            } else {
                Q.wxy_[3*i+0] = pts[2*i+0];
                Q.wxy_[3*i+1] = pts[2*i+1];
                Q.wxy_[3*i+2] = 0.;
            }
        }
        return Q;
    }

    /****************************************************************/
    /* Implementation of ElementMapping */
    /****************************************************************/
    ElementMapping::ElementMapping( const Mesh& M, bool border, int i )
        : border_( border )
    {
        // vertice_ = vecteur de vertex
        //std::cout << "[ElementMapping] constructor for element " << i << std::endl;
        if ( border ){
            //std::cout << "border"<<std::endl;
            for (int ind_local = 0; ind_local < 2; ind_local ++){
                vertices_.push_back(M.get_edge_vertex(i,ind_local));
                //std::cout<<"Vertex numéro "<<ind_local<<"\nX : " << vertices_[ind_local].x << "\nY : " << vertices_[ind_local].y << std::endl;;
            }
        }
        if (not border){
            //std::cout<<"triangle"<<std::endl;
            for (int ind_local = 0; ind_local < 3; ind_local ++){
                vertices_.push_back(M.get_triangle_vertex(i, ind_local));
                //std::cout<<"Vertex numéro "<<ind_local<<"\nX : " << vertices_[ind_local].x << "\nY : "<<vertices_[ind_local].y << std::endl;
            }
        }
    }

    vertex ElementMapping::transform( vertex x_r ) const
    {
        //std::cout << "[ElementMapping] transform reference to world space" << '\n';
        
        vertex r ;
        if ( border_ ){
            double x = 0;
            x = vertices_[0].x * (1 - x_r.x ) + vertices_[1].x * (x_r.x);
       	    double y = 0;
       	    y = vertices_[0].y * (1 - x_r.x) + vertices_[1].y * (x_r.x);
            r.x = x;
            r.y = y;
        }
        if (not border_){
	    double x = 0;
            x = vertices_[0].x * (1 - x_r.x - x_r.y) + vertices_[1].x * (x_r.x) + vertices_[2].x * (x_r.y);
       	    double y = 0;
       	    y = vertices_[0].y * (1 - x_r.x - x_r.y) + vertices_[1].y * (x_r.x) + vertices_[2].y * (x_r.y);
            r.x = x;
            r.y = y;
            
        	}
        return r ;
        }



    DenseMatrix ElementMapping::jacobian_matrix( vertex x_r ) const
    {
        //std::cout << "[ElementMapping] compute jacobian matrix" << '\n';
        
        DenseMatrix J ;
        if (border_){
        	J.set_size(2,1);
        	J.set(0, 0, vertices_[1].x - vertices_[0].x);
        	J.set(1, 0, vertices_[1].y - vertices_[0].y);        	
        	
        }
	if (not border_){
        	J.set_size(2,2);
        	J.set(0, 0, vertices_[1].x - vertices_[0].x);
        	J.set(0, 1, vertices_[2].x - vertices_[0].x);        	
        	J.set(1, 0, vertices_[1].y - vertices_[0].y);        	
        	J.set(1, 1, vertices_[2].y - vertices_[0].y);        	
        }
        return J ;
    	
    }

    double ElementMapping::jacobian( vertex x_r ) const
    {
        //std::cout << "[ElementMapping] compute jacobian determinant" << '\n';
        DenseMatrix J = jacobian_matrix( x_r );
        double det;
        if (border_){
		det = std::sqrt(J.get(0,0)*J.get(0,0) + J.get(1,0)*J.get(1,0));
	}
	if (not border_){
		det = J.det_2x2();
	}
        return det ;
    }

    /****************************************************************/
    /* Implementation of ShapeFunctions */
    /****************************************************************/
    ShapeFunctions::ShapeFunctions( int dim, int order )
        : dim_( dim ), order_( order )
    {
        //std::cout << "[ShapeFunctions] constructor in dimension " << dim << '\n';
    }

    int ShapeFunctions::nb_functions() const
    {
        //std::cout << "[ShapeFunctions] number of functions" << '\n';
        int nb = 0;
        if (dim_ == 1){nb = 2;}
        if (dim_ == 2){nb = 3;}
        return nb;
    }

    double ShapeFunctions::evaluate( int i, vertex x_r ) const
    {
        //std::cout << "[ShapeFunctions] evaluate shape function " << i << '\n';
        double x;
        if (dim_ == 1){
        	if ( i == 0 ){x = 1 - x_r.x;}
        	if ( i == 1 ){x = x_r.x;}
        }
        if (dim_ == 2){
        	if (i == 0){x = 1 - x_r.x - x_r.y;}
        	if (i == 1){x = x_r.x;}
        	if (i == 2){x = x_r.y;}
        }
        return x;
    }

    vec2 ShapeFunctions::evaluate_grad( int i, vertex x_r ) const
    {
        //std::cout << "[ShapeFunctions] evaluate gradient shape function " << i << '\n';
	
        vec2 g ;

        if (dim_ == 1){
        	if ( i == 0 ){g.x = -1; g.y = 0 ;}
        	if ( i == 1 ){g.x = 1; g.y = 0 ;}
        }
        if (dim_ == 2){
        	if (i == 0){g.x = -1; g.y = -1 ;}
        	if (i == 1){g.x = 1; g.y = 0 ;}
        	if (i == 2){g.x = 0; g.y = 1 ;}
        }
        return g ;
    }

    /****************************************************************/
    /* Implementation of Finite Element functions */
    /****************************************************************/
    void assemble_elementary_matrix(
        const ElementMapping& elt_mapping,
        const ShapeFunctions& reference_functions,
        const Quadrature& quadrature,
        double (*coefficient)(vertex),//toute fonction (k) qui prend un vertex \
        en entré et renvoi un double.
        DenseMatrix& Ke )
    {
        //std::cout << "Compute elementary matrix" << '\n';
        int p = reference_functions.nb_functions();// = nombre de fonctions d'interpolation

	for (int i = 0; i < p; i ++){//ligne i
		for (int j = 0; j < p; j ++){//colonne j
			double sum = 0;			
			for (int q = 0; q < p; q ++){ // point d'interpolation q
				double wq = quadrature.weight(q); // poids au point de Gauss
				vertex x_r = quadrature.point(q); // point de Gauss
				
				DenseMatrix J = elt_mapping.jacobian_matrix(x_r);
				J = J.invert_2x2();
				J = J.transpose();
				double det_J = elt_mapping.jacobian(x_r); //déterminant du Jacobien
				vertex x = elt_mapping.transform(x_r); //transformé du point de Gauss
				vec2 gradJ_i = reference_functions.evaluate_grad(i, x_r);
				vec2 gradJ_j = reference_functions.evaluate_grad(j, x_r);
				gradJ_i = J.mult_2x2_2(gradJ_i);
				gradJ_j = J.mult_2x2_2(gradJ_j);
				double dot_grad = dot(gradJ_i,gradJ_j);
				sum = sum + (wq * coefficient(x) * dot_grad * det_J);
				}
			Ke.set(i,j,sum);
			 
		}		
	}
    }

    void local_to_global_matrix(
        const Mesh& M,
        int t,
        const DenseMatrix& Ke,
        SparseMatrix& K )
    {
        //std::cout << "Ke -> K" << '\n';
        
	for (int ind_local_i = 0; ind_local_i < 3; ind_local_i ++){
		int ind_glob_i = M.get_triangle_vertex_index(t, ind_local_i);
		for (int ind_local_j = 0; ind_local_j < 3; ind_local_j ++){
			int ind_glob_j = M.get_triangle_vertex_index(t, ind_local_j);
			K.add(ind_glob_i, ind_glob_j, Ke.get(ind_local_i, ind_local_j));
		}
	}
    }

    void assemble_elementary_vector(
        const ElementMapping& elt_mapping,
        const ShapeFunctions& reference_functions,
        const Quadrature& quadrature,
        double (*source)(vertex),
        std::vector< double >& Fe )
    {
        //std::cout << "compute elementary vector (source term)" << '\n';
        int p = reference_functions.nb_functions();// = nombre de fonctions
        for (int i = 0; i < p; i ++){
		double sum = 0;			
		for (int q = 0; q < p; q ++){
			vertex x_r = quadrature.point(q); // point de Gauss 
			double wq = quadrature.weight(q); // poids au point de Gauss
			double val_i = reference_functions.evaluate( q, x_r);
			double det_J = elt_mapping.jacobian(x_r); //déterminant du Jacobien
			sum += wq * val_i * source(elt_mapping.transform(x_r)) * det_J;
		}
		Fe.push_back(sum);
	}
    }

    void assemble_elementary_neumann_vector(
        const ElementMapping& elt_mapping_1D,
        const ShapeFunctions& reference_functions_1D,
        const Quadrature& quadrature_1D,
        double (*neumann)(vertex),
        std::vector< double >& Fe )
    {
        std::cout << "compute elementary vector (neumann condition)" << '\n';
        int p = reference_functions_1D.nb_functions();// = nombre de fonctions
        std::cout<<p<<std::endl;
        for (int i = 0; i < p; i ++){
		double sum = 0;			
		for (int q = 0; q < p; q ++){
			vertex x_r = quadrature_1D.point(q); // point de Gauss 
			double wq = quadrature_1D.weight(q); // poids au point de Gauss
			double val_i = reference_functions_1D.evaluate( q, x_r);
			double det_J = elt_mapping_1D.jacobian(x_r); //déterminant du Jacobien
			sum += wq * val_i * neumann(elt_mapping_1D.transform(x_r)) * det_J;
		}
		Fe.push_back(sum);
        
	}
    }

    void local_to_global_vector(
        const Mesh& M,
        bool border,
        int i,
        std::vector< double >& Fe,
        std::vector< double >& F )
    {        
        if (not border){
		for (int ind_local = 0; ind_local < 3; ind_local ++){
			int ind_glob = M.get_triangle_vertex_index(i, ind_local);
			F[ind_glob]+=Fe[ind_local];
			//std::cout<<Fe[ind_local];
		}
	}
    }

    void apply_dirichlet_boundary_conditions(
        const Mesh& M,
        const std::vector< bool >& attribute_is_dirichlet, /* size: nb of attributes */
        const std::vector< double >& values, /* size: nb of DOFs */
        SparseMatrix& K,
        std::vector< double >& F )
    {
        
        double P = 10000;

        std::vector< bool > vertices_traite(values.size(), false);//contiendra des true au indices de vertices qui seront traités pour ne pas traiter deux fois un même vertice.
        
        for (int edge_ind = 0; edge_ind < M.nb_edges(); edge_ind ++){	
        	int edge_att = M.get_edge_attribute(edge_ind);
        	if (attribute_is_dirichlet[edge_att]){
        		for (int loc_ind = 0; loc_ind < 2; loc_ind ++){
	        		int vertice_ind = M.get_edge_vertex_index(edge_ind,loc_ind);
	        		// attention, comme plusieurs bords partagent un même vertice, il faut vérifier que le vertice n'a pas déjà été traité :
	        		if( not vertices_traite[vertice_ind] ) {
		                        vertices_traite[vertice_ind] = true;
        				K.add(vertice_ind,vertice_ind,P);
        				F[vertice_ind] += P * values[vertice_ind];
        			}        			
        		}
        	}
        }
    }

    void solve_poisson_problem(
            const Mesh& M,
            double (*diffusion_coef)(vertex),
            double (*source_term)(vertex),
            double (*dirichlet_fct)(vertex),
            double (*neumann_fct)(vertex),
            std::vector<double>& solution,
            int verbose )
    {
        std::cout << "solve poisson problem" << '\n';
        // TODO
    }

}
