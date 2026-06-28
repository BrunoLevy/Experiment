/*
 *  OGF/Graphite: Geometry and Graphics Programming Library + Utilities
 *  Copyright (C) 2000-2015 INRIA - Project ALICE
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact for Graphite: Bruno Levy - Bruno.Levy@inria.fr
 *  Contact for this Plugin: Bruno Levy - Bruno.Levy@inria.fr
 *
 *     Project ALICE
 *     LORIA, INRIA Lorraine,
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX
 *     FRANCE
 *
 *  Note that the GNU General Public License does not permit incorporating
 *  the Software into proprietary programs.
 *
 * As an exception to the GPL, Graphite can be linked with the following
 * (non-GPL) libraries:
 *     Qt, tetgen, SuperLU, WildMagic and CGAL
 */

#include <OGF/Experiment/commands/mesh_grob_calabi_yau_commands.h>

#ifdef EXPERIMENT_WITH_CGAL
#include <OGF/Experiment/CGAL/delaunay_nd_cgal.h>
#endif

#include <OGF/Experiment/algo/tiny_cas.h>
#include <complex>

namespace {
    using namespace OGF;
    typedef vecng<6,double> vec6;
    typedef Matrix<6,double> mat6;

    class BiQuadric {
    public:
	// The complex equation, simplest form, rational
	std::complex<double> eval0(const vec6& p) const {
	    double x1 = p[0]; double y1 = p[1]; double z1 = p[2];
	    double x2 = p[3]; double y2 = p[4]; double z2 = p[5];
	    std::complex Z1 = std::complex(x1,y1) / (1.0-z1);
	    std::complex Z2 = std::complex(x2,y2) / (1.0-z2);
	    return 1.0 + Z2*Z2 + Z1*Z1*Z2;
	}

	// The complex equation, polynomial version
	std::complex<double> eval(const vec6& p) const {
	    double x1 = p[0]; double y1 = p[1]; double z1 = p[2];
	    double x2 = p[3]; double y2 = p[4]; double z2 = p[5];
	    return
		geo_sqr(1 - z1) * geo_sqr(1 - z2)              +
		geo_sqr(1 - z1) * geo_sqr(std::complex(x2,y2)) +
		(1 - z2) * geo_sqr(std::complex(x1,y1)) * std::complex(x2,y2);
	}

	// The implicit equations of the manifold,
	// there are four of them (defines a 2-dim submanifold in IR^6)
	double F(const vec6& p, index_t i) const {
	    switch(i) {
	    case 0:
		return geo_sqr(p[0]) + geo_sqr(p[1]) + geo_sqr(p[2]) - 1.0;
	    case 1:
		return geo_sqr(p[3]) + geo_sqr(p[4]) + geo_sqr(p[5]) - 1.0;
	    case 2:
		return real(eval(p));
	    case 3:
		return imag(eval(p));
	    }
	    geo_assert_not_reached;
	}

	// The gradients of the 4 implicit equations of the manifold
	vec6 gradF(const vec6& p, index_t i) const {
	    double x1 = p[0]; double y1 = p[1]; double z1 = p[2];
	    double x2 = p[3]; double y2 = p[4]; double z2 = p[5];
	    switch(i) {
	    case 0:
		return grad_eval_real(p);
	    case 1:
		return grad_eval_imag(p);
	    case 2: {
		double g[6] = {2.0*x1, 2.0*y1, 2.0*z1, 0.0, 0.0, 0.0};
		return vec6(g);
	    }
	    case 3: {
		double g[6] = {0.0, 0.0, 0.0, 2.0*x2, 2.0*y2, 2.0*z2};
		return vec6(g);
	    }
	    }
	    geo_assert_not_reached;
	}

	// A basis of the 2D tangent space to the manifold
	std::pair<vec6, vec6> tangent_space(const vec6& p) {
	    mat6 M;
	    M.load_zero();
	    for(index_t i=0; i<4; ++i) {
		vec6 g = gradF(p,i);
		for(index_t j=0; j<6; ++j) {
		    M(i,j) = g[j];
		}
	    }

	    // minimum pivot for matrix inversion
	    static constexpr double eps = 1e-6;

	    // Compute U: orthogonal to g1,g2,g3,g4
	    vec6 U;
	    mat6 M1 = M;
	    // Loop until matrix is invertible
	    for(;;) {
		// Generate random coefficients in last two rows
		for(index_t i=4; i<6; ++i) {
		    for(index_t j=0; j<6; ++j) {
			M1(i,j) = 2.0*(Numeric::random_float64() - 0.5);
		    }
		}
		mat6 M1_inv;
		if(M1.compute_inverse(M1_inv, eps)) {
		    U = M1_inv * vec6({0.0, 0.0, 0.0, 0.0, 1.0, 1.0});
		    break;
		}
	    }

	    // Compute V: orthogonal to g1,g2,g3,g4,U
	    vec6 V;
	    mat6 M2 = M;
	    for(index_t j=0; j<6; ++j) {
		M2(4,j) = U[j];
	    }
	    // Loop until matrix is invertible
	    for(;;) {
		// Generate random coefficients in last row
		for(index_t j=0; j<6; ++j) {
		    M2(5,j) = 2.0*(Numeric::random_float64() - 0.5);
		}
		mat6 M2_inv;
		if(M2.compute_inverse(M2_inv, eps)) {
		    V = M2_inv * vec6({0.0, 0.0, 0.0, 0.0, 0.0, 1.0});
		    break;
		}
	    }

	    return std::make_pair(normalize(U),normalize(V));
	}

    protected:
	// gradient of the real part of the polynomial complex eqn
	// Generated using TinyCAS
	vec6 grad_eval_real(const vec6& p) const {
	    double x1 = p[0]; double y1 = p[1]; double z1 = p[2];
	    double x2 = p[3]; double y2 = p[4]; double z2 = p[5];
	    return vec6({
		+2*x1*x2 -2*x1*x2*z2 -2*y1*y2 +2*y1*y2*z2,
		-2*y1*x2 +2*y1*x2*z2 -2*x1*y2 +2*x1*y2*z2,
		-2 +4*z2 -2*z2*z2 +2*z1 -4*z1*z2 +2*z1*z2*z2
		                   -2*x2*x2 +2*y2*y2 +2*z1*x2*x2 -2*z1*y2*y2,
		+2*x2 -4*z1*x2 +2*z1*z1*x2 +x1*x1 -y1*y1 -x1*x1*z2 +y1*y1*z2,
		-2*y2 +4*z1*y2 -2*z1*z1*y2 -2*x1*y1 +2*x1*y1*z2,
		-2 +2*z2 +4*z1 -4*z1*z2 -2*z1*z1 +2*z1*z1*z2
		                   -x1*x1*x2 +y1*y1*x2 +2*x1*y1*y2
	    });
	}

	// gradient of the imaginary part of the polynomial complex eqn
	// Generated using TinyCAS
	vec6 grad_eval_imag(const vec6& p) const {
	    double x1 = p[0]; double y1 = p[1]; double z1 = p[2];
	    double x2 = p[3]; double y2 = p[4]; double z2 = p[5];
	    return vec6({
		+2*x1*y2 -2*x1*y2*z2 +2*y1*x2 -2*y1*x2*z2,
		-2*y1*y2 +2*y1*y2*z2 +2*x1*x2 -2*x1*x2*z2,
		-4*x2*y2 +4*z1*x2*y2,
		-4*z1*y2 +2*z1*z1*y2 +2*x1*y1 -2*x1*y1*z2,
		+2*x2 -4*z1*x2 +2*z1*z1*x2 +x1*x1 -y1*y1 -x1*x1*z2 +y1*y1*z2,
		-x1*x1*y2 +y1*y1*y2 -2*x1*y1*x2
	    });
	}

    };
}


namespace OGF {

    MeshGrobCalabiYauCommands::MeshGrobCalabiYauCommands() {
    }

    MeshGrobCalabiYauCommands::~MeshGrobCalabiYauCommands() {
    }

    void MeshGrobCalabiYauCommands::remove_points(index_t nb_points_to_keep) {
	vector<index_t> remove_points(mesh_grob()->vertices.nb());
	for(index_t v: mesh_grob()->vertices) {
	    remove_points[v] = index_t(v >= nb_points_to_keep);
	}
	mesh_grob()->vertices.delete_elements(remove_points);
	mesh_grob()->update();
    }

    void MeshGrobCalabiYauCommands::show_tangent_basis(double size) {
	BiQuadric Q;
	index_t nb_v = mesh_grob()->vertices.nb();
	for(index_t v=0; v<nb_v; ++v) {
	    vec6 p = mesh_grob()->vertices.point<6>(v);
	    auto [U,V] = Q.tangent_space(p);
	    index_t first_v = mesh_grob()->vertices.create_vertices(4);
	    index_t v1 = first_v;
	    index_t v2 = first_v+1;
	    index_t v3 = first_v+2;
	    index_t v4 = first_v+3;
	    mesh_grob()->vertices.point<6>(v1) = p + size * U;
	    mesh_grob()->vertices.point<6>(v2) = p - size * U;
	    mesh_grob()->vertices.point<6>(v3) = p + size * V;
	    mesh_grob()->vertices.point<6>(v4) = p - size * V;
	    mesh_grob()->edges.create_edge(v,v1);
	    mesh_grob()->edges.create_edge(v,v2);
	    mesh_grob()->edges.create_edge(v,v3);
	    mesh_grob()->edges.create_edge(v,v4);
	}
	mesh_grob()->update();
    }

    void MeshGrobCalabiYauCommands::compute_Delaunay_highdim(
	index_t nb_pts, bool build_graph
    ) {
#ifdef EXPERIMENT_WITH_CGAL
	index_t dim = mesh_grob()->vertices.dimension();
	Delaunay_var del = new DelaunayNd_CGAL(coord_index_t(dim));
	del->set_stores_cicl(build_graph);
	if(nb_pts == 0) {
	    nb_pts = mesh_grob()->vertices.nb();
	}
	Logger::out("CGAL") << "Computing Delaunay" << std::endl;
	del->set_vertices(
	    nb_pts, mesh_grob()->vertices.point_ptr(0)
	);
	Logger::out("CGAL") << "Delaunay done" << std::endl;
	if(build_graph) {
	    vector<index_t> N;
	    double avg_arity = 0.0;
	    for(index_t v=0; v<del->nb_vertices(); ++v) {
		del->get_neighbors(v,N);
		avg_arity += double(N.size());
	    }
	    avg_arity /= double(del->nb_vertices());
	    Logger::out("CGAL") << "Average arity:" << avg_arity << std::endl;
	}
#else
	geo_argused(nb_pts);
	geo_argused(build_graph);
	Logger::err("CGAL") << "Not compiled with CGAL support" << std::endl;
#endif
    }


    void MeshGrobCalabiYauCommands::equation_check() {
	BiQuadric Q;

	for(index_t v: mesh_grob()->vertices) {
	    const vec6& p = mesh_grob()->vertices.point<6>(v);
	    std::complex f1 = Q.eval0(p);
	    std::complex f2 = Q.eval(p);
	    std::cerr << f1 << " " << f2 << std::endl;
	}
    }

    // Z1 = (x1 + iy1) / (1-z1)
    // Z2 = (x2 + iy2) / (1-z2)
    //
    // F1(x1,y1,z1,x2,y2,z2) = || [x1,y1,z1] ||^2 - 1 -> GF3 = 2 [x1,y1,z1,0,0,0]
    // F2(x1,y1,z1,x2,y2,z2) = || [x2,y2,z2] ||^2 - 1 -> GF4 = 2 [0,0,0,x2,y2,z2]
    // P = 1 + Z2^2 + Z1^2*Z2
    // consider Q = (1-z1)^2(1-z2)^2 P (same denominator for the terms of P)
    // Q = (1-z1)^2(1-z2)^2 + (1-z1)^2(x2+iy2)^2 + (1-z2)(x1+iy1)^2(x2+iy2)
    // F3(x1,y1,z1,x2,y2,z2) = Re(Q(Z1,Z2))
    // F4(x1,y1,z1,x2,y2,z2) = Im(Q(Z1,Z2))

    // TODO
    // DONE: Check equation of Q on data
    // DONE: Functions to compute derivatives of F1,F2,F3,F4
    // DONE: Function to compute tangent space basis
    // DONE: Visualize tangent space vectors
    // Sutherland-Hogdman re-entrant clipping in 6D

    void MeshGrobCalabiYauCommands::equation_gradients() {
	using ::OGF::TinyCAS::Polynom;

	typedef Polynom R;
	typedef std::complex<R> C;

	std::vector<std::string> var_names = {
	    "x1","y1","z1",
	    "x2","y2","z2"
	};

	index_t nb_var = index_t(var_names.size());
	R Rx1(nb_var,0u);
	R Ry1(nb_var,1u);
	R Rz1(nb_var,2u);
	R Rx2(nb_var,3u);
	R Ry2(nb_var,4u);
	R Rz2(nb_var,5u);

	R Rzero(nb_var,0.0);
	R Rone(nb_var,1.0);

	C one(Rone, Rzero);
	C i(Rzero, Rone);
	C x1(Rx1, Rzero);
	C y1(Ry1, Rzero);
	C z1(Rz1, Rzero);
	C x2(Rx2, Rzero);
	C y2(Ry2, Rzero);
	C z2(Rz2, Rzero);

	auto sq = [](const C& Z)->C {
	    return Z*Z;
	};

	C Q =
	    sq(one - z1) * sq(one - z2)              +
	    sq(one - z1) * sq(x2 + i*y2)             +
	      (one - z2) * sq(x1 + i*y1) * (x2+i*y2) ;

	std::cerr
	    << "F3=" << real(Q).to_string(var_names) << std::endl
	    << "F4=" << imag(Q).to_string(var_names) << std::endl;

	for(index_t i=0; i<6; ++i) {
	    std::cerr << "der F3_" << var_names[i] << "="
		      << real(Q).derivative(i).to_string(var_names)
		      << std::endl;
	}

	for(index_t i=0; i<6; ++i) {
	    std::cerr << "der F4_" << var_names[i] << "="
		      << imag(Q).derivative(i).to_string(var_names)
		      << std::endl;
	}


    }

}
