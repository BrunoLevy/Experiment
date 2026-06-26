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

namespace OGF {

    MeshGrobCalabiYauCommands::MeshGrobCalabiYauCommands() {
    }

    MeshGrobCalabiYauCommands::~MeshGrobCalabiYauCommands() {
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
	for(index_t v: mesh_grob()->vertices) {
	    const double* p =  mesh_grob()->vertices.point_ptr(v);
	    double x1 = p[0];
	    double y1 = p[1];
	    double z1 = p[2];
	    double x2 = p[3];
	    double y2 = p[4];
	    double z2 = p[5];

            //<- /!\ it is (1+z) and not (1-z) in the denominator
	    std::complex Z1 = std::complex(x1,y1) / (1.0-z1);
	    std::complex Z2 = std::complex(x2,y2) / (1.0-z2);
	    std::complex P = 1.0 + Z2*Z2 + Z1*Z1*Z2;
	    std::cerr << P << std::endl;
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
	    << real(Q).to_string(var_names) << std::endl
	    << imag(Q).to_string(var_names) << std::endl;

    }

}
