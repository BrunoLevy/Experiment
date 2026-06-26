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

#include <complex>

namespace OGF {

    MeshGrobCalabiYauCommands::MeshGrobCalabiYauCommands() {
    }

    MeshGrobCalabiYauCommands::~MeshGrobCalabiYauCommands() {
    }

    void MeshGrobCalabiYauCommands::compute_Delaunay_highdim(index_t nb_pts) {
#ifdef EXPERIMENT_WITH_CGAL
	index_t dim = mesh_grob()->vertices.dimension();
	Delaunay_var del = new DelaunayNd_CGAL(coord_index_t(dim));
	if(nb_pts == 0) {
	    nb_pts = mesh_grob()->vertices.nb();
	}
	Logger::out("CGAL") << "Computing Delaunay" << std::endl;
	del->set_vertices(
	    nb_pts, mesh_grob()->vertices.point_ptr(0)
	);
	Logger::out("CGAL") << "Delaunay done" << std::endl;
#else
	geo_argused(nb_pts);
	Logger::err("CGAL") << "Not compiled with CGAL support" << std::endl;
#endif
    }

    void MeshGrobCalabiYauCommands::check_equation() {
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
}
