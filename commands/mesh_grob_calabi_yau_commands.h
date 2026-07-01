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

#ifndef H__OGF_EXPERIMENT_COMMANDS_MESH_GROB_CALABI_YAU_COMMANDS__H
#define H__OGF_EXPERIMENT_COMMANDS_MESH_GROB_CALABI_YAU_COMMANDS__H

#include <OGF/Experiment/common/common.h>
#include <OGF/mesh/commands/mesh_grob_commands.h>

namespace OGF {

    gom_class Experiment_API MeshGrobCalabiYauCommands :
        public MeshGrobCommands {
    public:
        MeshGrobCalabiYauCommands();
        ~MeshGrobCalabiYauCommands() override;

    gom_slots:
	void remove_points(index_t nb_points_to_keep);
	void show_tangent_basis(double size = 1e-3);
	void tangent_RVD(
	    double r=1.0, bool three_only = false, index_t nb_CVT=0
	);
	void compute_Delaunay_highdim(index_t nb_pts=0, bool build_graph=false);
	void equation_check();
	void equation_gradients();
	void export_tet6(const std::string& filename);
	void project_onto_manifold();
    };
}

#endif
