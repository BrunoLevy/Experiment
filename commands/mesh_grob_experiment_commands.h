
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
 

#ifndef H__OGF_EXPERIMENT_COMMANDS_MESH_GROB_EXPERIMENT_COMMANDS__H
#define H__OGF_EXPERIMENT_COMMANDS_MESH_GROB_EXPERIMENT_COMMANDS__H

#include <OGF/Experiment/common/common.h>
#include <OGF/Experiment/algo/mesh_surface_intersection.h>
#include <OGF/mesh/commands/mesh_grob_commands.h>

namespace OGF {

    gom_class Experiment_API MeshGrobExperimentCommands :
        public MeshGrobCommands {
    public:
        MeshGrobExperimentCommands();
        ~MeshGrobExperimentCommands() override;

    gom_slots:
        /**
         * \brief Computes and remeshes intersection in a surface mesh.
         * \param[in] merge_vertices_and_facets pre-check duplicated
         *  vertices and facets
         * \param[in] check_neighboring_triangles check intersections
         *  also between triangles that share a vertex
         * \param[in] post_connect_facets connect the facets of the 
         *  result, and identify/disconnect non-manifold edges
         * \param[in] order_facets spatially order the facets for the AABB
         * \param[in] FPE check floating point exceptions
         * \param[in] check_constraints double-check the constraints for
         *  intersections
         * \param[in] barycentric use barycentric coordinate for local
         *  triangle geometry
         * \param[in] per_component_ids assign an operand id to each connected
         *  component
         */
        void intersect_surface(
            bool merge_vertices_and_facets = true,
            bool check_neighboring_triangles = true,
            bool post_connect_facets=true,
            bool order_facets=true,
            bool FPE=true,
            bool check_constraints=false,
            bool barycentric=false,
            bool per_component_ids=true,
            bool use_halfedges=false,
            bool verbose=false
        );

        /**
         * \brief Classifies the charts from a computed mesh intersection
         * \param[in] expr a boolean expression. 
         * \param[in] dry_run if set, just visualize the result
         */
        void classify_intersections(
            const std::string& expr = "union",
            bool dry_run=true,
            bool reorder=false
        );

        /**
         * Sort facets with AABB.
         */
        void sort_facets();
        
        /**
         * Computes a constrained Delaunay triangulation using the vertices
         *  and the edges of the mesh as constraints.
         */
        void constrained_delaunay_2d(
            bool use_my_code=true,
            bool Delaunay=true
        );


        /**
         * \brief Snap all point coordinates to 32-bit.
         */
        void floatify();
    } ;
}

#endif

