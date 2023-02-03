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
 */
 
#include <OGF/Experiment/commands/mesh_grob_experiment_commands.h>
#include <OGF/Experiment/algo/mesh_surface_intersection.h>

#include <geogram/mesh/mesh_reorder.h>
#include <geogram/mesh/mesh_repair.h>
#include <geogram/delaunay/delaunay.h>

namespace OGF {

    MeshGrobExperimentCommands::MeshGrobExperimentCommands() { 
    }
        
    MeshGrobExperimentCommands::~MeshGrobExperimentCommands() { 
    }        

    void MeshGrobExperimentCommands::intersect_surface(
        bool merge_vertices_and_facets,
        bool test_neighboring_triangles,
        bool post_connect_facets,
        bool order_facets,
        bool FPE,
        bool check_constraints,
        bool barycentric
    ) {
        MeshSurfaceIntersectionParams params;

        params.pre_detect_duplicated_vertices = merge_vertices_and_facets;
        params.pre_detect_duplicated_facets = merge_vertices_and_facets;
        params.detect_intersecting_neighbors = test_neighboring_triangles;
        params.post_connect_facets = post_connect_facets;
        params.debug_do_not_order_facets = !order_facets;
        params.debug_enable_FPE = FPE;
        params.debug_check_constraints = check_constraints;
        params.barycentric = barycentric;
        
        GEO::mesh_intersect_surface(*mesh_grob(),params);
        
        show_mesh();
        mesh_grob()->update();
    }

    void MeshGrobExperimentCommands::classify_intersections(
        const std::string& expr, const std::string& attribute
    ) {
        mesh_classify_intersections(*mesh_grob(),expr,attribute);
        if(attribute == "" || attribute == "filter") {
            Shader* shd = mesh_grob()->get_shader();
            if(shd != nullptr) {
                if(shd->has_property("facets_filter")) {
                    if(attribute == "") {
                        shd->set_property("facets_filter", "false");
                    } else {
                        shd->set_property("facets_filter", "true");                        
                    }
                }
            }
        } else {
            show_attribute("facets." + attribute);
        }
        mesh_grob()->update();
    }

    void MeshGrobExperimentCommands::sort_facets() {
        mesh_reorder(*mesh_grob(), MESH_ORDER_MORTON);
        mesh_grob()->update();
    }
    
    void MeshGrobExperimentCommands::constrained_delaunay_2d() {
        if(mesh_grob()->vertices.dimension() != 2) {
            mesh_grob()->vertices.set_dimension(2);
        }

        Delaunay_var del = Delaunay::create(2, "triangle");
        del->set_constraints(mesh_grob());
        del->set_vertices(0, nullptr);

        // For each triangle of the Delaunay triangulation,
        // create a triangle in the mesh (and flip if need be).
        for(index_t t=0; t<del->nb_cells(); ++t) {
            index_t i = index_t(del->cell_vertex(t,0));
            index_t j = index_t(del->cell_vertex(t,1));
            index_t k = index_t(del->cell_vertex(t,2));
            mesh_grob()->facets.create_triangle(i,j,k);
        }

        mesh_grob()->facets.connect();
        show_mesh();
        mesh_grob()->update();
    }

    void MeshGrobExperimentCommands::floatify() {
        index_t N =
            mesh_grob()->vertices.nb() * mesh_grob()->vertices.dimension();
        double* p = mesh_grob()->vertices.point_ptr(0);
        for(index_t i=0; i<N; ++i) {
            p[i] = double(float(p[i]));
        }
        mesh_grob()->update();
    }
}
