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
#include <OGF/Experiment/algo/CDT.h>

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
        bool barycentric,
        bool per_component_ids,
        bool use_halfedges,
        bool verbose
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
        params.per_component_ids = per_component_ids;
        params.use_halfedges = use_halfedges;
        params.verbose = verbose;
        
        GEO::mesh_intersect_surface(*mesh_grob(),params);
        
        show_mesh();
        mesh_grob()->update();
    }

    void MeshGrobExperimentCommands::classify_intersections(
        const std::string& expr, bool dry_run, bool reorder
    ) {
        std::string attribute = dry_run ? "filter" : "";
        mesh_classify_intersections(*mesh_grob(),expr,attribute,reorder);
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
    
    void MeshGrobExperimentCommands::constrained_delaunay_2d(
        bool use_my_code,
        bool Delaunay,
        bool constrained,
        bool quad,
        bool remove_external_triangles
    ) {
        if(mesh_grob()->vertices.dimension() != 2) {
            mesh_grob()->vertices.set_dimension(2);
        }

        if(use_my_code) {
            CDT cdt;
            cdt.set_delaunay(Delaunay);

            index_t n=0;
            if(quad) {
                n=4;
                vec2 p0(mesh_grob()->vertices.point_ptr(0));
                vec2 p1(mesh_grob()->vertices.point_ptr(1));
                vec2 p2(mesh_grob()->vertices.point_ptr(2));
                vec2 p3(mesh_grob()->vertices.point_ptr(3));                
                cdt.create_enclosing_quad(p0,p1,p2,p3);
            } else {
                n=3;
                vec2 p0(mesh_grob()->vertices.point_ptr(0));
                vec2 p1(mesh_grob()->vertices.point_ptr(1));
                vec2 p2(mesh_grob()->vertices.point_ptr(2));
                cdt.create_enclosing_triangle(p0,p1,p2);
            }
            
            cdt.insert(
                mesh_grob()->vertices.nb()-n,
                mesh_grob()->vertices.point_ptr(n)
            );
            
            if(constrained) {
                for(index_t e: mesh_grob()->edges) {
                    cdt.insert_constraint(
                        mesh_grob()->edges.vertex(e,0),
                        mesh_grob()->edges.vertex(e,1)                    
                    );
                }
                // Create the vertices coming from constraint intersections
                for(index_t v=mesh_grob()->vertices.nb(); v<cdt.nv(); ++v) {
                    mesh_grob()->vertices.create_vertex(
                        cdt.point(v).data()
                    );
                }
            }

            // Check whether all edges are Delaunay
#ifdef GEO_DEBUG            
            {
                if(Delaunay) {
                    for(index_t t=0; t<cdt.nT(); ++t) {
                        for(index_t le=0; le<3; ++le) {
                            if(!cdt.Tedge_is_Delaunay(t,le)) {
                                Logger::err("CDT")
                                    << "Triangle " << t << " edge " << le
                                    << " is not Delaunay!"
                                    << std::endl;
                            }
                        }
                    }
                }
                cdt.save("CDT_result.geogram");
            }
#endif
            
            if(remove_external_triangles) {
                cdt.remove_external_triangles();
            }
            
            for(index_t t=0; t<cdt.nT(); ++t) {
                index_t i = cdt.Tv(t,0);
                index_t j = cdt.Tv(t,1);
                index_t k = cdt.Tv(t,2);
                mesh_grob()->facets.create_triangle(i,j,k);
            }


            
        } else {
            Delaunay_var del = Delaunay::create(2, "triangle");
            del->set_constraints(mesh_grob());
            del->set_vertices(0, nullptr);
            for(index_t t=0; t<del->nb_cells(); ++t) {
                index_t i = index_t(del->cell_vertex(t,0));
                index_t j = index_t(del->cell_vertex(t,1));
                index_t k = index_t(del->cell_vertex(t,2));
                mesh_grob()->facets.create_triangle(i,j,k);
            }
        }
        mesh_grob()->facets.connect();
        mesh_grob()->vertices.set_dimension(3);
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
