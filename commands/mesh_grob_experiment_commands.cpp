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

#include <geogram/mesh/mesh_surface_intersection.h>
#include <geogram/mesh/mesh_reorder.h>
#include <geogram/mesh/mesh_repair.h>
#include <geogram/mesh/triangle_intersection.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/delaunay/delaunay.h>
#include <geogram/delaunay/CDT_2d.h>
#include <geogram/numerics/exact_geometry.h>

#include <OGF/Experiment/algo/CDT_with_interval.h>

#ifdef GEOGRAM_WITH_VORPALINE
#include <vorpalib/mesh/mesh_weiler_model.h>
#endif

namespace OGF {

    MeshGrobExperimentCommands::MeshGrobExperimentCommands() { 
    }
        
    MeshGrobExperimentCommands::~MeshGrobExperimentCommands() { 
    }        

    void MeshGrobExperimentCommands::intersect_surface(
        bool FPE,
        bool remove_external_shell,
        bool remove_internal_shells,
        bool simplify_coplanar_facets,
        bool detect_intersecting_neighbors,
        bool interpolate_attributes,
        bool delaunay,
        bool verbose,
        bool post_process,
        const std::string& expr,
        const NewMeshGrobName& skeleton
    ) {
        bool FPE_bkp = Process::FPE_enabled();
        Process::enable_FPE(FPE);

        MeshSurfaceIntersection intersection(*mesh_grob());
        intersection.set_delaunay(delaunay);
        intersection.set_detect_intersecting_neighbors(
            detect_intersecting_neighbors
        );
        intersection.set_verbose(verbose);
        intersection.set_radial_sort(
            remove_external_shell    ||
            remove_internal_shells   ||
            simplify_coplanar_facets ||
            expr != ""
        );
        intersection.set_interpolate_attributes(interpolate_attributes);
        if(skeleton != "") {
            intersection.set_build_skeleton(
                MeshGrob::find_or_create(scene_graph(),skeleton)
            );
        }
        intersection.intersect();
        Process::enable_FPE(FPE_bkp);

        if(expr != "") {
            intersection.classify(expr);
        } else {
            if(remove_external_shell) {
                intersection.remove_external_shell();
            }
            if(remove_internal_shells) {
                intersection.remove_internal_shells();
            }
        }

        if(simplify_coplanar_facets && !interpolate_attributes) {
            intersection.simplify_coplanar_facets();
        }
        
        // Still need to do that, because snap-rounding may have created
        // degeneracies
        if(remove_internal_shells || post_process) {
            mesh_repair(*mesh_grob());
        }
        
        show_mesh();
        mesh_grob()->update();
    }

    void MeshGrobExperimentCommands::build_Weiler_model(
        double expand_surfaces 
    ) {
#ifdef GEOGRAM_WITH_VORPALINE        
        WeilerModel weiler(*mesh_grob());
        weiler.set_verbose(true);
        weiler.set_delaunay(true);
        weiler.expand_surfaces(expand_surfaces);
        weiler.build();

        index_t nb_regions = weiler.nb_regions();
        Logger::out("Weiler") << "Found " << nb_regions
                              << " regions" << std::endl;
        for(index_t i=0; i<nb_regions; ++i) {
            MeshGrob* region = MeshGrob::find_or_create(
                scene_graph(), "region_" + String::to_string(i)
            );
            weiler.copy_region(i,*region);
            region->update();
        }
        
        show_mesh();
        mesh_grob()->update();
#else
        Logger::err("Weiler") << "Needs to be compiled with Vorpaline"
                              << std::endl;
#endif        
    }


    void MeshGrobExperimentCommands::sort_facets() {
        mesh_reorder(*mesh_grob(), MESH_ORDER_MORTON);
        mesh_grob()->update();
    }
    
    void MeshGrobExperimentCommands::constrained_delaunay_2d(
        bool use_my_code,
        bool use_intervals,
        bool Delaunay,
        bool constrained,
        bool quad,
        bool remove_external_triangles,
        bool FPE
    ) {
        bool FPE_bkp = Process::FPE_enabled();
        Process::enable_FPE(FPE);
        
        if(mesh_grob()->vertices.dimension() != 2) {
            mesh_grob()->vertices.set_dimension(2);
        }

        if(use_my_code) {
            CDT2d* cdt = (use_intervals) ? new CDT2d_with_interval : new CDT2d;
            cdt->set_delaunay(Delaunay);

            index_t n=0;
            if(quad) {
                n=4;
                vec2 p0(mesh_grob()->vertices.point_ptr(0));
                vec2 p1(mesh_grob()->vertices.point_ptr(1));
                vec2 p2(mesh_grob()->vertices.point_ptr(2));
                vec2 p3(mesh_grob()->vertices.point_ptr(3));                
                cdt->create_enclosing_quad(p0,p1,p2,p3);
            } else {
                n=3;
                vec2 p0(mesh_grob()->vertices.point_ptr(0));
                vec2 p1(mesh_grob()->vertices.point_ptr(1));
                vec2 p2(mesh_grob()->vertices.point_ptr(2));
                cdt->create_enclosing_triangle(p0,p1,p2);
            }

            vector<index_t> indices(mesh_grob()->vertices.nb()-n);
            cdt->insert(
                mesh_grob()->vertices.nb()-n,
                mesh_grob()->vertices.point_ptr(n),
                indices.data()
            );

            if(constrained) {
                for(index_t e: mesh_grob()->edges) {
                    index_t v1=mesh_grob()->edges.vertex(e,0);
                    index_t v2=mesh_grob()->edges.vertex(e,1);
                    if(v1 >= n) {
                        v1 = indices[v1-n];
                    }
                    if(v2 >= n) {
                        v2 = indices[v2-n];
                    }
                    cdt->insert_constraint(v1,v2);
                }
            }

            if(remove_external_triangles) {
                cdt->remove_external_triangles();
            }

            // Create vertices coming from constraint intersections
            for(index_t v = mesh_grob()->vertices.nb(); v < cdt->nv(); ++v) {
                mesh_grob()->vertices.create_vertex(cdt->point(v).data());
            }

            // Create triangles
            for(index_t t=0; t<cdt->nT(); ++t) {
                mesh_grob()->facets.create_triangle(
                    cdt->Tv(t,0), cdt->Tv(t,1), cdt->Tv(t,2)
                );
            }

            delete cdt;
            
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

        Process::enable_FPE(FPE_bkp);
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

    void MeshGrobExperimentCommands::show_triangle_triangle_intersections() {
        if(mesh_grob()->facets.nb() < 2) {
            Logger::out("TT") << "Need at least two triangles"
                              << std::endl;
            return;
        }
        if(!mesh_grob()->facets.are_simplices()) {
            Logger::out("TT") << "Mesh is not triangulated"
                              << std::endl;
            return;
        }
        for(index_t t1: mesh_grob()->facets) {
            index_t v1 = mesh_grob()->facets.vertex(t1,0);
            index_t v2 = mesh_grob()->facets.vertex(t1,1);
            index_t v3 = mesh_grob()->facets.vertex(t1,2);
            vec3 p1(mesh_grob()->vertices.point_ptr(v1));
            vec3 p2(mesh_grob()->vertices.point_ptr(v2));
            vec3 p3(mesh_grob()->vertices.point_ptr(v3));
            
            for(index_t t2: mesh_grob()->facets) {
                if(t1 == t2) {
                    continue;
                }

                index_t w1 = mesh_grob()->facets.vertex(t2,0);
                index_t w2 = mesh_grob()->facets.vertex(t2,1);
                index_t w3 = mesh_grob()->facets.vertex(t2,2);
                vec3 q1(mesh_grob()->vertices.point_ptr(w1));
                vec3 q2(mesh_grob()->vertices.point_ptr(w2));
                vec3 q3(mesh_grob()->vertices.point_ptr(w3));

                
                vector<TriangleIsect> II;
                
                triangles_intersections(
                    p1, p2, p3,
                    q1, q2, q3,
                    II
                );

                Logger::out("II") << t1 << " /\\ " << t2
                                  << "   " << II << std::endl;
            }
        }
    }

    void MeshGrobExperimentCommands::inflate(double howmuch) {
        compute_normals(*mesh_grob());
        vector<vec3> N(mesh_grob()->vertices.nb());
        for(index_t v: mesh_grob()->vertices) {
            N[v] = howmuch * normalize(
                Geom::mesh_vertex_normal(*mesh_grob(), v)
            );
        }
        for(index_t v: mesh_grob()->vertices) {
            double* p = mesh_grob()->vertices.point_ptr(v);
            p[0] += N[v].x;
            p[1] += N[v].y;
            p[2] += N[v].z;
        }
        mesh_grob()->update();
    }
    
}
