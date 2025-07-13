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
#include <OGF/scene_graph/types/scene_graph.h>

#include <geogram/mesh/mesh_surface_intersection.h>
#include <geogram/mesh/mesh_reorder.h>
#include <geogram/mesh/mesh_repair.h>
#include <geogram/mesh/triangle_intersection.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/delaunay/delaunay.h>
#include <geogram/delaunay/CDT_2d.h>
#include <geogram/delaunay/delaunay_2d.h>
#include <geogram/numerics/exact_geometry.h>
#include <geogram/basic/geometry.h>
#include <geogram/basic/stopwatch.h>
#include <geogram/mesh/mesh_AABB.h>

#include <OGF/Experiment/algo/CDT_with_interval.h>

#ifdef GEOGRAM_WITH_VORPALINE
#include <vorpalib/mesh/mesh_weiler_model.h>
#include <vorpalib/mesh/mesh_geomodel.h>
#endif

namespace {
    using namespace OGF;

    typedef vector<vec2> Polygon;

    static inline Sign point_is_in_half_plane(
        const vec2& p, const vec2& q1, const vec2& q2
    ) {
        return PCK::orient_2d(q1, q2, p);
    }

    static inline bool intersect_segments(
        const vec2& p1, const vec2& p2,
        const vec2& q1, const vec2& q2,
        vec2& result
    ) {

        vec2 Vp = p2 - p1;
        vec2 Vq = q2 - q1;
        vec2 pq = q1 - p1;

        double a =  Vp.x;
        double b = -Vq.x;
        double c =  Vp.y;
        double d = -Vq.y;

        double delta = a*d-b*c;
        if(delta == 0.0) {
            return false ;
        }

        double tp = (d * pq.x -b * pq.y) / delta;

        result = vec2(
            (1.0 - tp) * p1.x + tp * p2.x,
            (1.0 - tp) * p1.y + tp * p2.y
        );

        return true;
    }

    void clip_polygon_by_half_plane(
        const Polygon& P,
        const vec2& q1,
        const vec2& q2,
        Polygon& result
    ) {
        result.clear() ;

        if(P.size() == 0) {
            return ;
        }

        if(P.size() == 1) {
            if(point_is_in_half_plane(P[0], q1, q2)) {
                result.push_back(P[0]) ;
            }
            return ;
        }

        vec2 prev_p = P[P.size() - 1] ;
        Sign prev_status = point_is_in_half_plane(
            prev_p, q1, q2
        );

        for(unsigned int i=0; i<P.size(); i++) {
            vec2 p = P[i] ;
            Sign status = point_is_in_half_plane(
                p, q1, q2
            );
            if(
                status != prev_status &&
                status != ZERO &&
                prev_status != ZERO
            ) {
                vec2 intersect ;
                if(intersect_segments(prev_p, p, q1, q2, intersect)) {
                    result.push_back(intersect) ;
                }
            }

            switch(status) {
            case NEGATIVE:
                break ;
            case ZERO:
                result.push_back(p) ;
                break ;
            case POSITIVE:
                result.push_back(p) ;
                break ;
            }

            prev_p = p ;
            prev_status = status ;
        }
    }

    void convex_clip_polygon(
        const Polygon& P, const Polygon& clip, Polygon& result
    ) {
        Polygon tmp1 = P ;
        Polygon tmp2 ;
        Polygon* src = &tmp1 ;
        Polygon* dst = &tmp2 ;
        for(unsigned int i=0; i<clip.size(); i++) {
            unsigned int j = ((i+1) % clip.size()) ;
            const vec2& p1 = clip[i] ;
            const vec2& p2 = clip[j] ;
            clip_polygon_by_half_plane(*src, p1, p2, *dst);
            std::swap(src, dst) ;
        }
        result = *src ;
    }

    class VoronoiDiagram2d {
    public:
        VoronoiDiagram2d(
            Delaunay* delaunay, const Polygon& border
        ) : delaunay_(delaunay), border_(border) {
        }


        /**
         * \brief Gets a Voronoi cell of a vertex
         * \details The vertex is specified by a triangle and a local index in
         *  the triangle
         * \param[in] t0 the triangle
         * \param[in] lv the local index of the vertex in triangle \p t0
         * \param[out] cell a reference to the Voronoi cell
         */
        void get_Voronoi_cell(index_t t0, index_t lv, Polygon& cell) {
            cell.resize(0);
            index_t v = index_t(delaunay_->cell_to_v()[3*t0+lv]);
            bool on_border = false;
            index_t t = t0;

            // First, we turn around the vertex v. To do that, we compute
            // lv, the local index of v in the current triangle. Following
            // the standard numerotation of a triangle, edge number lv is
            // not incident to vertex v. The two other edges (lv+1)%3 and
            // (lv+2)%3 of the triangle are indicent to vertex v. By always
            // traversing (lv+1)%3, we turn around vertex v.
            do {
                index_t e = (lv+1)%3;
                index_t neigh_t = delaunay_->cell_to_cell()[3*t+e];
                if(neigh_t == NO_INDEX) {
                    on_border = true;
                    break;
                }
                cell.push_back(circumcenter(t));
                t = neigh_t;
                lv = find_vertex(t,v);
            } while(t != t0);


            // If one traversed edge is on the border of the convex hull, then
            // we empty the cell, and start turing around the vertex in
            // the other direction, i.e. by traversing this time
            // edge (lv+2)%3 until we reach the other edge on the border of
            // the convex hull that is incident to v.
            if(on_border) {
                cell.resize(0);
                cell.push_back(infinite_vertex(t,(lv + 1)%3));
                for(;;) {
                    cell.push_back(circumcenter(t));
                    index_t e = (lv+2)%3;
                    index_t neigh_t = delaunay_->cell_to_cell()[3*t+e];
                    if(neigh_t == NO_INDEX) {
                        cell.push_back(infinite_vertex(t, e));
                        break;
                    }
                    t = neigh_t;
                    lv = find_vertex(t,v);
                }
            }

            Polygon clipped;
            convex_clip_polygon(cell, border_, clipped);
            cell.swap(clipped);
        }

        /**
         * \brief Gets the circumcenter of a triangle.
         * \param[in] t the index of the triangle, in 0..delaunay->nb_cells()-1
         * \return the circumcenter of triangle \p t
         */
        vec2 circumcenter(index_t t) {
            index_t v1 = delaunay_->cell_to_v()[3*t];
            index_t v2 = delaunay_->cell_to_v()[3*t+1];
            index_t v3 = delaunay_->cell_to_v()[3*t+2];
            vec2 p1(delaunay_->vertex_ptr(v1));
            vec2 p2(delaunay_->vertex_ptr(v2));
            vec2 p3(delaunay_->vertex_ptr(v3));
            return GEO::Geom::triangle_circumcenter(p1,p2,p3);
        }

        /**
         * \brief Gets an infinite vertex in the direction normal to an
         *  edge on the boundary of the convex hull.
         * \param[in] t the index of the triangle, in 0..delaunay->nb_cells()-1
         * \param[in] e the local index of the edge, in {0,1,2}
         * \return a point located far away along the direction normal to the
         *  edge \p e of triangle \p t
         */
        vec2 infinite_vertex(index_t t, index_t e) {
            index_t lv1 = (e+1)%3;
            index_t lv2 = (e+2)%3;
            index_t v1 = index_t(delaunay_->cell_to_v()[3*t+lv1]);
            index_t v2 = index_t(delaunay_->cell_to_v()[3*t+lv2]);
            vec2 p1(delaunay_->vertex_ptr(v1));
            vec2 p2(delaunay_->vertex_ptr(v2));
            vec2 n = normalize(p2-p1);
            n = vec2(n.y, -n.x);
            return 0.5*(p1+p2)+100000.0*n;
        }

        /**
         * \brief Finds the local index of a vertex in a triangle.
         * \details Throws an assertion failure if the triangle \p t is
         *  not incident to vertex \p v
         * \param[in] t the triangle, in 0..delaunay->nb_cells()-1
         * \param[in] v the vertex, in 0..delaunay->nb_vertices()-1
         * \return the local index of v in t, in {0,1,2}
         */
        index_t find_vertex(index_t t, index_t v) {
            for(index_t lv=0; lv<3; ++lv) {
                if(index_t(delaunay_->cell_to_v()[3*t+lv]) == v) {
                    return lv;
                }
            }
            geo_assert_not_reached;
        }

    private:
        Delaunay* delaunay_;
        Polygon border_;
    };
}

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
        const NewMeshGrobName& skeleton,
        bool skeleton_trim_fins
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
            expr != ""               ||
            skeleton != ""
        );
        intersection.set_interpolate_attributes(interpolate_attributes);
        if(skeleton != "") {
            intersection.set_build_skeleton(
                MeshGrob::find_or_create(scene_graph(),skeleton),
                skeleton_trim_fins
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

    void MeshGrobExperimentCommands::add_2d_box(double expansion) {
           double xyzmin[3];
           double xyzmax[3];
           mesh_grob()->vertices.set_dimension(3);
           get_bbox(*mesh_grob(), xyzmin, xyzmax);
           mesh_grob()->vertices.set_dimension(2);

           double w = xyzmax[0] - xyzmin[0];
           double h = xyzmax[1] - xyzmin[1];
           double x1 = xyzmin[0] - expansion*w;
           double y1 = xyzmin[1] - expansion*h;
           double x2 = xyzmax[0] + expansion*w;
           double y2 = xyzmax[1] + expansion*h;

           mesh_grob()->vertices.create_vertex(vec2(x1,y1).data());
           mesh_grob()->vertices.create_vertex(vec2(x2,y1).data());
           mesh_grob()->vertices.create_vertex(vec2(x2,y2).data());
           mesh_grob()->vertices.create_vertex(vec2(x1,y2).data());

           index_t N = mesh_grob()->vertices.nb();
           vector<index_t> reorder;
           reorder.reserve(N);
           reorder.push_back(N-4);
           reorder.push_back(N-3);
           reorder.push_back(N-2);
           reorder.push_back(N-1);
           for(index_t i=0; i<N-4; ++i) {
               reorder.push_back(i);
           }

           mesh_grob()->vertices.permute_elements(reorder);

           mesh_grob()->update();
    }


    void MeshGrobExperimentCommands::compute_Voronoi_diagram(
        double expansion,
        const NewMeshGrobName& voronoi_mesh_name
    ) {

        MeshGrob* voronoi_mesh = MeshGrob::find_or_create(
            scene_graph(), voronoi_mesh_name
        );
        voronoi_mesh->clear();
        voronoi_mesh->vertices.set_dimension(2);


        double xyzmin[3];
        double xyzmax[3];
        mesh_grob()->vertices.set_dimension(3);
        get_bbox(*mesh_grob(), xyzmin, xyzmax);
        mesh_grob()->vertices.set_dimension(2);

        double w = xyzmax[0] - xyzmin[0];
        double h = xyzmax[1] - xyzmin[1];
        double x1 = xyzmin[0] - expansion*w;
        double y1 = xyzmin[1] - expansion*h;
        double x2 = xyzmax[0] + expansion*w;
        double y2 = xyzmax[1] + expansion*h;


        Delaunay_var delaunay = new Delaunay2d();
        delaunay->set_vertices(
            mesh_grob()->vertices.nb(),
            mesh_grob()->vertices.point_ptr(0)
        );

        Polygon square;
        square.push_back(vec2(x1,y1));
        square.push_back(vec2(x2,y1));
        square.push_back(vec2(x2,y2));
        square.push_back(vec2(x1,y2));

        VoronoiDiagram2d voronoi(delaunay,square);

        vector<bool> v_visited(mesh_grob()->vertices.nb(),false);
        for(index_t t=0; t<delaunay->nb_cells(); ++t) {
            for(index_t lv=0; lv<3; ++lv) {
                index_t v = index_t(delaunay->cell_to_v()[3*t+lv]);
                if(!v_visited[v]) {
                    Polygon C;
                    v_visited[v] = true;
                    voronoi.get_Voronoi_cell(t,lv,C);
                    index_t ofs = voronoi_mesh->vertices.nb();
                    for(const vec2& p: C) {
                        voronoi_mesh->vertices.create_vertex(p.data());
                    }
                    index_t f = voronoi_mesh->facets.create_polygon(C.size());
                    for(index_t i=0; i<C.size(); ++i) {
                        voronoi_mesh->facets.set_vertex(f, i, i+ofs);
                    }
                }
            }
        }

        voronoi_mesh->vertices.set_dimension(3);
        mesh_grob()->vertices.set_dimension(3);

        voronoi_mesh->update();
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
                index_t i = cdt->Tv(t,0);
                index_t j = cdt->Tv(t,1);
                index_t k = cdt->Tv(t,2);
                mesh_grob()->facets.create_triangle(i,j,k);
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

    void MeshGrobExperimentCommands::create_cube() {

	std::cerr << "create vertices" << std::endl;

	Mesh& mesh = *mesh_grob();

	std::cerr << mesh.vertices.dimension() << std::endl;

	// mesh.vertices.set_dimension(3);

	index_t cell_nu = 1464;
	index_t cell_nv = 744;
	index_t cell_nw = 376;
	index_t node_nu = cell_nu + 1;
	index_t node_nv = cell_nv + 1;
	index_t node_nw = cell_nw + 1;
	mesh.vertices.create_vertices(node_nu* node_nv* node_nw);

	std::cerr << "init vertices" << std::endl;

	FOR(k, node_nw) {
	    FOR(j, node_nv) {
		FOR (i, node_nu) {
		    mesh.vertices.point(node_nu* node_nv*k + j* node_nu + i) =
			vec3(
			    double(i), //  / double(node_nu),
			    double(j), //  / double(node_nv),
			    double(k)  //  / double(node_nw)
			);
		}
	    }
	}

	vector<index_t> to_kill(node_nu* node_nv* node_nw, 0);
	FOR(k, 32) {
	    FOR(j, 248) {
		FOR(i, 488) {
		    to_kill[node_nu * node_nv * k + j * node_nu + i] = 1;
		}
	    }
	}

	mesh.vertices.delete_elements(to_kill);

	mesh_grob()->update();
    }

/*****************************************************************/

    class Geometry {
    public:
	Geometry() : vertex_identity{mesh.vertices.attributes(), "identity"} {}

	void make_tetrahedron(double r) {
	    const double sq2 = std::sqrt(2.0);
	    const double sq3 = std::sqrt(3.0);
	    mesh.vertices.create_vertices(4);
	    mesh.vertices.point(0) =
		GEO::vec3{ 0,                    r,        0                  };
	    mesh.vertices.point(1) =
		GEO::vec3{ 0,                   -r / 3.0,  r * 2.0 * sq2 / 3.0};
	    mesh.vertices.point(2) =
		GEO::vec3{-r * sq3 * sq2 / 3.0, -r / 3.0, -r * sq2 / 3.0      };
	    mesh.vertices.point(3) =
		GEO::vec3{ r * sq3 * sq2 / 3.0, -r / 3.0, -r * sq2 / 3.0      };
	    mesh.facets.create_triangle(0, 1, 3);
	    mesh.facets.create_triangle(0, 2, 1);
	    mesh.facets.create_triangle(0, 3, 2);
	    mesh.facets.create_triangle(1, 2, 3);
	}

	void copy_from(const Geometry& source) {
	    vertex_identity.unbind();
	    mesh.copy(source.mesh, false);
	    vertex_identity.bind(mesh.vertices.attributes(), "identity");
	    for (const GEO::index_t vertex : source.mesh.vertices) {
		vertex_identity[vertex] = source.vertex_identity[vertex];
		// Error: Assertion failed: i < superclass::nb_elements()
	    }
	}

	void update_identity() {
	    for (GEO::index_t vertex : mesh.vertices) {
		vertex_identity[vertex] = vertex;
	    }
	}
	GEO::Mesh mesh;
	GEO::Attribute<GEO::index_t> vertex_identity;
    };

    void MeshGrobExperimentCommands::geogram_issue_227() {
	Geometry A;
	A.make_tetrahedron(1.0);
	A.update_identity();
	Geometry B;
	B.copy_from(A);
    }

    double mesh_area_with_syntaxic_sugar(const Mesh& M) {
	double result = 0.0;
	for(index_t f: M.facets) {
	    for(auto [ p1, p2, p3]: M.facets.triangle_points(f)) {
		result += Geom::triangle_area(p1,p2,p3);
	    }
	}
	return result;
    }

    double mesh_area_plain(const Mesh& M) {
	double result = 0.0;
	for(index_t f: M.facets) {
	    index_t v0 = M.facets.vertex(f,0);
	    const vec3& p0 = M.vertices.point(v0);
	    for(index_t le=0; le+1<M.facets.nb_vertices(f); ++le) {
		index_t v1 = M.facets.vertex(f,le);
		const vec3& p1 = M.vertices.point(v1);
		index_t v2 = M.facets.vertex(f,le+1);
		const vec3& p2 = M.vertices.point(v2);
		result += Geom::triangle_area(p0,p1,p2);
	    }
	}
	return result;
    }

    void MeshGrobExperimentCommands::test_syntaxic_sugar() {
	double a1 = 0.0;
	double a2 = 0.0;
	{
	    Stopwatch W1("plain");
	    a1 = mesh_area_plain(*mesh_grob());
	}
	{
	    Stopwatch W2("with sugar");
	    a2 = mesh_area_plain(*mesh_grob());
	}
	Logger::out("Sugar") << a1 << " " << a2 << std::endl;
    }

    void MeshGrobExperimentCommands::test_syntaxic_sugar_2() {
	for(vec3& p: mesh_grob()->vertices.points()) {
	    p += {0.0, 0.0, 1.0};
	}
	mesh_grob()->update();
    }

    void MeshGrobExperimentCommands::test_PR_252() {
	MeshGrob& m = *mesh_grob();
	m.vertices.create_vertices(8);
	m.vertices.point(0) = GEO::vec3{0.0, 0.0, 0.0};
	m.vertices.point(1) = GEO::vec3{0.0, 0.0, 1.0};
	m.vertices.point(2) = GEO::vec3{1.0, 0.0, 1.0};
	m.vertices.point(3) = GEO::vec3{1.0, 0.0, 0.0};
	m.vertices.point(4) = GEO::vec3{0.0, 1.0, 0.0};
	m.vertices.point(5) = GEO::vec3{0.0, 1.0, 1.0};
	m.vertices.point(6) = GEO::vec3{1.0, 1.0, 1.0};
	m.vertices.point(7) = GEO::vec3{1.0, 1.0, 0.0};
	m.facets.create_quad(0, 1, 2, 3); // facet 0 - (3, 2, 1, 0) if flipped
	m.facets.create_quad(0, 4, 5, 1); // facet 1 - (1, 5, 3, 0) if flipped
	m.facets.create_quad(0, 3, 7, 4); // facet 2 - (4, 7, 3, 0) if flipped
	m.facets.create_quad(1, 2, 6, 5);
	m.facets.create_quad(2, 3, 7, 6);
	m.facets.create_quad(4, 5, 6, 7);
	m.facets.connect();
	//GEO::mesh_reorient(m);
	GEO::mesh_repair(m,GEO::MeshRepairMode(0),0.0);
	m.update();
    }

/*****************************************************************/

}
