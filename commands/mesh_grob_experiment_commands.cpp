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
#include <OGF/Experiment/algo/tiny_cas.h>

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

#ifdef EXPERIMENT_WITH_CGAL
#include <OGF/Experiment/CGAL/delaunay_nd_cgal.h>
#endif

#ifdef GEOGRAM_WITH_VORPALINE
#include <vorpalib/mesh/mesh_weiler_model.h>
#include <vorpalib/mesh/mesh_geomodel.h>
#endif

#include <array>
#include <complex>

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


                TriangleIsects II;

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

    void MeshGrobExperimentCommands::test_discussion_364() {
        // Defines points
	std::vector<double> points_coordinates = {
	    0, 0,
	    1, 2,
	    3, 2,
	    2, 0,
	    2, 2,
	    1, 1,
	    2, 1
	};
	index_t n = static_cast<index_t>(points_coordinates.size() / 2);

	double eps = 0.1;
	double xmin = 0 - eps, ymin = 0 - eps, xmax = 3 + eps, ymax = 2 + eps;

	// Defines constraints
	std::vector<std::array<GEO::index_t, 2> > segments = {
	    {0, 1},
	    {1, 4},
	    {4, 2},
	    {2, 3},
	    {3, 0},
	    {4, 5},
	    {5, 6},
	    {6, 4}
	};

	// Computes triangulation
	GEO::CDT2d cdt;
	cdt.create_enclosing_rectangle(xmin, ymin, xmax, ymax);
	cdt.set_delaunay(true);
	cdt.insert(n, points_coordinates.data());
	for (const std::array<GEO::index_t, 2> & segment : segments) {
	    cdt.insert_constraint(segment.at(0)+4, segment.at(1)+4);
	}
	cdt.save("CDT.geogram");
    }

/*****************************************************************/
}

namespace OGF {

    void MeshGrobExperimentCommands::copy_edges() {
       for(index_t f: mesh_grob()->facets) {
	   index_t N = mesh_grob()->facets.nb_vertices(f);
	   for(index_t lv=0; lv<N; ++lv) {
	       index_t v1 = mesh_grob()->facets.vertex(f,lv);
	       index_t v2 = mesh_grob()->facets.vertex(f,(lv+1)%N);
	       mesh_grob()->edges.create_edge(v1,v2);
	   }
       }
       mesh_grob()->update();
    }

}


namespace OGF {

    void MeshGrobExperimentCommands::test_orient3d_SOS() {
	using OGF::TinyCAS::Polynom;

	std::vector<std::string> var_names = {
	    "x1","y1","z1",
	    "x2","y2","z2",
	    "x3","y3","z3",
	    "x4","y4","z4",
	    "eps"
	};

	index_t nb_var = index_t(var_names.size());

	Polynom x1(nb_var,0u);
	Polynom y1(nb_var,1u);
	Polynom z1(nb_var,2u);
	Polynom x2(nb_var,3u);
	Polynom y2(nb_var,4u);
	Polynom z2(nb_var,5u);
	Polynom x3(nb_var,6u);
	Polynom y3(nb_var,7u);
	Polynom z3(nb_var,8u);
	Polynom x4(nb_var,9u);
	Polynom y4(nb_var,10u);
	Polynom z4(nb_var,11u);

	Polynom eps(nb_var,12u);
	Polynom one(nb_var,1.0);

	Polynom D = det4x4(
	    x1+eps.pow(1),   y1+eps.pow(2),    z1+eps.pow(4),    one,
	    x2+eps.pow(8),   y2+eps.pow(16),   z2+eps.pow(32),   one,
	    x3+eps.pow(64),  y3+eps.pow(128),  z3+eps.pow(256),  one,
	    x4+eps.pow(512), y4+eps.pow(1024), z4+eps.pow(2048), one
	);

	Logger::out("o3dSOS") << "D= " << D.to_string(var_names) << std::endl;

	index_t max_eps_pow = 0;
	for(const Polynom::Monom& M: D.monoms_)	{
	    max_eps_pow = std::max(max_eps_pow, M.pow[12u]);
	}

	vector<Polynom> perturbation(max_eps_pow+1, Polynom(nb_var));
	for(Polynom::Monom M: D.monoms_) {
	    index_t eps_pow = M.pow[12u];
	    M.pow[12u] = 0;
	    perturbation[eps_pow] += M;
	}

	for(index_t p=0; p<=max_eps_pow; p++) {
	    if(perturbation[p].is_zero()) {
		continue;
	    }
	    Logger::out("o3dSOS") << "eps^" << p << " : "
				  << perturbation[p].to_string(var_names)
				  << std::endl;
	}

    }

    // Generated using ChatGPT:
    // Prompt: "Hello, can you please generate for me a C++ program that
    // creates a seven-pointed star using geogram ?" 04/13/2026
    // Fixes: had to refine, "Mesh facets set_nb_vertices() does not exist,
    // use create_polygon() instead".
    // Replaced ints with index_t (note: Claude directly uses intex_t it seems)
    // Created vec3's directly instead of setting x,y,z

    void MeshGrobExperimentCommands::seven_pointed_star() {
	Mesh& M = *mesh_grob();

	const int N = 7;
	const int V2D = 2 * N;          // 14 vertices in the star outline
	const double R_outer = 1.0;
	const double R_inner = 0.45;
	const double height   = 0.3;
	const double PI = 3.141592653589793;

	// Create vertices: bottom layer + top layer
	M.vertices.create_vertices(2 * V2D);

	for(index_t i = 0; i < N; ++i) {
	    double angle_outer = 2.0 * PI * double(i) / double(N);
	    double angle_inner = angle_outer + PI / double(N);

	    double ox = R_outer * std::cos(angle_outer);
	    double oy = R_outer * std::sin(angle_outer);

	    double ix = R_inner * std::cos(angle_inner);
	    double iy = R_inner * std::sin(angle_inner);

	    // Bottom layer (z = 0)
	    M.vertices.point(2*i) = vec3(ox,oy,0.0);
	    M.vertices.point(2*i+1) = vec3(ix,iy,0.0);

	    // Top layer (z = height)
	    M.vertices.point(V2D + 2*i) = vec3(ox,oy,height);
	    M.vertices.point(V2D + 2*i+1) = vec3(ix,iy,height);
	}

	// -------------------------------
	// Bottom face (reverse orientation)
	// -------------------------------
	{
	    index_t f = M.facets.create_polygon(V2D);
	    for(index_t k = 0; k < V2D; ++k) {
		M.facets.set_vertex(f, k, V2D - 1 - k);
	    }
	}

	// -------------------------------
	// Top face
	// -------------------------------
	{
	    index_t f = M.facets.create_polygon(V2D);
	    for(index_t k = 0; k < V2D; ++k) {
		M.facets.set_vertex(f, k, V2D + k);
	    }
	}

	// -------------------------------
	// Side faces (quads)
	// -------------------------------
	for(index_t i = 0; i < V2D; ++i) {
	    index_t next = (i + 1) % V2D;

	    index_t f = M.facets.create_polygon(4);
	    M.facets.set_vertex(f, 0, i);
	    M.facets.set_vertex(f, 1, next);
	    M.facets.set_vertex(f, 2, V2D + next);
	    M.facets.set_vertex(f, 3, V2D + i);
	}

	mesh_grob()->update();
    }

/* ==================================================================
   EGG SURFACE
   ================================================================== */

/**
 * @brief Evaluate the egg parametric surface at (u, v).
 *
 * The egg shape is obtained by deforming a sphere:
 *   - a, b  : radii in X and Y (typically equal -> round cross-section)
 *   - c     : half-height along Z
 *   - taper : asymmetry factor that makes one pole narrower (the "pointed" end)
 *
 * @param u      Polar angle in [0, PI]
 * @param v      Azimuth angle in [0, 2*PI)
 * @param a      Radius in X
 * @param b      Radius in Y
 * @param c      Radius in Z (base height half-extent)
 * @param taper  Taper strength (0 = sphere, 0.3 = typical egg)
 * @return vec3  3D point on the surface
 */
    vec3 egg_point(double u, double v,
		   double a = 1.0, double b = 1.0,
		   double c = 1.3, double taper = 0.35)
    {
	/* r(u) modulates radius along the Z axis to create the egg asymmetry.
	   r > 1 near the top (u near 0), r < 1 near the bottom (u near PI). */
	double r = 1.0 + taper * std::cos(u);
	return vec3(
	    a * r * std::sin(u) * std::cos(v),
	    b * r * std::sin(u) * std::sin(v),
	    c * std::cos(u)
	);
    }

/**
 * @brief Build a triangulated surface mesh of an egg shape.
 *
 * Topology: a UV sphere mesh with nu latitude bands and nv longitude slices.
 * The two pole vertices are shared, and each quad band is split into 2 triangles.
 *
 * @param[out] mesh   Output Geogram mesh (cleared on entry).
 * @param nu          Number of latitude subdivisions (>= 2, recommended >= 32)
 * @param nv          Number of longitude subdivisions (>= 3, recommended >= 48)
 * @param a           X radius
 * @param b           Y radius
 * @param c           Z half-height
 * @param taper       Asymmetry factor (0 = no taper / sphere, 0.35 = egg)
 */
    void create_egg_mesh(Mesh& mesh,
			 index_t nu    = 48,
			 index_t nv    = 64,
			 double  a     = 1.0,
			 double  b     = 1.0,
			 double  c     = 1.3,
			 double  taper = 0.35)
    {
	mesh.clear();
	mesh.vertices.set_dimension(3);

	/* ------------------------------------------------------------------
	   1. Generate vertices
	   Layout: 1 north pole + (nu-1) rings of nv vertices + 1 south pole
	   Total : 1 + (nu-1)*nv + 1 = (nu-1)*nv + 2
	   ------------------------------------------------------------------ */
	const double pi     = M_PI;
	const double two_pi = 2.0 * M_PI;

	/* North pole (u = 0) */
	{
	    vec3 p = egg_point(0.0, 0.0, a, b, c, taper);
	    index_t vid = mesh.vertices.create_vertex();
	    mesh.vertices.point(vid) = p;
	}

	/* Interior rings (u = i * PI / nu,  i = 1 ... nu-1) */
	for (index_t i = 1; i < nu; ++i) {
	    double u = static_cast<double>(i) * pi / static_cast<double>(nu);
	    for (index_t j = 0; j < nv; ++j) {
		double v_angle = static_cast<double>(j) * two_pi / static_cast<double>(nv);
		vec3 p = egg_point(u, v_angle, a, b, c, taper);
		index_t vid = mesh.vertices.create_vertex();
		mesh.vertices.point(vid) = p;
	    }
	}

	/* South pole (u = PI) */
	{
	    vec3 p = egg_point(pi, 0.0, a, b, c, taper);
	    index_t vid = mesh.vertices.create_vertex();
	    mesh.vertices.point(vid) = p;
	}

	/* ------------------------------------------------------------------
	   Helper lambda: vertex index in ring i (1-based), longitude slot j.
	   ------------------------------------------------------------------ */
	const index_t north_idx = 0;
	const index_t south_idx = 1 + (nu - 1) * nv;

	auto ring_vertex = [&](index_t i, index_t j) -> index_t {
	    /* i in [1, nu-1], j wrapped in [0, nv) */
	    return 1 + (i - 1) * nv + (j % nv);
	};

	/* ------------------------------------------------------------------
	   2. Generate triangles
	   ------------------------------------------------------------------ */

	/* North cap: triangles from north pole to first ring */
	for (index_t j = 0; j < nv; ++j) {
	    index_t v0 = north_idx;
	    index_t v1 = ring_vertex(1, j);
	    index_t v2 = ring_vertex(1, j + 1);
	    index_t f  = mesh.facets.create_triangle(v0, v1, v2);
	    (void)f;
	}

	/* Middle bands: each pair of adjacent rings forms quads split into 2 triangles */
	for (index_t i = 1; i + 1 < nu; ++i) {
	    for (index_t j = 0; j < nv; ++j) {
		index_t a0 = ring_vertex(i,     j);
		index_t a1 = ring_vertex(i,     j + 1);
		index_t b0 = ring_vertex(i + 1, j);
		index_t b1 = ring_vertex(i + 1, j + 1);

		mesh.facets.create_triangle(a0, b0, a1);  /* triangle 1 */
		mesh.facets.create_triangle(a1, b0, b1);  /* triangle 2 */
	    }
	}

	/* South cap: triangles from last ring to south pole */
	for (index_t j = 0; j < nv; ++j) {
	    index_t v0 = ring_vertex(nu - 1, j);
	    index_t v1 = south_idx;
	    index_t v2 = ring_vertex(nu - 1, j + 1);
	    mesh.facets.create_triangle(v0, v1, v2);
	}

	/* ------------------------------------------------------------------
	   3. Connect adjacency (needed by many Geogram algorithms)
	   ------------------------------------------------------------------ */
	mesh.facets.connect();
    }

    void MeshGrobExperimentCommands::egg() {
	mesh_grob()->clear();
	create_egg_mesh(*mesh_grob());
	mesh_grob()->update();
    }


/* ==================================================================
   CRACK SURFACE
   ================================================================== */

/**
 * @brief Build an independent crack surface mesh along the upper hemisphere
 *        of an egg.
 *
 * Design intent
 * -------------
 * The crack represents the fracture line along the egg equator.
 * It is intentionally NOT conforming to the egg surface:
 *
 *   - The outer boundary follows the egg surface at u = PI/2 (the equator),
 *     with a deterministic zigzag in the polar direction produced by a
 *     cosine wave.  This gives the characteristic jagged broken-edge look
 *     without any random numbers.
 *   - The inner boundary is a flat ellipse at z = crack_z (inside the egg).
 *     The crack surface therefore dips into the egg interior, which is
 *     intentional -- the interior does not need to match the egg geometry.
 *
 * Vertex layout
 * -------------
 *   outer ring  [0 .. nv-1]   -- on the egg surface at the perturbed equator
 *   inner ring  [nv .. 2nv-1] -- flat ellipse at z = crack_z
 *   center      [2*nv]        -- single point at (0, 0, crack_z)
 *
 * The crack has nv_teeth full zigzag teeth around the circumference.
 * Each tooth is a cosine bump in the polar angle u: at tooth peaks the
 * crack rises toward the top cap; at valleys it dips toward the bottom.
 *
 * @param[out] crack     Output Geogram mesh (cleared on entry).
 * @param nv             Azimuthal vertex count (>= 2 * nv_teeth).
 * @param a              Egg X radius (must match egg mesh params).
 * @param b              Egg Y radius.
 * @param c              Egg Z half-height.
 * @param taper          Egg taper (must match egg mesh params).
 * @param crack_z        Z height of the flat inner disk (0 = egg center).
 * @param tooth_amp      Amplitude of each zigzag tooth in polar-angle units
 *                       (radians).  0.3 gives a clearly visible crack.
 * @param nv_teeth       Number of crack teeth around the full circumference.
 */
    void create_egg_crack_mesh(Mesh& crack,
			       index_t nv         = 128,
			       double  a          = 1.0,
			       double  b          = 1.0,
			       double  c          = 1.3,
			       double  taper      = 0.35,
			       double  crack_z    = 0.0,
			       double  tooth_amp  = 0.25,
			       index_t nv_teeth   = 8)
    {
	crack.clear();
	crack.vertices.set_dimension(3);

	const double pi     = M_PI;
	const double two_pi = 2.0 * M_PI;

	/* ------------------------------------------------------------------
	   Tooth profile: deterministic zigzag in the polar angle u.

	   u_crack(v) = PI/2 + tooth_amp * cos(nv_teeth * v)

	   At cosine peaks  (cos = +1): u < PI/2  -> vertex is above the equator.
	   At cosine valleys(cos = -1): u > PI/2  -> vertex is below the equator.
	   This produces a clean, repeating sawtooth crack with nv_teeth bumps.
	   ------------------------------------------------------------------ */
	auto u_crack = [&](double v_angle) -> double {
	    return pi * 0.5 + tooth_amp * std::cos(static_cast<double>(nv_teeth) * v_angle);
	};

	/* ------------------------------------------------------------------
	   1. Outer ring vertices -- on the egg surface along the crack profile.
	   ------------------------------------------------------------------ */
	for (index_t j = 0; j < nv; ++j) {
	    double v_angle = static_cast<double>(j) * two_pi / static_cast<double>(nv);
	    double u       = u_crack(v_angle);
	    vec3   p       = egg_point(u, v_angle, a, b, c, taper);

	    index_t vid = crack.vertices.create_vertex();
	    crack.vertices.point(vid) = p;
	}

	/* ------------------------------------------------------------------
	   2. Inner ring vertices -- flat ellipse at z = crack_z.

	   The radius is 80% of the egg equatorial radius so the ribbon has
	   visible width everywhere.
	   ------------------------------------------------------------------ */
	const double inner_scale = 0.80;

	for (index_t j = 0; j < nv; ++j) {
	    double v_angle = static_cast<double>(j) * two_pi / static_cast<double>(nv);
	    vec3   p(a * inner_scale * std::cos(v_angle),
		     b * inner_scale * std::sin(v_angle),
		     crack_z);

	    index_t vid = crack.vertices.create_vertex();
	    crack.vertices.point(vid) = p;
	}

	/* ------------------------------------------------------------------
	   3. Center vertex -- closes the inner disk fan.
	   ------------------------------------------------------------------ */
	const index_t center_idx = 2 * nv;
	{
	    index_t vid = crack.vertices.create_vertex();
	    crack.vertices.point(vid) = vec3(0.0, 0.0, crack_z);
	}

	/* ------------------------------------------------------------------
	   Index helpers (wrap in [0, nv))
	   ------------------------------------------------------------------ */
	auto outer = [&](index_t j) -> index_t { return          j % nv; };
	auto inner = [&](index_t j) -> index_t { return nv + (j % nv); };

	/* ------------------------------------------------------------------
	   4. Ribbon triangles: outer ring -> inner ring.
	   Each quad is split into 2 triangles with consistent winding.
	   ------------------------------------------------------------------ */
	for (index_t j = 0; j < nv; ++j) {
	    index_t o0 = outer(j);
	    index_t o1 = outer(j + 1);
	    index_t i0 = inner(j);
	    index_t i1 = inner(j + 1);

	    crack.facets.create_triangle(o0, i0, o1);
	    crack.facets.create_triangle(o1, i0, i1);
	}

	/* ------------------------------------------------------------------
	   5. Inner disk fan: center -> inner ring.
	   ------------------------------------------------------------------ */
	for (index_t j = 0; j < nv; ++j) {
	    crack.facets.create_triangle(center_idx, inner(j + 1), inner(j));
	}

	/* ------------------------------------------------------------------
	   6. Connect adjacency.
	   ------------------------------------------------------------------ */
	crack.facets.connect();
    }


    void MeshGrobExperimentCommands::egg_crack() {
	mesh_grob()->clear();
	create_egg_crack_mesh(*mesh_grob());
	mesh_grob()->update();
    }


    void MeshGrobExperimentCommands::create_torus(
	double R,
	double r,
	index_t Nu,
	index_t Nv
    ) {
	// ------------------------------------------------------------------
	// Retrieve (or create) the target MeshGrob
	// ------------------------------------------------------------------
	MeshGrob& M = *mesh_grob();
	M.clear();

	// ------------------------------------------------------------------
	// Vertices
	// ------------------------------------------------------------------
	// We generate exactly Nu x Nv vertices (no duplicated seam).
	// Seam closure is handled by modular indexing in the face loop.
	// Parametric equations:
	//   x(u,v) = (R + r*cos(v)) * cos(u)
	//   y(u,v) = (R + r*cos(v)) * sin(u)
	//   z(u,v) =        r*sin(v)
	// with u in [0, 2pi) (around torus) and v in [0, 2pi) (around tube).

	M.vertices.set_dimension(3);

	for (index_t j = 0; j < Nv; ++j) {
	    double v  = 2.0 * M_PI * double(j) / double(Nv);
	    double cv = std::cos(v);
	    double sv = std::sin(v);

	    for (index_t i = 0; i < Nu; ++i) {
		double u  = 2.0 * M_PI * double(i) / double(Nu);
		double cu = std::cos(u);
		double su = std::sin(u);

		double x = (R + r * cv) * cu;
		double y = (R + r * cv) * su;
		double z =       r * sv;

		M.vertices.create_vertex(GEO::vec3(x, y, z).data());
	    }
	}

	// ------------------------------------------------------------------
	// Faces  (two triangles per quad, counter-clockwise winding)
	// ------------------------------------------------------------------
	// Vertex index helper: row j, column i  =>  j*Nu + i
	// Modular wrap: (i+1)%Nu and (j+1)%Nv close the seam exactly,
	// so the last ring of quads shares the very first row/column of
	// vertices — no gap, no duplicate geometry.
	auto vid = [&](index_t i, index_t j) -> index_t {
	    return (j % Nv) * Nu + (i % Nu);
	};

	for (index_t j = 0; j < Nv; ++j) {
	    for (index_t i = 0; i < Nu; ++i) {
		index_t v00 = vid(i,     j    );
		index_t v10 = vid(i + 1, j    );
		index_t v01 = vid(i,     j + 1);
		index_t v11 = vid(i + 1, j + 1);

		// Triangle 1: lower-left, lower-right, upper-right
		{
		    index_t f = M.facets.create_triangle(v00, v10, v11);
		    (void)f;
		}
		// Triangle 2: lower-left, upper-right, upper-left
		{
		    index_t f = M.facets.create_triangle(v00, v11, v01);
		    (void)f;
		}
	    }
	}

	// ------------------------------------------------------------------
	// Finalise
	// ------------------------------------------------------------------
	M.facets.connect();
	M.update();

	Logger::out("CreateTorus")
	    << "Torus created: "
	    << M.vertices.nb() << " vertices, "
	    << M.facets.nb()   << " triangles."
	    << std::endl;

    }
}
