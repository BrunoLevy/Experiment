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
 
#include <OGF/Experiment/commands/mesh_grob_experiment_commands.h>
#include <OGF/Experiment/commands/geometry.h>
#include <geogram/mesh/triangle_intersection.h>
#include <geogram/mesh/mesh_AABB.h>
#include <geogram/mesh/mesh_repair.h>
#include <geogram/mesh/index.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/delaunay/delaunay.h>
#include <geogram/numerics/predicates.h>
#include <geogram/numerics/expansion_nt.h>

#include <sstream>

namespace {
    using namespace OGF;

    /***********************************************************************/

    /**
     * \brief Meshes a single triangle with the constraints that come from
     *  the intersections with the other triangles.
     */
    class MeshInTriangle {
    public:
        
        class Vertex {
        public:
            enum { NO_INDEX = index_t(-1) };
            
            Vertex(
                Mesh& M,
                index_t f1, index_t f2,
                TriangleRegion R1, TriangleRegion R2
            ) {
                mesh = &M;
                init_sym(f1,f2,R1,R2);
                init_geometry();
                mesh_vertex_index = NO_INDEX;
            }


            Vertex(Mesh& M, index_t f, index_t lv) {
                mesh = &M;
                index_t v = M.facets.vertex(f,lv);
                set_indices(v);
                init_geometry();
                mesh_vertex_index = v;
                tt_sym.f1 = f;
                tt_sym.R1 = TriangleRegion(lv);
                tt_sym.f2 = index_t(-1);
            }

            Vertex(Mesh& M, index_t v) {
                mesh = &M;
                set_indices(v);
                init_geometry();
                mesh_vertex_index = v;
                tt_sym.f1 = index_t(-1);
                tt_sym.f2 = index_t(-1);
            }

            Vertex() {
                mesh = nullptr;                
                set_indices(NO_INDEX, NO_INDEX, NO_INDEX, NO_INDEX);
                mesh_vertex_index = NO_INDEX;
                tt_sym.f1 = index_t(-1);
                tt_sym.f2 = index_t(-1);
            }

            bool is_existing_vertex() const {
                return
                    sym.indices[1] == NO_INDEX &&
                    sym.indices[2] == NO_INDEX &&
                    sym.indices[3] == NO_INDEX ;
            }

            bool is_edge_edge_isect() const {
                return sym.indices[3] != NO_INDEX ;
            }

            bool is_edge_triangle_isect() const {
                return sym.indices[2] != NO_INDEX &&
                       sym.indices[3] == NO_INDEX ;
            }

            void print(std::ostream& out=std::cerr) const {
                out << "[ ";
                if(is_existing_vertex()) {
                    out << "v(" << sym.indices[0] << ")";
                } else if(is_edge_edge_isect()) {
                    out << "E(" << sym.indices[0] << "," << sym.indices[1]
                        << ") /\\ E("
                        << sym.indices[2] << "," << sym.indices[3]
                        << ")"; 
                } else if(is_edge_triangle_isect()) {
                    index_t f = sym.indices[2];
                    out << "E(" << sym.indices[0] << "," << sym.indices[1]
                        << ") /\\ F("
                        << sym.indices[2]
                        << ":"
                        << mesh->facets.vertex(f,0) << ","
                        << mesh->facets.vertex(f,1) << ","
                        << mesh->facets.vertex(f,2) 
                        << ")";
                }
                if(tt_sym.f1 != index_t(-1)) {
                    out << " (= ";
                    out << tt_sym.f1;
                    out << region_to_string(tt_sym.R1).substr(2);
                }
                if(tt_sym.f2 != index_t(-1)) {
                    out << " /\\ ";
                    out << tt_sym.f2;
                    out << region_to_string(tt_sym.R2).substr(2);
                }
                if(tt_sym.f1 != index_t(-1)) {
                    out << ") ";
                }
                out << " ]";
            }

            std::string to_string() const {
                std::ostringstream out;
                print(out);
                return out.str();
            }

            bool has_tt_sym() const {
                return tt_sym.f2 != index_t(-1);
            }
            
            
        protected:

            void init_sym(
                index_t f1, index_t f2,
                TriangleRegion R1, TriangleRegion R2
            ) {
                tt_sym.f1 = f1;
                tt_sym.f2 = f2;
                tt_sym.R1 = R1;
                tt_sym.R2 = R2;
                
                index_t v[6] = {
                    mesh->facets.vertex(f1,0),
                    mesh->facets.vertex(f1,1),
                    mesh->facets.vertex(f1,2),
                    mesh->facets.vertex(f2,0),
                    mesh->facets.vertex(f2,1),
                    mesh->facets.vertex(f2,2)
                };
            
                if(region_dim(R1) == 0) {
                    set_indices(v[R1]);
                    return;
                } 

                if(region_dim(R2) == 0) {
                    set_indices(v[R2]);
                    return;
                }
        
                if(region_dim(R1) == 1 && region_dim(R2) == 1) {
                    bindex E1 = get_edge_vertices(R1,v);
                    bindex E2 = get_edge_vertices(R2,v);
                    if(E1.indices[0] < E2.indices[0]) {
                        set_indices(
                            E1.indices[0], E1.indices[1],
                            E2.indices[0], E2.indices[1]
                        );
                    } else {
                        set_indices(
                            E2.indices[0], E2.indices[1],
                            E1.indices[0], E1.indices[1]
                        );
                    }
                    return;
                }
                
                if(region_dim(R1) == 1) {
                    bindex E = get_edge_vertices(R1,v);
                    set_indices(E.indices[0], E.indices[1], f2);
                    return;
                }
                
                if(region_dim(R2) == 1) {
                    bindex E = get_edge_vertices(R2,v);
                    set_indices(E.indices[0], E.indices[1], f1);
                    return;
                }
                geo_assert_not_reached;
            }


            template <class T>
            vecng<3,T> get_point() {
                if(is_existing_vertex()) {
                    return convert_vec3_generic<T>(
                        vec3(mesh->vertices.point_ptr(sym.indices[0]))
                    );
                }

                if(is_edge_edge_isect()) {

                    // We distinguish 3D edge /\ edge isect (then we use
                    // edge /\ facet isect) and 2D edge /\ edge isect (
                    // then develop special code). 
                    
                    geo_debug_assert(has_tt_sym());
                    
                    // The three vertices of f1, with the intersected
                    // edges as (p1,p2)
                    vec3 p1, p2, p3;
                    
                    // The intersected edge of f2
                    vec3 q1, q2;

                    {
                        TriangleRegion lvr1, lvr2;
                        index_t lv1, lv2, lv3;
                        // The two vertices of R1 first
                        GEO::get_edge_vertices(tt_sym.R1,lvr1, lvr2);
                        lv1 = index_t(lvr1);
                        lv2 = index_t(lvr2);
                        geo_debug_assert(lv1 < 3);
                        geo_debug_assert(lv2 < 3);                        
                        // The other vertex of f1
                        lv3 = index_t(tt_sym.R1 - 6);
                        geo_debug_assert(lv3 < 3);
                        
                        index_t v1 = mesh->facets.vertex(tt_sym.f1, lv1);
                        index_t v2 = mesh->facets.vertex(tt_sym.f1, lv2);
                        index_t v3 = mesh->facets.vertex(tt_sym.f1, lv3);

                        p1 = vec3(mesh->vertices.point_ptr(v1));
                        p2 = vec3(mesh->vertices.point_ptr(v2));
                        p3 = vec3(mesh->vertices.point_ptr(v3));
                    }

                    {
                        TriangleRegion lvr1, lvr2;
                        GEO::get_edge_vertices(tt_sym.R2, lvr1, lvr2);
                        index_t lv1 = index_t(lvr1);
                        index_t lv2 = index_t(lvr2);
                        geo_debug_assert(lv1 >= 3 && lv1 < 6);
                        geo_debug_assert(lv2 >= 3 && lv2 < 6);
                        lv1 -= 3;
                        lv2 -= 3;
                        index_t v1 = mesh->facets.vertex(tt_sym.f2,lv1);
                        index_t v2 = mesh->facets.vertex(tt_sym.f2,lv2);
                        q1 = vec3(mesh->vertices.point_ptr(v1));
                        q2 = vec3(mesh->vertices.point_ptr(v2));
                    }

                    if(
                        PCK::orient_3d(p1,p2,p3,q1) == ZERO &&
                        PCK::orient_3d(p1,p2,p3,q2) == ZERO
                    ) {
                        return get_segment_segment_intersection_2D<T>(
                            q1,q2,p1,p2,
                            ::GEO::Geom::triangle_normal_axis(p1,p2,p3)
                        );
                    } 

                    // We are in 3D, we can use segment /\ triangle
                    // intersection (even if we know that the intersection
                    // will be on on (p1,p2))
                    return get_segment_triangle_intersection<T>(
                        q1,q2,p1,p2,p3
                    );
                }
                
                geo_debug_assert(is_edge_triangle_isect());
                
                vec3 q1(mesh->vertices.point_ptr(sym.indices[0]));
                vec3 q2(mesh->vertices.point_ptr(sym.indices[1]));
                index_t v1 = mesh->facets.vertex(sym.indices[2],0);
                index_t v2 = mesh->facets.vertex(sym.indices[2],1);
                index_t v3 = mesh->facets.vertex(sym.indices[2],2);
                vec3 p1(mesh->vertices.point_ptr(v1));
                vec3 p2(mesh->vertices.point_ptr(v2));
                vec3 p3(mesh->vertices.point_ptr(v3));
                
                return get_segment_triangle_intersection<T>(q1,q2,p1,p2,p3);
            }

            void init_geometry() {
                point       = get_point<double>();
                point_exact = get_point<rational_nt>();
                // TODO:
                point.x = point_exact.x.estimate();
                point.y = point_exact.y.estimate();
                point.z = point_exact.z.estimate();                
            }
            
            void set_indices(
                index_t i0=NO_INDEX, index_t i1=NO_INDEX,
                index_t i2=NO_INDEX, index_t i3=NO_INDEX
            ) {
                sym.indices[0] = i0;
                sym.indices[1] = i1;
                sym.indices[2] = i2;
                sym.indices[3] = i3;                        
            }
            
            bindex get_edge_vertices(TriangleRegion R, index_t v[]) {
                TriangleRegion v1,v2;
                geo_debug_assert(region_dim(R) == 1);
                ::GEO::get_edge_vertices(R,v1,v2);
                return bindex(v[v1],v[v2]);
            }

        public:
            Mesh* mesh;
            vec3 point;
            vec3Q point_exact;

            // Symbolic information, global indices in mesh
            quadindex sym;

            // Global mesh vertex index once created
            index_t mesh_vertex_index;

            // Symbolic information - triangle-triangle isect
            struct {
                index_t f1,f2;        // global facet indices in mesh
                TriangleRegion R1,R2; // triangle regions
            } tt_sym; 
        };

        
        MeshInTriangle(Mesh& M) : mesh_(M), f1_(index_t(-1)) {
        }

        void begin_facet(index_t f) {
            f1_ = f;
            index_t v[3] = {
                mesh_.facets.vertex(f,0),
                mesh_.facets.vertex(f,1),
                mesh_.facets.vertex(f,2)
            };
            for(index_t lv=0; lv<3; ++lv) {
                add_vertex(Vertex(mesh_, f, lv));
            }
            f1_normal_axis_ = ::GEO::Geom::triangle_normal_axis(
                vec3(mesh_.vertices.point_ptr(v[0])),
                vec3(mesh_.vertices.point_ptr(v[1])),
                vec3(mesh_.vertices.point_ptr(v[2]))                
            );
            u_ = coord_index_t((f1_normal_axis_ + 1) % 3);
            v_ = coord_index_t((f1_normal_axis_ + 2) % 3);
            f1_orient_ = PCK::orient_2d(project(0), project(1), project(2));
        }
        
        index_t add_vertex(
            index_t f2,
            TriangleRegion R1, TriangleRegion R2
        ) {
            geo_debug_assert(f1_ != index_t(-1));
            Vertex V(mesh_, f1_, f2, R1, R2);
            index_t sz = vertex_.size();
            index_t v = add_vertex(V);
            // If it is a new vertex, and if it is on
            // an edge, add it to the list of edge's
            // vertices
            if(vertex_.size() > sz) {
                if(R1 == T1_RGN_E0) {
                    vertices_in_E_[0].push_back(v);
                } else if(R1 == T1_RGN_E1) {
                    vertices_in_E_[1].push_back(v);                
                } else if(R1 == T1_RGN_E2) {
                    vertices_in_E_[2].push_back(v);                
                }
            }
            return v;
        }

        void add_edge(
            index_t f2,
            TriangleRegion AR1, TriangleRegion AR2,
            TriangleRegion BR1, TriangleRegion BR2
        ) {
            index_t v1 = add_vertex(f2, AR1, AR2);
            index_t v2 = add_vertex(f2, BR1, BR2);
            edges_.push_back(std::make_pair(v1,v2));
        }

        void end_facet() {
            bool ok = true;

            Mesh debug_constraints;
            
            // Sanity check: the facet does not have zero area
            geo_assert(f1_orient_ != ZERO);

            // Sort vertices along triangle's edges
            {
                for(index_t e=0; e<3; ++e) {
                    index_t v_org = (e+1)%3; // e's origin
                    index_t v_opp = e;       // f1's vertex opposite to e
                    std::sort(
                        vertices_in_E_[e].begin(),
                        vertices_in_E_[e].end(),
                        [&](index_t v1, index_t v2)->bool{

                            Sign o12 = orient2d(v1,v2,v_opp,true);

                            // Supposed to be negative if
                            // v1 is between v_org and v2 (what we want)
                            Sign o_o12 = dot3d(v1,v_org,v2,true);
                            
                            if(o12 == ZERO || o_o12 == ZERO) {
                                log_err();
                                std::cerr << "Macro edge: " << e << std::endl;
                                std::cerr << "(found two coincident vertices)"
                                          << std::endl;
                                std::cerr <<vertex_[v1].to_string() <<std::endl;
                                std::cerr <<vertex_[v2].to_string() <<std::endl;
                                
                                get_constraints(debug_constraints);
                                ok = false;
                            }

                            bool result_1 = (o12 == f1_orient_);
                            bool result_2 = (o_o12 == NEGATIVE);

                            if(result_1 != result_2) {
                                log_err();
                                std::cerr << "Macro edge: " << e << std::endl;
                                std::cerr << "Different results in two tests: "
                                          << std::endl;
                                std::cerr << " 2d: " << o12
                                          << " ==? "  << f1_orient_
                                          << std::endl;
                                std::cerr << " 3d: " << o_o12
                                          << " negative?"
                                          << std::endl;

                                get_constraints(debug_constraints);
                                ok = false;
                            }
                            return result_1;
                        }
                    );
                }
            }

            if(!ok) {
                mesh_save(
                    debug_constraints,
                    "constraints_" + String::to_string(f1_) + ".geogram"
                );
            }
            
            // Create edges along f1's edges, in MeshInTriangle data structure
            {
                for(index_t e=0; e<3; ++e) {
                    index_t v0 = (e+1)%3;
                    index_t v1 = (e+2)%3;
                    index_t prev = v0;
                    for(index_t v: vertices_in_E_[e]) {
                        edges_.push_back(std::make_pair(prev, v));
                        prev = v;
                    }
                    edges_.push_back(std::make_pair(prev, v1));
                }
            }

            commit();
            clear();
        }

    protected:

        void commit() {

            bool OK_exact   = check_constraints(true);
            bool OK_inexact = check_constraints();
            
            // Create all vertices (or find them if they already exist)
            for(index_t i=0; i<vertex_.size(); ++i) {
                if(vertex_[i].is_existing_vertex()) {
                    vertex_[i].mesh_vertex_index = vertex_[i].sym.indices[0];
                } else {
                    auto it = v_table_.find(vertex_[i].sym);
                    if(it != v_table_.end()) {
                        vertex_[i].mesh_vertex_index = it->second;
                    } else {
                        vec3 p = vertex_[i].point;
                        index_t v = mesh_.vertices.create_vertex(p.data());
                        vertex_[i].mesh_vertex_index = v;
                        v_table_[vertex_[i].sym] = v;
                    }
                }
            }

            if(true) {
                // Create a 2D constrained Delaunay triangulation
                Mesh constraints;
                constraints.vertices.set_dimension(2);
                for(index_t i=0; i<vertex_.size(); ++i) {
                    constraints.vertices.create_vertex(
                        project(i).data()
                    );
                }
                for(std::pair<index_t,index_t>& E: edges_) {
                    constraints.edges.create_edge(E.first, E.second);
                }

                if(!OK_exact || !OK_inexact) {
                    std::cerr << "==============================>>>>"
                              << "There were errors, saving constraints..."
                              << std::endl;
                    mesh_save(
                        constraints,
                        "constraints_" + String::to_string(f1_) + ".geogram"
                    );
                    return;
                    // abort();
                }
                
                
                Delaunay_var del = Delaunay::create(2, "triangle");
                del->set_constraints(&constraints);
                del->set_vertices(0, nullptr);

                // There were intersections between the constraints,
                //    Triangle created steiner points
                // TODO: generate their correct symbolic info,
                //   or maybe generate them before, without letting
                //   Triangle do the job
                if(del->nb_vertices() > constraints.vertices.nb()) {
                    vec3 p0 = vertex_[0].point;
                    vec3 p1 = vertex_[1].point;
                    vec3 p2 = vertex_[2].point;
                    vec2 q0 = project(0);
                    vec2 q1 = project(1);
                    vec2 q2 = project(2);
                    double A = ::GEO::Geom::triangle_area(q0,q1,q2);
                    for(
                        index_t v = constraints.vertices.nb();
                        v < del->nb_vertices(); ++v
                    ) {
                        vec2 q(del->vertex_ptr(v));
                        double l0 = ::GEO::Geom::triangle_area(q,q1,q2);
                        double l1 = ::GEO::Geom::triangle_area(q0,q,q2);
                        double l2 = ::GEO::Geom::triangle_area(q0,q1,q);
                        vec3 p = (1.0/A)*(l0*p0 + l1*p1 + l2*p2);
                        index_t vv = mesh_.vertices.create_vertex(p.data());
                        Vertex V(mesh_,vv);
                        V.mesh_vertex_index = vv;
                        vertex_.push_back(V);
                    }
                }
                
                // If orientation in Delaunay and orientation in
                // facet do not match, we will need to flip
                // triangles.
                vec2 p1(del->vertex_ptr(index_t(del->cell_vertex(0,0))));
                vec2 p2(del->vertex_ptr(index_t(del->cell_vertex(0,1))));
                vec2 p3(del->vertex_ptr(index_t(del->cell_vertex(0,2))));
                
                bool flip = (PCK::orient_2d(p1,p2,p3) != f1_orient_);

                // For each triangle of the Delaunay triangulation,
                // create a triangle in the mesh (and flip if need be).
                for(index_t t=0; t<del->nb_cells(); ++t) {
                    index_t i = index_t(del->cell_vertex(t,0));
                    index_t j = index_t(del->cell_vertex(t,1));
                    index_t k = index_t(del->cell_vertex(t,2));
                    i = vertex_[i].mesh_vertex_index;
                    j = vertex_[j].mesh_vertex_index;
                    k = vertex_[k].mesh_vertex_index;                    
                    if(flip) {
                        std::swap(i,j);
                    }
                    mesh_.facets.create_triangle(i,j,k);
                }
            }
            
        }
        
        void get_constraints(Mesh& M, bool with_edges=true) const {
            if(M.vertices.nb() == 0) {
                M.vertices.set_dimension(2);
                for(index_t v=0; v<vertex_.size(); ++v) {
                    M.vertices.create_vertex(project(v).data());
                }
            }
            if(with_edges && M.edges.nb() == 0) {
                for(std::pair<index_t, index_t> E: edges_) {
                    M.edges.create_edge(E.first, E.second);
                }
            }
        }

        bool check_constraints(bool exact = false) const {
            Mesh debug_constraints;
            Attribute<bool> debug_vertex_show(
                debug_constraints.vertices.attributes(), "selection"
            );
            
            bool result = true;
            // Test duplicated vertices
            for(index_t i=0; i<vertex_.size(); ++i) {
                for(index_t j=i+1; j<vertex_.size(); ++j) {
                    if(same_point(i,j,exact)) {
                        result = false;
                        log_err(exact);
                        std::cerr << "Same point ("<<i<<","<<j<<")"<<std::endl;
                        std::cerr << "   " <<vertex_[i].to_string()<<std::endl;
                        std::cerr << "   " <<vertex_[j].to_string()<<std::endl;
                        get_constraints(debug_constraints, false);
                        debug_vertex_show[i] = true;
                        debug_vertex_show[j] = true;                        
                    }
                }
            }

            // Test vertices location, in triangle, on edge
            {
                vector<int> on_macro_edge(vertex_.size(), -1);
                for(index_t e=0; e<3; ++e) {
                    for(index_t v: vertices_in_E_[e]) {
                        on_macro_edge[v] = int(e);
                    }
                }
                for(index_t v=3; v<vertex_.size(); ++v) {
                    Sign o1 = orient2d(1,2,v,exact);
                    Sign o2 = orient2d(2,0,v,exact);
                    Sign o3 = orient2d(0,1,v,exact);
                    o1 = Sign(o1*f1_orient_);
                    o2 = Sign(o2*f1_orient_);
                    o3 = Sign(o3*f1_orient_);                    
                    
                    if(o1 == NEGATIVE || o2 == NEGATIVE || o3 == NEGATIVE) {
                        result = false;
                        log_err(exact);
                        std::cerr << "Point outside tri ("<<v<<")" <<std::endl;
                        std::cerr << o1 << " " << o2 << " " << o3  <<std::endl;
                        std::cerr << "on macro edge:" << on_macro_edge[v]
                                  << std::endl;
                        std::cerr << "   " <<vertex_[v].to_string()<<std::endl;
                        get_constraints(debug_constraints, false);
                        debug_vertex_show[v] = true;

                        const Vertex& V = vertex_[v];
                        
                        // Try to understand what's going on 
                        if(
                            (on_macro_edge[v] != -1) &&
                            (region_dim(V.tt_sym.R1) == 1) &&
                            (region_dim(V.tt_sym.R2) == 2)
                          ) {
                            index_t e = index_t(V.tt_sym.R1);
                            index_t w1 = mesh_.facets.vertex(
                                V.tt_sym.f1, (e+1)%3
                            );
                            index_t w2 = mesh_.facets.vertex(
                                V.tt_sym.f1, (e+2)%3
                            );
                            // index_t w3 = mesh_.facets.vertex(
                            //    V.tt_sym.f1, e
                            //);
                            index_t v1 = mesh_.facets.vertex(V.tt_sym.f2, 0);
                            index_t v2 = mesh_.facets.vertex(V.tt_sym.f2, 1);
                            index_t v3 = mesh_.facets.vertex(V.tt_sym.f2, 2);

                            vec3 q1(mesh_.vertices.point_ptr(w1));
                            vec3 q2(mesh_.vertices.point_ptr(w2));
                            // vec3 q3(mesh_.vertices.point_ptr(w3));

                            vec3 p1(mesh_.vertices.point_ptr(v1));
                            vec3 p2(mesh_.vertices.point_ptr(v2));
                            vec3 p3(mesh_.vertices.point_ptr(v3));

                            vec3Q P1 = convert_vec3_generic<rational_nt>(vertex_[0].point);
                            vec3Q P2 = convert_vec3_generic<rational_nt>(vertex_[1].point);
                            vec3Q P3 = convert_vec3_generic<rational_nt>(vertex_[2].point);
                            
                            vec3Q I =
                            get_segment_triangle_intersection<rational_nt>(
                                q1, q2, p1, p2, p3
                            );


                            Sign o1 = PCK::orient_2d_projected(
                                P2,P3,I,f1_normal_axis_
                            );

                            Sign o1_2 = PCK::orient_2d_projected(
                                P3,I,P2,f1_normal_axis_
                            );

                            Sign o1_3 = PCK::orient_2d_projected(
                                I,P2,P3,f1_normal_axis_
                            );

                            Sign o1_4 = PCK::orient_2d_projected(
                                P3,P2,I,f1_normal_axis_
                            );

                            Sign o1_5 = PCK::orient_2d_projected(
                                P2,I,P3,f1_normal_axis_
                            );

                            Sign o1_6 = PCK::orient_2d_projected(
                                I,P3,P2,f1_normal_axis_
                            );
                            
                            
                            Sign o2 = PCK::orient_2d_projected(
                                P3,P1,I,f1_normal_axis_
                            );
                            Sign o3 = PCK::orient_2d_projected(
                                P1,P2,I,f1_normal_axis_
                            );

                            o1 = Sign(o1*f1_orient_);
                            o2 = Sign(o2*f1_orient_);
                            o3 = Sign(o3*f1_orient_);                            
                            
                            std::cerr << "   Locate isect: "
                                      << o1 << " " << o2 << " " << o3
                                      << std::endl;

                            std::cerr << "   o1..........: "
                                      << o1   << " " << o1_2 << " " << o1_3 << " "
                                      << o1_4 << " " << o1_5 << " " << o1_6                                
                                      << std::endl;
                            

                            if(
                                o1 != o1_2 ||
                                o1 != o1_3 ||
                                int(o1) != -int(o1_4) ||
                                int(o1) != -int(o1_5) ||
                                int(o1) != -int(o1_6)
                            ) {
                                std::cerr << "=========> Incoherent o1! (rational_nt is broken ?)"
                                          << std::endl;
                                // This *should* not happen in exact mode
                                // (except when there is an *underflow* in expansion code...)
                                if(exact) {
                                    abort();
                                }
                            }

                            
                            Sign oo1 = orient2d(v,1,2,true);
                            Sign oo2 = orient2d(0,v,2,true);
                            Sign oo3 = orient2d(0,1,v,true);
                            oo1 = Sign(oo1*f1_orient_);
                            oo2 = Sign(oo2*f1_orient_);
                            oo3 = Sign(oo3*f1_orient_);                    


                            std::cerr << "Re-Locate isect: "
                                      << oo1 << " " << oo2 << " " << oo3
                                      << std::endl;
                            
                        }
                        
                    }
                    int nb_zeros = (o1==0)+(o2==0)+(o3==0);
                    geo_assert(nb_zeros < 2);

                    if(nb_zeros == 1 && (on_macro_edge[v] == -1)) {
                        result = false;
                        log_err(exact);
                        std::cerr << "Inner point on macro edge (" << v << ")"
                                  << std::endl;
                        std::cerr << "   " <<vertex_[v].to_string()<<std::endl;
                        get_constraints(debug_constraints, false);
                        debug_vertex_show[v] = true;
                    }
                }
            }

            
            // Test vertices on constrained edges
            for(index_t i=0; i<vertex_.size(); ++i) {
                for(std::pair<index_t,index_t> E: edges_) {
                    if(
                        i != E.first && i != E.second && 
                        vertex_on_edge(i, E.first, E.second, exact)
                    ) {
                        result = false;
                        log_err(exact);
                        std::cerr << "Point on edge !" << std::endl;
                        std::cerr << "    Pt:"<< vertex_[i].to_string()
                                  << std::endl;
                        std::cerr << "   Edge:"<< vertex_[E.first].to_string()
                                  << " --- "
                                  << vertex_[E.second].to_string() << std::endl;
                        get_constraints(debug_constraints, false);
                        debug_vertex_show[i] = true;
                        debug_vertex_show[E.first] = true;
                        debug_vertex_show[E.second] = true;
                        debug_constraints.edges.create_edge(E.first, E.second);
                    }
                }
            }
            
            // Test constrained edges intersetions
            for(index_t e1=0; e1<edges_.size(); ++e1) {
                std::pair<index_t,index_t> E1 = edges_[e1];                
                for(index_t e2=e1+1; e2 < edges_.size(); ++e2) {
                    std::pair<index_t,index_t> E2 = edges_[e2];
                    if(
                        E1.first == E2.first ||
                        E1.first == E2.second ||
                        E1.second == E2.first ||
                        E1.second == E2.second
                    ) {
                        continue;
                    }
                    if(edge_edge_intersect(
                           E1.first, E1.second, E2.first, E2.second, exact
                    )) {
                        result = false;
                        log_err(exact);
                        std::cerr << "Intersecting edges !" << std::endl;
                        std::cerr << "   Edge1: " 
                                  << vertex_[E1.first].to_string()
                                  << " --- "
                                  << vertex_[E1.second].to_string()
                                  << std::endl;
                        std::cerr << "   Edge2: " 
                                  << vertex_[E2.first].to_string()
                                  << " --- "
                                  << vertex_[E2.second].to_string()
                                  << std::endl;
                        get_constraints(debug_constraints, false);
                        debug_constraints.edges.create_edge(
                            E1.first, E1.second
                        );
                        debug_constraints.edges.create_edge(
                            E2.first, E2.second
                        );
                        debug_vertex_show[E1.first] = true;
                        debug_vertex_show[E1.second] = true;
                        debug_vertex_show[E2.first] = true;
                        debug_vertex_show[E2.second] = true;


                        if(exact) {
                            Mesh debug_show;
                            debug_show.vertices.set_dimension(2);
                            vec2 p;
                            p.x = vertex_[E1.first].point_exact[u_].estimate();
                            p.y = vertex_[E1.first].point_exact[v_].estimate();
                            debug_show.vertices.create_vertex(p.data());
                            p.x = vertex_[E1.second].point_exact[u_].estimate();
                            p.y = vertex_[E1.second].point_exact[v_].estimate();
                            debug_show.vertices.create_vertex(p.data());
                            p.x = vertex_[E2.first].point_exact[u_].estimate();
                            p.y = vertex_[E2.first].point_exact[v_].estimate();
                            debug_show.vertices.create_vertex(p.data());
                            p.x = vertex_[E2.second].point_exact[u_].estimate();
                            p.y = vertex_[E2.second].point_exact[v_].estimate();
                            debug_show.vertices.create_vertex(p.data());
                            debug_show.edges.create_edge(0,1);
                            debug_show.edges.create_edge(2,3);
                            mesh_save(debug_show,"debug_show_exact.geogram");
                        }
                    }
                }
            }
            
            if(!result) {
                mesh_save(
                    debug_constraints,
                    "debug_constraints_" + String::to_string(f1_) +
                    (exact ? "_exact" : "_inexact") + ".geogram"
                );
            }

            if(!result) {
                for(index_t e=0; e<3; ++e) {
                    std::cerr << "MacroEdge " << e << " : ";
                    for(index_t v: vertices_in_E_[e]) {
                        std::cerr << v << " ";
                    }
                    std::cerr << std::endl;
                }
            }
            
            return result;
        }

        bool same_point(index_t i, index_t j, bool exact=false) const {
            if(exact) {
                const vec3Q& p1 = vertex_[i].point_exact;
                const vec3Q& p2 = vertex_[j].point_exact;
                return (p1[0] == p2[0]) &&
                       (p1[1] == p2[1]) &&
                       (p1[2] == p2[2]) ;
            }
            
            const vec3& p1 = vertex_[i].point;
            const vec3& p2 = vertex_[j].point;            
            return (p1[0] == p2[0]) &&
                   (p1[1] == p2[1]) &&
                   (p1[2] == p2[2]) ;
        }

        bool vertex_on_edge(
            index_t i, index_t j, index_t k, bool exact=false
        ) const {
            if(orient2d(i,j,k,exact) != ZERO) {
                return false;
            }
            Sign o = dot3d(i,j,k,exact);
            return (int(o) <= 0);
        }

        bool edge_edge_intersect(
            index_t i, index_t j,
            index_t k, index_t l,
            bool exact = false
        ) const {
            
            Sign o1 = orient2d(i,j,k,exact);
            Sign o2 = orient2d(i,j,l,exact);

            if(o1*o2 > 0) {
                return false;
            }

            // Particular case: 1D
            if(o1 == 0 && o2 == 0) {
                Sign d1 = dot3d(k,i,j,exact);
                Sign d2 = dot3d(l,i,j,exact);
                return (d1 <= 0 || d2 <= 0);
            }
            
            Sign o3 = orient2d(k,l,i,exact);
            Sign o4 = orient2d(k,l,j,exact);
            return (o3*o4 <= 0);
        } 
        

        Sign orient2d(index_t v1,index_t v2,index_t v3,bool exact=false) const {
            if(exact) {
                return PCK::orient_2d_projected(
                    vertex_[v1].point_exact,
                    vertex_[v2].point_exact,
                    vertex_[v3].point_exact,
                    f1_normal_axis_
                );
            }
            return PCK::orient_2d(project(v1), project(v2), project(v3));
        }

        Sign dot3d(index_t v1, index_t v2, index_t v3, bool exact=false) const {
            if(exact) {
                return PCK::dot_3d(
                    vertex_[v1].point_exact,
                    vertex_[v2].point_exact,
                    vertex_[v3].point_exact                    
                );
            }
            return PCK::dot_3d(
                vertex_[v1].point,
                vertex_[v2].point,
                vertex_[v3].point
            );
        }
        
        vec2 project(vec3 p) const {
            return vec2(p[u_], p[v_]);
        }
        
        vec2 project(index_t v) const {
            return project(vertex_[v].point);
        }

        index_t add_vertex(const Vertex& V) {
            for(index_t i=0; i<vertex_.size(); ++i) {
                if(vertex_[i].sym == V.sym) {
                    return i;
                }
            }
            vertex_.push_back(V);
            return vertex_.size()-1;
        }
        
        void clear() {
            vertex_.resize(0);
            edges_.resize(0);
            for(index_t e=0; e<3; ++e) {
                vertices_in_E_[e].resize(0);
            }
            f1_ = index_t(-1);
        }

        void show_intersected_facets() const {
            std::vector<index_t> isect_facets;
            for(const Vertex& V: vertex_) {
                if(V.tt_sym.f2 != index_t(-1)) {
                    isect_facets.push_back(V.tt_sym.f2);
                }
            }
            sort_unique(isect_facets);
            std::cerr << "isect facets: ";
            std::cerr << f1_ << ";";
            for(index_t f: isect_facets) {
                std::cerr << f << ";";
            }
            std::cerr << std::endl;
        }
        
        void log_err(bool exact=false) const {
            if(exact) {
                std::cerr << "[Exact] ";
            } else {
                std::cerr << "[Inexact] ";                
            }
            std::cerr << "Houston, we got a problem (while remeshing facet "
                      << f1_ << "):" << std::endl;
        }
        
    private:
        Mesh& mesh_;
        index_t f1_;
        coord_index_t f1_normal_axis_;
        coord_index_t u_; // = (f1_normal_axis_ + 1)%3
        coord_index_t v_; // = (f1_normal_axis_ + 2)%3
        Sign f1_orient_;
        vector<Vertex> vertex_;
        vector<std::pair<index_t, index_t> > edges_;
        vector<index_t> vertices_in_E_[3];
        std::map<quadindex, index_t> v_table_;
    };

    /***********************************************************************/
    
    /**
     * \brief Stores information about a triangle-triangle intersection.
     * \details The intersection is a segment. Its extremities are indicated
     *  by the regions in f1 and f2 that created the intersection. If the
     *  intersection is just a point, then A and B regions are the same.
     */
    struct IsectInfo {
    public:

        /**
         * Swaps the two facets and updates the combinatorial
         * information accordingly.
         */
        void flip() {
            std::swap(f1,f2);
            A_rgn_f1 = swap_T1_T2(A_rgn_f1);
            A_rgn_f2 = swap_T1_T2(A_rgn_f2);
            std::swap(A_rgn_f1, A_rgn_f2);
            B_rgn_f1 = swap_T1_T2(B_rgn_f1);
            B_rgn_f2 = swap_T1_T2(B_rgn_f2);
            std::swap(B_rgn_f1, B_rgn_f2);
        }

        /**
         * \brief Tests whether intersection is just a point.
         * \details Points are encoded as segments with the
         *  same symbolic information for both vertices.
         */
        bool is_point() const {
            return
                A_rgn_f1 == B_rgn_f1 &&
                A_rgn_f2 == B_rgn_f2 ;
        }
        
        index_t f1; 
        index_t f2;
        TriangleRegion A_rgn_f1;
        TriangleRegion A_rgn_f2;
        TriangleRegion B_rgn_f1;
        TriangleRegion B_rgn_f2;
    };

    bool mesh_facets_intersect(
        Mesh& M, index_t f1, index_t f2, vector<TriangleIsect>& I
    ) {
        // cerr << "Intersect facets " << f1 << " /\\ " << f2 << std::endl;
        geo_debug_assert(M.facets.nb_vertices(f1) == 3);
        geo_debug_assert(M.facets.nb_vertices(f2) == 3);        
        vec3 p1(M.vertices.point_ptr(M.facets.vertex(f1,0)));
        vec3 p2(M.vertices.point_ptr(M.facets.vertex(f1,1)));
        vec3 p3(M.vertices.point_ptr(M.facets.vertex(f1,2)));
        vec3 q1(M.vertices.point_ptr(M.facets.vertex(f2,0)));
        vec3 q2(M.vertices.point_ptr(M.facets.vertex(f2,1)));
        vec3 q3(M.vertices.point_ptr(M.facets.vertex(f2,2)));
        return triangles_intersections(p1,p2,p3,q1,q2,q3,I);
    }
    
    void mesh_intersections(
        MeshInTriangle& TM, vector<IsectInfo>& intersections
    ) {
        // Sort intersections by f1, so that all intersections between f1
        // and another facet appear as a contiguous sequence.
        std::sort(
            intersections.begin(), intersections.end(),
            [](const IsectInfo& a, const IsectInfo& b) -> bool {
                return (a.f1 < b.f1) ? true  :
                       (a.f1 > b.f1) ? false :
                       (a.f2 < b.f2) ;
            }
        );

        index_t b=0;
        while(b < intersections.size()) {
            index_t e = b;
            while(
                e < intersections.size() &&
                intersections[e].f1 == intersections[b].f1
            ) {
                ++e;
            }
            
            TM.begin_facet(intersections[b].f1);
            
            for(index_t i=b; i<e; ++i) {

                const IsectInfo& II = intersections[i];
                
                // If intersection is a single point that
                // already exist in f1, skip it.
                if(II.is_point() && region_dim(II.A_rgn_f1) == 0) {
                    continue;
                }

                if(II.is_point()) {
                    TM.add_vertex(
                        II.f2,
                        II.A_rgn_f1, II.A_rgn_f2
                    );
                } else {
                    // If both ends in same edge of f1,
                    // then add individual vertices (because
                    // edges will be created by sorting
                    // vertices along triangle edges).
                    if(
                        II.A_rgn_f1 == II.B_rgn_f1 &&
                        region_dim(II.A_rgn_f1) == 1
                    ) {
                        TM.add_vertex(
                            II.f2,
                            II.A_rgn_f1, II.A_rgn_f2                            
                        );
                        TM.add_vertex(
                            II.f2,
                            II.B_rgn_f1, II.B_rgn_f2                            
                        );
                    } else {
                        TM.add_edge(
                            II.f2, 
                            II.A_rgn_f1, II.A_rgn_f2,
                            II.B_rgn_f1, II.B_rgn_f2
                        );
                    }
                }
            }
            TM.end_facet();
            b = e;
        }
        
    }

    bool same_region(TriangleRegion R1, TriangleRegion R2) {
        if(R1 == R2) {
            return true;
        }
        if(region_dim(R1) == 1 && region_dim(R2) == 0) {
            TriangleRegion v1,v2;
            get_edge_vertices(R1,v1,v2);
            if(R2 == v1 || R2 == v2) {
                return true;
            }
        }
        if(region_dim(R2) == 1 && region_dim(R1) == 0) {
            TriangleRegion v1,v2;
            get_edge_vertices(R2,v1,v2);
            if(R1 == v1 || R2 == v2) {
                return true;
            }
        }
        return false;
    }
    
    bool planar_intersection_is_valid(const IsectInfo& II) {
        return
            same_region(II.A_rgn_f1, II.B_rgn_f1) ||
            same_region(II.A_rgn_f2, II.B_rgn_f2) ;
    }
    
    void intersect_surface_impl(
        Mesh& M, bool test_adjacent_facets, bool order_facets
    ) {

        // Set symbolic perturbation mode to lexicographic order
        // on point coordinates instead of point indices only,
        // Needed to get compatible triangulations on coplanar faces
        // (example, cubes that touch on a facet).
        PCK::SOSMode SOS_bkp = PCK::get_SOS_mode();
        PCK::set_SOS_mode(PCK::SOS_LEXICO);
        
        vector<IsectInfo> intersections;
        MeshFacetsAABB AABB(M,order_facets);
        vector<TriangleIsect> I;
        AABB.compute_facet_bbox_intersections(
            [&](index_t f1, index_t f2) {

                if(f1 == f2) {
                    return;
                }
                
                if(
                    !test_adjacent_facets && (
                        (M.facets.find_adjacent(f1,f2)      != index_t(-1)) ||
                        (M.facets.find_common_vertex(f1,f2) != index_t(-1))
                    )
                ) {
                    return;
                }
                
                if(mesh_facets_intersect(M, f1, f2, I)) {

                    // Coplanar intersection, find all segments
                    if(I.size() > 2) {
                        for(index_t i1=0; i1< I.size(); ++i1) {
                            for(index_t i2=0; i2<i1; ++i2) {
                                IsectInfo II = {
                                    f1, f2,
                                    I[i1].first, I[i1].second,
                                    I[i2].first, I[i2].second
                                };
                                if(planar_intersection_is_valid(II)) {
                                    intersections.push_back(II);
                                    II.flip();
                                    intersections.push_back(II);
                                } 
                            }
                        }
                        return;
                    }

                    geo_debug_assert(I.size() != 0);
                    
                    TriangleRegion A_rgn_f1 = I[0].first;
                    TriangleRegion A_rgn_f2 = I[0].second;

                    TriangleRegion B_rgn_f1 = A_rgn_f1;
                    TriangleRegion B_rgn_f2 = A_rgn_f2;

                    if(I.size() == 2) {
                        B_rgn_f1 = I[1].first;
                        B_rgn_f2 = I[1].second;
                    }
                    
                    IsectInfo II = {
                        f1, f2,
                        A_rgn_f1, A_rgn_f2,
                        B_rgn_f1, B_rgn_f2
                    };
                    intersections.push_back(II);
                    II.flip();
                    intersections.push_back(II);                   
                }
            }
        );

        MeshInTriangle TM(M);
        mesh_intersections(TM, intersections);
        
        vector<index_t> has_intersections(M.facets.nb(), 0);
        for(const IsectInfo& II: intersections) {
            has_intersections[II.f1] = 1;
            has_intersections[II.f2] = 1;
        }
        
        M.facets.delete_elements(has_intersections);
        M.facets.connect();

        PCK::set_SOS_mode(SOS_bkp);
    }
}

namespace OGF {

    MeshGrobExperimentCommands::MeshGrobExperimentCommands() { 
    }
        
    MeshGrobExperimentCommands::~MeshGrobExperimentCommands() { 
    }        

    void MeshGrobExperimentCommands::intersect_surface(
        bool merge_vertices_and_facets,
        bool test_neighboring_triangles,
        bool post_connect_facets,
        bool order_facets
    ) {
        if(!mesh_grob()->facets.are_simplices()) {
            Logger::err("Intersect") << "Mesh is not triangulated"
                                     << std::endl;
            return;
        }

        if(merge_vertices_and_facets) {
            mesh_colocate_vertices_no_check(*mesh_grob());
            mesh_remove_bad_facets_no_check(*mesh_grob());
        }
         
        intersect_surface_impl(
            *mesh_grob(), test_neighboring_triangles, order_facets
        );

        if(post_connect_facets) {
            mesh_repair(
                *mesh_grob(),
                GEO::MeshRepairMode(
                    GEO::MESH_REPAIR_COLOCATE | GEO::MESH_REPAIR_DUP_F
                ),
                0.0
            );
        }
        
        
        show_mesh();
        mesh_grob()->update();
    }

    void MeshGrobExperimentCommands::sort_facets() {
        MeshFacetsAABB AABB(*mesh_grob()); // This sorts the facets
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
