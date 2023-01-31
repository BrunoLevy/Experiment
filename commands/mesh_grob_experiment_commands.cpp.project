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
            enum Type {
                UNINITIALIZED, MESH_VERTEX, PRIMARY_ISECT, SECONDARY_ISECT
            };

            Vertex(
                MeshInTriangle* M,
                index_t f1, index_t f2,
                TriangleRegion R1, TriangleRegion R2
            ) {
                geo_assert(f1 == M->f1_);                
                type = PRIMARY_ISECT;
                mesh_in_triangle = M;
                init_sym(f1,f2,R1,R2);
                init_geometry();
            }

            Vertex(MeshInTriangle* M, index_t f, index_t lv) {
                geo_assert(f == M->f1_);
                type = MESH_VERTEX;
                mesh_in_triangle = M;
                init_sym(f, NO_INDEX, TriangleRegion(lv), T2_RGN_T);
                init_geometry();
            }

            Vertex(
                MeshInTriangle* M,
                const vec3Q& point_exact_in,
                const vec2Q& UV_exact_in
            ) : point_exact(point_exact_in),
                UV_exact(UV_exact_in)
            {
                type = SECONDARY_ISECT;                
                mesh_in_triangle = M;
                init_sym(NO_INDEX, NO_INDEX, T1_RGN_T, T2_RGN_T);
            }
            
            Vertex() {
                type = UNINITIALIZED;                
                mesh_in_triangle = nullptr;
                init_sym(NO_INDEX, NO_INDEX, T1_RGN_T, T2_RGN_T);
                mesh_vertex_index = NO_INDEX;
            }

            Mesh& mesh() {
                return mesh_in_triangle->mesh();
            }
            
            void print(std::ostream& out=std::cerr) const {
                if(sym.f1 != index_t(-1)) {
                    out << " ( ";
                    out << sym.f1;
                    out << region_to_string(sym.R1).substr(2);
                }
                if(sym.f2 != index_t(-1)) {
                    out << " /\\ ";
                    out << sym.f2;
                    out << region_to_string(sym.R2).substr(2);
                }
                if(sym.f1 != index_t(-1)) {
                    out << " ) ";
                }
            }

            std::string to_string() const {
                std::ostringstream out;
                print(out);
                return out.str();
            }

            bool is_on_facet(index_t f) const {
                for(index_t g: facets) {
                    if(f == g) {
                        return true;
                    }
                }
                return false;
            }
            
        protected:
            
            vec3 mesh_facet_vertex(index_t f, index_t lv) const {
                return mesh_in_triangle->mesh_facet_vertex(f,lv);
            }
            
            void init_sym(
                index_t f1, index_t f2,
                TriangleRegion R1, TriangleRegion R2
            ) {
                sym.f1 = f1;
                sym.f2 = f2;
                sym.R1 = R1;
                sym.R2 = R2;
                mesh_vertex_index = NO_INDEX;
                if(region_dim(sym.R1) == 0) {
                    mesh_vertex_index = mesh().facets.vertex(
                        sym.f1, index_t(sym.R1)
                    );
                }
                if(region_dim(sym.R2) == 0) {
                    mesh_vertex_index = mesh().facets.vertex(
                        sym.f2, index_t(sym.R2)-3
                    );
                }
                if(f2 != NO_INDEX) {
                    facets.push_back(f2);
                }
            }

            void init_geometry() {

                if(mesh_vertex_index != NO_INDEX) {
                    vec3 point(mesh().vertices.point_ptr(mesh_vertex_index));
                    point_exact = convert_vec3_generic<rational_nt>(point);
                    UV_exact = vec2Q(
                        point_exact[mesh_in_triangle->u_],
                        point_exact[mesh_in_triangle->v_]
                    );
                    return;
                }
                
                geo_assert(sym.f1 != NO_INDEX && sym.f2 != NO_INDEX);
                
                // The three vertices of f1. If intersection is on an edge
                // of f1, then [p1,p2] are the vertices of that edge.
                vec3 p1,p2,p3;
                get_vertices(sym.R1,p1,p2,p3);

                // The three vertices of f2. If intersection is on an edge
                // of f2, then [q1,q2] are the vertices of that edge.
                vec3 q1,q2,q3;
                get_vertices(sym.R2,q1,q2,q3);                

                // Segment-segment intersection
                if(region_dim(sym.R1) == 1 && region_dim(sym.R2) == 1) {

                    // We distinguish 3D edge /\ edge isect (then we use
                    // edge /\ facet isect) and 2D edge /\ edge isect (
                    // then develop special code). 
                    
                    geo_debug_assert(sym.f2 != NO_INDEX);

                    // If [q1,q2] is in the supporting plane of (p1,p2,p3),
                    // compute intersection in 2D.
                    if(
                        PCK::orient_3d(p1,p2,p3,q1) == ZERO &&
                        PCK::orient_3d(p1,p2,p3,q2) == ZERO
                    ) {
                        point_exact =
                            get_segment_segment_intersection_2D<rational_nt>(
                                q1,q2,p1,p2,
                                ::GEO::Geom::triangle_normal_axis(p1,p2,p3)
                            );
                    } else {
                        // We are in 3D, we can use segment /\ triangle
                        // intersection (even if we know that the intersection
                        // will be on on (p1,p2))
                        point_exact =
                            get_segment_triangle_intersection<rational_nt>(
                                q1,q2,p1,p2,p3
                            );
                    }
                } else {
                    // Edge-triangle intersection
                    geo_assert(
                        (region_dim(sym.R1) == 1 && region_dim(sym.R2) == 2) ||
                        (region_dim(sym.R1) == 2 && region_dim(sym.R2) == 1)
                    );
                
                    if(region_dim(sym.R1) == 1) {
                        point_exact =
                            get_segment_triangle_intersection<rational_nt>(
                                p1,p2,q1,q2,q3
                            );
                    } else {
                        point_exact =
                            get_segment_triangle_intersection<rational_nt>(
                                q1,q2,p1,p2,p3
                            ) ;
                    }
                }

                UV_exact = vec2Q(
                    point_exact[mesh_in_triangle->u_],
                    point_exact[mesh_in_triangle->v_]
                );
            }

            /**
             * \brief Gets the three vertices of a facet from the
             *  TriangleRegion
             * \details If the TriangleRegion is one of the edges,
             *  gets the vertices of this edge as the first two 
             *  vertices \p p1 and \p p2
             */
            void get_vertices(TriangleRegion R, vec3& p1, vec3& p2, vec3& p3) {
                geo_debug_assert(region_dim(R) == 1 || region_dim(R) == 2);
                index_t f = is_in_T1(R) ? sym.f1 : sym.f2;
                index_t lv1 = 0;
                index_t lv2 = 1;
                index_t lv3 = 2;
                if(region_dim(R) == 1) {
                    index_t e = is_in_T1(R) ? index_t(R)-6 : index_t(R)-9;
                    lv1 = (lv1+e+1)%3;
                    lv2 = (lv2+e+1)%3;
                    lv3 = (lv3+e+1)%3;                    
                }
                p1 = mesh_facet_vertex(f,lv1);
                p2 = mesh_facet_vertex(f,lv2);
                p3 = mesh_facet_vertex(f,lv3);                
            }
            
        public:
            MeshInTriangle* mesh_in_triangle;

            vec3Q point_exact;
            vec2Q UV_exact;

            Type type;
            
            // Global mesh vertex index once created
            index_t mesh_vertex_index;

            // Symbolic information - triangle-triangle isect
            struct {
                index_t f1,f2;        // global facet indices in mesh
                TriangleRegion R1,R2; // triangle regions
            } sym;
            
            vector<index_t> facets;   // all the facets that touch this vertex
        };

        
        MeshInTriangle(Mesh& M, bool check_cnstr) :
            mesh_(M),
            f1_(index_t(-1)),
            check_cnstr_(check_cnstr)
        {
        }

        Mesh& mesh() {
            return mesh_;
        }
        
        void begin_facet(index_t f) {
            f1_ = f;
            f1_normal_axis_ = ::GEO::Geom::triangle_normal_axis(
                mesh_facet_vertex(f,0),
                mesh_facet_vertex(f,1),
                mesh_facet_vertex(f,2)                
            );
            u_ = coord_index_t((f1_normal_axis_ + 1) % 3);
            v_ = coord_index_t((f1_normal_axis_ + 2) % 3);
            for(index_t lv=0; lv<3; ++lv) {
                add_vertex(Vertex(this, f, lv));
            }
            f1_orient_ = PCK::orient_2d(project(0), project(1), project(2));
        }
        
        index_t add_vertex(
            index_t f2,
            TriangleRegion R1, TriangleRegion R2
        ) {
            geo_debug_assert(f1_ != index_t(-1));
            Vertex V(this, f1_, f2, R1, R2);
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
            // Sanity check: the facet does not have zero area
            geo_assert(f1_orient_ != ZERO);

            compute_constraints_intersections();
            
            // Sort vertices along macro-triangle's edges
            for(index_t e=0; e<3; ++e) {
                std::pair<index_t, index_t> E = std::make_pair(
                    (e+1)%3, (e+2)%3
                );
                sort_vertices_along_edge(E, vertices_in_E_[e]);
            }

            // Create edges along f1's edges, in MeshInTriangle data structure
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

            commit();
            clear();
        }

    protected:

        index_t add_vertex(const Vertex& V) {
            // Test whether there is already a vertex at the same position
            for(index_t i=0; i<vertex_.size(); ++i) {
                if(PCK::same_point(vertex_[i].UV_exact,V.UV_exact)) {
                    // If the vertex already exist, update its list of
                    // incident facets
                    if(
                        V.sym.f2 != index_t(-1) &&
                        !vertex_[i].is_on_facet(V.sym.f2)) {
                        vertex_[i].facets.push_back(V.sym.f2);
                    }
                    // Return found vertex
                    return i;
                }
            }
            // There was no vertex at the same position, so add it
            // to the list of vertices
            vertex_.push_back(V);
            return vertex_.size()-1;
        }
        
        void commit() {

            // Sanity check
            bool OK = check_cnstr_ ? check_constraints() : true;
            
            // Create all vertices (or find them if they already exist)
            for(index_t i=0; i<vertex_.size(); ++i) {

                // Vertex already exists in this MeshInTriangle
                if(vertex_[i].mesh_vertex_index != index_t(-1)) {
                    continue;
                }

                // Use exact geometry as key
                const vec3Q& K = vertex_[i].point_exact;
                auto it = g_v_table_.find(K);
                if(it != g_v_table_.end()) {
                    // Vertex alreay exists in target mesh
                    vertex_[i].mesh_vertex_index = it->second;
                } else {
                    // Vertex does not exist in target mesh,
                    // create it and update table
                    vec3 p(
                        vertex_[i].point_exact.x.estimate(),
                        vertex_[i].point_exact.y.estimate(),
                        vertex_[i].point_exact.z.estimate()
                    );
                    index_t v = mesh_.vertices.create_vertex(p.data());
                    vertex_[i].mesh_vertex_index = v;
                    g_v_table_[K] = v;
                }
            }

            // Create a 2D constrained Delaunay triangulation
            Mesh constraints;
            get_constraints(constraints);

            if(!OK) {
                std::cerr << "==============================>>>>"
                          << "There were errors, saving constraints..."
                          << std::endl;
                mesh_save(
                    constraints,
                    "constraints_" + String::to_string(f1_) + ".geogram"
                );
                abort();
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
                geo_assert_not_reached; // for now...
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

        static void print(const rational_nt& x) {
            std::cerr << "num: " << x.num().rep().length() << " : ";
            for(index_t i=0; i<x.num().rep().length(); ++i) {
                std::cerr << x.num().rep()[i] << " ";
            }
            std::cerr << std::endl << std::endl ;
            std::cerr << "den: " << x.denom().rep().length() << " : ";
            
            for(index_t i=0; i<x.denom().rep().length(); ++i) {
                std::cerr << x.denom().rep()[i] << " ";
            }
            std::cerr << std::endl << std::endl << std::endl;
        }

        index_t common_facet(index_t v1, index_t v2) const {
            const Vertex& V1 = vertex_[v1];
            const Vertex& V2 = vertex_[v2];
            index_t result = index_t(-1);
            index_t nb = 0;
            for(index_t f: V1.facets) {
                if(V2.is_on_facet(f)) {
                    result = f;
                    ++nb;
                }
            }
            if(nb > 1) {
                std::cerr
                    << "====> more than one facet common to edge extremities"
                    << std::endl;
            }
            return result;
        }
        
        void compute_constraints_intersections() {

            // I don't like to use vectors of vectors, but well...
            vector<vector<index_t> > new_vertices_in_edge(edges_.size());

            // Step 1: Compute all intersections
            for(index_t e1=0; e1<edges_.size(); ++e1) {
                index_t v1 = edges_[e1].first;
                index_t v2 = edges_[e1].second;
                for(index_t e2=e1+1; e2<edges_.size(); ++e2) {
                    index_t w1 = edges_[e2].first;
                    index_t w2 = edges_[e2].second;
                    if(v1 == w1 || v1 == w2 || v2 == w1 || v2 == w2) {
                        // NOTE/TODO: we could have two edges
                        // coming from two
                        // different triangles and having an intersection
                        // along an edge with an existing vertex.
                        continue;
                    }
                    // TODO: check whether intersection is one of the
                    // vertices.
                    // TODO: check co-linear edges

                    if(edge_edge_intersect(v1,v2,w1,w2)) {

                        index_t f1 = f1_;
                        index_t f2 = common_facet(v1,v2);
                        index_t f3 = common_facet(w1,w2);

                        geo_assert(f1 != index_t(-1));
                        geo_assert(f2 != index_t(-1));
                        geo_assert(f3 != index_t(-1));                        

                        vec3 P[9] = {
                            mesh_facet_vertex(f1,0),
                            mesh_facet_vertex(f1,1),
                            mesh_facet_vertex(f1,2),
                            mesh_facet_vertex(f2,0),
                            mesh_facet_vertex(f2,1),
                            mesh_facet_vertex(f2,2),
                            mesh_facet_vertex(f3,0),
                            mesh_facet_vertex(f3,1),
                            mesh_facet_vertex(f3,2)                            
                        };

                        vec3Q I;
                        if(
                          !get_three_planes_intersection(
                            I,
                            P[0], P[1], P[2],
                            P[3], P[4], P[5],
                            P[6], P[7], P[8]
                           )
                        ) {
                            std::cerr << " 2D cnstr isect" << std::endl;
                            const vec3Q& P1 = vertex_[v1].point_exact;
                            const vec3Q& P2 = vertex_[v2].point_exact;
                            const vec3Q& Q1 = vertex_[w1].point_exact;
                            const vec3Q& Q2 = vertex_[w2].point_exact;
                            if(
                                !get_segment_segment_intersection_2D_bis(
                                    I,
                                    P1,P2,Q1,Q2,
                                    f1_normal_axis_
                                )
                            ) {
                                std::cerr << " 1D cnstr isect" << std::endl;
                                geo_assert_not_reached;
                            }
                        }
                        
                        vec2Q I_UV(I[u_],I[v_]);
                        index_t x = add_vertex(Vertex(this,I,I_UV));
                        vertex_[x].facets.push_back(f2);
                        vertex_[x].facets.push_back(f3);

                        new_vertices_in_edge[e1].push_back(x);
                        new_vertices_in_edge[e2].push_back(x);
                    }
                }
            }

            // Step 2: sort intersections in edges
            for(index_t e=0; e<edges_.size(); ++e) {
                sort_vertices_along_edge(edges_[e], new_vertices_in_edge[e]);
            }

            // Step 3: "remesh" the edges that have intersections
            for(index_t e=0; e<new_vertices_in_edge.size(); ++e) {
                index_t prev_v = edges_[e].first;
                index_t last_v = edges_[e].second;
                new_vertices_in_edge[e].push_back(last_v);
                for(index_t i=0; i<new_vertices_in_edge[e].size(); ++i) {
                    index_t cur_v = new_vertices_in_edge[e][i];
                    if(i == 0) {
                        edges_[e].second = cur_v;
                    } else {
                        edges_.push_back(std::make_pair(prev_v, cur_v));
                    }
                    prev_v = cur_v;
                }
            }
            
        }
        
        bool check_constraints() const {
            Mesh debug_constraints;
            Attribute<bool> debug_vertex_show(
                debug_constraints.vertices.attributes(), "selection"
            );
            
            bool result = true;
            // Test duplicated vertices
            for(index_t i=0; i<vertex_.size(); ++i) {
                for(index_t j=i+1; j<vertex_.size(); ++j) {
                    if(same_point(i,j)) {
                        result = false;
                        log_err();
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
                    Sign o1 = orient2d(1,2,v);
                    Sign o2 = orient2d(2,0,v);
                    Sign o3 = orient2d(0,1,v);
                    o1 = Sign(o1*f1_orient_);
                    o2 = Sign(o2*f1_orient_);
                    o3 = Sign(o3*f1_orient_);                    
                    if(o1 == NEGATIVE || o2 == NEGATIVE || o3 == NEGATIVE) {
                        result = false;
                        log_err();
                        std::cerr << "Point outside tri ("<<v<<")" <<std::endl;
                        std::cerr << o1 << " " << o2 << " " << o3  <<std::endl;
                        std::cerr << "on macro edge:" << on_macro_edge[v]
                                  << std::endl;
                        std::cerr << "   " <<vertex_[v].to_string()<<std::endl;
                        get_constraints(debug_constraints, false);
                        debug_vertex_show[v] = true;
                    }
                    int nb_zeros = (o1==0)+(o2==0)+(o3==0);
                    geo_assert(nb_zeros < 2);
                    if(nb_zeros == 1 && (on_macro_edge[v] == -1)) {
                        result = false;
                        log_err();
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
                        vertex_on_edge(i, E.first, E.second)
                    ) {
                        result = false;
                        log_err();
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
                        // NOTE: we could have two edges coming from two
                        // different triangles and having an intersection
                        // along an edge with an existing vertex.
                        continue;
                    }
                    if(edge_edge_intersect(
                           E1.first, E1.second, E2.first, E2.second
                    )) {
                        result = false;
                        log_err();
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
                    }
                }
            }
            
            if(!result) {
                mesh_save(
                    debug_constraints,
                    "debug_constraints_" + String::to_string(f1_) +
                    ".geogram"
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


        void sort_vertices_along_edge(
            std::pair<index_t, index_t> E,
            vector<index_t>& V
        ) {
            index_t v_org = E.first;
            std::sort(
                V.begin(), V.end(),
                [&](index_t v1, index_t v2)->bool {
                    // Supposed to be negative if
                    // v1 is between v_org and v2 (what we want)
                    Sign o_o12 = dot2d(v1,v_org,v2);
                    if(o_o12 == ZERO) {
                        std::cerr << "Duplicated vertex in edge"
                                  << std::endl;
                    }
                    return (o_o12 == NEGATIVE);
                }
            );
            // There can be degenerate vertices if there is
            // a triple point in the triangle (quadruple point
            // in 3D) ... shit happens ! [Forrest Gump]
            V.erase(
                std::unique(V.begin(), V.end()), V.end()
            );
            // TODO: check intersections on edge extremities...
        }
        
        vec3 mesh_facet_vertex(index_t f, index_t lv) const {
            index_t v = mesh_.facets.vertex(f,lv);
            return vec3(mesh_.vertices.point_ptr(v));
        }
        
        bool same_point(index_t i, index_t j) const {
            return PCK::same_point(
                vertex_[i].UV_exact,
                vertex_[j].UV_exact
            );
        }

        bool vertex_on_edge(index_t i, index_t j, index_t k) const {
            if(orient2d(i,j,k) != ZERO) {
                return false;
            }
            Sign o = dot2d(i,j,k);
            return (int(o) <= 0);
        }

        bool edge_edge_intersect(
            index_t i, index_t j,
            index_t k, index_t l
        ) const {
            
            Sign o1 = orient2d(i,j,k);
            Sign o2 = orient2d(i,j,l);

            if(o1*o2 > 0) {
                return false;
            }

            // Particular case: 1D
            if(o1 == 0 && o2 == 0) {
                Sign d1 = dot2d(k,i,j);
                Sign d2 = dot2d(l,i,j);
                bool result = (d1 <= 0 || d2 <= 0);
                if(result) {
                    std::cerr << "Isect in co-linear edges" << std::endl;
                }
                return result;
            }
            
            Sign o3 = orient2d(k,l,i);
            Sign o4 = orient2d(k,l,j);
            return (o3*o4 <= 0);
        } 
        

        Sign orient2d(index_t v1,index_t v2,index_t v3) const {
            return PCK::orient_2d(
                vertex_[v1].UV_exact,
                vertex_[v2].UV_exact,
                vertex_[v3].UV_exact
            );
        }

        Sign dot2d(index_t v1, index_t v2, index_t v3) const {
            return PCK::dot_2d(
                vertex_[v1].UV_exact,
                vertex_[v2].UV_exact,
                vertex_[v3].UV_exact
            );
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
                if(V.sym.f2 != index_t(-1)) {
                    isect_facets.push_back(V.sym.f2);
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
        
        void log_err() const {
            std::cerr << "Houston, we got a problem (while remeshing facet "
                      << f1_ << "):" << std::endl;
        }

        vec2 project(index_t i) const {
            return vec2(
                vertex_[i].point_exact[u_].estimate(),
                vertex_[i].point_exact[v_].estimate()
            );
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
        std::map<vec3Q, index_t, vec3QLexicoCompare> g_v_table_;
        std::map<trindex, index_t> t_v_table_;
        bool check_cnstr_;
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
    
    void remesh_intersected_triangles(
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
        Mesh& M, bool test_adjacent_facets, bool order_facets,
        bool FPE, bool check_constraints
    ) {

        // Set symbolic perturbation mode to lexicographic order
        // on point coordinates instead of point indices only,
        // Needed to get compatible triangulations on coplanar faces
        // (example, cubes that touch on a facet).
        PCK::SOSMode SOS_bkp = PCK::get_SOS_mode();
        PCK::set_SOS_mode(PCK::SOS_LEXICO);

        // Exact arithmetics is exact ... until we encounter
        // underflows/overflows (and underflows can happen quite
        // often !!) -> I want to detect them.
        bool FPE_bkp = Process::FPE_enabled();
        Process::enable_FPE(FPE);
        
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

        MeshInTriangle TM(M,check_constraints);
        remesh_intersected_triangles(TM, intersections);
        
        vector<index_t> has_intersections(M.facets.nb(), 0);
        for(const IsectInfo& II: intersections) {
            has_intersections[II.f1] = 1;
            has_intersections[II.f2] = 1;
        }
        
        M.facets.delete_elements(has_intersections);
        M.facets.connect();

        if(!FPE_bkp) {
            Process::enable_FPE(false);
        }
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
        bool order_facets,
        bool FPE,
        bool check_constraints
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

        const double SCALING = double(1ull << 20);
        const double INV_SCALING = 1.0/SCALING;
        
        {
            double* p = mesh_grob()->vertices.point_ptr(0);
            index_t N = mesh_grob()->vertices.nb() *
                        mesh_grob()->vertices.dimension();
            for(index_t i=0; i<N; ++i) {
                p[i] *= SCALING;
            }
        }
        
        intersect_surface_impl(
            *mesh_grob(), test_neighboring_triangles,
            order_facets, FPE, check_constraints
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

        {
            double* p = mesh_grob()->vertices.point_ptr(0);
            index_t N = mesh_grob()->vertices.nb() *
                        mesh_grob()->vertices.dimension();
            for(index_t i=0; i<N; ++i) {
                p[i] *= INV_SCALING;
            }
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
