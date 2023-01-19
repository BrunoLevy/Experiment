
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
#include <geogram/mesh/triangle_intersection.h>
#include <geogram/mesh/mesh_AABB.h>
#include <geogram/mesh/index.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/delaunay/delaunay.h>
#include <geogram/numerics/predicates.h>

#include <geogram/numerics/expansion_nt.h>

namespace {
    using namespace OGF;

    /***********************************************************************/
    
    /**
     * \brief Computes the intersection between a segment and a triangle
     * \pre The intersection exists
     * \param[in) q1 , q2 the two vertices of the segment
     * \param[in] p1 , p2 , p3 the three vertices of the triangle
     * \return the intersection between the segment and the triangle
     */
    vec3 segment_triangle_isect(
        vec3 q1, vec3 q2, vec3 p1, vec3 p2, vec3 p3
    ) {
        vec3 D = q2-q1;
        vec3 E1 = p3-p2;
        vec3 E2 = p1-p3;
        vec3 AO = q1-p1;
        vec3 N  = cross(E1,E2);
        double d =  dot(D,N);
        double t = -dot(AO,N)/d;
        return t*q2 + (1.0-t)*q1;
    }

    /**
     * \brief Compute the intersection between two segments
     * \pre The intersection exists
     * \param[in] P0 , P1 the two extremities of the first segment
     * \param[in] Q0 , Q1 the two extremities of the second segment
     * \return the intersection between the two segments
     */
    vec3 segment_segment_isect(vec3 P0, vec3 P1, vec3 Q0, vec3 Q1) {

        // Taken from https://github.com/davideberly/GeometricTools/blob/
        // master/GTE/Samples/Distance/DistanceSegments3/
        // DistanceSegments3Console.cpp
        
        double sqrDistance,s,t;
        vec3 closest[2];

        geo_argused(sqrDistance);
        geo_argused(s);
        geo_argused(t);        
        
        double const SMALL_NUM = 1e-8;
        vec3   u = P1 - P0;
        vec3   v = Q1 - Q0;
        vec3   w = P0 - Q0;
        double a = dot(u, u);         // always >= 0
        double b = dot(u, v);
        double c = dot(v, v);         // always >= 0
        double d = dot(u, w);
        double e = dot(v, w);
        double D = a*c - b*b;        // always >= 0
        double sc, sN, sD = D;       // sc = sN / sD, default sD = D >= 0
        double tc, tN, tD = D;       // tc = tN / tD, default tD = D >= 0
        
        // compute the line parameters of the two closest points
        if (D < SMALL_NUM) { // the lines are almost parallel
            sN = 0.0;        // force using point P0 on segment S1
            sD = 1.0;        // to prevent possible division by 0.0 later
            tN = e;
            tD = c;
        } else {             // get the closest points on the infinite lines
            sN = (b*e - c*d);
            tN = (a*e - b*d);
            if (sN < 0.0) {        // sc < 0 => the s=0 edge is visible
                sN = 0.0;
                tN = e;
                tD = c;
            } else if (sN > sD) {  // sc > 1  => the s=1 edge is visible
                sN = sD;
                tN = e + b;
                tD = c;
            }
        }
        
        if (tN < 0.0) {            // tc < 0 => the t=0 edge is visible
            tN = 0.0;
            // recompute sc for this edge
            if (-d < 0.0) {
                sN = 0.0;
            } else if (-d > a) {
                sN = sD;
            } else {
                sN = -d;
                sD = a;
            }
        }
        else if (tN > tD) {      // tc > 1  => the t=1 edge is visible
            tN = tD;
            // recompute sc for this edge
            if ((-d + b) < 0.0) {
                sN = 0;
            } else if ((-d + b) > a) {
                sN = sD;
            } else {
                sN = (-d + b);
                sD = a;
            }
        }
        // finally do the division to get sc and tc
        sc = (std::fabs(sN) < SMALL_NUM ? 0.0 : sN / sD);
        tc = (std::fabs(tN) < SMALL_NUM ? 0.0 : tN / tD);
        
        // get the difference of the two closest points
        s = sc;
        t = tc;
        closest[0] = (1.0 - sc) * P0 + sc * P1;
        closest[1] = (1.0 - tc) * Q0 + tc * Q1;
        vec3 diff = closest[0] - closest[1];
        sqrDistance = dot(diff, diff);
        return closest[1];
    }    

    void print(std::ostream& out, const GEO::expansion_nt& x) {
        out << "expansion_nt(estimate="
            << x.estimate();
        out << ", components=[";
        for(GEO::index_t i=0; i<x.length(); ++i) {
            out << x.component(i);
            if(i != x.length()-1) {
                out << " ";
            }
        }
        out << "]";
        out << ")";
    }

    vec3 segment_segment_isect_2D(
        vec3 p1, vec3 p2, vec3 q1, vec3 q2, coord_index_t ax
    ) {
        coord_index_t u = coord_index_t((ax + 1)%3);
        coord_index_t v = coord_index_t((ax + 2)%3);

        // [ a b ] [ l1 ]   [ e ]
        // [ c d ] [ l2 ] = [ f ]
        
        double a = p2[u]-p1[u];
        double b = q1[u]-q2[u];
        double c = p2[v]-p1[v];
        double d = q1[v]-q2[v];

        double e = q1[u]-p1[u];
        double f = q1[v]-p1[v];

        // [ a b ] [ s ]   [ e ]
        // [ c d ] [ t ] = [ f ]

        // [ s ]               [ d -b] [ e ]
        // [ t ] = 1/(ad-bc) * [-c  a] [ f ]

        double det = a*d-b*c;
        
        double s = ( d*e-b*f)/det;
        double t = (-c*e+a*f)/det;

        return s*p2+(1.0-s)*p1;
    }
    
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
                init_sym(     M, f1, f2, R1, R2);
                init_geometry(M);
                mesh_vertex_index = NO_INDEX;
            }


            Vertex(Mesh& M, index_t v) {
                set_indices(v);
                init_geometry(M);
                mesh_vertex_index = v;
                has_tt_sym_ = false;
                f1_ = index_t(-1);
                f2_ = index_t(-1);
            }

            
            Vertex() {
                set_indices(NO_INDEX, NO_INDEX, NO_INDEX, NO_INDEX);
                mesh_vertex_index = NO_INDEX;
                has_tt_sym_ = false;                
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

            void debug_show() const {
                std::cerr << "SYM:" << sym;
                if(has_tt_sym_) {
                    std::cerr
                        << "   TTSYM: "
                        << f1_ << " " << f2_ << " " << std::make_pair(R1_,R2_);
                }
                std::cerr << std::endl;
            }
            
        protected:
            void init_sym(
                Mesh& M,
                index_t f1, index_t f2,
                TriangleRegion R1, TriangleRegion R2
            ) {
                has_tt_sym_ = true;
                f1_ = f1;
                f2_ = f2;
                R1_ = R1;
                R2_ = R2;
                
                index_t v[6] = {
                    M.facets.vertex(f1,0),
                    M.facets.vertex(f1,1),
                    M.facets.vertex(f1,2),
                    M.facets.vertex(f2,0),
                    M.facets.vertex(f2,1),
                    M.facets.vertex(f2,2)
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


            void init_geometry(const Mesh& M) {
                if(is_existing_vertex()) {
                    point = vec3(M.vertices.point_ptr(sym.indices[0]));
                    return;
                }

                if(is_edge_edge_isect()) {

                    // We distinguish 3D edge /\ edge isect (then we use
                    // edge /\ facet isect) and 2D edge /\ edge isect (
                    // then develop special code). 
                    
                    geo_debug_assert(has_tt_sym_);
                    
                    // The three vertices of f1, with the intersected
                    // edges as (p1,p2)
                    vec3 p1, p2, p3;
                    
                    // The intersected edge of f2
                    vec3 q1, q2;

                    {
                        TriangleRegion lvr1, lvr2;
                        index_t lv1, lv2, lv3;
                        // The two vertices of R1 first
                        GEO::get_edge_vertices(R1_,lvr1, lvr2);
                        lv1 = index_t(lvr1);
                        lv2 = index_t(lvr2);
                        geo_debug_assert(lv1 < 3);
                        geo_debug_assert(lv2 < 3);                        
                        // The other vertex of f1
                        lv3 = index_t(R1_ - 6);
                        
                        index_t v1 = M.facets.vertex(f1_, lv1);
                        index_t v2 = M.facets.vertex(f1_, lv2);
                        index_t v3 = M.facets.vertex(f1_, lv3);

                        p1 = vec3(M.vertices.point_ptr(v1));
                        p2 = vec3(M.vertices.point_ptr(v2));
                        p3 = vec3(M.vertices.point_ptr(v3));
                    }

                    {
                        TriangleRegion lvr1, lvr2;
                        GEO::get_edge_vertices(R2_, lvr1, lvr2);
                        index_t lv1 = index_t(lvr1);
                        index_t lv2 = index_t(lvr2);
                        geo_debug_assert(lv1 >= 3 && lv1 < 6);
                        geo_debug_assert(lv2 >= 3 && lv2 < 6);
                        lv1 -= 3;
                        lv2 -= 3;
                        index_t v1 = M.facets.vertex(f2_,lv1);
                        index_t v2 = M.facets.vertex(f2_,lv2);
                        q1 = vec3(M.vertices.point_ptr(v1));
                        q2 = vec3(M.vertices.point_ptr(v2));
                    }

                    if(
                        PCK::orient_3d(p1,p2,p3,q1) == ZERO &&
                        PCK::orient_3d(p1,p2,p3,q2) == ZERO
                    ) {
                        // point = segment_segment_isect(q1,q2,p1,p2);
                        point = segment_segment_isect_2D(
                            q1,q2,p1,p2,
                            ::GEO::Geom::triangle_normal_axis(p1,p2,p3)
                        );
                    } else {
                        // We are in 3D, we can use segment /\ triangle
                        // intersection (even if we know that the intersection
                        // will be on on (p1,p2))
                        point = segment_triangle_isect(q1,q2,p1,p2,p3);
                    }
                    return;
                }
                
                geo_debug_assert(is_edge_triangle_isect());
                vec3 q1(M.vertices.point_ptr(sym.indices[0]));
                vec3 q2(M.vertices.point_ptr(sym.indices[1]));
                index_t v1 = M.facets.vertex(sym.indices[2],0);
                index_t v2 = M.facets.vertex(sym.indices[2],1);
                index_t v3 = M.facets.vertex(sym.indices[2],2);
                vec3 p1(M.vertices.point_ptr(v1));
                vec3 p2(M.vertices.point_ptr(v2));
                vec3 p3(M.vertices.point_ptr(v3));
                point = segment_triangle_isect(q1,q2,p1,p2,p3);
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
            vec3 point;
            quadindex sym;
            index_t mesh_vertex_index;

            // Initial combinatorial information,
            // Needed for:
            //  - edge-edge intersection, to distinguish 2D case
            //  - debugging
            bool has_tt_sym_;
            index_t f1_,f2_;
            TriangleRegion R1_,R2_;
            
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
                add_vertex(Vertex(mesh_, v[lv]));
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

            // DEBUG: show list of intersected facets
            if(false) {
                std::vector<index_t> isect_facets;
                for(const Vertex& V: vertex_) {
                    if(V.has_tt_sym_) {
                        isect_facets.push_back(V.f2_);
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
            
            bool ok = true;

            if(f1_orient_ == ZERO) {
                std::cerr << ">>>>> Facet " << f1_ << " has ZERO orient2d"
                          << std::endl;
            }
            
            // Sort vertices along triangle's edges
            {
//                std::cerr << "Remesh " << f1_ << std::endl;
                for(index_t e=0; e<3; ++e) {

/*                    
                    std::cerr << "Macro edge " << e << " has "
                              << vertices_in_E_[e].size()
                              << " additional vert(ex)(ices)"
                              << std::endl;

                    for(index_t v: vertices_in_E_[e]) {
                        vertex_[v].debug_show(); 
                    }
*/
                    index_t v_org = (e+1)%3; // e's origin
                    index_t v_opp = e;       // f1's vertex opposite to e
                    std::sort(
                        vertices_in_E_[e].begin(),
                        vertices_in_E_[e].end(),
                        [&](index_t v1, index_t v2)->bool{

                            Sign o12 = orient2d(v1,v2,v_opp);

                            Sign o_o12 = PCK::dot_3d(
                                vertex_[v1].point.data(),
                                vertex_[v_org].point.data(),
                                vertex_[v2].point.data()
                            );
                            
                            if(o12 == ZERO || o_o12 == ZERO) {
                                std::cerr << "Houston, we got a problem !"
                                          << std::endl;
                                std::cerr << "While remeshing "
                                          << f1_
                                          << std::endl;
                                std::cerr << "Macro edge: " << e
                                          << std::endl;
                                std::cerr << "(found two coincident vertices)"
                                          << std::endl;
                                vertex_[v1].debug_show();
                                vertex_[v2].debug_show();
                                ok = false;
                            }

                            bool result_1 = (o12 == f1_orient_);
                            bool result_2 = (o_o12 == NEGATIVE);

                            if(result_1 != result_2) {
                                std::cerr << "Houston, we got a problem !"
                                          << std::endl;
                                std::cerr << "While remeshing "
                                          << f1_
                                          << std::endl;
                                std::cerr << "Macro edge: " << e
                                          << std::endl;
                                std::cerr << "Different results in two tests: "
                                          << std::endl;
                                std::cerr << " 2d: " << o12
                                          << " ==? "  << f1_orient_
                                          << std::endl;
                                std::cerr << " 3d: " << o_o12
                                          << " negative?"
                                          << std::endl;
                            }
                            
                            return result_1;
                        }
                    );
                }
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

            // if(ok) {
                commit();
            //}
            clear();
        }

    protected:

        void commit() {
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

                mesh_save(
                    constraints,
                    "constraints_" + String::to_string(f1_) + ".geogram"
                );
                
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

            if(false) {
                // Create edges in Mesh (for visualization purposes)
                for(std::pair<index_t, index_t>& E: edges_) {
                    mesh_.edges.create_edge(
                        vertex_[E.first].mesh_vertex_index,
                        vertex_[E.second].mesh_vertex_index
                    );
                }
            }
        }

        Sign orient2d(index_t v1, index_t v2, index_t v3) {
            return PCK::orient_2d(project(v1), project(v2), project(v3));
        }
        
        vec2 project(index_t v) {
            const vec3& p = vertex_[v].point;
            return vec2(p[u_], p[v_]);
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
    
    void intersect_surface_impl(Mesh& M, bool test_adjacent_facets) {

        // Set symbolic perturbation mode to lexicographic order
        // on point coordinates instead of point indices only,
        // Needed to get compatible triangulations on coplanar faces
        // (example, cubes that touch on a facet).
        PCK::SOSMode SOS_bkp = PCK::get_SOS_mode();
        PCK::set_SOS_mode(PCK::SOS_LEXICO);
        
        vector<IsectInfo> intersections;
        MeshFacetsAABB AABB(M);
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
        bool test_neighboring_triangles 
    ) {
        intersect_surface_impl(*mesh_grob(), test_neighboring_triangles);
        show_mesh();
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

    static void show_o3(vec3 p[], index_t i, index_t j, index_t k, index_t l) {
        Logger::out("Expe") <<
            "o3(" << i << "," << j << "," << k << "," << l << ")="
                  << PCK::orient_3d(p[i], p[j], p[k], p[l])
                  << std::endl;
    }

    static void show_o3(
        vec3 p[], index_t i, index_t j, index_t k, index_t l, index_t m
    ) {
        Sign o1 = PCK::orient_3d(p[i], p[j], p[l], p[m]);
        Sign o2 = PCK::orient_3d(p[i], p[j], p[m], p[k]);
        Sign o3 = PCK::orient_3d(p[i], p[j], p[k], p[l]);
        Logger::out("Expe") <<
            "(" << i << "," << j <<") / (" << k << "," << l << "," << m << ")"
                << "   "
                << o1 << "," << o2 << "," << o3
                << std::endl;
        show_o3(p, i, k, l, m);
        show_o3(p, j, k, l, m);
    }

    
    void MeshGrobExperimentCommands::debug_ze_case() {
        if(
            mesh_grob()->vertices.nb() != 11 ||
            mesh_grob()->facets.nb()   != 8
        ) {
            Logger::err("Expe") << "This is not ze case"
                                << std::endl;
            return;
        }

        vec3 P[11];
        for(index_t i=0; i<11; ++i) {
            P[i] = vec3(mesh_grob()->vertices.point_ptr(i));
        }

        show_o3(P,0,   8,9,10);
        show_o3(P,0,1, 8,9,10);
        show_o3(P,0,2, 8,9,10);
        show_o3(P,0,3, 8,9,10);
        show_o3(P,0,4, 8,9,10);
        show_o3(P,0,5, 8,9,10);
        show_o3(P,0,6, 8,9,10);

        Logger::out("Expe") << std::endl;
        
        show_o3(P,0,   7,8,10);
        show_o3(P,0,1, 7,8,10);
        show_o3(P,0,2, 7,8,10);
        show_o3(P,0,3, 7,8,10);
        show_o3(P,0,4, 7,8,10);
        show_o3(P,0,5, 7,8,10);
        show_o3(P,0,6, 7,8,10);                
        
    }
    
}
