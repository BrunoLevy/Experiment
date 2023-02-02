
/*
 *
 */

// BUG: PR10: generates several triangles several times (?) 
// BUG: PR3:  generates several triangles several times (?)
//            / misses isects ?
//
// TODO: can we treat the three edges of the macro-triangle and the
// other edges the same way ? (that is, starting with the three edges,
// and adding the intersections only when detecting them). It would be
// a bit stupid, because we know already that they are in the three
// edges...
// TODO: each time a point is added, debug-check whether it is in
// one of the three edges of the triangle
//
// TODO: figure out what happens for ~/three_cubes.obj
//
// TODO: sort_uniq all edges to detect duplicated edges


#include <OGF/Experiment/algo/mesh_surface_intersection.h>
#include <OGF/Experiment/algo/exact_geometry.h>
#include <geogram/mesh/triangle_intersection.h>
#include <geogram/mesh/mesh_AABB.h>
#include <geogram/mesh/mesh_repair.h>
#include <geogram/mesh/mesh_fill_holes.h>
#include <geogram/mesh/index.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/delaunay/delaunay.h>
#include <geogram/numerics/predicates.h>
#include <geogram/numerics/expansion_nt.h>

#include <sstream>

namespace {
    using namespace GEO;

    /***********************************************************************/


    /**
     * \brief Meshes a single triangle with the constraints that come from
     *  the intersections with the other triangles.
     */
    class MeshInTriangle {
    public:

        /***************************************************************/
        
        class Edge {
        public:
            Edge(
                index_t v1_in = index_t(-1),
                index_t v2_in = index_t(-1),
                index_t f2    = index_t(-1),
                TriangleRegion R2 = T2_RGN_T
            ) :
                v1(v1_in),
                v2(v2_in)
            {
                sym.f2 = f2;
                sym.R2 = R2;
            }
            index_t v1;
            index_t v2;
            struct {
                index_t        f2;
                TriangleRegion R2;
            } sym;
            vector<index_t> vertices_in_edge;
        };

        /***************************************************************/
        
        class Vertex {
        public:
            enum { NO_INDEX = index_t(-1) };
            enum Type {
                UNINITIALIZED, MESH_VERTEX, PRIMARY_ISECT, SECONDARY_ISECT
            };

            /**
             * \brief Constructor for macro-triangle vertices.
             */
            Vertex(MeshInTriangle* M, index_t f, index_t lv) {
                geo_assert(f == M->f1_);
                type = MESH_VERTEX;
                mesh_in_triangle = M;
                init_sym(f, NO_INDEX, TriangleRegion(lv), T2_RGN_T);
                init_geometry();
            }

            /**
             * \brief Constructor for intersections with other facets.
             */
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

            /**
             * \brief Constructor for intersections between constraints.
             */
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

            Mesh& mesh() const {
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
            
        protected:

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
                        sym.f1, index_t(sym.R1)-index_t(T1_RGN_P0)
                    );
                }
                if(region_dim(sym.R2) == 0) {
                    mesh_vertex_index = mesh().facets.vertex(
                        sym.f2, index_t(sym.R2)-index_t(T2_RGN_P0)
                    );
                }
            }

            void init_geometry() {

                static vec2 uv[3] = {
                    vec2(0.0, 0.0),
                    vec2(1.0, 0.0),
                    vec2(0.0, 1.0)
                };

                mesh_vertex_index = NO_INDEX;
                
                // Case 1: f1 vertex
                if(region_dim(sym.R1) == 0) {
                    index_t lv = index_t(sym.R1);
                    geo_assert(lv < 3);
                    mesh_vertex_index = mesh().facets.vertex(sym.f1,lv);
                    point_exact = mesh_vertex<vec3Q>(mesh_vertex_index);
                    if(barycentric()) {
                        UV_exact.x = rational_nt(uv[lv].x);
                        UV_exact.y = rational_nt(uv[lv].y);
                    } else {
                        UV_exact = mesh_vertex_project<vec2Q>(
                            mesh_vertex_index
                        );
                    }
                    return;
                }

                geo_assert(sym.f1 != NO_INDEX && sym.f2 != NO_INDEX);

                // Case 2: f2 vertex
                if(region_dim(sym.R2) == 0) {
                    index_t lv = index_t(sym.R2)-3;
                    geo_assert(lv < 3);
                    mesh_vertex_index = mesh().facets.vertex(sym.f2, lv);
                    point_exact = mesh_vertex<vec3Q>(mesh_vertex_index);
                    if(barycentric()) {
                        vec2E q  = mesh_facet_vertex_project<vec2E>(
                            sym.f2, lv
                        );
                        vec2E p1 = mesh_facet_vertex_project<vec2E>(
                            sym.f1, 0
                        );
                        vec2E p2 = mesh_facet_vertex_project<vec2E>(
                            sym.f1, 1
                        );
                        vec2E p3 = mesh_facet_vertex_project<vec2E>(
                            sym.f1, 2
                        );
                        expansion_nt D = det(p2-p1,p3-p1);
                        UV_exact = vec2Q(
                            rational_nt(det(q-p1,p3-p1),D),
                            rational_nt(det(p2-p1,q-p1),D)
                        );
                    } else {
                        UV_exact = mesh_vertex_project<vec2Q>(
                            mesh_vertex_index
                        );
                    }
                    return;
                }

                // case 3: f1 /\ f2 edge or f1 edge /\ f2 edge in 3D
                if(
                    (region_dim(sym.R1) == 2 || region_dim(sym.R1) == 1) &&
                     region_dim(sym.R2) == 1
                ) {
                    
                    // Compute u,v,t using Moller & Trumbore's algorithm
                    // see: https://stackoverflow.com/questions/42740765/
                    //  intersection-between-line-and-triangle-in-3d

                    vec3 p1 = mesh_facet_vertex(sym.f1, 0);
                    vec3 p2 = mesh_facet_vertex(sym.f1, 1);
                    vec3 p3 = mesh_facet_vertex(sym.f1, 2);
                    
                    index_t e = index_t(sym.R2)-index_t(T2_RGN_E0);
                    geo_debug_assert(e<3);
                    vec3 q1 = mesh_facet_vertex(sym.f2, (e+1)%3);
                    vec3 q2 = mesh_facet_vertex(sym.f2, (e+2)%3);

                    bool seg_seg_two_D = (
                        region_dim(sym.R1) == 1 &&
                        PCK::orient_3d(p1,p2,p3,q1) == ZERO &&
                        PCK::orient_3d(p1,p2,p3,q2) == ZERO) ;

                    if(!seg_seg_two_D) {
                        vec3E D   = make_vec3<vec3E>(q1,q2);
                        vec3E E1  = make_vec3<vec3E>(p1,p2);
                        vec3E E2  = make_vec3<vec3E>(p1,p3);
                        vec3E AO  = make_vec3<vec3E>(p1,q1);
                        vec3E N   = cross(E1,E2);
                    
                        expansion_nt d = -dot(D,N);
                        geo_assert(d.sign() != ZERO);
                        rational_nt t(dot(AO,N),d);
                        point_exact = mix(t,q1,q2);

                        if(barycentric()) {
                            vec3E DAO = cross(AO,D);
                            UV_exact.x = rational_nt( dot(E2,DAO),d);
                            UV_exact.y = rational_nt(-dot(E1,DAO),d);
                        } else {
                            UV_exact.x = point_exact[u_coord()];
                            UV_exact.y = point_exact[v_coord()];
                        }
                        return;
                    }
                }

                // case 4: f1 edge /\ f2
                if(region_dim(sym.R1) == 1 && region_dim(sym.R2) == 2) {
                    index_t e = index_t(sym.R1)-index_t(T1_RGN_E0);
                    geo_debug_assert(e<3);
                    vec3 p1 = mesh_facet_vertex(sym.f1, (e+1)%3);
                    vec3 p2 = mesh_facet_vertex(sym.f1, (e+2)%3);

                    vec3 q1 = mesh_facet_vertex(sym.f2,0);
                    vec3 q2 = mesh_facet_vertex(sym.f2,1);
                    vec3 q3 = mesh_facet_vertex(sym.f2,2);
                    
                    vec3E E1  = make_vec3<vec3E>(q1,q2);
                    vec3E E2  = make_vec3<vec3E>(q1,q3);
                    vec3E N   = cross(E1,E2);

                    vec3E D   = make_vec3<vec3E>(p1,p2);
                    vec3E AO  = make_vec3<vec3E>(q1,p1);

                    expansion_nt d = -dot(D,N);
                    rational_nt t(dot(AO,N),d);

                    point_exact = mix(t, p1, p2);
                    if(barycentric()) {
                        UV_exact    = mix(t, uv[(e+1)%3], uv[(e+2)%3]);
                    } else {
                        UV_exact.x = point_exact[u_coord()];
                        UV_exact.y = point_exact[v_coord()];
                    }
                    return;
                }

                // case 5: f1 edge /\ f2 edge in 2D
                if(region_dim(sym.R1) == 1 && region_dim(sym.R2) == 1) {

                    index_t e1 = index_t(sym.R1) - index_t(T1_RGN_E0);
                    geo_debug_assert(e1 < 3);
                    index_t e2 = index_t(sym.R2) - index_t(T2_RGN_E0);
                    geo_debug_assert(e2 < 3);

                    vec2 p1 = mesh_facet_vertex_project<vec2>(
                        sym.f1, (e1+1)%3
                    );
                    vec2 p2 = mesh_facet_vertex_project<vec2>(
                        sym.f1, (e1+2)%3
                    );

                    vec2 q1 = mesh_facet_vertex_project<vec2>(
                        sym.f2, (e2+1)%3
                    );
                    vec2 q2 = mesh_facet_vertex_project<vec2>(
                        sym.f2, (e2+2)%3
                    );

                    vec2E D1 = make_vec2<vec2E>(p1,p2);
                    vec2E D2 = make_vec2<vec2E>(q1,q2);

                    expansion_nt d = det(D1,D2);
                    geo_assert(d.sign() != ZERO);
                    
                    vec2E AO = make_vec2<vec2E>(p1,q1);
                    vec3 P1 = mesh_facet_vertex(sym.f1, (e1+1)%3);
                    vec3 P2 = mesh_facet_vertex(sym.f1, (e1+2)%3);
                    rational_nt t(det(AO,D2),d);
                    point_exact = mix(t,P1,P2);
                    if(barycentric()) {
                        UV_exact = mix(t, uv[(e1+1)%3], uv[(e1+2)%3]);
                    } else {
                        UV_exact.x = point_exact[u_coord()];
                        UV_exact.y = point_exact[v_coord()];
                    }
                    return;
                }

                // Normally we enumerated all possible cases
                geo_assert_not_reached;
            }

            bool barycentric() const {
                return mesh_in_triangle->barycentric_;
            }

            coord_index_t u_coord() const {
                return mesh_in_triangle->u_;
            }

            coord_index_t v_coord() const {
                return mesh_in_triangle->v_;                
            }
            
            template <class VEC = vec3> VEC
            mesh_vertex(index_t v) const {
                typedef typename VEC::value_type value_type;
                const double* p = mesh().vertices.point_ptr(v);
                return VEC(
                    value_type(p[0]),
                    value_type(p[1]),
                    value_type(p[2])                    
                );
            }
            
            template <class VEC = vec3> VEC
            mesh_facet_vertex(index_t f, index_t lv) const {
                index_t v = mesh().facets.vertex(f,lv);
                return mesh_vertex<VEC>(v);
            }

            template <class VEC = vec2> VEC
            mesh_vertex_project(index_t v) const {
                typedef typename VEC::value_type value_type;
                const double* p = mesh().vertices.point_ptr(v);
                return VEC(
                    value_type(p[u_coord()]),
                    value_type(p[v_coord()])
                );
            }
            
            template <class VEC = vec2> VEC
            mesh_facet_vertex_project(index_t f, index_t lv) const {
                index_t v = mesh().facets.vertex(f,lv);
                return mesh_vertex_project<VEC>(v);
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
        };

        /***************************************************************/
        
        MeshInTriangle(Mesh& M) :
            mesh_(M),
            f1_(index_t(-1)),
            check_cnstr_(true),
            barycentric_(true)
        {
        }

        void set_check_constraints(bool x) {
            check_cnstr_ = x;
        }

        void set_barycentric(bool x) {
            barycentric_ = x;
        }

        Mesh& mesh() {
            return mesh_;
        }
        
        void begin_facet(index_t f) {
            f1_ = f;
            latest_f2_ = index_t(-1);
            latest_f2_count_ = 0;
            
            vec3 p1 = mesh_facet_vertex(f,0);
            vec3 p2 = mesh_facet_vertex(f,1);
            vec3 p3 = mesh_facet_vertex(f,2);
            
            f1_normal_axis_ = ::GEO::Geom::triangle_normal_axis(
                p1,p2,p3
            );
            
            u_ = coord_index_t((f1_normal_axis_ + 1) % 3);
            v_ = coord_index_t((f1_normal_axis_ + 2) % 3);
            
            for(index_t lv=0; lv<3; ++lv) {
                add_vertex(Vertex(this, f, lv));
            }

            edges_.push_back(Edge(1,2));
            edges_.push_back(Edge(2,0));
            edges_.push_back(Edge(0,1));
            
            if(barycentric_) {
                f1_orient_ = PCK::orient_2d(
                    vec2(0.0, 0.0),
                    vec2(1.0, 0.0),
                    vec2(0.0, 1.0)
                );
            } else {
                f1_orient_ = PCK::orient_2d(
                    vec2(p1[u_],p1[v_]),
                    vec2(p2[u_],p2[v_]),
                    vec2(p3[u_],p3[v_])                    
                );
            }
            // Sanity check: the facet does not have zero area
            geo_assert(f1_orient_ != ZERO);
            has_planar_isect_ = false;
        }
        
        index_t add_vertex(
            index_t f2,
            TriangleRegion R1, TriangleRegion R2
        ) {
            geo_debug_assert(f1_ != index_t(-1));

            // If the same f2 comes more than twice, then
            // we got a planar facet /\ facet intersection
            // (and it is good to know it, see get_constraints())
            if(f2 != index_t(-1) && f2 == latest_f2_) {
                ++latest_f2_count_;
                if(latest_f2_count_ > 2) {
                    if(!has_planar_isect_) {
                        std::cerr << f1_
                                  << ":switching to planar isect mode"
                                  << std::endl;
                    }
                    has_planar_isect_ = true;
                }
            } else {
                latest_f2_ = f2;
                latest_f2_count_ = 0;
            }
            
            Vertex V(this, f1_, f2, R1, R2);
            index_t sz = vertex_.size();
            index_t v = add_vertex(V);
            // If it is a new vertex, and if it is on
            // an edge, add it to the list of edge's
            // vertices
            if(vertex_.size() > sz) {
                if(region_dim(R1) == 1) {
                    index_t e = index_t(R1) - index_t(T1_RGN_E0);
                    edges_[e].vertices_in_edge.push_back(v);
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

            TriangleRegion R = T2_RGN_T;
            
            if(region_dim(AR2) == 1 && region_dim(BR2) == 0) {
                index_t e = index_t(AR2) - index_t(T2_RGN_E0);
                index_t v = index_t(BR2) - index_t(T2_RGN_P0);
                if(e != v) {
                    R = AR2;
                }
            } else if(region_dim(BR2) == 1 && region_dim(AR2) == 0) {
                index_t e = index_t(BR2) - index_t(T2_RGN_E0);
                index_t v = index_t(AR2) - index_t(T2_RGN_P0);
                if(e != v) {
                    R = BR2;
                }
            }
            
            if(region_dim(R) == 1) {
                edges_.push_back(Edge(v1,v2,f2,R));
            } else {
                edges_.push_back(Edge(v1,v2,f2));                
            }
        }

        void end_facet() {
            compute_constraints_intersections();
            commit();
            clear();
        }

    protected:

        index_t add_vertex(const Vertex& V) {
            // Test whether there is already a vertex at the same position
            for(index_t i=0; i<vertex_.size(); ++i) {
                if(PCK::same_point(vertex_[i].UV_exact,V.UV_exact)) {
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

            // std::cerr << " Nb edges=" << edges_.size() << std::endl;
            
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
            //    Triangle created steiner points.
            // It is not normal because we have generated all
            // intersections before. 
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
		   vec2 p;

                   // If barycentric coodinates are used and the facet
                   // has a planar intersection, then we need to make sure
                   // that the facet will be meshed the same way in other
                   // places (that is, in the other facet with which it has
                   // a planar intersection), thus we compute the (projected)
                   // coordinates in the triangle (and pass them to Triangle).
		   if(barycentric_ && has_planar_isect_) {
                       const vec2Q& uv = vertex_[v].UV_exact;
                       vec2 p1 = mesh_facet_vertex_project(f1_,0);
                       vec2 p2 = mesh_facet_vertex_project(f1_,1);
                       vec2 p3 = mesh_facet_vertex_project(f1_,2);
                       vec2Q p_exact = u_P1P2_plus_v_P1P3(uv.x, uv.y, p1, p2, p3);
                       p.x = p_exact.x.estimate();
                       p.y = p_exact.y.estimate();
                   } else {
		       p = vec2(
                           vertex_[v].UV_exact.x.estimate(),
                           vertex_[v].UV_exact.y.estimate()
                       );
		    }
                    M.vertices.create_vertex(p.data());
                }
            }
            if(with_edges && M.edges.nb() == 0) {
                index_t i=0;
                for(const Edge& E: edges_) {
                    M.edges.create_edge(E.v1, E.v2); 
                    ++i;
                }
            }
        }

        void compute_constraints_intersections() {

            // Step 1: Compute all intersections
            // (starting at 3, ignoring the 3 macro-edges)
            for(index_t e1=3; e1<edges_.size(); ++e1) {
                index_t v1 = edges_[e1].v1;
                index_t v2 = edges_[e1].v2;
                for(index_t e2=e1+1; e2<edges_.size(); ++e2) {
                    index_t w1 = edges_[e2].v1;
                    index_t w2 = edges_[e2].v2;
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
                        index_t f2 = edges_[e1].sym.f2; 
                        index_t f3 = edges_[e2].sym.f2; 

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
                        vec2Q I_uv;
                        
                        if(
                            get_three_planes_intersection(
                                I,
                                P[0], P[1], P[2],
                                P[3], P[4], P[5],
                                P[6], P[7], P[8]
                            )
                        ) {
                            // Barycentric mode:
                            // compute (u,v) from 3D geometry (direct from
                            // data, reduces degree of used expansions)
                            if(barycentric_) {
                                
                                vec3E E1 = make_vec3<vec3E>(P[0],P[1]); 
                                vec3E E2 = make_vec3<vec3E>(P[0],P[2]); 

                                vec3E N1 = triangle_normal<vec3E>(P[3],P[4],P[5]);
                                vec3E N2 = triangle_normal<vec3E>(P[6],P[7],P[8]);
                                
                                vec2E C1(dot(E1,N1),dot(E1,N2));
                                vec2E C2(dot(E2,N1),dot(E2,N2));
                            
                                vec2E B(
                                    dot(make_vec3<vec3E>(P[0],P[3]),N1),
                                    dot(make_vec3<vec3E>(P[0],P[6]),N2)
                                );

                                expansion_nt d = det(C1,C2);
                                geo_assert(d.sign() != ZERO);

                                I_uv.x = rational_nt(det(B,C2),d);
                                I_uv.y = rational_nt(det(C1,B),d);
                            } else {
                                I_uv.x = I[u_];
                                I_uv.y = I[v_];
                            }
                            
                        } else {
                            std::cerr << "2D constraints intersection" << std::endl;

                            std::cerr << edges_[e1].sym.f2 << "." << region_to_string(edges_[e1].sym.R2)
                                      << " /\\ "
                                      << edges_[e2].sym.f2 << "." << region_to_string(edges_[e2].sym.R2)
                                      << std::endl;
                            
                            geo_assert_not_reached;
                        }
                        
                        index_t x = add_vertex(Vertex(this,I,I_uv));
                        //vertex_[x].facets.push_back(f2);
                        //vertex_[x].facets.push_back(f3);

                        edges_[e1].vertices_in_edge.push_back(x);
                        edges_[e2].vertices_in_edge.push_back(x);
                    }
                }
            }

            // Step 2: sort intersections in edges
            for(Edge& E: edges_) {
                vector<index_t>& V = E.vertices_in_edge;
                index_t v_org = E.v1;
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
            }

            // Step 3: "remesh" the edges that have intersections
            index_t ne = edges_.size();
            for(index_t e=0; e<ne; ++e) {
                index_t prev_v = edges_[e].v1;
                index_t last_v = edges_[e].v2;
                edges_[e].vertices_in_edge.push_back(last_v);
                for(index_t i=0; i<edges_[e].vertices_in_edge.size(); ++i) {
                    index_t cur_v = edges_[e].vertices_in_edge[i];
                    if(i == 0) {
                        edges_[e].v2 = cur_v;
                    } else {
                        edges_.push_back(Edge(prev_v, cur_v));
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
                    for(index_t v: edges_[e].vertices_in_edge) {
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
            if(false) // TO BE ADAPTED
            for(index_t i=0; i<vertex_.size(); ++i) {
                for(const Edge& E: edges_) {
                    if(
                        i != E.v1 && i != E.v2 && 
                        vertex_on_edge(i, E.v1, E.v2)
                    ) {
                        result = false;
                        log_err();
                        std::cerr << "Point on edge !" << std::endl;
                        std::cerr << "    Pt:"<< vertex_[i].to_string()
                                  << std::endl;
                        std::cerr << "   Edge:"<< vertex_[E.v1].to_string()
                                  << " --- "
                                  << vertex_[E.v2].to_string() << std::endl;
                        get_constraints(debug_constraints, false);
                        debug_vertex_show[i] = true;
                        debug_vertex_show[E.v1] = true;
                        debug_vertex_show[E.v2] = true;
                        debug_constraints.edges.create_edge(E.v1, E.v2);
                    }
                }
            }
            
            // Test constrained edges intersetions
            for(index_t e1=0; e1<edges_.size(); ++e1) {
                const Edge& E1 = edges_[e1];                
                for(index_t e2=e1+1; e2 < edges_.size(); ++e2) {
                    const Edge& E2 = edges_[e2];
                    if(
                        E1.v1 == E2.v1 ||
                        E1.v1 == E2.v2 ||
                        E1.v2 == E2.v1 ||
                        E1.v2 == E2.v2
                    ) {
                        // NOTE: we could have two edges coming from two
                        // different triangles and having an intersection
                        // along an edge with an existing vertex.
                        continue;
                    }
                    if(edge_edge_intersect(E1.v1, E1.v2, E2.v1, E2.v2)) {
                        result = false;
                        log_err();
                        std::cerr << "Intersecting edges !" << std::endl;
                        std::cerr << "   Edge1: " 
                                  << vertex_[E1.v1].to_string()
                                  << " --- "
                                  << vertex_[E1.v2].to_string()
                                  << std::endl;
                        std::cerr << "   Edge2: " 
                                  << vertex_[E2.v1].to_string()
                                  << " --- "
                                  << vertex_[E2.v2].to_string()
                                  << std::endl;
                        get_constraints(debug_constraints, false);
                        debug_constraints.edges.create_edge(
                            E1.v1, E1.v2
                        );
                        debug_constraints.edges.create_edge(
                            E2.v1, E2.v2
                        );
                        debug_vertex_show[E1.v1] = true;
                        debug_vertex_show[E1.v2] = true;
                        debug_vertex_show[E2.v1] = true;
                        debug_vertex_show[E2.v2] = true;
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

            return result;
        }

        vec3 mesh_facet_vertex(index_t f, index_t lv) const {
            index_t v = mesh_.facets.vertex(f,lv);
            return vec3(mesh_.vertices.point_ptr(v));
        }

        vec2 mesh_facet_vertex_project(index_t f, index_t lv) const {
            index_t v = mesh_.facets.vertex(f,lv);
	    const double* p = mesh_.vertices.point_ptr(v);
            return vec2(p[u_],p[v_]);
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
            f1_ = index_t(-1);
        }

        void log_err() const {
            std::cerr << "Houston, we got a problem (while remeshing facet "
                      << f1_ << "):" << std::endl;
        }

    private:
        Mesh& mesh_;
        index_t f1_;
        index_t latest_f2_;
        index_t latest_f2_count_;
        coord_index_t f1_normal_axis_;
        coord_index_t u_; // = (f1_normal_axis_ + 1)%3
        coord_index_t v_; // = (f1_normal_axis_ + 2)%3
        Sign f1_orient_;
        vector<Vertex> vertex_;
        vector<Edge> edges_;
        std::map<vec3Q, index_t, vec3QLexicoCompare> g_v_table_;
        std::map<trindex, index_t> t_v_table_;
        bool check_cnstr_;
        bool barycentric_;
        bool has_planar_isect_;
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
                    // then add individual vertices and
                    // do not add edge (because
                    // edges in macro edges will be
                    // created by sorting vertices
                    // along triangle edges).
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
            if(R1 == v1 || R1 == v2) {
                return true;
            }
        }
        return false;
    }
    
    bool planar_intersection_is_valid(const IsectInfo& II) {
        bool result = 
            same_region(II.A_rgn_f1, II.B_rgn_f1) ||
            same_region(II.A_rgn_f2, II.B_rgn_f2) ;
        return result;
    }
    
    void intersect_surface_impl(Mesh& M, const MeshIntersectionParams& params) {

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
        Process::enable_FPE(params.debug_enable_FPE);
        
        vector<IsectInfo> intersections;
        MeshFacetsAABB AABB(M,!params.debug_do_not_order_facets);
        vector<TriangleIsect> I;
        AABB.compute_facet_bbox_intersections(
            [&](index_t f1, index_t f2) {

                if(f1 == f2) {
                    return;
                }
                
                if(
                    !params.detect_intersecting_neighbors && (
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
        TM.set_check_constraints(params.debug_check_constraints);
        TM.set_barycentric(params.barycentric);
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

namespace GEO {
    
    void mesh_intersect_surface(Mesh& M, const MeshIntersectionParams& params) {
        if(!M.facets.are_simplices()) {
            tessellate_facets(M,3);
        }

        if(params.pre_detect_duplicated_vertices) {
            mesh_colocate_vertices_no_check(M);
        }
        
        if(params.pre_detect_duplicated_facets) {        
            mesh_remove_bad_facets_no_check(M);
        }

        const double SCALING = double(1ull << 20);
        const double INV_SCALING = 1.0/SCALING;

        // Pre-scale everything by 2^20 to avoid underflows
        // (note: this just adds 20 to the exponents of all
        //  coordinates).
        {
            double* p = M.vertices.point_ptr(0);
            index_t N = M.vertices.nb() *
                        M.vertices.dimension();
            for(index_t i=0; i<N; ++i) {
                p[i] *= SCALING;
            }
        }
        
        intersect_surface_impl(M, params);

        if(params.post_connect_facets) {
            mesh_repair(
                M,
                GEO::MeshRepairMode(
                    GEO::MESH_REPAIR_COLOCATE | GEO::MESH_REPAIR_DUP_F
                ),
                0.0
            );
        }

        // Scale-back everything
        {
            double* p = M.vertices.point_ptr(0);
            index_t N = M.vertices.nb() *
                        M.vertices.dimension();
            for(index_t i=0; i<N; ++i) {
                p[i] *= INV_SCALING;
            }
        }
    }
}
