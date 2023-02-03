
/*
 *
 */

// TODO: co-linear constraints intersections

#include <OGF/Experiment/algo/mesh_surface_intersection.h>
#include <OGF/Experiment/algo/exact_geometry.h>
#include <geogram/mesh/triangle_intersection.h>
#include <geogram/mesh/mesh_AABB.h>
#include <geogram/mesh/mesh_repair.h>
#include <geogram/mesh/mesh_fill_holes.h>
#include <geogram/mesh/index.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/delaunay/delaunay.h>
#include <geogram/numerics/predicates.h>
#include <geogram/numerics/expansion_nt.h>

#include <sstream>
#include <stack>

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

        /**
         * \brief An edge of the mesh.
         * \details It represents the constraints to be used by the
         *  constrained triangulation to remesh the facet. It contains
         *  a list of vertices coming from the intersection with other
         *  constrained edge.
         */
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
            
            void add_vertex_in_edge(index_t v) {
                if(v == v1 || v == v2) {
                    return;
                }
                for(index_t w: vertices_in_edge) {
                    if(w == v) {
                        return;
                    }
                }
                vertices_in_edge.push_back(v);
            }

            // The two extremities of the constraint.
            index_t v1;
            index_t v2;

            // Symbolic information: this edge was generated
            // by the intersection with f2's region R2.
            struct {
                index_t        f2;
                TriangleRegion R2;
            } sym;

            // The indices of all vertices that fall in this edge.
            vector<index_t> vertices_in_edge;
        };

        /***************************************************************/

        /**
         * \brief A vertex of the triangulation
         * \details Stores geometric information in exact precision, both
         *  in 3D and in local 2D coordinates. Local 2D coordinates can 
         *  be either projected or barycentric, depending on MeshInTriangle's
         *  barycentric_ flag. It also stores symbolic information, that is,
         *  facet indices and regions that generated the vertex.
         */
        class Vertex {
        public:
            enum { NO_INDEX = index_t(-1) };
            enum Type {
                UNINITIALIZED, MESH_VERTEX, PRIMARY_ISECT, SECONDARY_ISECT
            };

            /**
             * \brief Constructor for macro-triangle vertices.
             * \param[in] f facet index, supposed to correspond to
             *  MeshInTriangle's current facet
             * \param[in] lv local vertex index in \p f
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
             * \param[in] f1 , f2 the two facets. \p f1 is suposed to
             *  correspond to MeshInTriangle's current facet
             * \param[in] R1 , R2 the two facet regions.
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
             * \param[in] point_exact_in exact 3D coordinates 
             *   of the intersection
             * \param[in] UV_exact_in exact 2D coordinates. They can be either
             *   projected or barycentric, depending on MeshInTriangle's 
             *   barycentric_ flag.
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

            /**
             * \brief Default constructor
             */
            Vertex() {
                type = UNINITIALIZED;                
                mesh_in_triangle = nullptr;
                init_sym(NO_INDEX, NO_INDEX, T1_RGN_T, T2_RGN_T);
                mesh_vertex_index = NO_INDEX;
            }

            /**
             * \brief Gets the mesh
             * \return a reference to the mesh
             */
            Mesh& mesh() const {
                return mesh_in_triangle->mesh();
            }

            /**
             * \brief Prints this vertex
             * \details Displays the combinatorial information
             * \param[out] out an optional stream where to print
             */
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

            /**
             * \brief Gets a string representation of this Vertex
             * \return a string with the combinatorial information 
             *  of this Vertex
             */
            std::string to_string() const {
                std::ostringstream out;
                print(out);
                return out.str();
            }
            
        protected:

            /**
             * \brief Initializes the symbolic information of this Vertex
             * \param[in] f1 , f2 the two facets. \p f1 is suposed to
             *  correspond to MeshInTriangle's current facet
             * \param[in] R1 , R2 the two facet regions.
             */
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

            /**
             * \brief Initializes the geometry of this vertex
             * \details Computes the exact 3D and 2D position of this vertex
             *  based on the mesh and the combinatorial information
             */
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
                        geo_debug_assert(D.sign() != ZERO);
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
                        geo_debug_assert(d.sign() != ZERO);
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
                    geo_debug_assert(d.sign() != ZERO);
                    
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

            /**
             * \brief Indicates whether barycentric or projected coordinates
             *  are used.
             * \retval true if barycentric coordinates are used
             * \retval false otherwise
             */
            bool barycentric() const {
                return mesh_in_triangle->barycentric_;
            }

            /**
             * \brief Gets the coordinate used for the U axis
             * \return one of 0,1,2
             */
            coord_index_t u_coord() const {
                return mesh_in_triangle->u_;
            }

            /**
             * \brief Gets the coordinate used for the V axis
             * \return one of 0,1,2
             */
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

        Mesh& mesh() const {
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
                    /*
                    if(!has_planar_isect_) {
                        std::cerr << f1_
                                  << ":switching to planar isect mode"
                                  << std::endl;
                    }
                    */
                    has_planar_isect_ = true;
                }
            } else {
                latest_f2_ = f2;
                latest_f2_count_ = 0;
            }

            // If vertex is a macro-vertex, return it directly.
            if(region_dim(R1) == 0) {
                return index_t(R1);
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
                    edges_[e].add_vertex_in_edge(v);
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

            // If both extremities are on the same edge of f1,
            // we do not add the edge, because it will be generating
            // when remeshing the edge of f1
            if(region_dim(regions_convex_hull(AR1,BR1)) == 1) {
                return;
            }

            // The combinatorial information of the edge indicates whether
            // both extremities are on the same edge of f2
            edges_.push_back(Edge(v1,v2,f2,regions_convex_hull(AR2,BR2)));
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
                index_t new_t = mesh_.facets.create_triangle(i,j,k);
                // Copy all attributes from initial facet
                mesh_.facets.attributes().copy_item(new_t, f1_);
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
                       vec2Q p_exact = u_P1P2_plus_v_P1P3(uv.x,uv.y,p1,p2,p3);
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

            // Note: intersections that fall exactly
            // on vertices already present in this TriangleInMesh
            // are properly handled by two mechanisms:
            // - MeshInTriangle::add_vertex() compares the new
            //   vertex with all previous ones (and returns the
            //   existing vertex if one was found at the same position)
            // - Edge::add_vertex_in_edge(v) compares v with edge extremities
            //   and vertices already present in edge.
            // These two mechanisms (and the two nested edge-edge loops) make
            // the overall algorithm O(n^2). Not a big drama if what happens
            // within a single triangle stays reasonable...
            //
            // (we could capture some of them earlier, using the orient2d()
            //  computed by segment-segment intersection predicate).
            
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
                    
                    // TODO: check and handle co-linear edges
                    if(edge_edge_intersect(v1,v2,w1,w2)) {
                        vec3Q I;
                        vec2Q UV;
                        get_edge_edge_intersection(e1,e2,I,UV);
                        index_t x = add_vertex(Vertex(this,I,UV));
                        edges_[e1].add_vertex_in_edge(x);
                        edges_[e2].add_vertex_in_edge(x);
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
                //   Well, normally we already detected them when
                // inserting the points in the edges
                // (Edge::add_to_vertices_in_edge())
                V.erase(
                    std::unique(V.begin(), V.end()), V.end()
                );
            }

            // Step 3: "remesh" the edges that have intersections
            index_t nE = edges_.size();
            for(index_t e=0; e<nE; ++e) {
                if(edges_[e].vertices_in_edge.size() == 0) {
                    continue;
                }
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

        void get_edge_edge_intersection(
            index_t e1, index_t e2, vec3Q& I, vec2Q& UV
        ) const {
            index_t f1 = f1_;
            index_t f2 = edges_[e1].sym.f2; 
            index_t f3 = edges_[e2].sym.f2; 
            
            geo_assert(f1 != index_t(-1));
            geo_assert(f2 != index_t(-1));
            geo_assert(f3 != index_t(-1));                        
            
            vec3 P[9] = {
                mesh_facet_vertex(f1,0), mesh_facet_vertex(f1,1),
                mesh_facet_vertex(f1,2),
                mesh_facet_vertex(f2,0), mesh_facet_vertex(f2,1),
                mesh_facet_vertex(f2,2),
                mesh_facet_vertex(f3,0), mesh_facet_vertex(f3,1),
                mesh_facet_vertex(f3,2)
            };

            if(!get_three_planes_intersection(
                    I,
                    P[0], P[1], P[2],
                    P[3], P[4], P[5],
                    P[6], P[7], P[8]
            )) {
                get_edge_edge_intersection_2D(e1,e2,I,UV);
                return;
            }
            
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
                geo_debug_assert(d.sign() != ZERO);
                UV.x = rational_nt(det(B,C2),d);
                UV.y = rational_nt(det(C1,B),d);
            } else {
                UV.x = I[u_];
                UV.y = I[v_];
            }
        }             

        void get_edge_edge_intersection_2D(
            index_t e1, index_t e2, vec3Q& I, vec2Q& UV
        ) const {
            const Edge& E1 = edges_[e1];
            const Edge& E2 = edges_[e2];
            geo_assert(region_dim(E1.sym.R2) == 1);
            geo_assert(region_dim(E2.sym.R2) == 1);
            index_t le1 = index_t(E1.sym.R2)-index_t(T2_RGN_E0);
            index_t le2 = index_t(E2.sym.R2)-index_t(T2_RGN_E0);
            geo_assert(le1 < 3);
            geo_assert(le2 < 3);

            vec2 p1_uv = mesh_facet_vertex_project(E1.sym.f2, (le1+1)%3);
            vec2 p2_uv = mesh_facet_vertex_project(E1.sym.f2, (le1+2)%3);
            vec2 q1_uv = mesh_facet_vertex_project(E2.sym.f2, (le2+1)%3);
            vec2 q2_uv = mesh_facet_vertex_project(E2.sym.f2, (le2+2)%3);
            vec2E C1 = make_vec2<vec2E>(p1_uv, p2_uv);
            vec2E C2 = make_vec2<vec2E>(q2_uv, q1_uv);
            vec2E B  = make_vec2<vec2E>(p1_uv, q1_uv);
            
            expansion_nt d = det(C1,C2);
            geo_debug_assert(d.sign() != ZERO);
            rational_nt t(det(B,C2),d);
            I = mix(
                t,
                mesh_facet_vertex(E1.sym.f2,(le1+1)%3),
                mesh_facet_vertex(E1.sym.f2,(le1+2)%3)
            );
            
            if(barycentric_) {
                vec2 a1_uv = mesh_facet_vertex_project(f1_,0);
                vec2 a2_uv = mesh_facet_vertex_project(f1_,1);
                vec2 a3_uv = mesh_facet_vertex_project(f1_,2);
                vec2E a12_uv = make_vec2<vec2E>(a1_uv, a2_uv);
                vec2E a13_uv = make_vec2<vec2E>(a1_uv, a3_uv);
                vec2E p12_uv = make_vec2<vec2E>(p1_uv, p2_uv);
                vec2E q12_uv = make_vec2<vec2E>(q1_uv, q2_uv);
                vec2E C1(det(a12_uv, p12_uv), det(a12_uv, q12_uv));
                vec2E C2(det(a13_uv, p12_uv), det(a13_uv, q12_uv));
                vec2E B(
                    det(make_vec2<vec2E>(a1_uv, p1_uv), p12_uv),
                    det(make_vec2<vec2E>(a1_uv, q1_uv), q12_uv)
                );
                expansion_nt d = det(C1,C2);
                geo_debug_assert(d.sign() != ZERO);
                UV.x = rational_nt(det(B,C2),d);
                UV.y = rational_nt(det(C1,B),d);                
            } else {
                UV.x = I[u_];
                UV.y = I[v_];
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
                value_type(p[u_]),
                value_type(p[v_])
            );
        }
            
        template <class VEC = vec2> VEC
        mesh_facet_vertex_project(index_t f, index_t lv) const {
            index_t v = mesh().facets.vertex(f,lv);
            return mesh_vertex_project<VEC>(v);
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

    /**
     * \brief Computes the intersection between two mesh triangular facets
     * \details This function is just a wrapper around triangles_intersections()
     *  for Mesh facets.
     * \param[in] M the mesh
     * \param[in] f1 , f2 the two facets
     * \param[out] I a vector of triangle intersections
     * \retval true if there was an intersection
     * \retval false otherwise
     */
    bool mesh_facets_intersect(
        Mesh& M, index_t f1, index_t f2, vector<TriangleIsect>& I
    ) {
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

    void mesh_intersect_surface_compute_arrangement(
        Mesh& M, const MeshSurfaceIntersectionParams& params
    ) {

        // Step 1: Preparation
        // -------------------
        
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

        // Step 2: Get intersections
        // -------------------------
        
        vector<IsectInfo> intersections;
        MeshFacetsAABB AABB(M,!params.debug_do_not_order_facets);
        vector<TriangleIsect> I;
        AABB.compute_facet_bbox_intersections(
            [&](index_t f1, index_t f2) {

                // Needed (maybe I should change that in AABB class)
                if(f1 == f2) {
                    return;
                }

                // Optionally skip facet pairs that share a vertex or an edge
                if(
                    !params.detect_intersecting_neighbors && (
                        (M.facets.find_adjacent(f1,f2)      != index_t(-1)) ||
                        (M.facets.find_common_vertex(f1,f2) != index_t(-1))
                    )
                ) {
                    return;
                }
                
                if(mesh_facets_intersect(M, f1, f2, I)) {

                    if(I.size() > 2) {
                        // Coplanar intersection: to generate the edges,
                        // test validity of all possible pairs of vertices.
                        for(index_t i1=0; i1< I.size(); ++i1) {
                            for(index_t i2=0; i2<i1; ++i2) {
                                IsectInfo II = {
                                    f1, f2,
                                    I[i1].first, I[i1].second,
                                    I[i2].first, I[i2].second
                                };

                                // Valid edges are the ones where both
                                // extremities are on the same edge of f1
                                // or on the same edge of f2
                                
                                TriangleRegion AB1 =
                                  regions_convex_hull(II.A_rgn_f1, II.B_rgn_f1);
                                
                                TriangleRegion AB2 =
                                  regions_convex_hull(II.A_rgn_f2, II.B_rgn_f2);
                                
                                if(
                                    region_dim(AB1) == 1 ||
                                    region_dim(AB2) == 1
                                ) {
                                    intersections.push_back(II);
                                    II.flip();
                                    intersections.push_back(II);
                                } 
                            }
                        }
                    } else {
                        // Intersection is either a segment
                        // or a vertex of f2.
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
            }
        );

        // Step 3: Remesh intersected triangles
        // ------------------------------------

        // The MeshInTriangle, that implements local triangulation
        // in each intersected triangular facet. It also keeps a global map of 
        // vertices, indexed by their exact geometry.         
        MeshInTriangle TM(M);
        TM.set_check_constraints(params.debug_check_constraints);
        TM.set_barycentric(params.barycentric);

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

        // Now iterate on all intersections, and identify the [b,e[ intervals
        // that correspond to the same f1 facet.
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

                // Each IsectIfo is either an individual vertex
                // or a segment with two vertices.
                // Each vertex is represented combinatorially.
                // The MeshInTriangle knows how to compute the
                // geometry from the combinatorial information.
                
                if(II.is_point()) {
                    TM.add_vertex(
                        II.f2,
                        II.A_rgn_f1, II.A_rgn_f2
                    );
                } else {
                    TM.add_edge(
                        II.f2, 
                        II.A_rgn_f1, II.A_rgn_f2,
                        II.B_rgn_f1, II.B_rgn_f2
                    );
                }
            }
            TM.end_facet();
            b = e;
        }

        // Step 4: Epilogue
        // ----------------
        
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

    /*****************************************************************/

    index_t compute_charts(Mesh& M, const std::string& attribute = "chart") {
        Attribute<index_t> chart(M.facets.attributes(), attribute);
        for(index_t f: M.facets) {
            chart[f] = index_t(-1);
        }
        std::stack<index_t> S;
        index_t cur_chart = 0;
        for(index_t f: M.facets) {
            if(chart[f] == index_t(-1)) {
                chart[f] = cur_chart;
                S.push(f);
                while(!S.empty()) {
                    index_t g = S.top();
                    S.pop();
                    for(
                        index_t le=0;
                        le<M.facets.nb_vertices(g); ++le
                    ) {
                        index_t h = M.facets.adjacent(g,le);
                        if(h != index_t(-1) && chart[h] == index_t(-1)) {
                            chart[h] = cur_chart;
                            S.push(h);
                        }
                    }
                }
                ++cur_chart;
            }
        }
        return cur_chart;
    }
}




namespace GEO {
    
    void mesh_intersect_surface(
        Mesh& M, const MeshSurfaceIntersectionParams& params
    ) {
        if(!M.facets.are_simplices()) {
            tessellate_facets(M,3);
        }

        Attribute<index_t> operand_bit;
        operand_bit.bind_if_is_defined(M.facets.attributes(), "operand_bit");
        if(!operand_bit.is_bound()) {
            compute_charts(M,"operand_bit");
            operand_bit.bind(M.facets.attributes(), "operand_bit");
            for(index_t f: M.facets) {
                operand_bit[f] = (1 << operand_bit[f]);
            }
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
        
        mesh_intersect_surface_compute_arrangement(M, params);

        if(params.post_connect_facets) {
            /*
            mesh_colocate_vertices_no_check(M);
            mesh_remove_bad_facets_no_check(M);
            mesh_connect_and_reorient_facets_no_check(M);
            */
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

namespace {
    using namespace GEO;

    /**
     * \brief A simple parser for boolean expressions
     * \details 
     *  - Variables: A..Z or x0..x31 
     *  - and:        '&' or '*'
     *  - or:         '|' or '+'
     *  - xor:        '^'
     *  - difference: '-'
     */
    class BooleanExprParser {
    public:
        BooleanExprParser(
            const std::string& expr
        ) : expr_(expr) {
        }

        bool eval(index_t x) {
            x_   = x;
            ptr_ = expr_.begin();
            return parse_or();
        }

    protected:

        bool parse_or() {
            bool left = parse_and();
            while(
                cur_char() == '|' ||
                cur_char() == '^' ||
                cur_char() == '+' ||
                cur_char() == '-'
            ) {
                char op = cur_char();
                next_char();
                bool right = parse_and();
                left = (op == '-') ? (left && !right) :
                       (op == '^') ? (left ^   right) :
                                     (left ||  right) ;
            }
            return left;
        }

        bool parse_and() {
            bool left = parse_factor();
            while(cur_char() == '&' || cur_char() == '*') {
                next_char();
                bool right = parse_factor();
                left = left && right;
            }
            return left;
        }

        bool parse_factor() {
            if(cur_char() == '!' || cur_char() == '~' || cur_char() == '-') {
                next_char();
                return !parse_factor();
            }
            if(cur_char() == '(') {
                next_char();
                bool result = parse_or();
                if(cur_char() != ')') {
                    throw std::logic_error(
                        std::string("Unmatched parenthesis: ")+cur_char()
                    );
                }
                next_char();
                return result;
            }
            if((cur_char() >= 'A' && cur_char() <= 'Z') || cur_char() == 'x') {
                return parse_variable();
            }
            throw std::logic_error("Syntax error");
        }

        bool parse_variable() {
            int bit = 0;
            if(cur_char() >= 'A' && cur_char() <= 'Z') {
                bit = int(cur_char()) - int('A');
                next_char();
            } else {
                if(cur_char() != 'x') {
                    throw std::logic_error("Syntax error in variable");
                }
                next_char();
                while(cur_char() >= '0' && cur_char() <= '9') {
                    bit = bit * 10 + (int(cur_char()) - '0');
                    next_char();
                }
            }
            if(bit > 31) {
                throw std::logic_error("Bit larger than 31");
            }
            return ((x_ & (1 << bit)) != 0);
        }

        char cur_char() const {
            return *ptr_;
        }
        
        void next_char() {
            if(ptr_ == expr_.end()) {
                throw std::logic_error("Unexpected end of string");
            }
            ptr_++;
        }
        
    private:
        std::string expr_;
        std::string::iterator ptr_;
        index_t x_;
    };
}

namespace GEO {
    
    void mesh_classify_intersections(
        Mesh& M, std::function<bool(index_t)> eqn, const std::string& attribute
    ) {
        MeshFacetsAABB AABB(M);
        index_t nb_charts = compute_charts(M);
        Attribute<index_t> chart(M.facets.attributes(), "chart");
        Attribute<index_t> operand_bit(M.facets.attributes(), "operand_bit");
        Attribute<bool> selection(M.facets.attributes(), attribute);
        vector<index_t> chart_facet(nb_charts, index_t(-1));
        for(index_t f: M.facets) {
            index_t c = chart[f];
            if(chart_facet[c] == index_t(-1)) {
                chart_facet[c] = f;
                index_t parity = 0;
                vec3 D(
                    Numeric::random_float64(),
                    Numeric::random_float64(),
                    Numeric::random_float64()
                );
                vec3 g = Geom::mesh_facet_center(M,f);
                AABB.ray_all_intersections(
                    Ray(g,D),
                    [&](const MeshFacetsAABB::Intersection & I) {
                        if(I.f != f) {
                            parity = parity ^ operand_bit[I.f];
                        }
                    }
                );
                try {
                    selection[f] =
                        eqn(parity | operand_bit[f]) !=
                        eqn(parity & ~operand_bit[f] ) ;
                } catch(const std::logic_error& e) {
                    Logger::err("Classify") << "Error while parsing expression:"
                                            << e.what()
                                            << std::endl;
                    return;
                }
            } else {
                selection[f] = selection[chart_facet[c]];
            }
        }
    }

    void mesh_classify_intersections(
        Mesh& M, const std::string& expr, const std::string& attribute
    ) {
        BooleanExprParser eqn(expr);
        index_t operand_all_bits;
        {
            index_t max_operand_bit = 0;
            Attribute<index_t> operand_bit;
            operand_bit.bind_if_is_defined(M.facets.attributes(),"operand_bit");
            if(!operand_bit.is_bound()) {
                Logger::err("Classify")
                    << "operand_bit: no such facet attribute"
                    << std::endl;
                return;
            }
            for(index_t f: M.facets) {
                max_operand_bit = std::max(max_operand_bit, operand_bit[f]);
            }
            operand_all_bits = (max_operand_bit << 1)-1;
        }
        
        
        try {
            mesh_classify_intersections(
                M,
                [&](index_t x)->bool {
                    return
                        (expr == "union")        ? (x != 0)                :
                        (expr == "intersection") ? (x == operand_all_bits) : 
                        eqn.eval(x);
                },
                attribute
            );
        } catch(const std::logic_error& e) {
            Logger::err("Classify") << "Error while parsing expression:"
                                    << e.what()
                                    << std::endl;
            return;
        }
    }    
}
