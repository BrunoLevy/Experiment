
/*
 *
 */

// TODO: exact classification
// TODO: expansion::optimize(), when do we need to call it ?
//       Would be better to call it internally when needed.
//       Influence on performance is important.
// TODO: would be good to be able to do computations with
//       points that have their coordinates stored in the
//       stack.



#include <OGF/Experiment/algo/mesh_surface_intersection.h>
#include <OGF/Experiment/algo/exact_geometry.h>

#include <geogram/mesh/triangle_intersection.h>
#include <geogram/mesh/mesh_AABB.h>
#include <geogram/mesh/mesh_repair.h>
#include <geogram/mesh/mesh_fill_holes.h>
#include <geogram/mesh/index.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/mesh/mesh_geometry.h>
#include <geogram/delaunay/CDT_2d.h>
#include <geogram/numerics/predicates.h>
#include <geogram/numerics/expansion_nt.h>
#include <geogram/basic/stopwatch.h>

#include <sstream>
#include <stack>

namespace GEO {
    void GEOGRAM_API SOS_sort(const double** begin, const double** end, GEO::index_t dim);
}

namespace {
    using namespace GEO;

    /***********************************************************************/

    /**
     * \brief Removes all the triangles with their three vertices aligned
     */
    void remove_linear_triangles(Mesh& M) {
        vector<index_t> remove_f(M.facets.nb());
        for(index_t f: M.facets) {
            index_t v1 = M.facets.vertex(f,0);
            index_t v2 = M.facets.vertex(f,1);
            index_t v3 = M.facets.vertex(f,2);            
            const double* p1 = M.vertices.point_ptr(v1);
            const double* p2 = M.vertices.point_ptr(v2);
            const double* p3 = M.vertices.point_ptr(v3);
            remove_f[f] = PCK::aligned_3d(p1,p2,p3);
        }
        M.facets.delete_elements(remove_f);
    }

    /***********************************************************************/
    
    /**
     * \brief Meshes a single triangle with the constraints that come from
     *  the intersections with the other triangles.
     */
    class MeshInTriangle : public CDTBase2d {
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
            
            // The two extremities of the constraint.
            index_t v1;
            index_t v2;

            // Symbolic information: this edge was generated
            // by the intersection with f2's region R2.
            struct {
                index_t        f2;
                TriangleRegion R2;
            } sym;
        };

        /***************************************************************/

        /**
         * \brief A vertex of the triangulation
         * \details Stores geometric information in exact precision, both
         *  in 3D and in local 2D coordinates. It also stores symbolic information, 
         *  that is, facet indices and regions that generated the vertex.
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
            Vertex(
                MeshInTriangle* M, index_t f, index_t lv
            ) :
                point_exact(vec3HE_noinit()),
                UV_exact(vec2HE_noinit())
            {
                geo_assert(f == M->f1_);
                type = MESH_VERTEX;
                mesh_in_triangle = M;
                init_sym(f, NO_INDEX, TriangleRegion(lv), T2_RGN_T);
                init_geometry();
                post_init_geometry();
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
            ) :
                point_exact(vec3HE_noinit()),
                UV_exact(vec2HE_noinit())
            {
                geo_assert(f1 == M->f1_);                
                type = PRIMARY_ISECT;
                mesh_in_triangle = M;
                init_sym(f1,f2,R1,R2);
                init_geometry();
                post_init_geometry();
            }

            /**
             * \brief Constructor for intersections between constraints.
             * \param[in] point_exact_in exact 3D coordinates 
             *   of the intersection
             * \param[in] UV_exact_in exact 2D coordinates. 
             */
            Vertex(
                MeshInTriangle* M,
                const vec3HE& point_exact_in,
                const vec2HE& UV_exact_in
            ) : point_exact(point_exact_in),
                UV_exact(UV_exact_in)
            {
                type = SECONDARY_ISECT;                
                mesh_in_triangle = M;
                init_sym(NO_INDEX, NO_INDEX, T1_RGN_T, T2_RGN_T);
                post_init_geometry();
            }

            /**
             * \brief Default constructor
             */
            Vertex() :
                point_exact(vec3HE_noinit()),
                UV_exact(vec2HE_noinit())
            {
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

                mesh_vertex_index = NO_INDEX;
                
                // Case 1: f1 vertex
                if(region_dim(sym.R1) == 0) {
                    index_t lv = index_t(sym.R1);
                    geo_assert(lv < 3);
                    mesh_vertex_index = mesh().facets.vertex(sym.f1,lv);
                    point_exact = mesh_vertex_vec3HE(mesh_vertex_index);
                    UV_exact = mesh_vertex_project_vec2HE(mesh_vertex_index);
                    return;
                }

                geo_assert(sym.f1 != NO_INDEX && sym.f2 != NO_INDEX);

                // Case 2: f2 vertex
                if(region_dim(sym.R2) == 0) {
                    index_t lv = index_t(sym.R2)-3;
                    geo_assert(lv < 3);
                    mesh_vertex_index = mesh().facets.vertex(sym.f2, lv);
                    point_exact = mesh_vertex_vec3HE(mesh_vertex_index);
                    UV_exact = mesh_vertex_project_vec2HE(mesh_vertex_index);
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
                        UV_exact.x = point_exact[u_coord()];
                        UV_exact.y = point_exact[v_coord()];
                        UV_exact.w = point_exact.w;
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
                    UV_exact.x = point_exact[u_coord()];
                    UV_exact.y = point_exact[v_coord()];
                    UV_exact.w = point_exact.w;
                    return;
                }

                // case 5: f1 edge /\ f2 edge in 2D
                if(region_dim(sym.R1) == 1 && region_dim(sym.R2) == 1) {
                    index_t e1 = index_t(sym.R1) - index_t(T1_RGN_E0);
                    geo_debug_assert(e1 < 3);
                    index_t e2 = index_t(sym.R2) - index_t(T2_RGN_E0);
                    geo_debug_assert(e2 < 3);

                    vec2 p1 = mesh_facet_vertex_project(
                        sym.f1, (e1+1)%3
                    );
                    vec2 p2 = mesh_facet_vertex_project(
                        sym.f1, (e1+2)%3
                    );

                    vec2 q1 = mesh_facet_vertex_project(
                        sym.f2, (e2+1)%3
                    );
                    vec2 q2 = mesh_facet_vertex_project(
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
                    UV_exact.x = point_exact[u_coord()];
                    UV_exact.y = point_exact[v_coord()];
                    UV_exact.w = point_exact.w;
                    return;
                }

                // Normally we enumerated all possible cases
                geo_assert_not_reached;
            }

            /**
             * \brief Optimizes exact numbers in generated
             *  points and computes approximate coordinates.
             */
            void post_init_geometry() {
                point_exact.x.optimize();
                point_exact.y.optimize();
                point_exact.z.optimize();
                point_exact.w.optimize();
                UV_exact.x.optimize();
                UV_exact.y.optimize();
                UV_exact.w.optimize();
                if(mesh_in_triangle->delaunay_) {
                    double w = UV_exact.w.estimate();
                    UV_approx.x = UV_exact.x.estimate()/w;
                    UV_approx.y = UV_exact.y.estimate()/w;                    
                }
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
            
            vec3HE mesh_vertex_vec3HE(index_t v) const {
                const double* p = mesh().vertices.point_ptr(v);
                return vec3HE(
                    expansion_nt(p[0]),
                    expansion_nt(p[1]),
                    expansion_nt(p[2]),
                    expansion_nt(1.0)
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

            vec2HE mesh_vertex_project_vec2HE(index_t v) const {
                const double* p = mesh().vertices.point_ptr(v);
                return vec2HE(
                    expansion_nt(p[u_coord()]),
                    expansion_nt(p[v_coord()]),
                    expansion_nt(1.0)
                );
            }
            
            template <class VEC = vec2> VEC
            mesh_facet_vertex_project(index_t f, index_t lv) const {
                index_t v = mesh().facets.vertex(f,lv);
                return mesh_vertex_project<VEC>(v);
            }
            
        public:
            MeshInTriangle* mesh_in_triangle;

            vec3HE point_exact;
            vec2HE UV_exact;
            vec2   UV_approx;

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
            f1_(index_t(-1))
        {
        }

        Mesh& mesh() const {
            return mesh_;
        }

        void save_constraints(const std::string& filename) {
            Mesh M;
            get_constraints(M);
            mesh_save(M,filename);
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
                vertex_.push_back(Vertex(this, f, lv));
            }

            CDTBase2d::create_enclosing_triangle(0,1,2);

            edges_.push_back(Edge(1,2));
            edges_.push_back(Edge(2,0));
            edges_.push_back(Edge(0,1));
            
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

            // Create the vertex
            vertex_.push_back(Vertex(this, f1_, f2, R1, R2));
            // Insert it into the triangulation

            index_t v = CDTBase2d::insert(vertex_.size()-1);
            
            // If it was an existing vertex, return the existing vertex
            if(vertex_.size() > CDTBase2d::nv()) {
                vertex_.pop_back();
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
            // we do not add the edge, because it will be generated
            // when remeshing the edge of f1
            if(region_dim(regions_convex_hull(AR1,BR1)) == 1) {
                return;
            }

            // Generate also the combinatorial information of the edge,
            // that indicates whether both extremities are on the same
            // edge of f2 (useful later to compute the intersections)
            edges_.push_back(Edge(v1,v2,f2,regions_convex_hull(AR2,BR2)));

            // Constraints will be added to the triangulation during commit()
        }

        void end_facet() {
            commit();
            if(false && edges_.size() > 1000) {
                save_constraints(
                    "constraints_" + String::to_string(f1_) + ".geogram"
                );
                save("triangles_" + String::to_string(f1_) + ".geogram");
            }
            clear();
        }

    protected:

        void commit() {
            if(false) {
                Mesh constraints;
                get_constraints(constraints);
                mesh_save(constraints, "constraints.geogram");
            }

            for(const Edge& E: edges_) {
                CDTBase2d::insert_constraint(E.v1, E.v2);
            }
            
            // Create vertices and facets in target mesh
            create_vertices();
            for(index_t t=0; t<CDTBase2d::nT(); ++t) {
                index_t i = CDTBase2d::Tv(t,0);
                index_t j = CDTBase2d::Tv(t,1);
                index_t k = CDTBase2d::Tv(t,2);
                i = vertex_[i].mesh_vertex_index;
                j = vertex_[j].mesh_vertex_index;
                k = vertex_[k].mesh_vertex_index;                    
                index_t new_t = mesh_.facets.create_triangle(i,j,k);
                // Copy all attributes from initial facet
                mesh_.facets.attributes().copy_item(new_t, f1_);
            }
        }
        

        /**
         * \brief Creates the vertices in the target mesh or
         *  finds them if they already exist.
         */
        void create_vertices() {
            for(index_t i=0; i<vertex_.size(); ++i) {
                // Vertex already exists in this MeshInTriangle
                if(vertex_[i].mesh_vertex_index != index_t(-1)) {
                    continue;
                }

                // Use exact geometry as key
                const vec3HE& K = vertex_[i].point_exact;
                auto it = g_v_table_.find(K);
                if(it != g_v_table_.end()) {
                    // Vertex alreay exists in target mesh
                    vertex_[i].mesh_vertex_index = it->second;
                } else {
                    // Vertex does not exist in target mesh,
                    // create it and update table
                    double w = vertex_[i].point_exact.w.estimate();
                    vec3 p(
                        vertex_[i].point_exact.x.estimate() / w,
                        vertex_[i].point_exact.y.estimate() / w,
                        vertex_[i].point_exact.z.estimate() / w
                    );
                    index_t v = mesh_.vertices.create_vertex(p.data());
                    vertex_[i].mesh_vertex_index = v;
                    g_v_table_[K] = v;
                }
            }
        }
        
        void get_constraints(Mesh& M, bool with_edges=true) const {
            if(M.vertices.nb() == 0) {
                M.vertices.set_dimension(2);
                for(index_t v=0; v<vertex_.size(); ++v) {
		   vec2 p;
                   double w = vertex_[v].UV_exact.w.estimate();
                   p = vec2(
                       vertex_[v].UV_exact.x.estimate() / w,
                       vertex_[v].UV_exact.y.estimate() / w
                   );
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
        
        void get_edge_edge_intersection(
            index_t e1, index_t e2, vec3HE& I, vec2HE& UV
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
            
            UV.x = I[u_];
            UV.y = I[v_];
            UV.w = I.w;
        }             

        void get_edge_edge_intersection_2D(
            index_t e1, index_t e2, vec3HE& I, vec2HE& UV
        ) const {
            const Edge& E1 = edges_[e1];
            const Edge& E2 = edges_[e2];

            if(
                region_dim(E1.sym.R2) == 1 &&
                region_dim(E2.sym.R2) == 1
            ) {
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
            
                UV.x = I[u_];
                UV.y = I[v_];
                UV.w = I.w;
            } else {
                geo_assert(region_dim(E1.sym.R2) == 1 || region_dim(E2.sym.R2) == 1);
                index_t f1 = E1.sym.f2;
                TriangleRegion R1 = E1.sym.R2;
                index_t f2 = E2.sym.f2;
                TriangleRegion R2 = E2.sym.R2;
                if(region_dim(R1) == 1) {
                    std::swap(f1,f2);
                    std::swap(R1,R2);
                }

                vec3 p1 = mesh_facet_vertex(f1,0);
                vec3 p2 = mesh_facet_vertex(f1,1);
                vec3 p3 = mesh_facet_vertex(f1,2);                
                
                index_t e = index_t(R2) - index_t(T2_RGN_E0);
                geo_assert(e < 3);

                vec3 q1 = mesh_facet_vertex(f2,(e+1)%3);
                vec3 q2 = mesh_facet_vertex(f2,(e+2)%3);

                vec3E D   = make_vec3<vec3E>(q1,q2);
                vec3E E1  = make_vec3<vec3E>(p1,p2);
                vec3E E2  = make_vec3<vec3E>(p1,p3);
                vec3E AO  = make_vec3<vec3E>(p1,q1);
                vec3E N   = cross(E1,E2);
                
                expansion_nt d = -dot(D,N);
                geo_debug_assert(d.sign() != ZERO);
                rational_nt t(dot(AO,N),d);
                I = mix(t,q1,q2);
                UV.x = I[u_];
                UV.y = I[v_];
                UV.w = I.w;
                return;
            }
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
        
        void clear() override {
            vertex_.resize(0);
            edges_.resize(0);
            f1_ = index_t(-1);
            CDTBase2d::clear();
        }

        void log_err() const {
            std::cerr << "Houston, we got a problem (while remeshing facet "
                      << f1_ << "):" << std::endl;
        }
        
    protected:

        /********************** CDTBase2d overrides ***********************/
        
        Sign orient2d(index_t v1,index_t v2,index_t v3) const override {
            return PCK::orient_2d(
                vertex_[v1].UV_exact,
                vertex_[v2].UV_exact,
                vertex_[v3].UV_exact
            );
        }

        /**
         * \brief Tests the relative position of a point with respect
         *  to the circumscribed circle of a triangle
         * \param[in] v1 , v2 , v3 the three vertices of the triangle
         *  oriented anticlockwise
         * \param[in] v4 the point to be tested
         * \retval POSITIVE if the point is inside the circle
         * \retval ZERO if the point is on the circle
         * \retval NEGATIVE if the point is outside the circle
         */
        Sign incircle(
            index_t v1,index_t v2,index_t v3,index_t v4
        ) const override {

            /*
            return PCK::in_circle_2d_SOS(
                vertex_[v1].UV_approx.data(),
                vertex_[v2].UV_approx.data(),
                vertex_[v3].UV_approx.data(),
                vertex_[v4].UV_approx.data()
            );
            */


            /*
            {
                const vec2& p1 = vertex_[v1].UV_approx;
                const vec2& p2 = vertex_[v2].UV_approx;
                const vec2& p3 = vertex_[v3].UV_approx;
                const vec2& p4 = vertex_[v4].UV_approx;            
                return PCK::orient_2dlifted_SOS(
                    p1.data(), p2.data(), p3.data(), p4.data(),
                    length2(p1), length2(p2), length2(p3), length2(p4)
                );
            }
            */


            if(false) {
                const vec2& p0 = vertex_[v1].UV_approx;
                const vec2& p1 = vertex_[v2].UV_approx;
                const vec2& p2 = vertex_[v3].UV_approx;
                const vec2& p3 = vertex_[v4].UV_approx;            
                
                double h0 = length2(vertex_[v1].UV_approx);
                double h1 = length2(vertex_[v2].UV_approx);
                double h2 = length2(vertex_[v3].UV_approx);
                double h3 = length2(vertex_[v4].UV_approx);            
                
                expansion_nt a11(expansion_nt::DIFF, p1.x, p0.x);
                expansion_nt a12(expansion_nt::DIFF, p1.y, p0.y);
                expansion_nt a13(expansion_nt::DIFF, h0, h1);
                
                expansion_nt a21(expansion_nt::DIFF, p2.x, p0.x);
                expansion_nt a22(expansion_nt::DIFF, p2.y, p0.y);
                expansion_nt a23(expansion_nt::DIFF, h0, h2);
                
                expansion_nt a31(expansion_nt::DIFF, p3.x, p0.x);
                expansion_nt a32(expansion_nt::DIFF, p3.y, p0.y);
                expansion_nt a33(expansion_nt::DIFF, h0, h3);

                expansion_nt Delta1 = det2x2(a21, a22, a31, a32);
                expansion_nt Delta2 = det2x2(a11, a12, a31, a32);
                expansion_nt Delta3 = det2x2(a11, a12, a21, a22);
                
                Sign Delta3_sign = Delta3.sign();
                geo_assert(Delta3_sign != ZERO);
                
                expansion_nt r = a13*Delta1-a23*Delta2+a33*Delta3;
                Sign r_sign = r.sign();
                
                if(r_sign == ZERO) {
                    const double* p_sort[4];
                    p_sort[0] = p0.data();
                    p_sort[1] = p1.data();
                    p_sort[2] = p2.data();
                    p_sort[3] = p3.data();
                    GEO::SOS_sort(p_sort, p_sort + 4, 3);
                    
                    for(index_t i = 0; i < 4; ++i) {
                        if(p_sort[i] == p0.data()) {
                            expansion_nt z = Delta2-Delta1+Delta3;
                            Sign z_sign = z.sign();
                            if(z_sign != ZERO) {
                                return Sign(Delta3_sign*z_sign);
                            }
                        } else if(p_sort[i] == p1.data()) {
                            Sign Delta1_sign = Delta1.sign();
                            if(Delta1_sign != ZERO) {
                                return Sign(Delta3_sign * Delta1_sign);
                            }
                        } else if(p_sort[i] == p2.data()) {
                            Sign Delta2_sign = Delta2.sign();
                            if(Delta2_sign != ZERO) {
                                return Sign(-Delta3_sign * Delta2_sign);
                            }
                        } else if(p_sort[i] == p3.data()) {
                            return NEGATIVE;
                        }
                    }
                }
                return Sign(Delta3_sign * r_sign);
            }

            if(true) {
                const vec2HE& p0 = vertex_[v1].UV_exact;
                const vec2HE& p1 = vertex_[v2].UV_exact;
                const vec2HE& p2 = vertex_[v3].UV_exact;
                const vec2HE& p3 = vertex_[v4].UV_exact;
                
                double h0 = length2(vertex_[v1].UV_approx);
                double h1 = length2(vertex_[v2].UV_approx);
                double h2 = length2(vertex_[v3].UV_approx);
                double h3 = length2(vertex_[v4].UV_approx);            

                expansion_nt a13(expansion_nt::DIFF, h0, h1);
                expansion_nt a23(expansion_nt::DIFF, h0, h2);
                expansion_nt a33(expansion_nt::DIFF, h0, h3);                
                
                vec2HE U1 = p1-p0;
                const expansion_nt& w1 = U1.w;
                Sign sw1 = w1.sign();
                
                vec2HE U2 = p2-p0;
                const expansion_nt& w2 = U2.w;
                Sign sw2 = w2.sign();

                vec2HE U3 = p3-p0;
                const expansion_nt& w3 = U3.w;
                Sign sw3 = w3.sign();                

                expansion_nt w2w3Delta1 = det2x2(U2.x, U2.y, U3.x, U3.y);
                expansion_nt w1w3Delta2 = det2x2(U1.x, U1.y, U3.x, U3.y);
                expansion_nt w1w2Delta3 = det2x2(U1.x, U1.y, U2.x, U2.y);
                
                Sign Delta3_sign = Sign(w1w2Delta3.sign()*sw1*sw2);
                geo_assert(Delta3_sign != ZERO);
                
                expansion_nt w1w2w3R = a13*w1*w2w3Delta1-a23*w2*w1w3Delta2+a33*w3*w1w2Delta3;
                Sign R_sign = Sign(w1w2w3R.sign()*sw1*sw2*sw3);

                // Simulation of simplicity
                if(R_sign == ZERO) {
                    const vec2HE* p_sort[4] = {
                        &vertex_[v1].UV_exact,
                        &vertex_[v2].UV_exact,
                        &vertex_[v3].UV_exact,
                        &vertex_[v4].UV_exact,
                    };
                    std::sort(
                        p_sort, p_sort+4,
                        [](const vec2HE* p1, const vec2HE* p2)->bool{
                            vec2HELexicoCompare cmp;
                            return cmp(*p1,*p2);
                        }
                    );
                    for(index_t i = 0; i < 4; ++i) {
                        if(p_sort[i] == &p0) {
                            expansion_nt w1w2w3Z = w2*w1w3Delta2-w1*w2w3Delta1+w3*w1w2Delta3;
                            Sign Z_sign = Sign(w1w2w3Z.sign()*sw1*sw2*sw3);
                            if(Z_sign != ZERO) {
                                return Sign(Delta3_sign*Z_sign);
                            }
                        } else if(p_sort[i] == &p1) {
                            Sign Delta1_sign = Sign(w2w3Delta1.sign()*sw2*sw3);
                            if(Delta1_sign != ZERO) {
                                return Sign(Delta3_sign * Delta1_sign);
                            }
                        } else if(p_sort[i] == &p2) {
                            Sign Delta2_sign = Sign(w1w3Delta2.sign()*sw1*sw3);
                            if(Delta2_sign != ZERO) {
                                return Sign(-Delta3_sign * Delta2_sign);
                            }
                        } else if(p_sort[i] == &p3) {
                            return NEGATIVE;
                        }
                    }
                }
                
                return Sign(Delta3_sign * R_sign);
            }

        }

        /**
         * \brief Given two segments that have an intersection, create the
         *  intersection
         * \details The intersection is given both as the indices of segment
         *  extremities (i,j) and (k,l), that one can use to retreive the 
         *  points in derived classes, and constraint indices E1 and E2, that
         *  derived classes may use to retreive symbolic information attached
         *  to the constraint
         * \param[in] E1 the index of the first edge, corresponding to the
         *  value of ncnstr() when insert_constraint() was called for
         *  that edge
         * \param[in] i , j the vertices of the first segment
         * \param[in] E2 the index of the second edge, corresponding to the
         *  value of ncnstr() when insert_constraint() was called for
         *  that edge
         * \param[in] k , l the vertices of the second segment
         * \return the index of a newly created vertex that corresponds to
         *  the intersection between [\p i , \p j] and [\p k , \p l]
         */
        index_t create_intersection(
            index_t e1, index_t i, index_t j,
            index_t e2, index_t k, index_t l
        ) override {
            geo_argused(i);
            geo_argused(j);
            geo_argused(k);
            geo_argused(l);

            vec3HE I  = vec3HE_noinit();
            vec2HE UV = vec2HE_noinit();
            get_edge_edge_intersection(e1,e2,I,UV);
            vertex_.push_back(Vertex(this,I,UV));
            index_t x = vertex_.size()-1;
            CDTBase2d::v2T_.push_back(index_t(-1));
            geo_debug_assert(x == CDTBase2d::nv_);
            ++CDTBase2d::nv_;
            return x;
        }

    public:
        void save(const std::string& filename) const override {
            Mesh M;
            M.vertices.set_dimension(2);
            for(index_t v=0; v<CDTBase2d::nv(); ++v) {
                vec2 p;
                double w = vertex_[v].UV_exact.w.estimate();
                p = vec2(
                    vertex_[v].UV_exact.x.estimate() / w,
                    vertex_[v].UV_exact.y.estimate() / w
                );
                M.vertices.create_vertex(p.data());
            }
            for(index_t t=0; t<CDTBase2d::nT(); ++t) {
                M.facets.create_triangle(
                    CDTBase2d::Tv(t,0),
                    CDTBase2d::Tv(t,1),
                    CDTBase2d::Tv(t,2)
                );
            }

            Attribute<double> tex_coord;
            tex_coord.create_vector_attribute(
                M.facet_corners.attributes(), "tex_coord", 2
            );
            static double triangle_tex[3][2] = {
                {0.0, 0.0},
                {1.0, 0.0},
                {0.0, 1.0}
            };
            for(index_t c: M.facet_corners) {
                tex_coord[2*c]   = triangle_tex[c%3][0];
                tex_coord[2*c+1] = triangle_tex[c%3][1];
            }
            
            mesh_save(M, filename);
        }
        
    private:
        Mesh& mesh_;
        index_t f1_;
        index_t latest_f2_;
        index_t latest_f2_count_;
        coord_index_t f1_normal_axis_;
        coord_index_t u_; // = (f1_normal_axis_ + 1)%3
        coord_index_t v_; // = (f1_normal_axis_ + 2)%3
        vector<Vertex> vertex_;
        vector<Edge> edges_;
        std::map<vec3HE, index_t, vec3HELexicoCompare> g_v_table_;
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
        {
            Stopwatch W("Detect isect");
            MeshFacetsAABB AABB(M,!params.debug_do_not_order_facets);
            vector<std::pair<index_t, index_t> > FF;

            // Get candidate pairs of intersecting facets
            AABB.compute_facet_bbox_intersections(
                [&](index_t f1, index_t f2) {
                    // Needed (maybe I should change that in AABB class)
                    if(f1 == f2) {
                        return;
                    }
                    geo_assert(f1 < f2);

                    // Optionally skip facet pairs that
                    // share a vertex or an edge
                    if(
                        !params.detect_intersecting_neighbors && (
                            (M.facets.find_adjacent(f1,f2)  != index_t(-1)) ||
                            (M.facets.find_common_vertex(f1,f2) != index_t(-1))
                        )
                    ) {
                        return;
                    }

                    FF.push_back(std::make_pair(f1,f2));
                }
            );

            // Compute facet-facet intersections in parallel
            Process::spinlock lock = GEOGRAM_SPINLOCK_INIT;
            parallel_for_slice(
                0,FF.size(), [&](index_t b, index_t e) {
                    vector<TriangleIsect> I;
                    for(index_t i=b; i<e; ++i){
                        index_t f1 = FF[i].first;
                        index_t f2 = FF[i].second;
                        
                        if(mesh_facets_intersect(M, f1, f2, I)) {

                            Process::acquire_spinlock(lock);                            
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
                                            regions_convex_hull(II.A_rgn_f1,II.B_rgn_f1);
                                        
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
                            Process::release_spinlock(lock);
                        }
                    }
                }
            );
        }

        // Step 3: Remesh intersected triangles
        // ------------------------------------
        {
            Stopwatch W("Remesh isect");
            
            // The MeshInTriangle, that implements local triangulation
            // in each intersected triangular facet.
            // It also keeps a global map of 
            // vertices, indexed by their exact geometry.         
            MeshInTriangle TM(M);
            // TM.set_check_constraints(params.debug_check_constraints);
            // TM.set_barycentric(params.barycentric);
            // TM.set_use_halfedges(params.use_halfedges);
            TM.set_delaunay(params.delaunay);
            
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

            index_t nf = M.facets.nb();
            
            // Now iterate on all intersections, and identify
            // the [b,e[ intervals that correspond to the same f1 facet.
            index_t b=0;
            while(b < intersections.size()) {
                index_t e = b;
                while(
                    e < intersections.size() &&
                    intersections[e].f1 == intersections[b].f1
                ) {
                    ++e;
                }

                if(params.verbose) {
                    std::cerr << "Isects in " << intersections[b].f1
                              << " / " << nf                    
                              << "    : " << (e-b)
                              << std::endl;
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

    /**
     * \brief Enumerates the connected components in a facet attribute
     * \param[in] M a reference to the mesh
     * \param[in] attribute the name of the facet attribute
     * \return the number of found connected components
     */
    index_t get_surface_connected_components(Mesh& M, const std::string& attribute = "chart") {
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
            get_surface_connected_components(M,"operand_bit");
            operand_bit.bind(M.facets.attributes(), "operand_bit");
            for(index_t f: M.facets) {
                operand_bit[f] =
                    params.per_component_ids ? (1 << operand_bit[f]) : 1;
            }
        }
        
        if(params.pre_detect_duplicated_vertices) {
            remove_linear_triangles(M);
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

    /**
     * \brief Gets the number of bits set in
     *  a 32 bits integer
     * \param[in] x the integer
     * \return the number of "ones" in the 
     *  binary form of x
     */
    inline index_t nb_bits_set(index_t x) {
        index_t result = 0;
        for(index_t i=0; i<32; ++i) {
            result += (x&1);
            x = x >> 1;
        }
        return result;
    }
}

namespace {
    using namespace GEO;

    // TODO: something less hacky...
    bool facet_on_border(MeshFacetsAABB& AABB, index_t f) {
        vec3 g = Geom::mesh_facet_center(*(AABB.mesh()),f);
        for(index_t i=0; i<1000; ++i) {
            vec3 D(
                Numeric::random_float64(),
                Numeric::random_float64(),
                Numeric::random_float64()
            );
            
            if(!AABB.ray_intersection(
                   Ray(g,D), Numeric::max_float64(), f
            )) {
                return true;
            }
            if(!AABB.ray_intersection(
                   Ray(g,-D), Numeric::max_float64(), f
            )) {
                return true;
            }
        }
        return false;
    }
    
    void mesh_classify_union(
        Mesh& M, 
        const std::string& attribute,
        bool reorder
    ) {
        MeshFacetsAABB AABB(M,reorder);
        index_t nb_charts = get_surface_connected_components(M);        
        Attribute<index_t> chart(M.facets.attributes(), "chart");
        Attribute<bool> selection;
        vector<index_t> delete_f;
        if(attribute != "") {
            selection.bind(M.facets.attributes(), attribute);            
        } else {
            delete_f.assign(M.facets.nb(), 0);
        }

        vector<index_t> chart_facet(nb_charts, index_t(-1));
        for(index_t f: M.facets) {
            index_t c = chart[f];
            if(chart_facet[c] == index_t(-1)) {
                bool f_is_selected = facet_on_border(AABB,f);
                chart_facet[c] = f;
                if(selection.is_bound()) {
                    selection[f] = f_is_selected;
                } else {
                    delete_f[f] = !f_is_selected;
                }
            } else {
                if(selection.is_bound()) {
                    selection[f] = selection[chart_facet[c]];
                } else {
                    delete_f[f] = delete_f[chart_facet[c]];
                }
            }
        } 
        if(!selection.is_bound()) {
            M.facets.delete_elements(delete_f);
        }
        mesh_repair(
            M,
            GEO::MeshRepairMode(
                GEO::MESH_REPAIR_COLOCATE | GEO::MESH_REPAIR_DUP_F
            ),
            0.0
        );
    }
    
}

namespace GEO {


    
    void mesh_classify_intersections(
        Mesh& M, std::function<bool(index_t)> eqn,
        const std::string& attribute,
        bool reorder
    ) {
        MeshFacetsAABB AABB(M,reorder);
        index_t nb_charts = get_surface_connected_components(M);
        Attribute<index_t> chart(M.facets.attributes(), "chart");
        Attribute<index_t> operand_bit(M.facets.attributes(), "operand_bit");
        Attribute<bool> selection;
        vector<index_t> delete_f;
        if(attribute != "") {
            selection.bind(M.facets.attributes(), attribute);            
        } else {
            delete_f.assign(M.facets.nb(), 0);
        }

        vector<index_t> chart_facet(nb_charts, index_t(-1));
        try {
            for(index_t f: M.facets) {
                index_t c = chart[f];
                if(chart_facet[c] == index_t(-1)) {
                    bool f_is_selected = false;
                    chart_facet[c] = f;
                    vec3 g = Geom::mesh_facet_center(M,f);
                    // Picking f's normal is not a good idea,
                    // because for industrial parts it will
                    // encounter many degenerate ray/triangle
                    // intersections.
                    // TODO: we need to detect them and launch
                    // another ray whenever the ray intersects
                    // the surface on a vertex or on an edge.
                    vec3 D(
                        Numeric::random_float64(),
                        Numeric::random_float64(),
                        Numeric::random_float64()
                    );
                    index_t parity = 0;
                    AABB.ray_all_intersections(
                        Ray(g,D),
                        [&](const MeshFacetsAABB::Intersection & I) {
                            if(I.f != f) {
                                parity = parity ^ operand_bit[I.f];
                            }
                        }
                    );
                    if(nb_bits_set(operand_bit[f]) == 1) {
                        // Facet f is on the boundary of the result if
                        // crossing f changes the result of eqn,
                        // in other words, if eqn gives a different
                        // result with and without f's object bit set
                        f_is_selected =
                            eqn(parity |  operand_bit[f]) !=
                            eqn(parity & ~operand_bit[f]) ;
                    } else {
                        // Now if f is on several objects (that is, has
                        // several bit sets), then we determine whether
                        // it is on the boundary of the result by raytracing
                        // in two different directions, and seeing if eqn
                        // gives a different result. 
                        index_t parity2 = 0;
                        AABB.ray_all_intersections(
                            Ray(g,-D),
                            [&](const MeshFacetsAABB::Intersection & I) {
                                if(I.f != f) {
                                    parity2 = parity2 ^ operand_bit[I.f];
                                }
                            }
                        );
                        f_is_selected = (eqn(parity) != eqn(parity2));
                    }
                    if(selection.is_bound()) {
                        selection[f] = f_is_selected;
                    } else {
                        delete_f[f] = !f_is_selected;
                    }
                } else {
                    if(selection.is_bound()) {
                        selection[f] = selection[chart_facet[c]];
                    } else {
                        delete_f[f] = delete_f[chart_facet[c]];
                    }
                }
            }
        } catch(const std::logic_error& e) {
            Logger::err("Classify") << "Error while parsing expression:"
                                    << e.what()
                                    << std::endl;
            return;
        }
        if(!selection.is_bound()) {
            M.facets.delete_elements(delete_f);
        }
        mesh_repair(
            M,
            GEO::MeshRepairMode(
                GEO::MESH_REPAIR_COLOCATE | GEO::MESH_REPAIR_DUP_F
            ),
            0.0
        );
    }

    void mesh_classify_intersections(
        Mesh& M, const std::string& expr,
        const std::string& attribute, bool reorder
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
        
        if(expr == "union") {
            mesh_classify_union(M, attribute, reorder);
            return;
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
                attribute,
                reorder
            );
        } catch(const std::logic_error& e) {
            Logger::err("Classify") << "Error while parsing expression:"
                                    << e.what()
                                    << std::endl;
            return;
        }
    }    
}
