/*
 *
 */
 
#ifndef H__OGF_EXPERIMENT_CDT__H
#define H__OGF_EXPERIMENT_CDT__H

#include <OGF/Experiment/common/common.h>
#include <geogram/basic/geometry.h>
#include <geogram/numerics/predicates.h>

#define CDT_DEBUG
#ifdef CDT_DEBUG
#define CDT_LOG(X) std::cerr << X << std::endl
#else
#define CDT_LOG(X)
#endif

namespace GEO {

    // Under development, unfinished !!
    class Experiment_API CDT {
    public:
        typedef std::pair<index_t, index_t> Edge;
    
        CDT() : delaunay_(true) {
        }
    
        index_t insert(const vec2& P);
        
        void insert_constraint(index_t i, index_t j);
        void insert_constraint_simple(index_t i, index_t j);
        void insert_constraint_optimized(index_t i, index_t j);

        /**
         * \brief Specifies whether a constrained Delaunay
         *  triangulation should be constructed, or just a
         *  plain constrained triangulation
         * \param[in] delaunay true if a Delaunay triangulation
         *  should be constructed, false otherwise.
         */
        void set_delaunay(bool delaunay) {
            delaunay_ = delaunay;
        }
        
        /**
         * \brief Gets the number of triangles
         */
        index_t nT() const {
            return T_.size()/3;
        }

        /**
         * \brief Gets the number of vertices
         */
        index_t nv() const {
            return point_.size();
        }
    
        /**
         * \brief Gets a vertex of a triangle
         * \param[in] t the triangle
         * \param[in] lv the local index of the vertex, in 0,1,2
         * \return the global index of the vertex
         */
        index_t Tv(index_t t, index_t lv) const {
            geo_debug_assert(t<nT());
            geo_debug_assert(lv<3);
            return T_[3*t+lv];
        }

        /**
         * \brief Finds the local index of a vertex in a triangle
         * \param[in] t the triangle
         * \param[in] v the vertex
         * \return lv such that Tv(t,lv) = v
         */
        index_t Tv_find(index_t t, index_t v) const {
            geo_debug_assert(t<nT());
            geo_debug_assert(v<nv());
            for(index_t lv=0; lv<3; ++lv) {
                if(T_[3*t+lv] == v) {
                    return lv;
                }
            }
            geo_assert_not_reached;
            return index_t(-1);
        }
        
        /**
         * \brief Gets a triangle adjacent to a triangle
         * \param[in] t the triangle
         * \param[in] le the local edge index, in 0,1,2
         * \return the triangle adjacent to \p t accross \p le,
         *  or index_t(-1) if there is no such triangle
         */
        index_t Tadj(index_t t, index_t le) const {
            geo_debug_assert(t<nT());
            geo_debug_assert(le<3);
            return Tadj_[3*t+le];
        }
        
        /**
         * \brief Finds the edge accross which a triangle is
         *  adjacent to another one
         * \param[in] t1 , t2 the two triangles
         * \return the local edge index le such that
         *  Tadj(t1,le) = t2
         */
        index_t Tadj_find(index_t t1, index_t t2) const {
            geo_debug_assert(t1<nT());
            geo_debug_assert(t2<nT());            
            for(index_t le=0; le<3; ++le) {
                if(Tadj_[3*t1+le] == t2) {
                    return le;
                }
            }
            geo_assert_not_reached;
            return index_t(-1);
        }

        /**
         * \brief Gets a triangle incident to a given vertex
         * \param[in] v a vertex
         * \return a triangle t such that there exists lv in 
         *  0,1,2 such that Tv(t,lv) = v
         */
        index_t vT(index_t v) const {
            geo_debug_assert(v < nv());
            return v2T_[v];
        }

        /**
         * \brief Gets a triangle incident to an edge
         * \param[in] v1 , v2 the two vertices of the edge
         * \return a triangle such that edge 0 corresponds to
         *  v1 , v2.
         * \details The function rotates the triangle if need
         *  be (modifies the triangulation).
         */
        index_t eT(index_t v1, index_t v2);


        /**
         * \brief Flips an edge
         * \details Unefficient, because it needs internally
         *  to retreive a triangle incident to the edge, by
         *  turning around one of the vertices.
         * \param[in,out] E the edge to be flipped. On exit the
         *  flipped edge.
         */
        void swap_edge(Edge& E);

        /**
         * \brief Saves this CDT to a geogram mesh file.
         * \param[in] filename where to save this CDT
         */
        void save(const std::string& filename) const;
    
    protected:
    
        /**
         * \brief Swaps an edge.
         * \details Swaps edge 0 of \p t1.
         *    Vertex 0 of \p t1 is vertex 0 of
         *    the two new triangles.
         * \param[in] t1 a triangle index. Its edge
         *    opposite to vertex 0 is swapped
         */
        void swap_edge(index_t t1);
    
        /**
         * \brief Sets all the combinatorial information
         *  of a triangle
         * \param[in] t the triangle
         * \param[in] v1 , v2 , v3 the three vertices 
         * \param[in] adj1 , adj2 , adj3 the three triangles
         *  adjacent to \p t
         */
        void Tset(
            index_t t,
            index_t v1,   index_t v2,   index_t v3,
            index_t adj1, index_t adj2, index_t adj3
        ) {
            geo_debug_assert(t < nT());
            geo_debug_assert(v1 < nv());
            geo_debug_assert(v2 < nv());
            geo_debug_assert(v3 < nv());                        
            geo_debug_assert(adj1 < nT() || adj1 == index_t(-1));
            geo_debug_assert(adj2 < nT() || adj2 == index_t(-1));
            geo_debug_assert(adj3 < nT() || adj3 == index_t(-1));
            geo_debug_assert(v1 != v2);
            geo_debug_assert(v2 != v3);
            geo_debug_assert(v3 != v1);            
            geo_debug_assert(adj1 != adj2 || adj1 == index_t(-1));
            geo_debug_assert(adj2 != adj3 || adj2 == index_t(-1));
            geo_debug_assert(adj3 != adj1 || adj3 == index_t(-1));            
            T_[3*t  ]    = v1;
            T_[3*t+1]    = v2;
            T_[3*t+2]    = v3;                        
            Tadj_[3*t  ] = adj1;
            Tadj_[3*t+1] = adj2;
            Tadj_[3*t+2] = adj3;
            v2T_[v1] = t;
            v2T_[v2] = t;
            v2T_[v3] = t;        
        }

        /**
         * \brief Rotates indices in triangle t in such a way
         *  that a given vertex becomes vertex 0
         * \details On exit, vertex \p lv of \p t becomes vertex 0
         * \param[in] t a triangle index
         * \param[in] lv local vertex index in 0,1,2
         */
        void Trot(index_t t, index_t lv) {
            geo_debug_assert(t < nT());
            index_t i = 3*t+(lv%3);
            index_t j = 3*t+((lv+1)%3);
            index_t k = 3*t+((lv+2)%3);        
            Tset(
                t,
                T_[i],    T_[j],    T_[k],
                Tadj_[i], Tadj_[j], Tadj_[k]
            );
        }

        /**
         * \brief Consistency check for a triangle
         * \details in debug mode, aborts if inconsistency is detected
         * \param[in] t the triangle to be tested
         */
        void Tcheck(index_t t) {
            geo_argused(t);
#ifdef GEO_DEBUG
            if(t == index_t(-1)) {
                return;
            }
            for(index_t e=0; e<3; ++e) {
                geo_debug_assert(Tv(t,e) != Tv(t,(e+1)%3));
                if(Tadj(t,e) == index_t(-1)) {
                    continue;
                }
                geo_debug_assert(Tadj(t,e) != Tadj(t,(e+1)%3));
                index_t t2 = Tadj(t,e);
                index_t e2 = Tadj_find(t2,t);
                geo_debug_assert(Tadj(t2,e2) == t);
            }
#endif
        }

        void check_consistency() {
#ifdef GEO_DEBUG            
            for(index_t t=0; t<nT(); ++t) {
                Tcheck(t);
            }
#endif            
        }

        
        /**
         * \brief Sets a triangle adjacency relation
         * \param[in] t a triangle
         * \param[in] le local edge index, in 0,1,2
         * \param[in] adj the triangle adjacent to \p t 
         *  accross \p le
         */
        void Tadj_set(index_t t, index_t le, index_t adj) {
            geo_debug_assert(t < nT());
            geo_debug_assert(adj < nT());
            geo_debug_assert(le < 3);
            Tadj_[3*t+le] = adj;
        }

        /**
         * \brief Replaces a triangle adjacency relation
         * \param[in] t a triangle
         * \param[in] t1 a triangle adjacent to \p t
         * \param[in] t2 the triangle that will replace \p t1
         */
        void Tadj_replace(index_t t, index_t t1, index_t t2) {
            if(t == index_t(-1)) {
                return;
            }
            geo_debug_assert(t < nT());
            geo_debug_assert(t1 < nT());
            geo_debug_assert(t2 < nT());
            index_t le = Tadj_find(t,t1);
            Tadj_set(t, le, t2);
            Tcheck(t);            
        }

        /**
         * \brief Creates a new triangle
         * \return the index of the new triange
         */
        index_t Tnew() {
            index_t t = nT();
            T_.resize((t+1)*3, index_t(-1));
            Tadj_.resize((t+1)*3, index_t(-1));
            Tflags_.resize(t+1,0);
            return t;
        }

        void Tmark(index_t t) {
            Tflags_[t] = 1;
        }

        void Tunmark(index_t t) {
            Tflags_[t] = 0;
        }

        bool Tis_marked(index_t t) {
            return Tflags_[t] == 1;
        }
    
        /**
         * \brief Calls a user-defined function for each triangle 
         * around a vertex
         * \param[in] v the vertex
         * \param[in] doit the function, that takes as argument the 
         *  current triangle t and the local index lv of \p v in t. 
         *  The function returns true if  iteration is finished and can be 
         *  exited, false otherwise.
         */
        void for_each_T_around_v(
            index_t v, std::function<bool(index_t t, index_t lv)> doit
        ) {
            index_t t = vT(v);
            index_t lv = index_t(-1);
            do {
                lv = Tv_find(t,v);
                if(doit(t,lv)) {
                    return;
                }
                t = Tadj(t, (lv+1)%3);
            } while(t != vT(v) && t != index_t(-1));
            
            // We are done, this was an interior vertex
            if(t != index_t(-1)) {
                return;
            }
            
            // It was a vertex on the border, so we need
            // to traverse the triangle fan in the other
            // direction until we reach the border again
            t = vT(v);
            lv = Tv_find(t,v);
            t = Tadj(t, (lv+2)%3);
            while(t != index_t(-1)) {
                lv = Tv_find(t,v);
                if(doit(t,lv)) {
                    return;
                }
                t = Tadj(t, (lv+2)%3);
            }
        }
        
        /**
         * \brief Locates a vertex
         * \param[in] v the vertex index
         * \param[out] o1 , o2 , o3 the three orientation predicates
         * \return a triangle that contains \p v
         */
        index_t locate(index_t v, Sign& o1, Sign& o2, Sign& o3) const;
        
        /**
         * \brief Orientation predicate
         * \param[in] i , j , k three vertices
         * \return the sign of det(pj-pi,pk-pi)
         */
        Sign orient2d(index_t i, index_t j, index_t k) const;

        /**
         * \brief Incircle predicate
         * \param[in] i , j , k the three vertices of a triangle
         * \param[in] l another vertex
         * \retval POSITIVE if \p l is inside the circumscribed circle of
         *  the triangle
         * \retval ZERO if \p l is on the circumscribed circle of
         *  the triangle
         * \retval NEGATIVE if \p l is outside the circumscribed circle of
         *  the triangle
         */
        Sign incircle(index_t i, index_t j, index_t k, index_t l) const;

        /**
         * \brief Tests whether two segments intersect
         * \param[in] i , j the extremities of the first segment
         * \param[in] k , l the extremities of the second segment
         * \retval true if the two segments have an intersection (if they just
         *  touch it does not count)
         * \retval false otherwise
         */
        bool seg_seg_intersect(index_t i,index_t j,index_t k,index_t l) const;

        /**
         * \brief Tests whether triange t and its neighbor accross edge 0 form 
         *  a strictly convex quad
         * \retval true if triange \p t and its neighbor accross edge 0 form
         *  a strictly convex quad
         * \retval false otherwise
         */
        bool is_convex_quad(index_t t) const;

        template <class LIST> void debug_show_triangles(const LIST& L, const std::string& name) {
            geo_argused(L);
            geo_argused(name);
#ifdef CDT_DEBUG            
            std::cerr << name << "= ";
            for(auto t: L) {
                std::cerr << t
                          << (Tis_marked(t) ? "" : "*")
                          << "->(" << Tv(t,1) << "," << Tv(t,2) << ")"
                          << " ";
            }
            std::cerr << std::endl;
        }
#endif
        
    protected:
        vector<vec2> point_;
        vector<index_t> T_;
        vector<index_t> Tadj_;
        vector<index_t> v2T_;
        vector<bool>    Tflags_;
        bool delaunay_; 
    };

}

#endif