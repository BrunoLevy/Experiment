/*
 *
 */
 
#ifndef H__OGF_EXPERIMENT_CDT__H
#define H__OGF_EXPERIMENT_CDT__H

#include <OGF/Experiment/common/common.h>
#include <geogram/basic/geometry.h>
#include <geogram/numerics/predicates.h>

namespace GEO {

    /**
     * \brief Base class for constrained Delaunay triangulation
     * \details Manages the combinatorics of the constrained Delaunay
     *  triangulation. The points need to be stored elsewhere, and manipulated 
     *  through indices, with two predicates:
     *  - orient2d(i,j,k)
     *  - incircle(i,j,k,l)
     *  and one construction:
     *  - create_intersection(i,j,k,l)
     *  See \ref CDT for an example of implementation
     */
    class Experiment_API CDTBase {
    public:
        enum { NO_INDEX = index_t(-1) };

        /**
         * \brief CDTBase constructor
         */
        CDTBase();

        /**
         * \brief CDTBase destructor
         */
        virtual ~CDTBase();

        /**
         * \brief Removes everything from this triangulation
         */
        virtual void clear();
        
        /**
         * \brief Inserts a constraint
         * \param[in] i , j the indices of the two vertices
         *  of the constrained segment
         */
        void insert_constraint(index_t i, index_t j);

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
            return nv_;
        }

        /**
         * \brief Gets the number of constraints
         */
        index_t ncnstr() const {
            return ncnstr_;
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
         * \brief Saves this CDT to a geogram mesh file.
         * \param[in] filename where to save this CDT
         */
        virtual void save(const std::string& filename) const = 0;


        /**
         * \brief Tests whether a triangle edge is Delaunay
         * \details returns true also for constrained edges and edges on borders
         */
        bool Tedge_is_Delaunay(index_t t, index_t le) const;

    protected:
        /**
         * \brief Inserts a new point 
         * \param[in] v the index of the new point, supposed to be
         *  equal to nv()
         * \return the index of the created point. May be different
         *  from v if the point already existed in the triangulation
         */
        index_t insert(index_t v);


        /**
         * \brief Creates the combinatorics for a first large enclosing
         *  triangle
         * \param[in] v1 , v2 , v3 the three vertices of the first triangle,
         *  in 0,1,2
         * \details create_enclosing_triangle() or create_enclosing_quad() 
         *  need to be called before anything else
         */
        void create_enclosing_triangle(index_t v1, index_t v2, index_t v3);

        /**
         * \brief Creates the combinatorics for a first large enclosing
         *  quad
         * \param[in] v1 , v2 , v3 , v4 the four vertices of the quad,
         *  in 0,1,2,3
         * \details create_enclosing_triangle() or create_enclosing_quad() 
         *  need to be called before anything else
         */
        void create_enclosing_quad(
            index_t v1, index_t v2, index_t v3, index_t v4
        );

        /**
         * \brief Constants for triangle flags
         */
        enum {T_MARKED_MASK = 1, T_IN_LIST_MASK = 2};

        index_t Tnext(index_t t) {
            geo_debug_assert(t < nT());
            geo_debug_assert(Tis_in_list(t));
            return Tnext_[t];
        }

        index_t Tprev(index_t t) {
            geo_debug_assert(t < nT());
            geo_debug_assert(Tis_in_list(t));            
            return Tprev_[t];
        }
        
        /**
         * \brief Doubly connected triangle list
         */
        struct DList {
            DList() : back(NO_INDEX), front(NO_INDEX) {
            }
            bool empty() {
                geo_debug_assert(
                    (back==NO_INDEX)==(front==NO_INDEX)
                );
                return (back==NO_INDEX);
            }
            index_t back;
            index_t front;
        };

        void DList_clear(DList& L) {
            for(index_t t=L.front; t!=NO_INDEX; t = Tnext_[t]) {
                Tunflag_as_in_list(t);
            }
            L.back = NO_INDEX;
            L.front = NO_INDEX;
        }
        
        void DList_push_back(DList& L, index_t t) {
            geo_debug_assert(!Tis_in_list(t));
            Tflag_as_in_list(t);
            if(L.empty()) {
                L.back = t;
                L.front = t;
                Tnext_[t] = NO_INDEX;
                Tprev_[t] = NO_INDEX;
            } else {
                Tnext_[t] = NO_INDEX;
                Tnext_[L.back] = t;
                Tprev_[t] = L.back;
                L.back = t;
            }
        }

        index_t DList_pop_back(DList& L) {
            geo_debug_assert(!L.empty());
            index_t t = L.back;
            L.back = Tprev_[L.back];
            if(L.back == NO_INDEX) {
                geo_debug_assert(L.front == t);
                L.front = NO_INDEX;
            } else {
                Tnext_[L.back] = NO_INDEX;
            }
            geo_debug_assert(Tis_in_list(t));
            Tunflag_as_in_list(t);
            return t;
        }

        void DList_push_front(DList& L, index_t t) {
            geo_debug_assert(!Tis_in_list(t));
            Tflag_as_in_list(t);
            if(L.empty()) {
                L.back = t;
                L.front = t;
                Tnext_[t] = NO_INDEX;
                Tprev_[t] = NO_INDEX;
            } else {
                Tprev_[t] = NO_INDEX;
                Tprev_[L.front] = t;
                Tnext_[t] = L.front;
                L.front = t;
            }
        }

        index_t DList_pop_front(DList& L) {
            geo_debug_assert(!L.empty());
            index_t t = L.front;
            L.front = Tnext_[L.front];
            if(L.front == NO_INDEX) {
                geo_debug_assert(L.back == t);
                L.back = NO_INDEX;
            } else {
                Tprev_[L.front] = NO_INDEX;
            }
            geo_debug_assert(Tis_in_list(t));
            Tunflag_as_in_list(t);
            return t;
        }


        void DList_remove(DList& L, index_t t) {
            if(t == L.front) {
                DList_pop_front(L);
            } else if(t == L.back) {
                DList_pop_back(L);
            } else {
                geo_debug_assert(Tis_in_list(t));
                index_t t_prev = Tprev(t);
                index_t t_next = Tnext(t);
                Tprev_[t_next] = t_prev;
                Tnext_[t_prev] = t_next;
                Tunflag_as_in_list(t);
            }
        }
        
        /**
         * \brief Inserts a vertex in an edge
         * \param[in] v the vertex to be inserted
         * \param[in] t a triangle incident to the edge
         * \param[in] le the local index of the edge in \p t
         * \param[out] S optional DList of created triangles
         */
        void insert_vertex_in_edge(
            index_t v, index_t t, index_t le,
            DList* S = nullptr
        );

        /**
         * \brief Inserts a vertex in a triangle
         * \param[in] v the vertex to be inserted
         * \param[in] t the triangle
         * \param[out] S optional DList of created triangles
         */
        void insert_vertex_in_triangle(
            index_t v, index_t t,
            DList* S = nullptr
        );
        
        /**
         * \brief Finds the edges intersected by a constraint
         * \param[in] i , j the two vertices of the constraint
         * \param[out] Q for each intersected edge, a triangle t
         *  will be pushed-back to Q, such that vT(t,1) and
         *  vT(t,2) are the extremities of the intersected edge.
         *  In addition, each triangle t is marked.
         * \details If a vertex k that is exactly on the constraint
         *  is found, then traversal stops there and k is returned.
         *  One can find the remaining intersections by continuing
         *  to call the function until \p j is returned.
         * \return the first vertex on [i,j] encountered when
         *  traversing the segment [i,j]. 
         */
        index_t find_intersected_edges(index_t i, index_t j, DList& Q);

        /**
         * \brief Constrains an edge by iteratively flipping
         *  the intersected edges.
         * \param[in] i , j the extremities of the edge
         * \param[in] Q the list of intersected edges, computed
         * \param[out] N the new edges, that need to be re-Delaunized
         *  by find_intersected_edges()
         */
        void constrain_edges(index_t i, index_t j, DList& Q, DList& N);

        /**
         * \brief Restore Delaunay condition starting from the
         *  triangles incident to a given vertex.
         * \param[in] v the vertex
         * \param[in] S a stack of triangles, initialized with 
         *  the triangles incident to the vertex. Each triangle t
         *  are Trot()-ed in such a way that the vertex 
         *  corresponds to Vt(t,0)
         * \details Each time a triangle edge is swapped, the
         *  two new neighbors are recursively examined.
         */
        void Delaunayize_vertex_neighbors(index_t v, DList& S);
        
        /**
         * \brief Restore Delaunay condition for a set of
         *  edges after inserting a constrained edge.
         * \param[in] i , j the extremities of the constraint
         * \param[in] N the edges for which Delaunay condition
         *  should be restored.
         */
        void Delaunayize_new_edges(DList& N);

        
        /**
         * \brief Sets all the combinatorial information
         *  of a triangle and clears edge flags
         * \param[in] t the triangle
         * \param[in] v1 , v2 , v3 the three vertices 
         * \param[in] adj1 , adj2 , adj3 the three triangles
         *  adjacent to \p t
         * \param[in] e1cnstr , e2cnstr , e3cnstr optional
         *  edge constraints
         */
        void Tset(
            index_t t,
            index_t v1,   index_t v2,   index_t v3,
            index_t adj1, index_t adj2, index_t adj3,
            index_t e1cnstr = NO_INDEX,
            index_t e2cnstr = NO_INDEX,
            index_t e3cnstr = NO_INDEX            
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
            Tecnstr_[3*t]   = e1cnstr;
            Tecnstr_[3*t+1] = e2cnstr;
            Tecnstr_[3*t+2] = e3cnstr;
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
            geo_debug_assert(lv < 3);
            index_t i = 3*t+lv;
            index_t j = 3*t+((lv+1)%3);
            index_t k = 3*t+((lv+2)%3);
            Tset(
                t,
                T_[i], T_[j], T_[k],
                Tadj_[i], Tadj_[j], Tadj_[k],
                Tecnstr_[i], Tecnstr_[j], Tecnstr_[k]
            );
        }

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

        /**
         * \brief Consistency check for all the triangles
         * \details in debug mode, aborts if inconsistency is detected
         */        
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
         * \brief After having changed connections from triangle
         *  to a neighbor, creates connections from neighbor
         *  to triangle.
         * \details edge flags are copied from the neighbor to \p t1.
         *  If there is no triangle accross \p le1, then
         *  nothing is done
         * \param[in] t1 a triangle
         * \param[in] le1 a local edge of \p t1, in 0,1,2
         * \param[in] prev_t2_adj_e2 the triangle adjacent to t2 that
         *  \p t1 will replace, where t2 = Tadj(t1,le1)
         */
        void Tadj_back_connect(
            index_t t1, index_t le1, index_t prev_t2_adj_e2
        ) {
            geo_debug_assert(t1 < nT());
            geo_debug_assert(le1 < 3);
            index_t t2 = Tadj(t1,le1);
            if(t2 == NO_INDEX) {
                return;
            }
            index_t le2 = Tadj_find(t2,prev_t2_adj_e2);
            Tadj_set(t2,le2,t1);
            Tset_edge_cnstr(t1,le1,Tedge_cnstr(t2,le2));
        }
        
        /**
         * \brief Creates a new triangle
         * \return the index of the new triange
         */
        index_t Tnew() {
            index_t t = nT();
            index_t nc = (t+1)*3; // new number of corners
            T_.resize(nc, NO_INDEX);
            Tadj_.resize(nc, NO_INDEX);
            Tecnstr_.resize(nc, NO_INDEX);
            Tflags_.resize(t+1,0);
            Tnext_.resize(t+1,NO_INDEX);
            Tprev_.resize(t+1,NO_INDEX);
            return t;
        }

        
        /**
         * \brief Marks a triangle as intersected by the constraint
         */
        void Tmark(index_t t) {
            geo_debug_assert(t < nT());            
            Tflags_[t] |= T_MARKED_MASK;
        }

        /**
         * \brief Unmarks a triangle as intersected by the constraint
         */
        void Tunmark(index_t t) {
            geo_debug_assert(t < nT());            
            Tflags_[t] &= Numeric::uint8(~T_MARKED_MASK);
        }

        /**
         * \brief Tests whether a triangle is marked as intersected
         */
        bool Tis_marked(index_t t) {
            geo_debug_assert(t < nT());            
            return ((Tflags_[t] & T_MARKED_MASK) != 0);
        }

        /**
         * \brief Marks a triangle as in a DList
         */
        void Tflag_as_in_list(index_t t) {
            geo_debug_assert(t < nT());            
            Tflags_[t] |= T_IN_LIST_MASK;
        }

        /**
         * \brief Unmarks a triangle as in a DList
         */
        void Tunflag_as_in_list(index_t t) {
            geo_debug_assert(t < nT());            
            Tflags_[t] &= Numeric::uint8(~T_IN_LIST_MASK);
        }
        
        /**
         * \brief Tests whether a triangle is in a DList
         */
        bool Tis_in_list(index_t t) {
            geo_debug_assert(t < nT());            
            return ((Tflags_[t] & T_IN_LIST_MASK) != 0);
        }

        
        /**
         * \brief Tests whether an edge is constrained
         * \param[in] t a triangle
         * \param[in] le local edge index, in 0,1,2
         * \retval true if the edge is constrained
         * \retval false otherwise
         */
        index_t Tedge_cnstr(index_t t, index_t le) const {
            geo_debug_assert(t < nT());
            geo_debug_assert(le < 3);
            return Tecnstr_[3*t+le];
        }
        
        /**
         * \brief Sets an edge as constrained
         * \param[in] t a triangle
         * \param[in] le local edge index, in 0,1,2
         * \param[in] cnstr_id identifier of the constrained edge
         */
        void Tset_edge_cnstr(
            index_t t, index_t le, index_t cnstr_id 
        ) {
            geo_debug_assert(t < nT());
            geo_debug_assert(le < 3);
            Tecnstr_[3*t+le] = cnstr_id;
        }

        /**
         * \brief Sets an edge as constrained in triangle and in neighbor
         * \param[in] t a triangle
         * \param[in] le local edge index, in 0,1,2
         * \param[in] cnstr_id identifier of the constrained edge
         */
        void Tset_edge_cnstr_with_neighbor(
            index_t t, index_t le, index_t cnstr_id 
        ) {
            geo_debug_assert(t < nT());
            geo_debug_assert(le < 3);
            Tset_edge_cnstr(t, le, cnstr_id);
            index_t t2 = Tadj(t,le);
            if(t2 != NO_INDEX) {
                index_t le2 = Tadj_find(t2,t);
                Tset_edge_cnstr(t2,le2,cnstr_id);
            }
        }
        
        /**
         * \brief Tests whether an edge is constrained
         * \param[in] t a triangle
         * \param[in] le local edge index, in 0,1,2
         * \retval true if the edge is constrained
         * \retval false otherwise
         */
        bool Tedge_is_constrained(index_t t, index_t le) const {
            return (Tedge_cnstr(t,le) != NO_INDEX);
        }

        /**
         * \brief Calls a user-defined function for each triangle 
         * around a vertex
         * \param[in] v the vertex
         * \param[in] doit the function, that takes as argument the 
         *  current triangle t and the local index lv of \p v in t. 
         *  The function returns true if iteration is finished and can be 
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
         * \brief Tests whether triange t and its neighbor accross edge 0 form 
         *  a strictly convex quad
         * \retval true if triange \p t and its neighbor accross edge 0 form
         *  a strictly convex quad
         * \retval false otherwise
         */
        bool is_convex_quad(index_t t) const;

        /**
         * \brief Orientation predicate
         * \param[in] i , j , k three vertices
         * \return the sign of det(pj-pi,pk-pi)
         */
        virtual Sign orient2d(index_t i,index_t j,index_t k) const=0;

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
        virtual Sign incircle(index_t i,index_t j,index_t k,index_t l) const=0;

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
        virtual index_t create_intersection(
            index_t E1, index_t i, index_t j,
            index_t E2, index_t k, index_t l
        ) = 0;
        
    protected:
        index_t nv_;
        index_t ncnstr_;
        vector<index_t> T_;       /**< triangles vertices array              */
        vector<index_t> Tadj_;    /**< triangles adjacency array             */
        vector<index_t> v2T_;     /**< vertex to triangle back pointer       */
        vector<uint8_t> Tflags_;  /**< triangle flags                        */
        vector<index_t> Tecnstr_; /**< triangle edge constraint              */
        vector<index_t> Tnext_;   /**< doubly connected triangle list        */
        vector<index_t> Tprev_;   /**< doubly connected triangle list        */
        bool delaunay_;           /**< if set, compute a CDT, else just a CT */
        Sign orient_012_;         /**< global triangles orientation          */
    };

    /*****************************************************************/
    
    /**
     * \brief Constrained Delaunay triangulation
     * \details
     *   Example:
     *   \code
     *    CDT cdt;
     *    vec2 p1 = ..., p2 = ..., p3 = ...;
     *    cdt.create_enclosing_triangle(p1,p2,p3); // or enclosing_quad
     *    // insert points
     *    for(...) {
     *      vec2 p = ...;
     *      index_t v = cdt.insert(p);
     *      ...
     *    }
     *    // insert constraints
     *    for(...) {
     *       index_t v1=..., v2=...;
     *       cdt.insert_constraint(v1,v2);   
     *    }
     *    // get triangles
     *    for(index_t t=0; t<cdt.nT(); ++t) {
     *       index_t v1 = cdt.Tv(t,0); 
     *       index_t v2 = cdt.Tv(t,1); 
     *       index_t v3 = cdt.Tv(t,2); 
     *       ... do something with v1,v2,v3
     *    }
     *   \endcode   
     *   If some constraints are intersecting, new vertices are generated. They
     *   can be accessed using the function vec2 CDT::point(index_t v). Vertices
     *   coming from an intersection are between indices nv1 and CDT::nv(), 
     *   where nv1 is the value of CDT::nv() before inserting the constraints
     *   (nv1 corresponds to the number of times CDT::insert() was called plus
     *   the number of points in the enclosing polygon).
     *   If one only wants a constrained triangulation (not Delaunay), 
     *   one can call CDT::set_Delaunay(first) before inserting the points.
     */
    class Experiment_API CDT: public CDTBase {
    public:

        /**
         * \copydoc CDTBase::clear()
         */
        void clear() override;

        /**
         * \brief Creates a first large enclosing triangle
         * \param[in] p1 , p2 , p3 the three vertices of the first triangle
         * \details create_enclosing_triangle() or create_enclosing_quad() 
         *  need to be called before anything else
         */
        void create_enclosing_triangle(
            const vec2& p1, const vec2& p2, const vec2& p3
        );

        /**
         * \brief Creates a first large enclosing quad
         * \param[in] p1 , p2 , p3 , p4 the four vertices of the quad
         * \details The quad needs to be convex. 
         *  create_enclosing_triangle() or create_enclosing_quad() 
         *  need to be called before anything else. 
         */
        void create_enclosing_quad(
            const vec2& p1, const vec2& p2, const vec2& p3, const vec2& p4
        );
        
        /**
         * \brief Inserts a point
         * \details Start by inserting three points forming
         *  a large triangle around all the other points
         * \return the index of the created point
         */
        index_t insert(const vec2& p) {
            point_.push_back(p);
            index_t v = CDTBase::insert(point_.size()-1);
            // If inserted point already existed in
            // triangulation, then nv() did not increase
            if(point_.size() > nv()) {
                point_.pop_back();
            }
            return v;
        }

        /**
         * \copydoc CDTBase::save()
         */
        void save(const std::string& filename) const override;


        vec2 point(index_t v) const {
            geo_debug_assert(v < nv());
            return point_[v];
        }
        
    protected:
        /**
         * \copydoc CDTBase::orient_2d()
         */
        Sign orient2d(index_t i, index_t j, index_t k) const override;

        /**
         * \copydoc CDTBase::incircle()
         */
        Sign incircle(index_t i,index_t j,index_t k,index_t l) const override;

        /**
         * \copydoc CDTBase::create_intersection()
         */
        index_t create_intersection(
            index_t E1, index_t i, index_t j,
            index_t E2, index_t k, index_t l
        ) override;
        
    protected:
        vector<vec2> point_;
    };

    /*****************************************************************/    
}

#endif
