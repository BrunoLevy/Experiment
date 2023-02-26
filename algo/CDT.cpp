/**
 *
 */

// Reference: S. W. Sloan, a fast algorithm for generating
// constrained Delaunay triangulation, 1992,
// Computers and Structures
//
// Specificities of this implementation:
//
// - Edges are systematically manipulated through triangles,
// and these triangles are rotated in-place in the mesh,
// in such a way that the edge we are talking about is
// systematically edge 0 (with vertices 1 and 2).
//
// - The constraint-enforcing step manipulates a queue Q
// of edges encoded this way. It examines pairs of triangles
// t1,t2=Tadj(t1,0), decides whether to swap their common
// edge (based on convexity test and intersection of
// t1's edge 0 with the constraint). In fact, this intersection
// test only depends on the combinatorics of (t1,t2) (two cases)
// and the position of t1' vertex 0 relative to the constraint
// (two cases), that makes 4 cases in total. In these cases,
//    - t1 can either leave Q or be enqueued again
//    - t2 was always in Q already (because it has an edge
//      that has an intersection with the constraint), but
//      there is one case where it leaves Q
// We need to mark/unmark triangles for one case: when t2 is unmarked,
// it means there is no intersection.


// TODO:
// 1) can we avoid computing o ? (stored in Tflags)
// 2) predicate cache:
//     - test whether needed.
//     - find a "noalloc" implementation (std::make_heap ?)
// 3) faster locate()
// 4) management of boundary: can we have "vertex at infinity" like in CGAL ?
// 5) tag/remove external triangles

// NOTE - TOREAD:
// https://www.sciencedirect.com/science/article/pii/S0890540112000752
// (other ways of doing exact computations using FP)

#include <OGF/Experiment/algo/CDT.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/basic/algorithm.h>
#include <stack>

#ifdef GEO_DEBUG
#define CDT_DEBUG
#endif

#ifdef CDT_DEBUG
#define CDT_LOG(X) std::cerr << X << std::endl
#else
#define CDT_LOG(X)
#endif

namespace GEO {

    CDTBase::CDTBase() : nv_(0), ncnstr_(0), delaunay_(true) {
    }

    CDTBase::~CDTBase() {
    }

    void CDTBase::clear() {
        nv_ = 0;
        ncnstr_ = 0;
        T_.resize(0);
        Tadj_.resize(0);
        v2T_.resize(0);
        Tflags_.resize(0);
        Tecnstr_.resize(0);
        Tnext_.resize(0);
        Tprev_.resize(0);
    }

    void CDTBase::create_enclosing_triangle(
        index_t v0, index_t v1, index_t v2
    ) {
        nv_ = 3;
        v2T_.resize(3);
        geo_debug_assert(v0 <= 3);
        geo_debug_assert(v1 <= 3);
        geo_debug_assert(v2 <= 3);
        index_t t0 = Tnew();
        Tset(t0, v0, v1, v2, index_t(-1), index_t(-1), index_t(-1));
        orient_012_ = orient2d(0,1,2);
    }

    void CDTBase::create_enclosing_quad(
        index_t v0, index_t v1, index_t v2, index_t v3
    ) {
        nv_ = 4;
        v2T_.resize(4);
        geo_debug_assert(v0 <= 4);
        geo_debug_assert(v1 <= 4);
        geo_debug_assert(v2 <= 4);
        geo_debug_assert(v3 <= 4);        
        index_t t0 = Tnew();
        index_t t1 = Tnew();        
        Tset(t0, v0, v1, v3, t1, NO_INDEX, NO_INDEX);
        Tset(t1, v3, v1, v2, NO_INDEX, NO_INDEX, t0);
        orient_012_ = orient2d(0,1,2);
        geo_debug_assert(is_convex_quad(t0));
        if(incircle(v0,v1,v2,v3) == POSITIVE) {
            swap_edge(t0);
        }
    }
    
    index_t CDTBase::insert(index_t v) {
        geo_debug_assert(v == nv());
        v2T_.push_back(index_t(-1));
        ++nv_;
        
        // Phase 1: find triangle that contains vertex i
        Sign o1,o2,o3;
        index_t t = locate(v,o1,o2,o3);
        int nb_z = (o1 == ZERO) + (o2 == ZERO) + (o3 == ZERO);
        geo_debug_assert(nb_z != 3);

        // Duplicated vertex
        if(nb_z == 2) {
            CDT_LOG("duplicated vertex");
            v = (o1 != ZERO) ? Tv(t,0) :
                (o2 != ZERO) ? Tv(t,1) :
                               Tv(t,2) ;
            v2T_.pop_back();
            --nv_;
            return v;
        }

        // Stack of triangle edges to examine for flipping
        // Note: it is always edge 0 that we examine, since new
        // triangles are always created with v as vertex 0.
        DList S;

        // Phase 2: split triangle
        // Particular case: v is on edge
        if(nb_z == 1) {
            CDT_LOG("vertex on edge");            
            index_t le = (o1 == ZERO) ? 0 :
                         (o2 == ZERO) ? 1 :
                          2 ;
            insert_vertex_in_edge(v,t,le,&S);
        } else {
            insert_vertex_in_triangle(v,t,&S);
        }

        check_consistency();
        
        if(!delaunay_) {
            return v;
        }

        Delaunayize_vertex_neighbors(v,S);
        check_consistency();
        
        return v;
    }

    
    void CDTBase::insert_constraint(index_t i, index_t j) {

        CDT_LOG("insert constraint: " << i << "-" << j);

        ++ncnstr_;

        // First vertex coming from constraints intersection
        index_t first_v_isect = nv_;
        
        DList Q; // Queue of edges to constrain
        DList N; // New edges to re-Delaunayize

        while(i != j) {
            // Step 1: find all the edges that have an intersection 
            // with the constraint [i,j], enqueue them in Q.
            // Stop at vertex on constraint or constraint intersection
            // if any (returned in k)
            index_t k = find_intersected_edges(i,j,Q);
            // Step 2: constrain edges
            constrain_edges(i,k,Q,N);
            check_consistency();
            // Step 3: restore Delaunay condition
            if(delaunay_) {
                Delaunayize_new_edges(N);
            }
            check_consistency();            
            i = k;
        }

        if(!delaunay_) {
            return;
        }

        // Delaunayize triangles around vertices coming from
        // constraint intersections
        DList S;        
        for(index_t v=first_v_isect; v<nv(); ++v) {
            // We cannot use for_each_triangle_around_vertex()
            // because we need to Trot() t during traveral,
            // to have v has t's vertex 0
            // But the good news is that v is never on the border,
            // so traversal is easier.
            index_t t0 = vT(v); // Need to store it, because we Trot()
            index_t t = t0;
            do {
                index_t lv = Tv_find(t,v);
                Trot(t,lv);
                geo_debug_assert(Tv(t,0) == v);
                DList_push_back(S,t);
                t = Tadj(t, 1);
                geo_debug_assert(t != NO_INDEX);
            } while(t != t0);            
            Delaunayize_vertex_neighbors(v,S);
        }
    }

    index_t CDTBase::find_intersected_edges(
        index_t i, index_t j, DList& Q
    ) {

        CDT_LOG("Find intersected edges: " << i << "-" << j);
        
        // Walk from i to j, detect intersected edges and push
        // them to Q. At each step, we are either on a triangle
        // (generic case) or on a vertex (when we start from i, or
        // when there is a vertex that is exactly on [i,k]).
        // We keep track of the previous triangle or previous vertex
        // to make sure we don't go backwards.
        
        index_t t_prev = NO_INDEX;
        index_t v_prev = NO_INDEX;
        index_t t      = NO_INDEX;
        index_t v      = i;
        index_t t_next = NO_INDEX;
        index_t v_next = NO_INDEX;

        // We stop at the first encountered vertex
        // The function could be also used to traverse
        // the whole intersected segment, but doing it
        // like that is better to schedule constraint
        // enforcement / re-Delaunay, especially when
        // new vertices are inserted.
        while(v == NO_INDEX || v == i) {
            CDT_LOG(
                "   t=" << int(t) << " v=" << int(v) << "   "
                "t_prev=" << int(t_prev) << " v_prev=" << int(v_prev) << "   "
            );

            if(v != NO_INDEX) {
                // We are on a vertex (when we start from i, or when there
                // is a vertex exactly on [i,j], or when there was a constraints
                // intersection right before
                geo_debug_assert(t == NO_INDEX);
                // Turn around the triangles incident to v
                for_each_T_around_v(
                    v, [&](index_t t_around_v, index_t le) {
                        // If triangle around vertex is the triangle
                        // we came from, continue iteration around vertex
                        if(t_around_v == t_prev) {
                            return false;
                        }
                        index_t v1 = Tv(t_around_v, (le + 1)%3);
                        index_t v2 = Tv(t_around_v, (le + 2)%3);

                        // Are we arrived at j ? 
                        if(v1 == j || v2 == j) {
                            v_next = j;
                            t_next = NO_INDEX;
                            // Edge is flagged as constraint here, because
                            // it will not be seen by constraint enforcement.
                            index_t le_cnstr_edge = (v1 == j) ? le : (le+1)%3;
                            Tset_edge_cnstr_with_neighbor(
                                t_around_v, le_cnstr_edge, ncnstr_-1
                            );
                            return true;
                        }
                        
                        Sign o1 = orient2d(i,j,v1);
                        Sign o2 = orient2d(i,j,v2);                        
                        Sign o3 = orient2d(v1,v2,j);
                        Sign o4 = orient_012_; //equivalent to orient2d(v1,v2,i)
                        if(o1*o2 < 0 && o3*o4 < 0) {
                            Trot(t_around_v,le); // so that le becomes edge 0
                            t_next = t_around_v; // added to Q during next round
                            v_next = NO_INDEX;
                            return true;
                        } else {
                            // Special case: v1 or v2 is exactly on [i,j]
                            // Edge is flagged as constraint here, because
                            // it will not be seen by constraint enforcement.
                            geo_debug_assert(o1 != ZERO || o2 != ZERO);
                            if(o1 == ZERO && o3*o4 < 0 && v1 != v_prev) {
                                t_next = NO_INDEX;
                                v_next = v1;
                                Tset_edge_cnstr_with_neighbor(
                                    t_around_v, (le + 2)%3, ncnstr_-1
                                );
                                return true;
                            } else if(o2 == ZERO && o3*o4 < 0 && v2 != v_prev) {
                                t_next = NO_INDEX;
                                v_next = v2;
                                Tset_edge_cnstr_with_neighbor(
                                    t_around_v, (le + 1)%3, ncnstr_-1
                                );                                
                                return true;
                            }
                        }
                        return false;
                    }
                );
            } else {
                // Generic case: we are on a triangle
                geo_debug_assert(t != NO_INDEX);
                // Are we arrived at j ? 
                if(Tv(t,0) == j || Tv(t,1) == j || Tv(t,2) == j) {
                    v_next = j;
                    t_next = NO_INDEX;
                } else {
                    // Test the three edges of the triangle
                    for(index_t le = 0; le<3; ++le) {
                        // Skip the edge we are coming from
                        if(Tadj(t,le) == t_prev) {
                            continue;
                        }
                        // Test whether [v1,v2] intersects the
                        // support line of (i,j). No need to test
                        // the *segment* [i,j]: we know the line enters
                        // the triangle, it is how we came here,
                        // and we know it leaves it, else j would
                        // have been one of the triangle's vertices.
                        index_t v1 = Tv(t, (le + 1)%3);
                        index_t v2 = Tv(t, (le + 2)%3);
                        Sign o1 = orient2d(i,j,v1);
                        Sign o2 = orient2d(i,j,v2);                        
                        if(o1*o2 < 0) {
                            // [v1,v2] has a frank intersection with [i,j]
                            Trot(t,le); // So that edge 0 is intersected edge
                            if(Tedge_is_constrained(t,0)) {
                                CDT_LOG("   ====> Constraints intersection");
                                v_next = create_intersection(
                                    ncnstr()-1, i, j,
                                    Tedge_cnstr(t,0), v1, v2
                                );
                                insert_vertex_in_edge(v_next,t,0);
                                t_next = NO_INDEX;
                            } else {
                                CDT_LOG("   Intersection: t="
                                        << t << " E=" << v1 << "-" << v2
                                       );
                                DList_push_back(Q,t);
                                Tmark(t);
                                t_next = Tadj(t,0);
                                v_next = NO_INDEX;
                            }
                            break;
                        } else {
                            // Special case: v1 or v2 is exactly on [i,j]
                            geo_debug_assert(o1 != ZERO || o2 != ZERO);
                            if(o1 == ZERO) {
                                t_next = NO_INDEX;
                                v_next = v1;
                                break;
                            } else if(o2 == ZERO) {
                                t_next = NO_INDEX;
                                v_next = v2;
                                break;
                            }
                        }
                    }
                }
            }
            t_prev = t;
            v_prev = v;
            t = t_next;            
            v = v_next;
            if(v != NO_INDEX) {
                return v;
            }
        }
        return NO_INDEX;
    }
    
    void CDTBase::constrain_edges(index_t i, index_t j, DList& Q, DList& N) {

        // Edge le of triangle t has no isect with cnstr, it is a new edge
        auto new_edge = [&](index_t t,index_t le) {
            Trot(t,le);
            Tunmark(t);
            if(
                (Tv(t,1) == i && Tv(t,2) == j) ||
                (Tv(t,1) == j && Tv(t,2) == i)
            ) {
                Tset_edge_cnstr_with_neighbor(t,0,ncnstr_-1);
            } else {
                DList_push_back(N,t);
            }
        };

        // Edge le of triangle t still has an isect with cnstr
        auto isect_edge = [&](index_t t, index_t le) {
            Trot(t,le);
            DList_push_front(Q,t);
        };

        while(!Q.empty()) {
            index_t t1 = DList_pop_back(Q);
            if(!Tis_marked(t1)) {
                continue;
            }
            if(!is_convex_quad(t1)) {
                // If the only remaining edge to flip does not form a convex
                // quad, it means we are going to flip forever ! (shoud not
                // happen)
                geo_assert(!Q.empty());
                DList_push_front(Q,t1);
            } else {
                index_t t2 = Tadj(t1,0);
                bool no_isect  = !Tis_marked(t2);
                index_t v0     = Tv(t1,0);
                bool t2v0_t1v2 = (Tis_marked(t2) && Tv(t2,0) == Tv(t1,2));
                bool t2v0_t1v1 = (Tis_marked(t2) && Tv(t2,0) == Tv(t1,1));
                geo_argused(t2v0_t1v1);
                
                swap_edge(t1);
                if(no_isect) {
                    new_edge(t1,2);
                } else {
                    // See comment at beginning of file (a variation in Sloan's
                    // method that makes better use of the combinatorics)
                    Sign o = orient2d(i,j,v0);
                    if(t2v0_t1v2) {
                        if(o >= 0) {
                            new_edge(t1,2);
                        } else {
                            isect_edge(t1,2);
                        }
                    } else {
                        geo_debug_assert(t2v0_t1v1);
                        if(o > 0) {
                            Trot(t2,1); 
                            isect_edge(t1,0);
                        } else {
                            DList_remove(Q,t2);
                            new_edge(t2,1);
                            isect_edge(t1,0);
                        }
                    }
                }
            }
        }
    }

    void CDTBase::Delaunayize_vertex_neighbors(index_t v, DList& S) {
        while(!S.empty()) {
            index_t t1 = DList_pop_back(S);
            geo_debug_assert(Tv(t1,0) == v);
            if(Tedge_is_constrained(t1,0)) {
                continue;
            }
            index_t t2 = Tadj(t1,0);
            if(t2 == index_t(-1)) {
                continue;
            }
            index_t v1 = Tv(t2,0);
            index_t v2 = Tv(t2,1);
            index_t v3 = Tv(t2,2);
            if(incircle(v1,v2,v3,v) == POSITIVE) {
                swap_edge(t1);
                DList_push_back(S,t1);
                DList_push_back(S,t2);                    
            }
        }
    }
    
    void CDTBase::Delaunayize_new_edges(DList& N) {
        bool swap_occured = true;
        while(swap_occured) {
            swap_occured = false;
            for(index_t t1 = N.front; t1 != NO_INDEX; t1 = Tnext(t1)) {
                if(Tedge_is_constrained(t1,0)) {
                    continue;
                }
                index_t v1 = Tv(t1,1);
                index_t v2 = Tv(t1,2);
                index_t v0 = Tv(t1,0);
                index_t t2 = Tadj(t1,0);
                if(t2 == NO_INDEX) {
                    continue;
                }
                index_t e2 = Tadj_find(t2,t1);
                index_t v3 = Tv(t2,e2);
                if(incircle(v0,v1,v2,v3) == POSITIVE) {
                    swap_edge(t1);
                    swap_occured = true;
                } 
            }
        }
        DList_clear(N);
    }
    
    index_t CDTBase::locate(index_t v, Sign& o1, Sign& o2, Sign& o3) const {
        // This version: linear scan
        // TODO: faster version
        for(index_t t=0; t<nT(); ++t) {
            index_t i = Tv(t,0);
            index_t j = Tv(t,1);
            index_t k = Tv(t,2);
            o1 = orient2d(v,j,k);
            o2 = orient2d(v,k,i);
            o3 = orient2d(v,i,j);
            if(o1*o2 >= 0 && o2*o3 >= 0 && o3*o1 >= 0) {
                return t;
            }
        }
        geo_assert_not_reached;
    }

    /***************** Triangulation surgery (boring code ahead) *********/
    
    void CDTBase::insert_vertex_in_edge(
        index_t v, index_t t, index_t le1,
        DList* S
    ) {
        index_t cnstr = Tedge_cnstr(t,le1);
        index_t t1 = t;
        index_t t2 = Tadj(t1,le1);
        index_t v1 = Tv(t1,le1);
        index_t v2 = Tv(t1,(le1+1)%3);
        index_t v3 = Tv(t1,(le1+2)%3);
        index_t t1_adj2 = Tadj(t1,(le1+1)%3);
        index_t t1_adj3 = Tadj(t1,(le1+2)%3);
        if(t2 != NO_INDEX) {
            // New vertex is on an edge of t1 and t1 has a neighbor
            // accross that edge. Discard the two triangles t1 and t2
            // adjacent to the edge, and create four new triangles
            // (t1 and t2 are recycled).
            index_t le2 = Tadj_find(t2,t1);
            geo_debug_assert(Tv(t2, (le2+1)%3) == v3);
            geo_debug_assert(Tv(t2, (le2+2)%3) == v2);
            index_t v4 = Tv(t2,le2);
            index_t t2_adj2 = Tadj(t2,(le2+1)%3);
            index_t t2_adj3 = Tadj(t2,(le2+2)%3);
            index_t t3 = Tnew();
            index_t t4 = Tnew();
            Tset(t1,v,v1,v2,t1_adj3,t2,t4);
            Tset(t2,v,v2,v4,t2_adj2,t3,t1);
            Tset(t3,v,v4,v3,t2_adj3,t4,t2);
            Tset(t4,v,v3,v1,t1_adj2,t1,t3);
            Tadj_back_connect(t1,0,t1);
            Tadj_back_connect(t2,0,t2);
            Tadj_back_connect(t3,0,t2);
            Tadj_back_connect(t4,0,t1);
            Tset_edge_cnstr(t1,1,cnstr);
            Tset_edge_cnstr(t2,2,cnstr);
            Tset_edge_cnstr(t3,1,cnstr);
            Tset_edge_cnstr(t4,2,cnstr);
            if(S != nullptr) {
                DList_push_back(*S, t1);
                DList_push_back(*S, t2);
                DList_push_back(*S, t3);
                DList_push_back(*S, t4);                
            }
        } else {
            // New vertex is on an edge of t1 and t1 has no neighbor
            // accross that edge. Discard t1 and replace it with two
            // new triangles (recycle t1).
            t2 = Tnew();
            Tset(t1,v,v1,v2,t1_adj3,NO_INDEX,t2);
            Tset(t2,v,v3,v1,t1_adj2,t1,NO_INDEX);
            Tadj_back_connect(t1,0,t1);
            Tadj_back_connect(t2,0,t1);            
            Tset_edge_cnstr(t1,1,cnstr);
            Tset_edge_cnstr(t2,2,cnstr);
            if(S != nullptr) {
                DList_push_back(*S, t1);
                DList_push_back(*S, t2);
            }
        }
    }

    void CDTBase::insert_vertex_in_triangle(
        index_t v, index_t t, DList* S
    ) {
        // New vertex is in t1. Discard t1 and replace it with three
        // new triangles (recycle t1).
        index_t t1 = t;
        index_t v1 = Tv(t1,0);
        index_t v2 = Tv(t1,1);
        index_t v3 = Tv(t1,2);
        index_t adj1 = Tadj(t1,0);
        index_t adj2 = Tadj(t1,1);
        index_t adj3 = Tadj(t1,2);
        index_t t2 = Tnew();
        index_t t3 = Tnew();
        Tset(t1,v,v2,v3,adj1,t2,t3);
        Tset(t2,v,v3,v1,adj2,t3,t1);
        Tset(t3,v,v1,v2,adj3,t1,t2);
        Tadj_back_connect(t1,0,t1);
        Tadj_back_connect(t2,0,t1);
        Tadj_back_connect(t3,0,t1);
        if(S != nullptr) {
            DList_push_back(*S,t1);
            DList_push_back(*S,t2);
            DList_push_back(*S,t3);            
        }
    }
    
    void CDTBase::swap_edge(index_t t1) {
        geo_debug_assert(!Tedge_is_constrained(t1,0));
        index_t v1 = Tv(t1,0);
        index_t v2 = Tv(t1,1);
        index_t v3 = Tv(t1,2);                        
        index_t t1_adj2 = Tadj(t1,1);
        index_t t1_adj3 = Tadj(t1,2);
        index_t t2 = Tadj(t1,0);
        index_t le2 = Tadj_find(t2,t1);
        index_t v4 = Tv(t2,le2);
        geo_debug_assert(Tv(t2,(le2+1)%3) == v3);
        geo_debug_assert(Tv(t2,(le2+2)%3) == v2);
        Tcheck(t1);
        Tcheck(t2);
        index_t t2_adj2 = Tadj(t2,(le2+1)%3);
        index_t t2_adj3 = Tadj(t2,(le2+2)%3);
        Tset(t1,v1,v4,v3,t2_adj3,t1_adj2,t2);
        Tset(t2,v1,v2,v4,t2_adj2,t1,t1_adj3);
        Tadj_back_connect(t1,0,t2);
        Tadj_back_connect(t1,1,t1);
        Tadj_back_connect(t2,0,t2);
        Tadj_back_connect(t2,2,t1);        
        Tcheck(t1);
        Tcheck(t2);        
    }
    
    bool CDTBase::is_convex_quad(index_t t) const {
        index_t v1 = Tv(t,0);
        index_t v2 = Tv(t,1);
        index_t v3 = Tv(t,2);        
        index_t t2 = Tadj(t,0);
        index_t le2 = Tadj_find(t2,t);
        index_t v4 = Tv(t2,le2);

        Sign o1 = orient2d(v2,v1,v4);
        Sign o2 = orient2d(v4,v2,v3);
        Sign o3 = orient2d(v3,v4,v1);
        Sign o4 = orient2d(v1,v3,v2);

        bool result = (
            o1*o2 > 0 &&
            o2*o3 > 0 &&
            o3*o4 > 0 &&
            o4*o1 > 0
        );
        
        return result;
    }

    bool CDTBase::Tedge_is_Delaunay(index_t t1, index_t le1) const {
        if(Tedge_is_constrained(t1,le1)) {
            return true;
        }
        index_t t2 = Tadj(t1,le1);
        if(t2 == NO_INDEX) {
            return true;
        }
        index_t le2 = Tadj_find(t2,t1);
        index_t v1 = Tv(t1,le1);
        index_t v2 = Tv(t1,(le1+1)%3);
        index_t v3 = Tv(t1,(le1+2)%3);
        index_t v4 = Tv(t2,le2);
        return incircle(v1,v2,v3,v4) <= 0;
    }
    
    /********************************************************************/

    void CDT::clear() {
        CDTBase::clear();
        point_.resize(0);
    }
    
    void CDT::create_enclosing_triangle(
        const vec2& p1, const vec2& p2, const vec2& p3
    ) {
        geo_assert(nv() == 0);
        geo_assert(nT() == 0);        
        point_.push_back(p1);
        point_.push_back(p2);
        point_.push_back(p3);
        CDTBase::create_enclosing_triangle(0,1,2);
    }

    void CDT::create_enclosing_quad(
        const vec2& p1, const vec2& p2, const vec2& p3, const vec2& p4
    ) {
        geo_assert(nv() == 0);
        geo_assert(nT() == 0);        
        point_.push_back(p1);
        point_.push_back(p2);
        point_.push_back(p3);        
        point_.push_back(p4);
        CDTBase::create_enclosing_quad(0,1,2,3);
    }
    
    Sign CDT::orient2d(index_t i, index_t j, index_t k) const {
        geo_debug_assert(i < nv());
        geo_debug_assert(j < nv());
        geo_debug_assert(k < nv());        
        return PCK::orient_2d(point_[i], point_[j], point_[k]);
    }

    Sign CDT::incircle(index_t i, index_t j, index_t k, index_t l) const {
        geo_debug_assert(i < nv());
        geo_debug_assert(j < nv());
        geo_debug_assert(k < nv());
        geo_debug_assert(l < nv());                
        return PCK::in_circle_2d_SOS(
            point_[i].data(), point_[j].data(), point_[k].data(),
            point_[l].data()
        );
    }

    index_t CDT::create_intersection(
        index_t E1, index_t i, index_t j,
        index_t E2, index_t k, index_t l
    ) {
        geo_argused(E1);
        geo_argused(E2);        
        geo_debug_assert(i < nv());
        geo_debug_assert(j < nv());
        geo_debug_assert(k < nv());
        geo_debug_assert(l < nv());
        geo_debug_assert(E1 < ncnstr());
        geo_debug_assert(E2 < ncnstr());
        vec2 U = point_[j] - point_[i];
        vec2 V = point_[l] - point_[k];
        vec2 D = point_[k] - point_[i];
        double delta = det(U,V);
        double t = det(D,V)/delta;
        vec2 P = point_[i] + t*U;
        point_.push_back(P);
        v2T_.push_back(NO_INDEX);
        index_t v = nv_;
        ++nv_;
        return v;
    }
    
    void CDT::save(const std::string& filename) const {
        Mesh M;
        M.vertices.set_dimension(2);
        for(const vec2& P: point_) {
            M.vertices.create_vertex(P.data());
        }
        for(index_t t=0; t<nT(); ++t) {
            index_t i = Tv(t,0);
            index_t j = Tv(t,1);
            index_t k = Tv(t,2);
            M.facets.create_triangle(i,j,k);
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

        Attribute<bool> constraint(M.facet_corners.attributes(), "constraint");
        for(index_t c: M.facet_corners) {
            index_t t  = c/3;
            index_t lv = c%3;
            constraint[c] =
                Tedge_is_constrained(t, (lv+1)%3) ||
                Tedge_is_constrained(t, (lv+2)%3) ; 
        }
        
        mesh_save(M, filename);
    }
}

