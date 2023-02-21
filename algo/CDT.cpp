/**
 *
 */

#include <OGF/Experiment/algo/CDT.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>
#include <geogram/basic/algorithm.h>

#include <stack>
#include <deque>

// TODO:
// 1) restore Delaunay after constraint forcing
// 2) constaints walking: manage points on constraint
// 3) manage intersecting constraints on the fly
// 4) doubly connected triangle list for S,Q,N
// 5) parameterization by predicates orient2d and incircle
// 6) can we avoid computing o ? (stored in Tflags)


namespace GEO {

    index_t CDT::insert(const vec2& P) {
        CDT_LOG("insert");
        index_t v = point_.size();
        point_.push_back(P);
        v2T_.push_back(index_t(-1));

        // The big triangle comes first
        if(point_.size() == 3) {
            index_t t0 = Tnew();
            Tset(t0, 0, 1, 2, index_t(-1), index_t(-1), index_t(-1));
        }
        if(point_.size() <= 3) {
            return v;
        }
        
        // Phase 1: find triangle that contains vertex i
        Sign o1,o2,o3;
        index_t t1 = locate(v,o1,o2,o3);
        int nb_z = (o1 == ZERO) + (o2 == ZERO) + (o3 == ZERO);
        geo_debug_assert(nb_z != 3);

        // Duplicated vertex
        if(nb_z == 2) {
            CDT_LOG("======> duplicated vertex");
            if(o1 != ZERO) {
                v = Tv(t1,0);
            } else if(o2 != ZERO) {
                v = Tv(t1,1);
            } else {
                geo_debug_assert(o3 !=ZERO);
                v = Tv(t1,2);
            }
            return v;
        }

        // Stack of triangle edges to examine for flipping
        // Note: it is always edge 0 that we examine, since new
        // triangles are always created with v as vertex 0.
        std::stack<index_t> S;

        // Phase 2: split triangle

        // Particular case: v is on edge
        if(nb_z == 1) {
            index_t le1 = (o1 == ZERO) ? 0 :
                          (o2 == ZERO) ? 1 :
                           2 ;
            index_t t2 = Tadj(t1,le1);

            index_t v1 = Tv(t1,le1);
            index_t v2 = Tv(t1,(le1+1)%3);
            index_t v3 = Tv(t1,(le1+2)%3);
            index_t t1_adj2 = Tadj(t1,(le1+1)%3);
            index_t t1_adj3 = Tadj(t1,(le1+2)%3);
                
            if(t2 != index_t(-1)) {
                CDT_LOG("on edge");
                    
                // New vertex is on an edge of t1 and t1 has a neighbor
                // accross that edge. Discard the two triangles t1 and t2
                // adjacent to the edge, and create four new triangles
                // (t1 and t2 are recycled).
                    
                index_t le2 = Tadj_find(t2,t1);

                // I'm always terrified when I write this type of code,
                // so add sanity checks.
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

                Tadj_replace(t1_adj2, t1, t4);
                Tadj_replace(t2_adj3, t2, t3);
                    
                S.push(t1);
                S.push(t2);
                S.push(t3);
                S.push(t4);
                    
            } else {
                CDT_LOG("on border edge");
                    
                // New vertex is on an edge of t1 and t1 has no neighbor
                // accross that edge. Discard t1 and replace it with two
                // new triangles (recycle t1).
                index_t t2 = Tnew();
                    
                Tset(t1,v,v1,v2,t1_adj3,index_t(-1),t2);
                Tset(t2,v,v3,v1,t1_adj2,t1,index_t(-1));

                Tadj_replace(t1_adj2, t1, t2);                    
                    
                S.push(t1);
                S.push(t2);
            }
        } else {
            // New vertex is in t1. Discard t1 and replace it with three
            // new triangles (recycle t1).
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

            Tadj_replace(adj2, t1, t2);
            Tadj_replace(adj3, t1, t3);                
                
            S.push(t1);
            S.push(t2);
            S.push(t3);
        }

        CDT_LOG("inserted, check consistency");
        check_consistency();
        
        if(!delaunay_) {
            return v;
        }

        CDT_LOG("Delaunay-ize");
        
        // Phase 3: Delaunay (by edge flipping)
        while(S.size() != 0) {
            index_t t1 = S.top();
            S.pop();
            index_t t2 = Tadj(t1,0);
            if(t2 == index_t(-1)) {
                continue;
            }
            index_t v1 = Tv(t2,0);
            index_t v2 = Tv(t2,1);
            index_t v3 = Tv(t2,2);
            if(incircle(v1,v2,v3,v) == POSITIVE) {
                swap_edge(t1);
                S.push(t1);
                S.push(t2);                    
            }
        }

        CDT_LOG("check consistency");
        
        check_consistency();

        CDT_LOG("insert OK");
        
        return v;
    }

    void CDT::swap_edge(index_t t1) {
        CDT_LOG("swap edge");
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
        Tadj_replace(t1_adj3, t1, t2);
        Tadj_replace(t2_adj3, t2, t1);
        Tcheck(t1);
        Tcheck(t2);        
        CDT_LOG("swapped edge");        
    }

    void CDT::insert_constraint(index_t i, index_t j) {
        // insert_constraint_simple(i,j);
        insert_constraint_optimized(i,j);
    }
        
    void CDT::insert_constraint_simple(index_t i, index_t j) {
        std::deque< Edge > Q; // Queue of edges to constrain        
        vector< Edge > N;     // New edges to re-Delaunize

        // Step 1: detect intersected edges
        CDT_LOG("Detect intersected edges");
        {
            index_t t_prev = index_t(-1);
            index_t t = vT(i);
            index_t le = index_t(-1);
            bool finished = false;
            for_each_T_around_v(
                i, [&](index_t t_in, index_t le_in) {
                    CDT_LOG("Testing " << t_in);
                    index_t v1 = Tv(t_in, (le_in + 1)%3);
                    index_t v2 = Tv(t_in, (le_in + 2)%3);
                    if(v1 == j || v2 == j) {
                        finished = true;
                        return true;
                    }
                    if(seg_seg_intersect(i,j,v1,v2)) {
                        CDT_LOG("Found intersection : " << v1 << " " << v2);
                        t  = t_in;
                        le = le_in;
                        Trot(t,le);
                        Q.push_back(std::make_pair(Tv(t,1), Tv(t,2)));
                        return true;
                    }
                    return false;
                }
            );
            
            while(!finished) {
                t_prev = t;
                t = Tadj(t,0);
                CDT_LOG("Testing " << t);
                if(Tv(t,0) == j || Tv(t,1) == j || Tv(t,2) == j) {
                    CDT_LOG("Finished");
                    finished = true;
                } else {
                    for(le = 0; le<3; ++le) {
                        if(Tadj(t,le) == t_prev) {
                            continue;
                        }
                        index_t v1 = Tv(t, (le + 1)%3);
                        index_t v2 = Tv(t, (le + 2)%3);
                        if(seg_seg_intersect(i,j,v1,v2)) {
                            CDT_LOG("Found intersection : " << v1 << " " << v2);
                            Trot(t,le);
                            Q.push_back(std::make_pair(Tv(t,1), Tv(t,2)));                            
                            break;
                        }
                    }
                }
            }
        }

        // Step 2: constrain edges
#ifdef CDT_DEBUG        
        save("CDT_0.geogram");
#endif        
        index_t step = 1;
        while(Q.size() != 0) {
            Edge E = Q.back();
            Q.pop_back();
            if(!is_convex_quad(eT(E.first, E.second))) {
                if(Q.size() == 0) {
                    CDT_LOG("... infinite iteration");
#ifdef CDT_DEBUG                    
                    save("last_CDT.geogram");
#endif                    
                    abort();
                }
                Q.push_front(E);
            } else {
                CDT_LOG("Swap");
                swap_edge(E);
                if(seg_seg_intersect(i,j,E.first,E.second)) {
                    Q.push_front(E);
                } else {
                    N.push_back(E);
                }
            }
#ifdef CDT_DEBUG            
            save("CDT_"+String::to_string(step)+".geogram");
#endif            
            step += 1;
        }

        if(!delaunay_) {
            return;
        }

        CDT_LOG("check consistency before re-Delaunize swapped edges");
        check_consistency();
        
        step = 0;
        save("CDT2_"+String::to_string(step)+".geogram");        

        std::cerr << "N=";
        for(auto E: N) {
            std::cerr << "(" << E.first << "," << E.second << ") ";
        }
        std::cerr << std::endl;

        // Step 3: restore Delaunay, version using edges
        bool swap_occured = true;
        while(swap_occured) {
            swap_occured = false;
            for(Edge& E: N) {
                index_t v1 = E.first;
                index_t v2 = E.second;
                CDT_LOG("Test edge " << v1 << "-" << v2);                
                if(
                    (v1 == i && v2 == j) ||
                    (v2 == i && v1 == j)
                ) {
                    CDT_LOG("  ... is the constrained edge (skipping)");
                    continue;
                }
                index_t t1 = eT(v1,v2);
                index_t v0 = Tv(t1,0);
                index_t t2 = Tadj(t1,0);
                index_t e2 = Tadj_find(t2,t1);
                index_t v3 = Tv(t2,e2);
                CDT_LOG("Test trgls " << t1 << "," << t2);
                if(
                    (incircle(v0,v1,v2,v3) || incircle(v3,v2,v1,v0)) &&
                    is_convex_quad(t1)
                ) {
                    std::cerr << "Swap" << std::endl;
                    swap_edge(t1);
                    E.first  = Tv(t1,1);
                    E.second = Tv(t1,2);                    
                    ++step;
                    save("CDT2_"+String::to_string(step)+".geogram");                            
                    swap_occured = true;
                }
            }
        }
    }
    
    void CDT::insert_constraint_optimized(index_t i, index_t j) {
        CDT_LOG("Insert constraint " << i << " " << j);

        std::deque<index_t> Q; // Queue of edges to constrain
        vector<index_t> N; // New triangles to re-Delaunize

        // Step 1: find all the edges that have an intersection 
        // with the constraint [i,j], enqueue them in Q.

        // TODO: -optimize predicate (do not need to test 4 orients
        //        per segment...)
        //       -handle corner-cases (constraint passes through existing vrtx)

        CDT_LOG("Detect intersected edges");

        {
            index_t t_prev = index_t(-1);
            index_t t = vT(i);
            index_t le = index_t(-1);
            bool finished = false;
            for_each_T_around_v(
                i, [&](index_t t_in, index_t le_in) {
                    CDT_LOG("Testing " << t_in);
                    index_t v1 = Tv(t_in, (le_in + 1)%3);
                    index_t v2 = Tv(t_in, (le_in + 2)%3);
                    if(v1 == j || v2 == j) {
                        finished = true;
                        return true;
                    }
                    if(seg_seg_intersect(i,j,v1,v2)) {
                        CDT_LOG("Found intersection : " << v1 << " " << v2);
                        t  = t_in;
                        le = le_in;
                        Trot(t,le);
                        Q.push_back(t);
                        Tmark(t);
                        return true;
                    }
                    return false;
                }
            );
            
            while(!finished) {
                t_prev = t;
                t = Tadj(t,0);
                CDT_LOG("Testing " << t);
                if(Tv(t,0) == j || Tv(t,1) == j || Tv(t,2) == j) {
                    CDT_LOG("Finished");
                    finished = true;
                } else {
                    for(le = 0; le<3; ++le) {
                        if(Tadj(t,le) == t_prev) {
                            continue;
                        }
                        index_t v1 = Tv(t, (le + 1)%3);
                        index_t v2 = Tv(t, (le + 2)%3);
                        if(seg_seg_intersect(i,j,v1,v2)) {
                            CDT_LOG("Found intersection : " << v1 << " " << v2);
                            Trot(t,le);
                            Q.push_back(t);
                            Tmark(t);
                            break;
                        }
                    }
                }
            }
        }

                
        // Step 2: constrain edges
        // Version using triangle-as-edges
#ifdef CDT_DEBUG        
        save("CDT_0.geogram");
#endif        
        index_t step = 1;
        while(Q.size() != 0) {
            CDT_LOG("========== STEP " << step);
            debug_show_triangles(Q,"Q");
            index_t t1 = Q.back();
            Q.pop_back();
            CDT_LOG("Examining " << t1 << "," << Tadj(t1,0));
            if(!Tis_marked(t1)) {
                CDT_LOG("   Not marked");
                continue;
            }
            if(!is_convex_quad(t1)) {
                CDT_LOG("   Not convex");
                if(Q.size() == 0) {
                    CDT_LOG("... infinite iteration");
#ifdef CDT_DEBUG                    
                    save("last_CDT.geogram");
#endif                    
                    abort();
                }
                Q.push_front(t1);
            } else {
                CDT_LOG("Swap");
                index_t t2 = Tadj(t1,0);

                bool no_isect  = !Tis_marked(t2);
                index_t v0     = Tv(t1,0);
                bool t2v0_t1v2 = (Tis_marked(t2) && Tv(t2,0) == Tv(t1,2));
                bool t2v0_t1v1 = (Tis_marked(t2) && Tv(t2,0) == Tv(t1,1));
                
                swap_edge(t1);
                if(no_isect) {
                    CDT_LOG("   OK (no isect)");
                    Tunmark(t1);
                    Trot(t1,2); // so that new edge is edge 0
                    N.push_back(t1);
                    CDT_LOG("New edge: " << t1 << "->(" << Tv(t1,1) << "," << Tv(t1,2) << ")");
                } else {
                    Sign o = orient2d(i,j,v0);
                    if(t2v0_t1v2) {
                        if(o >= 0) {
                            CDT_LOG(t1 << ":no isect  " << t2 << ":isect");
                            Tunmark(t1);
                            Trot(t1,2); // so that new edge is edge 0                            
                            N.push_back(t1);
                            CDT_LOG("New edge: " << t1 << "->(" << Tv(t1,1) << "," << Tv(t1,2) << ")");                            
                        } else {
                            CDT_LOG(t1 << ":isect  " << t2 << ":isect");
                            Trot(t1,2);                            
                            Q.push_front(t1);
                        }
                    } else {
                        geo_debug_assert(t2v0_t1v1);
                        if(o > 0) {
                            Trot(t2,1);                            
                            CDT_LOG(t1 << ":isect  " << t2 << ":isect");
                            Q.push_front(t1);                            
                        } else {
                            CDT_LOG(t1 << ":isect  " << t2 << ":no isect");
                            Tunmark(t2);
                            Trot(t2,1); // so that new edge is edge 0
                            N.push_back(t2);
                            Q.push_front(t1);
                            CDT_LOG("New edge: " << t2 << "->(" << Tv(t2,1) << "," << Tv(t2,2) << ")");                            
                        }
                    }
                }
            }
#ifdef CDT_DEBUG            
            save("CDT_"+String::to_string(step)+".geogram");
#endif            
            step += 1;
        }

        if(true || !delaunay_) {
            return;
        }

        CDT_LOG("check consistency before re-Delaunize swapped edges");
        check_consistency();
        
        step = 0;
        save("CDT2_"+String::to_string(step)+".geogram");        

        debug_show_triangles(N,"N");


        // Restore CDT condition
        bool swap_occured = true;
        while(swap_occured) {
            swap_occured = false;
            for(index_t t1: N) {
                index_t v1 = Tv(t1,1);
                index_t v2 = Tv(t1,2); 
                if(
                    (v1 == i && v2 == j) ||
                    (v2 == i && v1 == j)
                ) {
                    continue;
                }
                index_t v0 = Tv(t1,0);
                index_t t2 = Tadj(t1,0);
                index_t e2 = Tadj_find(t2,t1);
                index_t v3 = Tv(t2,e2);
                std::cerr << "Test " << t1 << "," << t2 << std::endl;
                if(incircle(v0,v1,v2,v3) || incircle(v3,v2,v1,v0)) {
                    std::cerr << "Swap" << std::endl;
                    swap_edge(t1);
                    ++step;
                    save("CDT2_"+String::to_string(step)+".geogram");                            
                    swap_occured = true;
                }
            }
        }
    }

    index_t CDT::eT(index_t v1, index_t v2) {
        index_t result = index_t(-1);
        for_each_T_around_v(
            v1,
            [&](index_t t, index_t lv)->bool {
                if(Tv(t, (lv+1)%3) == v2) {
                    if(Tv(t, (lv+2)%3) != v1) {
                        Trot(t, (lv+2)%3);
                    }
                    result = t;
                    return true;
                } else if(Tv(t, (lv+1)%3) == v1) {
                    if(Tv(t, (lv+2)%3) != v2) {
                        Trot(t, (lv+2)%3);
                    }
                    result = t;                    
                    return true;
                }
                return false;
            }
        );
        geo_debug_assert(result != index_t(-1));
        geo_debug_assert(
            (Tv(result,1) == v1 && Tv(result,2) == v2) ||
            (Tv(result,1) == v2 && Tv(result,2) == v1) 
        );
        return result;
    }

    void CDT::swap_edge(Edge& E) {
        index_t t = eT(E.first, E.second);
        swap_edge(t);
        E = std::make_pair(Tv(t,0), Tv(t,1));
    }
    
    index_t CDT::locate(index_t v, Sign& o1, Sign& o2, Sign& o3) const {
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

    Sign CDT::orient2d(index_t i, index_t j, index_t k) const {
        return PCK::orient_2d(point_[i], point_[j], point_[k]);
    }

    Sign CDT::incircle(index_t i, index_t j, index_t k, index_t l) const {
        return PCK::in_circle_2d_SOS(
            point_[i].data(), point_[j].data(), point_[k].data(),
            point_[l].data()
        );
    }

    bool CDT::seg_seg_intersect(
        index_t i, index_t j, index_t k, index_t l
    ) const  {
        Sign o1 = orient2d(i,j,k);
        Sign o2 = orient2d(i,j,l);
        if(o1*o2 >=0) {
            return false;
        }
        Sign o3 = orient2d(k,l,i);
        Sign o4 = orient2d(k,l,j);
        return (o3*o4 < 0);
    }

    bool CDT::is_convex_quad(index_t t) const {
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

        CDT_LOG(
             v1 << " " << v2 << " " <<
             v4 << " " << v3 << " " <<
             "convex?=" << result
         );
        
        return result;
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
        mesh_save(M, filename);
    }
}

