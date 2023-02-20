/**
 *
 */

#include <OGF/Experiment/algo/CDT.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>

#include <stack>
#include <deque>

namespace GEO {

    index_t CDT::insert(const vec2& P) {
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
        geo_assert(nb_z != 3);

        // Duplicated vertex
        if(nb_z == 2) {
            std::cerr << "======> duplicated vertex" << std::endl;
            if(o1 != ZERO) {
                v = Tv(t1,0);
            } else if(o2 != ZERO) {
                v = Tv(t1,1);
            } else {
                geo_assert(o3 !=ZERO);
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
                std::cerr << "on edge" << std::endl;
                    
                // New vertex is on an edge of t1 and t1 has a neighbor
                // accross that edge. Discard the two triangles t1 and t2
                // adjacent to the edge, and create four new triangles
                // (t1 and t2 are recycled).
                    
                index_t le2 = Tadj_find(t2,t1);

                // I'm always terrified when I write this type of code,
                // so add sanity checks.
                geo_assert(Tv(t2, (le2+1)%3) == v3);
                geo_assert(Tv(t2, (le2+2)%3) == v2);
                    
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
                // Tadj_replace(t1_adj3, t1, t1);
                // Tadj_replace(t2_adj2, t2, t2);
                Tadj_replace(t2_adj3, t2, t3);
                    
                S.push(t1);
                S.push(t2);
                S.push(t3);
                S.push(t4);
                    
            } else {
                std::cerr << "on border edge" << std::endl;
                    
                // New vertex is on an edge of t1 and t1 has no neighbor
                // accross that edge. Discard t1 and replace it with two
                // new triangles (recycle t1).
                index_t t2 = Tnew();
                    
                Tset(t1,v,v1,v2,t1_adj3,index_t(-1),t2);
                Tset(t2,v,v3,v1,t1_adj2,t1,index_t(-1));

                Tadj_replace(t1_adj2, t1, t2);                    
                // Tadj_replace(t1_adj3, t1, t1);
                    
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

            // Tadj_replace(adj1, t1, t1);
            Tadj_replace(adj2, t1, t2);
            Tadj_replace(adj3, t1, t3);                
                
            S.push(t1);
            S.push(t2);
            S.push(t3);
        }

        // Phase 3: Delaunay
        if(!delaunay_) {
            return v;
        }
        
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
        return v;
    }

    void CDT::swap_edge(index_t t1) {
        index_t v1 = Tv(t1,0);
        index_t v2 = Tv(t1,1);
        index_t v3 = Tv(t1,2);                        
        index_t t1_adj2 = Tadj(t1,1);
        index_t t1_adj3 = Tadj(t1,2);

        index_t t2 = Tadj(t1,0);
        index_t le2 = Tadj_find(t2,t1);
        index_t v4 = Tv(t2,le2);
        geo_assert(Tv(t2,(le2+1)%3) == v3);
        geo_assert(Tv(t2,(le2+2)%3) == v2);
        index_t t2_adj2 = Tadj(t2,(le2+1)%3);
        index_t t2_adj3 = Tadj(t2,(le2+2)%3);
        Tset(t1,v1,v4,v3,t2_adj3,t1_adj2,t2);
        Tset(t2,v1,v2,v4,t2_adj2,t1,t1_adj3);

        // Tadj_replace(t1_adj2, t1, t1);
        Tadj_replace(t1_adj3, t1, t2);
        // Tadj_replace(t2_adj2, t2, t2);
        Tadj_replace(t2_adj3, t2, t1);
    }

    void CDT::insert_constraint(index_t i, index_t j) {
        std::cerr << "Insert constraint " << i << " " << j << std::endl;

        std::deque<index_t> Q; // Queue of edges to constrain
        vector<index_t> N; // New triangles to re-Delaunay-ize

        // Step 1: find all the edges that have an intersection 
        // with the constraint [i,j], enqueue them in Q.

        // TODO: -optimize predicate (do not need to test 4 orients
        //        per segment...)
        //       -handle corner-cases (constraint passes through existing vrtx)

        std::cerr << "Detect intersected edges" << std::endl;

        {
            index_t t_prev = index_t(-1);
            index_t t = vT(i);
            index_t le = index_t(-1);
            bool finished = false;
            for_each_T_around_v(
                i, [&](index_t t_in, index_t le_in) {
                    std::cerr << "Testing " << t_in << std::endl;
                    index_t v1 = Tv(t_in, (le_in + 1)%3);
                    index_t v2 = Tv(t_in, (le_in + 2)%3);
                    if(v1 == j || v2 == j) {
                        finished = true;
                        return true;
                    }
                    if(seg_seg_intersect(i,j,v1,v2)) {
                        std::cerr << "Found intersection : "
                                  << v1 << " " << v2
                                  << std::endl;
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
                std::cerr << "Testing " << t << std::endl;
                if(Tv(t,0) == j || Tv(t,1) == j || Tv(t,2) == j) {
                    std::cerr << "Finished" << std::endl;
                    finished = true;
                } else {
                    for(le = 0; le<3; ++le) {
                        if(Tadj(t,le) == t_prev) {
                            continue;
                        }
                        index_t v1 = Tv(t, (le + 1)%3);
                        index_t v2 = Tv(t, (le + 2)%3);
                        if(seg_seg_intersect(i,j,v1,v2)) {
                            std::cerr << "Found intersection : "
                                      << v1 << " " << v2
                                      << std::endl;
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
        save("CDT_0.geogram");
        index_t step = 1;
        while(Q.size() != 0) {
            std::cerr << "========== STEP " << step << std::endl;
            std::cerr << "Q= ";
            for(index_t q : Q) {
                std::cerr << q
                          << (Tis_marked(q) ? ' ' : '*')
                          << " ";
            }
            std::cerr << std::endl;
            index_t t1 = Q.back();
            Q.pop_back();
            std::cerr << "Examining " << t1 << "," << Tadj(t1,0) << std::endl;
            if(!Tis_marked(t1)) {
                std::cerr << "   Not marked"
                          << std::endl;
                continue;
            }
            if(!is_convex_quad(t1)) {
                std::cerr << "Not convex" << std::endl;
                if(Q.size() == 0) {
                    std::cerr << "... infinite iteration" << std::endl;
                    save("last_CDT.geogram");
                    abort();
                }
                Q.push_front(t1);
            } else {
                std::cerr << "Swap" << std::endl;                
                index_t t2 = Tadj(t1,0);

                bool no_isect    = !Tis_marked(t2);
                index_t v0 = Tv(t1,0);
                bool t2v0_t1v2 = (Tis_marked(t2) && Tv(t2,0) == Tv(t1,2));
                bool t2v0_t1v1 = (Tis_marked(t2) && Tv(t2,0) == Tv(t1,1));
                
                swap_edge(t1);
                if(no_isect) {
                    std::cerr << "   OK (no isect)" << std::endl;
                    Tunmark(t1);
                    N.push_back(t1);
                } else {
                    Sign o = orient2d(i,j,v0);
                    // geo_assert(o != ZERO);
                    if(t2v0_t1v2) {
                        Trot(t1,2);
                        if(o == POSITIVE) {
                            std::cerr << t1 << ":no isect  "
                                      << t2 << ":isect"
                                      << std::endl;
                            Tunmark(t1);
                            N.push_back(t1);
                        } else {
                            std::cerr << t1 << ":isect  "
                                      << t2 << ":isect"
                                      << std::endl;
                            Q.push_front(t1);
                        }
                    } else {
                        geo_assert(t2v0_t1v1);
                        Trot(t2,1);
                        if(o == POSITIVE) {
                            std::cerr << t1 << ":isect  "
                                      << t2 << ":isect"
                                      << std::endl;
                            Q.push_front(t1);                            
                        } else {
                            std::cerr << t1 << ":isect  "
                                      << t2 << ":no isect"
                                      << std::endl;
                            Tunmark(t2);
                            N.push_back(t2);
                            Q.push_front(t1);
                        }
                    }
                }
            }
            save("CDT_"+String::to_string(step)+".geogram");
            step += 1;
        }

        // Step 3: restore Delaunay 
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

        std::cerr << v1 << " " << v2 << " "
                  << v4 << " " << v3
                  << "  convex?="
                  << result
                  << std::endl;
        
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

