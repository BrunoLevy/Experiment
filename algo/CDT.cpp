/**
 *
 */

#include <OGF/Experiment/algo/CDT.h>
#include <geogram/mesh/mesh.h>
#include <geogram/mesh/mesh_io.h>
#include <stack>

namespace GEO {

    index_t CDT::insert(const vec2& P) {
        index_t v = point_.size();
        point_.push_back(P);

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
        // TODO
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
        return PCK::orient_2d(
            point_[i], point_[j], point_[k]
        );
    }

    Sign CDT::incircle(index_t i, index_t j, index_t k, index_t l) const {
        return PCK::in_circle_2d_SOS(
            point_[i].data(), point_[j].data(), point_[k].data(),
            point_[l].data()
        );
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
        mesh_save(M, filename);
    }
}

