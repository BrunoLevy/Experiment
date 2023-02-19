
/*
 *
 */
 
#ifndef H__OGF_EXPERIMENT_CDT__H
#define H__OGF_EXPERIMENT_CDT__H

#include <OGF/Experiment/common/common.h>
#include <geogram/basic/geometry.h>
#include <geogram/numerics/predicates.h>

namespace GEO {

    // Under development, buggy !!
    class Experiment_API CDT {
    public:
    
    index_t insert(const vec2& P);
    void insert_constraint(index_t i, index_t j);

    index_t nT() const {
        return T_.size()/3;
    }

    index_t Tv(index_t t, index_t lv) const {
        geo_debug_assert(t<nT());
        geo_debug_assert(lv<3);
        return T_[3*t+lv];
    }
        
    index_t Tadj(index_t t, index_t le) const {
        geo_debug_assert(t<nT());
        geo_debug_assert(le<3);
        return Tadj_[3*t+le];
    }

    void save(const std::string& filename) const;
    
    protected:
    /**
     * \brief Swaps edge 0 of \p t1.
     *  Vertex 0 of \p t1 is vertex 0 of
     *  the two new triangles.
     */
    void swap_edge(index_t t1);

        
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

    void Tset(
        index_t t,
        index_t v1,   index_t v2,   index_t v3,
        index_t adj1, index_t adj2, index_t adj3
    ) {
        geo_assert(t < nT());
        T_[3*t  ]    = v1;
        T_[3*t+1]    = v2;
        T_[3*t+2]    = v3;                        
        Tadj_[3*t  ] = adj1;
        Tadj_[3*t+1] = adj2;
        Tadj_[3*t+2] = adj3;                        
    }

    void Tadj_set(index_t t, index_t le, index_t adj) {
        geo_assert(t < nT());
        geo_assert(adj < nT());
        geo_assert(le < 3);
        Tadj_[3*t+le] = adj;
    }

    void Tadj_replace(index_t t, index_t t1, index_t t2) {
        if(t == index_t(-1)) {
            return;
        }
        geo_assert(t < nT());
        geo_assert(t1 < nT());
        geo_assert(t2 < nT());
        index_t le = Tadj_find(t,t1);
        Tadj_set(t, le, t2);
    }
        
    index_t Tnew() {
        index_t t = nT();
        T_.resize((t+1)*3, index_t(-1));
        Tadj_.resize((t+1)*3, index_t(-1));
        return t;
    }

    index_t locate(index_t v, Sign& o1, Sign& o2, Sign& o3) const;

    Sign orient2d(index_t i, index_t j, index_t k) const;

    Sign incircle(index_t i, index_t j, index_t k, index_t l) const;
    
    protected:
    vector<vec2> point_;
    vector<index_t> T_;
    vector<index_t> Tadj_;
    };

}

#endif
