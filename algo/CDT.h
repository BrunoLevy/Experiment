
/*
 *
 */
 
#ifndef H__OGF_EXPERIMENT_CDT__H
#define H__OGF_EXPERIMENT_CDT__H

#include <OGF/Experiment/common/common.h>
#include <geogram/basic/geometry.h>
#include <geogram/numerics/predicates.h>

namespace GEO {

    // Under development, unfinished !!
    class Experiment_API CDT {
    public:
    
    index_t insert(const vec2& P);
    
    void insert_constraint(index_t i, index_t j);

    /**
     * \brief Gets the number of triangles
     */
    index_t nT() const {
        return T_.size()/3;
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
        geo_assert(t < nT());
        T_[3*t  ]    = v1;
        T_[3*t+1]    = v2;
        T_[3*t+2]    = v3;                        
        Tadj_[3*t  ] = adj1;
        Tadj_[3*t+1] = adj2;
        Tadj_[3*t+2] = adj3;                        
    }

    /**
     * \brief Sets a triangle adjacency relation
     * \param[in] t a triangle
     * \param[in] le local edge index, in 0,1,2
     * \param[in] adj the triangle adjacent to \p t 
     *  accross \p le
     */
    void Tadj_set(index_t t, index_t le, index_t adj) {
        geo_assert(t < nT());
        geo_assert(adj < nT());
        geo_assert(le < 3);
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
        geo_assert(t < nT());
        geo_assert(t1 < nT());
        geo_assert(t2 < nT());
        index_t le = Tadj_find(t,t1);
        Tadj_set(t, le, t2);
    }

    /**
     * \brief Creates a new triangle
     * \return the index of the new triange
     */
    index_t Tnew() {
        index_t t = nT();
        T_.resize((t+1)*3, index_t(-1));
        Tadj_.resize((t+1)*3, index_t(-1));
        return t;
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
    
    protected:
    vector<vec2> point_;
    vector<index_t> T_;
    vector<index_t> Tadj_;
    };

}

#endif
