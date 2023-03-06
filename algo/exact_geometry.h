#ifndef H__OGF_EXPERIMENT_COMMANDS_EXACT_GEOMETRY__H
#define H__OGF_EXPERIMENT_COMMANDS_EXACT_GEOMETRY__H

#include <OGF/Experiment/common/common.h>
#include <geogram/basic/geometry.h>
#include <geogram/numerics/expansion_nt.h>

namespace GEO {


    /**
     * \brief vec2 with coordinates as expansions
     * \details Coordinates support +,-,*
     */
    typedef vecng<2,expansion_nt> vec2E;

    /**
     * \brief vec3 with coordinates as expansions
     * \details Coordinates support +,-,*
     */
    typedef vecng<3,expansion_nt> vec3E;    

    /**
     * \brief 2D vector in homogeneous coordinates
     *  with coordinates as arithmetic expansions
     * \details Coordinates support +,-,* and / by
     *  multiplying w.
     */
    struct vec2HE {

        vec2HE() :
            x(expansion_nt::UNINITIALIZED),
            y(expansion_nt::UNINITIALIZED),
            w(expansion_nt::UNINITIALIZED)            
        {
        }

        vec2HE(
            const expansion_nt& x_in,
            const expansion_nt& y_in,
            const expansion_nt& w_in
        ) : x(x_in), y(y_in), w(w_in) {
        }

        vec2HE(
            expansion_nt&& x_in,
            expansion_nt&& y_in,
            expansion_nt&& w_in
        ) : x(x_in), y(y_in), w(w_in) {
        }
        
        vec2HE(const vec2HE& rhs) :
            x(rhs.x), y(rhs.y), w(rhs.w) {
        }

        vec2HE(vec2HE&& rhs) :
            x(rhs.x), y(rhs.y), w(rhs.w) {
        }

        vec2HE& operator=(const vec2HE& rhs) {
            if(&rhs != this) {
                x=rhs.x;
                y=rhs.y;
                w=rhs.w;
            }
            return *this;
        }

        vec2HE& operator=(vec2HE&& rhs) {
            if(&rhs != this) {
                x=rhs.x;
                y=rhs.y;
                w=rhs.w;
            }
            return *this;
        }

        expansion_nt* data() {
            return &x;
        }

        const expansion_nt* data() const {
            return &x;
        }
        
        expansion_nt& operator[](coord_index_t i) {
            geo_debug_assert(i < 2);
            return data()[i];
        }

        const expansion_nt& operator[](coord_index_t i) const {
            geo_debug_assert(i < 2);
            return data()[i];
        }

        void optimize() {
            x.optimize();
            y.optimize();
            w.optimize();
        }
        
        expansion_nt x;
        expansion_nt y;
        expansion_nt w;
    };

    /**
     * \brief 3D vector in homogeneous coordinates
     *  with coordinates as arithmetic expansions.
     * \details Coordinates support +,-,* and / by
     *  multiplying w.
     */
    struct vec3HE {
        vec3HE() :
            x(expansion_nt::UNINITIALIZED),
            y(expansion_nt::UNINITIALIZED),
            z(expansion_nt::UNINITIALIZED),
            w(expansion_nt::UNINITIALIZED)            
        {
        }

        vec3HE(
            const expansion_nt& x_in,
            const expansion_nt& y_in,
            const expansion_nt& z_in,
            const expansion_nt& w_in            
        ) : x(x_in), y(y_in), z(z_in), w(w_in) {
        }

        vec3HE(
            expansion_nt&& x_in,
            expansion_nt&& y_in,
            expansion_nt&& z_in,
            expansion_nt&& w_in
        ) : x(x_in), y(y_in), z(z_in), w(w_in) {
        }
        
        vec3HE(const vec3HE& rhs) :
            x(rhs.x), y(rhs.y), z(rhs.z), w(rhs.w) {
        }

        vec3HE(vec3HE&& rhs) :
            x(rhs.x), y(rhs.y), z(rhs.z), w(rhs.w) {
        }

        vec3HE& operator=(const vec3HE& rhs) {
            if(&rhs != this) {
                x=rhs.x;
                y=rhs.y;
                z=rhs.z;                
                w=rhs.w;
            }
            return *this;
        }

        vec3HE& operator=(vec3HE&& rhs) {
            if(&rhs != this) {
                x=rhs.x;
                y=rhs.y;
                z=rhs.z;                
                w=rhs.w;
            }
            return *this;
        }

        expansion_nt* data() {
            return &x;
        }

        const expansion_nt* data() const {
            return &x;
        }
        
        expansion_nt& operator[](coord_index_t i) {
            geo_debug_assert(i < 3);
            return data()[i];
        }

        const expansion_nt& operator[](coord_index_t i) const {
            geo_debug_assert(i < 3);
            return data()[i];
        }

        void optimize() {
            x.optimize();
            y.optimize();
            z.optimize();            
            w.optimize();
        }
        
        expansion_nt x;
        expansion_nt y;
        expansion_nt z;        
        expansion_nt w;
    };

    inline vec3HE vec3HE_noinit() {
        return vec3HE(
            expansion_nt(expansion_nt::UNINITIALIZED),
            expansion_nt(expansion_nt::UNINITIALIZED),
            expansion_nt(expansion_nt::UNINITIALIZED),
            expansion_nt(expansion_nt::UNINITIALIZED)                       
        );
    }

    inline vec2HE vec2HE_noinit() {
        return vec2HE(
            expansion_nt(expansion_nt::UNINITIALIZED),
            expansion_nt(expansion_nt::UNINITIALIZED),
            expansion_nt(expansion_nt::UNINITIALIZED)                        
        );
    }
    
    vec2HE Experiment_API operator-(const vec2HE& p1, const vec2HE& p2);
    
    vec3HE Experiment_API operator-(const vec3HE& p1, const vec3HE& p2);

    /**
     * \brief Comparator class for vec3HE
     * \detail Used to create maps indexed by vec3HE
     */
    class Experiment_API vec2HELexicoCompare {
    public:
       /**
        * \brief Compares two vec3HE
        * \retval true if \p v1 is before \p v2 in the lexicographic
        *  order
        * \retval false otherwise
        */
       bool operator()(const vec2HE& v1, const vec2HE& v2) const; 
    };
    
    /**
     * \brief Comparator class for vec3HE
     * \detail Used to create maps indexed by vec3HE
     */
    class Experiment_API vec3HELexicoCompare {
    public:
       /**
        * \brief Compares two vec3HE
        * \retval true if \p v1 is before \p v2 in the lexicographic
        *  order
        * \retval false otherwise
        */
       bool operator()(const vec3HE& v1, const vec3HE& v2) const; 
    };

    vec3HE Experiment_API mix(
        const rational_nt& t, const vec3& p1, const vec3& p2
    );

    vec2HE Experiment_API mix(
        const rational_nt& t, const vec2& p1, const vec2& p2
    );

    vec2HE Experiment_API mix(
        const rational_nt& t, const vec2HE& p1, const vec2HE& p2
    );

    vec3HE Experiment_API mix(
        const rational_nt& t, const vec3HE& p1, const vec3HE& p2
    );

    /**
     * \brief Specialization optimized using low-level API
     */
    template<> expansion_nt Experiment_API det(const vec2E& v1, const vec2E& v2);

    /**
     * \brief Specialization optimized using low-level API
     */
    template<> expansion_nt Experiment_API dot(const vec2E& v1, const vec2E& v2);

    /**
     * \brief Specialization optimized using low-level API
     */
    template<> expansion_nt Experiment_API dot(const vec3E& v1, const vec3E& v2);
    
    namespace PCK {

        template <class T> inline bool same_point(
            const vecng<3,T>& v1, const vecng<3,T>& v2
        ) {
            // operator== is well optimized for expansion_nt and rational_nt
            return (v1.x == v2.x) && (v1.y == v2.y) && (v2.z == v1.z);
        }

        template <class T> inline bool same_point(
            const vecng<2,T>& v1, const vecng<2,T>& v2
        ) {
            // operator== is well optimized for expansion_nt and rational_nt
            return (v1.x == v2.x) && (v1.y == v2.y);
        }

        bool Experiment_API same_point(const vec2HE& v1, const vec2HE& v2);
        bool Experiment_API same_point(const vec3HE& v1, const vec3HE& v2);  
        
        Sign Experiment_API orient_2d(
            const vec2HE& p0, const vec2HE& p1, const vec2HE& p2
        );
        
        Sign Experiment_API dot_2d(
            const vec2HE& p0, const vec2HE& p1, const vec2HE& p2
        );

        Sign orient_2dlifted_SOS(
            const vec2HE& p0, const vec2HE& p1,
            const vec2HE& p2, const vec2HE& p3,
            double h0, double h1, double h2, double h3
        );
    }

    /**
     * \brief Converts a 3d vector with double coordinates
     *  into a 3d vector with coordinates of arbitrary type
     * \param[in] p the vector to be converted
     * \return the converted vector
     * \tparam VEC3 the type of the returned vector
     */
    template <class VEC3>
    inline VEC3 make_vec3(const vec3& p) {
        typedef typename VEC3::value_type value_type;
        return VEC3(value_type(p.x),value_type(p.y),value_type(p.z));
    }

    /**
     * \brief Creates a vector with coordinates of arbitrary type
     *  from two points with double coordinates
     * \param[in] p1 , p2 the two vectors
     * \return The vector \p p2 - \p p1 
     * \tparam VEC3 the type of the returned vector
     */
    template <class VEC3 = vec3>
    inline VEC3 make_vec3(const vec3& p1, const vec3& p2) {
        typedef typename VEC3::value_type value_type;
        return VEC3(
            value_type(p2.x) - value_type(p1.x),
            value_type(p2.y) - value_type(p1.y),
            value_type(p2.z) - value_type(p1.z)
        );
    }

    /**
     * \brief Specialization for vec3E
     */
    template <>
    inline vec3E make_vec3<vec3E>(const vec3& p1, const vec3& p2) {
        return vec3E(
            expansion_nt(expansion_nt::DIFF, p2.x, p1.x),
            expansion_nt(expansion_nt::DIFF, p2.y, p1.y),
            expansion_nt(expansion_nt::DIFF, p2.z, p1.z)            
        );
    }
    
    /**
     * \brief Creates a vector with coordinates of arbitrary type
     *  from two points with double coordinates
     * \param[in] p1 , p2 the two vectors
     * \return The vector \p p2 - \p p1 
     * \tparam VEC2 the type of the returned vector
     */
    template <class VEC2>
    inline VEC2 make_vec2(
        const vec2& p1, const vec2& p2
    ) {
        typedef typename VEC2::value_type value_type;        
        return VEC2(
            value_type(p2.x) - value_type(p1.x),
            value_type(p2.y) - value_type(p1.y)
        );
    }

    /**
     * \brief Specialization for vec2E
     */
    template <>
    inline vec2E make_vec2<vec2E>(const vec2& p1, const vec2& p2) {
        return vec2E(
            expansion_nt(expansion_nt::DIFF, p2.x, p1.x),
            expansion_nt(expansion_nt::DIFF, p2.y, p1.y)
        );
    }
    
    /**
     * \brief Computes the normal to a triangle from its three
     *  vertices 
     * \param[in] p1 , p2 , p3 the three vertices of the triangle
     * \return the normal to the triangle with coordinates of 
     *  arbitrary type
     * \tparam VEC3 the type of the returned vector
     */
    template <class VEC3>
    inline VEC3 triangle_normal(
        const vec3& p1, const vec3& p2, const vec3& p3
    ) {
        return cross(
            make_vec3<VEC3>(p1,p2),
            make_vec3<VEC3>(p1,p3)
        );
    }

    /**
     * \brief Specialization for vec3E
     */
    template <> Experiment_API vec3E triangle_normal<vec3E>(
        const vec3& p1, const vec3& p2, const vec3& p3
    );

    /**
     * \brief Interpolates a point linearly in a triangle
     * \return \p p1 + \p u * (\p p2 - \p p1) + \p v * (\p p3 - \p p1) 
     */
    inline vec2HE u_P1P2_plus_v_P1P3(
        rational_nt u, rational_nt v,
        const vec2& p1, const vec2& p2, const vec2& p3
    ) {
        vec2E E1 = make_vec2<vec2E>(p1,p2);
        vec2E E2 = make_vec2<vec2E>(p1,p3);		      
        vec2E Pnum = (u.num()*v.denom())*E1 + (v.num()*u.denom())*E2;
        expansion_nt Pden = u.denom()*v.denom();
        return vec2HE(
            Pnum.x,
            Pnum.y,
            Pden
        );
    }
    
    bool Experiment_API get_three_planes_intersection(
        vec3HE& result,
        const vec3& p1, const vec3& p2, const vec3& p3,
        const vec3& q1, const vec3& q2, const vec3& q3,
        const vec3& r1, const vec3& r2, const vec3& r3
    );

    /************************************************************************/

    
}

#endif
