#ifndef H__OGF_EXPERIMENT_COMMANDS_EXACT_GEOMETRY__H
#define H__OGF_EXPERIMENT_COMMANDS_EXACT_GEOMETRY__H

#include <OGF/Experiment/common/common.h>
#include <geogram/basic/geometry.h>
#include <geogram/numerics/expansion_nt.h>

namespace GEO {

    typedef vecng<2,expansion_nt> vec2E;    
    typedef vecng<2,rational_nt>  vec2Q;
    
    typedef vecng<3,expansion_nt> vec3E;    
    typedef vecng<3,rational_nt>  vec3Q;

    inline vec3Q mix(rational_nt t, const vec3& p1, const vec3& p2) {
        rational_nt tt = rational_nt(1.0)-t;
        return vec3Q(
            tt*rational_nt(p1.x)+t*rational_nt(p2.x),
            tt*rational_nt(p1.y)+t*rational_nt(p2.y),
            tt*rational_nt(p1.z)+t*rational_nt(p2.z)
        );
    }

    inline vec2Q mix(rational_nt t, const vec2& p1, const vec2& p2) {
        rational_nt tt = rational_nt(1.0)-t;
        return vec2Q(
            tt*rational_nt(p1.x)+t*rational_nt(p2.x),
            tt*rational_nt(p1.y)+t*rational_nt(p2.y)
        );
    }

    inline vec2Q mix(rational_nt t, const vec2Q& p1, const vec2Q& p2) {
        rational_nt tt = rational_nt(1.0)-t;
        return vec2Q(
            tt*p1.x+t*p2.x,
            tt*p1.y+t*p2.y
        );
    }

    inline vec3Q mix(rational_nt t, const vec3Q& p1, const vec3Q& p2) {
        rational_nt tt = rational_nt(1.0)-t;
        return vec3Q(
            tt*p1.x+t*p2.x,
            tt*p1.y+t*p2.y,
            tt*p1.z+t*p2.z
        );
    }

    class vec3QLexicoCompare {
    public:
        bool operator()(const vec3Q& v1, const vec3Q& v2) const {
            Sign s = v2.x.compare(v1.x);
            if(s == POSITIVE) {
                return true;
            }
            if(s == NEGATIVE) {
                return false;
            }
            s = v2.y.compare(v1.y);
            if(s == POSITIVE) {
                return true;
            }
            if(s == NEGATIVE) {
                return false;
            }
            s = v2.z.compare(v1.z);
            return (s == POSITIVE);
        }
    };

    namespace PCK {

        template <class T> inline bool same_point(
            const vecng<3,T>& v1, const vecng<3,T>& v2
        ) {
            return (v1.x == v2.x) && (v1.y == v2.y) && (v2.z == v1.z);
        }

        template <class T> inline bool same_point(
            const vecng<2,T>& v1, const vecng<2,T>& v2
        ) {
            return (v1.x == v2.x) && (v1.y == v2.y);
        }

        inline Sign orient_2d(
            const vec2Q& p0, const vec2Q& p1, const vec2Q& p2
        ) {
            rational_nt a11 = p1.x - p0.x;
            rational_nt a12 = p1.y - p0.y;
            rational_nt a21 = p2.x - p0.x;
            rational_nt a22 = p2.y - p0.y;
            rational_nt Delta = det2x2(
                a11,a12,
                a21,a22
            );
            return Delta.sign();
        }
        
        inline Sign dot_2d(const vec2Q& p0, const vec2Q& p1, const vec2Q& p2) {
            return dot(p1-p0,p2-p0).sign();
        }
        
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

    inline vec2Q u_P1P2_plus_v_P1P3(
        rational_nt u, rational_nt v, const vec2& p1, const vec2& p2, const vec2& p3
    ) {
        vec2E E1 = make_vec2<vec2E>(p1,p2);
        vec2E E2 = make_vec2<vec2E>(p1,p3);		      
        vec2E Pnum = (u.num()*v.denom())*E1 + (v.num()*u.denom())*E2;
        expansion_nt Pden = u.denom()*v.denom();
        return vec2Q(
            rational_nt(Pnum.x, Pden),
            rational_nt(Pnum.y, Pden)
        );
    }
    
    bool Experiment_API get_three_planes_intersection(
        vec3Q& result,
        const vec3& p1, const vec3& p2, const vec3& p3,
        const vec3& q1, const vec3& q2, const vec3& q3,
        const vec3& r1, const vec3& r2, const vec3& r3
    );

    /*********************************************************************************/

    /**
     * \brief 2D vector in homogeneous coordinates
     *  with coordinates as arithmetic expansions.
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
        
        expansion_nt x;
        expansion_nt y;
        expansion_nt w;
    };

    /**
     * \brief 3D vector in homogeneous coordinates
     *  with coordinates as arithmetic expansions.
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
        
        expansion_nt x;
        expansion_nt y;
        expansion_nt z;        
        expansion_nt w;
    };
    
    class Experiment_API vec3HELexicoCompare {
    public:
       bool operator()(const vec3HE& v1, const vec3HE& v2) const; 
    };
    
}

#endif
