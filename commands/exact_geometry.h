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
            Sign s = (v2.x - v1.x).sign();
            if(s == POSITIVE) {
                return true;
            }
            if(s == NEGATIVE) {
                return false;
            }
            s = (v2.y - v1.y).sign();
            if(s == POSITIVE) {
                return true;
            }
            if(s == NEGATIVE) {
                return false;
            }
            s = (v2.z - v1.z).sign();
            return s == POSITIVE;
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
    
    inline bool get_three_planes_intersection(
        vec3Q& result,
        const vec3& p1, const vec3& p2, const vec3& p3,
        const vec3& q1, const vec3& q2, const vec3& q3,
        const vec3& r1, const vec3& r2, const vec3& r3
    ) {

        vec3E N1 = triangle_normal<vec3E>(p1,p2,p3);
        vec3E N2 = triangle_normal<vec3E>(q1,q2,q3);
        vec3E N3 = triangle_normal<vec3E>(r1,r2,r3);

        vec3E B(
            dot(N1,make_vec3<vec3E>(p1)),
            dot(N2,make_vec3<vec3E>(q1)),
            dot(N3,make_vec3<vec3E>(r1))
        );
        
        expansion_nt Delta = det3x3(
            N1.x, N1.y, N1.z,
            N2.x, N2.y, N2.z,
            N3.x, N3.y, N3.z
        );

        if(Delta.sign() == ZERO) {
            return false;
        }
        
        expansion_nt X = det3x3(
            B.x, N1.y, N1.z,
            B.y, N2.y, N2.z,
            B.z, N3.y, N3.z
        );

        expansion_nt Y = det3x3(
            N1.x, B.x, N1.z,
            N2.x, B.y, N2.z,
            N3.x, B.z, N3.z
        );

        expansion_nt Z = det3x3(
            N1.x, N1.y, B.x,
            N2.x, N2.y, B.y,
            N3.x, N3.y, B.z
        );

        result.x=rational_nt(X,Delta);
        result.y=rational_nt(Y,Delta);
        result.z=rational_nt(Z,Delta);

        return true;
    }
}

#endif
