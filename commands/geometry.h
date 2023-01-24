#include <OGF/Experiment/common/common.h>
#include <geogram/basic/geometry.h>
#include <geogram/numerics/expansion_nt.h>

namespace GEO {

    typedef vecng<3,rational_nt> vec3Q;

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
        
        /**
         * \brief Computes the orientation predicate in 3d using
         *  vertices with rational coordinates.
         * \details Computes the sign of the signed volume of
         *  the tetrahedron p0, p1, p2, p3.
         * \param[in] p0 , p1 , p2 , p3 vertices of the tetrahedron
         * \retval POSITIVE if the tetrahedron is oriented positively
         * \retval ZERO if the tetrahedron is flat
         * \retval NEGATIVE if the tetrahedron is oriented negatively
         * \todo check whether orientation is inverted as compared to 
         *   Shewchuk's version.
         */
        Sign orient_3d(
            const vec3Q& p0, const vec3Q& p1, const vec3Q& p2, const vec3Q& p3
        ) {
	    rational_nt a11 = p1[0] - p0[0] ;
	    rational_nt a12 = p1[1] - p0[1] ;
	    rational_nt a13 = p1[2] - p0[2] ;
	    
	    rational_nt a21 = p2[0] - p0[0] ;
	    rational_nt a22 = p2[1] - p0[1] ;
	    rational_nt a23 = p2[2] - p0[2] ;
	    
	    rational_nt a31 = p3[0] - p0[0] ;
	    rational_nt a32 = p3[1] - p0[1] ;
	    rational_nt a33 = p3[2] - p0[2] ;
	    
	    rational_nt Delta = det3x3(
		a11,a12,a13,
		a21,a22,a23,
		a31,a32,a33
	    );

	    return Delta.sign();
        }


        Sign orient_2d_projected_impl(
            const vec3Q& p0, const vec3Q& p1, const vec3Q& p2,
            coord_index_t axis
        ) {
        
            coord_index_t u = coord_index_t((axis + 1)%3);
            coord_index_t v = coord_index_t((axis + 2)%3);

            rational_nt a11 = p1[u] - p0[u];
            rational_nt a12 = p1[v] - p0[v];
            
            rational_nt a21 = p2[u] - p0[u] ;
            rational_nt a22 = p2[v] - p0[v] ;
            
            rational_nt Delta = det2x2(
                a11,a12,
                a21,a22
            );

            return Delta.sign();
        }


        Sign orient_2d_projected(
            const vec3Q& p0, const vec3Q& p1, const vec3Q& p2,
            coord_index_t axis
        ) {
            Sign o1 = orient_2d_projected_impl(p0,p1,p2,axis);
            
#ifdef GEO_DEBUG
            Sign o2 = orient_2d_projected_impl(p1,p2,p0,axis);
            Sign o3 = orient_2d_projected_impl(p2,p0,p1,axis);
            Sign o4 = orient_2d_projected_impl(p2,p1,p0,axis);
            Sign o5 = orient_2d_projected_impl(p0,p2,p1,axis);
            Sign o6 = orient_2d_projected_impl(p1,p0,p2,axis);
            geo_debug_assert(o1 == o2 && o1 == o3);
            geo_debug_assert(o4 == o5 && o4 == o6);
            geo_debug_assert(int(o1) == -int(o4));

#endif
            return o1;            
        }
        
        Sign dot_3d(const vec3Q& p0, const vec3Q& p1, const vec3Q& p2) {
            return dot(p1-p0,p2-p0).sign();
        }
        
    }
    
    /**
     * \brief Converts a 3d vector with double coordinates
     *  into a 3d vector with coordinates of type T.
     * \param[in] p the vector to be converted
     * \return the converted vector, with coordinates of type
     *  T
     * \tparam T the type of the coordinates
     */
    template <class T>
    inline vecng<3,T> convert_vec3_generic(const vec3& p) {
        return vecng<3,T>(
            T(p.x), T(p.y), T(p.z)
        );
    }

    /**
     * \brief Creates a vector with coordinates of type T
     *  from two points with double coordinates
     * \param[in] p1 , p2 the two vectors
     * \return The vector \p p2 - \p p1 
     *  with coordinates of type T
     * \tparam T the type of the coordinates
     */
    template <class T>
    inline vecng<3,T> make_vec3_generic(
        const vec3& p1, const vec3& p2
    ) {
        return vecng<3,T>(
            T(p2.x) - T(p1.x),
            T(p2.y) - T(p1.y),
            T(p2.z) - T(p1.z)
        );
    }
    
    /**
     * \brief Computes the intersection between a segment and a triangle
     * \pre The intersection exists
     * \tparam T type of the vector coordinate (double or rational_nt)
     * \param[in) q1 , q2 the two vertices of the segment
     * \param[in] p1 , p2 , p3 the three vertices of the triangle
     * \return the intersection between the segment and the triangle
     */
    template <class T>
    inline vecng<3,T> get_segment_triangle_intersection(
        const vec3& q1, const vec3& q2,
        const vec3& p1, const vec3& p2, const vec3& p3
    ) {
        vecng<3,T> D  = make_vec3_generic<T>(q1,q2);
        vecng<3,T> E1 = make_vec3_generic<T>(p2,p3);
        vecng<3,T> E2 = make_vec3_generic<T>(p3,p1);
        vecng<3,T> AO = make_vec3_generic<T>(p1,q1);
        vecng<3,T> N  = cross(E1,E2);
        T d =  dot(D,N);
        geo_assert(geo_sgn(d) != ZERO);
        T t = -dot(AO,N)/d;
        return
            t          * convert_vec3_generic<T>(q2) +
            (T(1.0)-t) * convert_vec3_generic<T>(q1) ;
    }

    /**
     * \brief Computes the intersection between two 3D segments
     *  based on their projection onto one of the axes
     * \pre the intersection exists in 3D
     * \tparam T the type of the coordinates
     * \param[in] p1 , p2 the extremities of the first segment
     * \param[in] q1 , q2 the extremities of the second segment
     * \param[in] ax one of 0,1,2, the axis along which to project
     * \return the 3D intersection of the two segments
     */
    template <class T>    
    inline vecng<3,T> get_segment_segment_intersection_2D(
        const vec3& p1, const vec3& p2,
        const vec3& q1, const vec3& q2, coord_index_t ax
    ) {
        coord_index_t u = coord_index_t((ax + 1)%3);
        coord_index_t v = coord_index_t((ax + 2)%3);

        // [ a b ] [ l1 ]   [ e ]
        // [ c d ] [ l2 ] = [ f ]
        
        T a = T(p2[u])-T(p1[u]);
        T b = T(q1[u])-T(q2[u]);
        T c = T(p2[v])-T(p1[v]);
        T d = T(q1[v])-T(q2[v]);

        T e = T(q1[u])-T(p1[u]);
        T f = T(q1[v])-T(p1[v]);

        // [ a b ] [ s ]   [ e ]
        // [ c d ] [ t ] = [ f ]

        // [ s ]               [ d -b] [ e ]
        // [ t ] = 1/(ad-bc) * [-c  a] [ f ]

        T det = a*d-b*c;

        geo_assert(geo_sgn(det) != ZERO);
        
        T s = ( d*e-b*f)/det;
        
        return
            s          * convert_vec3_generic<T>(p2) +
            (T(1.0)-s) * convert_vec3_generic<T>(p1) ;
    }

}
