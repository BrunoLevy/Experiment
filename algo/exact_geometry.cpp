/*
 *
 */

#include <OGF/Experiment/algo/exact_geometry.h>

namespace {
    using namespace GEO;
    /**
     * \brief Compares two rational numbers given as separate
     *   numerators and denominators.
     * \param[in] a_num , a_denom defines a = \p a_num / \p a_denom
     * \param[in] b_num , b_denom defines b = \p b_num / \p b_denom
     * \return the sign of a - b
     */
    inline Sign ratio_compare(
        const expansion_nt& a_num,
        const expansion_nt& a_denom,
        const expansion_nt& b_num,
        const expansion_nt& b_denom
    ) {
	if(a_denom == b_denom) {
	    const expansion& diff_num = expansion_diff(
		a_num.rep(), b_num.rep()
	    );
	    return Sign(diff_num.sign() * a_denom.sign());
	}
	const expansion& num_a = expansion_product(
	    a_num.rep(), b_denom.rep()
	);
	const expansion& num_b = expansion_product(
	    b_num.rep(), a_denom.rep()
	);
	const expansion& diff_num = expansion_diff(num_a, num_b);
	return Sign(
	    diff_num.sign() * a_denom.sign() * b_denom.sign()
	);
    }
}

namespace GEO {
    namespace PCK {
        Sign orient_2d(
            const vec2Q& p0, const vec2Q& p1, const vec2Q& p2
        ) {

            // Special case: per-point denominators are the same
            // (homogeneous coordinates encoded in rational_nt)
            if(
                p0.x.denom() == p0.y.denom() &&
                p1.x.denom() == p1.y.denom()
            ) {
                expansion& Delta = expansion_det3x3(
                    p0.x.num().rep(), p0.y.num().rep(), p0.x.denom().rep(),
                    p1.x.num().rep(), p1.y.num().rep(), p1.x.denom().rep(),
                    p2.x.num().rep(), p2.y.num().rep(), p2.x.denom().rep()
                );
                return Sign(
                    Delta.sign()*
                    p0.x.denom().rep().sign()*
                    p1.x.denom().rep().sign()*
                    p2.x.denom().rep().sign()
                );
            }

            // General case, using (slower) rational_nt type.
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
        
        Sign dot_2d(const vec2Q& p0, const vec2Q& p1, const vec2Q& p2) {
            return dot(p1-p0,p2-p0).sign();
        }
        
    }

    bool get_three_planes_intersection(
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

    bool vec3HELexicoCompare::operator()(const vec3HE& v1, const vec3HE& v2) const {
        Sign s = ratio_compare(v2.x, v2.w, v1.x, v1.w);
        if(s == POSITIVE) {
            return true;
        }
        if(s == NEGATIVE) {
            return false;
        }

        s = ratio_compare(v2.y, v2.w, v1.y, v1.w);
        if(s == POSITIVE) {
            return true;
        }
        if(s == NEGATIVE) {
            return false;
        }
        
        s = ratio_compare(v2.z, v2.w, v1.z, v1.w);
        return (s == POSITIVE);
    }
            
}

