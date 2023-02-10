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

    inline rational_nt one_minus_t(const rational_nt& t) {
        // Next line equivalent to: return rational_nt(1.0)-t;
        return rational_nt(
            t.denom()-t.num(),
            t.denom() 
        );
    }
}

namespace GEO {

    vec3Q mix(const rational_nt& t, const vec3& p1, const vec3& p2) {
        rational_nt tt = one_minus_t(t);
        return vec3Q(
            tt*rational_nt(p1.x)+t*rational_nt(p2.x),
            tt*rational_nt(p1.y)+t*rational_nt(p2.y),
            tt*rational_nt(p1.z)+t*rational_nt(p2.z)
        );
    }

    vec2Q mix(const rational_nt& t, const vec2& p1, const vec2& p2) {
        rational_nt tt = one_minus_t(t);
        return vec2Q(
            tt*rational_nt(p1.x)+t*rational_nt(p2.x),
            tt*rational_nt(p1.y)+t*rational_nt(p2.y)
        );
    }

    
    vec2Q mix(const rational_nt& t, const vec2Q& p1, const vec2Q& p2) {
        rational_nt tt = one_minus_t(t);
        return vec2Q(
            tt*p1.x+t*p2.x,
            tt*p1.y+t*p2.y
        );
    }

    vec3Q mix(const rational_nt& t, const vec3Q& p1, const vec3Q& p2) {
        rational_nt tt = one_minus_t(t);
        return vec3Q(
            tt*p1.x+t*p2.x,
            tt*p1.y+t*p2.y,
            tt*p1.z+t*p2.z
        );
    }


    template<> expansion_nt det(const vec2E& v1, const vec2E& v2) {
        expansion* result = expansion::new_expansion_on_heap(
            expansion::det2x2_capacity(
                v1.x.rep(), v1.y.rep(),
                v2.x.rep(), v2.y.rep()
            )
        );
        result->assign_det2x2(
            v1.x.rep(), v1.y.rep(),
            v2.x.rep(), v2.y.rep()
        );
        return expansion_nt(result);
    }

    template<> expansion_nt dot(const vec2E& v1, const vec2E& v2) {
        const expansion& m1 = expansion_product(v1.x.rep(), v2.x.rep());
        const expansion& m2 = expansion_product(v1.y.rep(), v2.y.rep());
        return expansion_nt(expansion_nt::SUM, m1, m2);
    }

    template<> expansion_nt dot(const vec3E& v1, const vec3E& v2) {
        const expansion& m1 = expansion_product(v1.x.rep(), v2.x.rep());
        const expansion& m2 = expansion_product(v1.y.rep(), v2.y.rep());
        const expansion& m3 = expansion_product(v1.z.rep(), v2.z.rep());
        return expansion_nt(expansion_nt::SUM,m1,m2,m3);
    }

    template <> vec3E triangle_normal<vec3E>(
        const vec3& p1, const vec3& p2, const vec3& p3
    ) {
        const expansion& Ux = expansion_diff(p2.x,p1.x);
        const expansion& Uy = expansion_diff(p2.y,p1.y);
        const expansion& Uz = expansion_diff(p2.z,p1.z);
        const expansion& Vx = expansion_diff(p3.x,p1.x);
        const expansion& Vy = expansion_diff(p3.y,p1.y);
        const expansion& Vz = expansion_diff(p3.z,p1.z);
        expansion* Nx = expansion::new_expansion_on_heap(
            expansion::det2x2_capacity(Uy,Uz,Vy,Vz)
        );
        Nx->assign_det2x2(Uy,Uz,Vy,Vz);
        expansion* Ny = expansion::new_expansion_on_heap(
            expansion::det2x2_capacity(Uz,Ux,Vz,Vx)
        );
        Ny->assign_det2x2(Uz,Ux,Vz,Vx);
        expansion* Nz = expansion::new_expansion_on_heap(
            expansion::det2x2_capacity(Ux,Uy,Vx,Vy)
        );
        Nz->assign_det2x2(Ux,Uy,Vx,Vy);
        return vec3E(expansion_nt(Nx), expansion_nt(Ny), expansion_nt(Nz));
    }
    
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

