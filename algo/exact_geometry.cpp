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
    
    vec2HE operator-(const vec2HE& p1, const vec2HE& p2) {
        if(p1.w == p2.w) {
            return vec2HE(
                expansion_nt(expansion_nt::DIFF, p2.x.rep(), p1.x.rep()),
                expansion_nt(expansion_nt::DIFF, p2.y.rep(), p1.y.rep()),
                p1.w
            );
        }
        expansion& x1 = expansion_product(p1.x.rep(), p2.w.rep());
        expansion& y1 = expansion_product(p1.y.rep(), p2.w.rep());        
        expansion& x2 = expansion_product(p2.x.rep(), p1.w.rep());
        expansion& y2 = expansion_product(p2.y.rep(), p1.w.rep());
        return vec2HE(
            expansion_nt(expansion_nt::DIFF, x2, x1),
            expansion_nt(expansion_nt::DIFF, y2, y1),
            expansion_nt(expansion_nt::PRODUCT, p1.w.rep(), p2.w.rep())
        );
    }
    
    vec3HE operator-(const vec3HE& p1, const vec3HE& p2) {
        if(p1.w == p2.w) {
            return vec3HE(
                expansion_nt(expansion_nt::DIFF, p2.x.rep(), p1.x.rep()),
                expansion_nt(expansion_nt::DIFF, p2.y.rep(), p1.y.rep()),
                expansion_nt(expansion_nt::DIFF, p2.z.rep(), p1.z.rep()),
                p1.w
            );
        }
        expansion& x1 = expansion_product(p1.x.rep(), p2.w.rep());
        expansion& y1 = expansion_product(p1.y.rep(), p2.w.rep());
        expansion& z1 = expansion_product(p1.z.rep(), p2.w.rep());   
        expansion& x2 = expansion_product(p2.x.rep(), p1.w.rep());
        expansion& y2 = expansion_product(p2.y.rep(), p1.w.rep());
        expansion& z2 = expansion_product(p2.z.rep(), p1.w.rep());
        return vec3HE(
            expansion_nt(expansion_nt::DIFF, x2, x1),
            expansion_nt(expansion_nt::DIFF, y2, y1),
            expansion_nt(expansion_nt::DIFF, z2, z1),            
            expansion_nt(expansion_nt::PRODUCT, p1.w.rep(), p2.w.rep())
        );
    }

    bool vec2HELexicoCompare::operator()(
        const vec2HE& v1, const vec2HE& v2
    ) const {
        Sign s = ratio_compare(v2.x, v2.w, v1.x, v1.w);
        if(s == POSITIVE) {
            return true;
        }
        if(s == NEGATIVE) {
            return false;
        }
        s = ratio_compare(v2.y, v2.w, v1.y, v1.w);
        return (s == POSITIVE);
    }
    
    bool vec3HELexicoCompare::operator()(
        const vec3HE& v1, const vec3HE& v2
    ) const {
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
            

    vec3HE mix(const rational_nt& t, const vec3& p1, const vec3& p2) {
        expansion& st_d = const_cast<expansion&>(t.denom().rep());
        st_d.optimize();
        expansion& t_n  = const_cast<expansion&>(t.num().rep());
        t_n.optimize();
        const expansion& s_n  = expansion_diff(st_d, t_n);
        const expansion& sx   = expansion_product(s_n, p1.x);
        const expansion& tx   = expansion_product(t_n, p2.x);
        const expansion& sy   = expansion_product(s_n, p1.y);
        const expansion& ty   = expansion_product(t_n, p2.y);
        const expansion& sz   = expansion_product(s_n, p1.z);
        const expansion& tz   = expansion_product(t_n, p2.z);
        return vec3HE(
            expansion_nt(expansion_nt::SUM, sx,tx),
            expansion_nt(expansion_nt::SUM, sy,ty),
            expansion_nt(expansion_nt::SUM, sz,tz),
            expansion_nt(st_d)
        );
    }

    vec2HE mix(const rational_nt& t, const vec2& p1, const vec2& p2) {
        expansion& st_d = const_cast<expansion&>(t.denom().rep());
        st_d.optimize();
        expansion& t_n  = const_cast<expansion&>(t.num().rep());
        t_n.optimize();
        const expansion& s_n  = expansion_diff(st_d, t_n);
        const expansion& sx   = expansion_product(s_n, p1.x);
        const expansion& tx   = expansion_product(t_n, p2.x);
        const expansion& sy   = expansion_product(s_n, p1.y);
        const expansion& ty   = expansion_product(t_n, p2.y);
        return vec2HE(
            expansion_nt(expansion_nt::SUM, sx,tx),
            expansion_nt(expansion_nt::SUM, sy,ty),
            expansion_nt(st_d)
        );
    }
    
    vec2HE mix(const rational_nt& t, const vec2HE& p1, const vec2HE& p2) {
        expansion& st_d = expansion_product(t.denom().rep(),p1.w.rep());
        st_d.optimize();
        expansion& t_n  = const_cast<expansion&>(t.num().rep());
        t_n.optimize();
        const expansion& s_n  = expansion_diff(st_d, t_n);
        if(p1.w == p2.w) {
            const expansion& sx   = expansion_product(s_n, p1.x.rep());
            const expansion& tx   = expansion_product(t_n, p2.x.rep());
            const expansion& sy   = expansion_product(s_n, p1.y.rep());
            const expansion& ty   = expansion_product(t_n, p2.y.rep());
            return vec2HE(
                expansion_nt(expansion_nt::SUM, sx,tx),
                expansion_nt(expansion_nt::SUM, sy,ty),
                expansion_nt(st_d)
            );
        }
        expansion& st_d_2 = expansion_product(st_d, p2.w.rep());
        st_d_2.optimize();
        const expansion& t_n_2  = expansion_product(t_n, p1.w.rep());
        const expansion& s_n_2  = expansion_product(s_n, p2.w.rep());
        const expansion& sx   = expansion_product(s_n_2, p1.x.rep());
        const expansion& tx   = expansion_product(t_n_2, p2.x.rep());
        const expansion& sy   = expansion_product(s_n_2, p1.y.rep());
        const expansion& ty   = expansion_product(t_n_2, p2.y.rep());
        return vec2HE(
            expansion_nt(expansion_nt::SUM, sx,tx),
            expansion_nt(expansion_nt::SUM, sy,ty),
            expansion_nt(st_d_2)
        );
    }

    vec3HE mix(const rational_nt& t, const vec3HE& p1, const vec3HE& p2) {
        expansion& st_d = expansion_product(t.denom().rep(),p1.w.rep());
        st_d.optimize();
        expansion& t_n  = const_cast<expansion&>(t.num().rep());
        t_n.optimize();
        const expansion& s_n  = expansion_diff(st_d, t_n);
        if(p1.w == p2.w) {
            const expansion& sx   = expansion_product(s_n, p1.x.rep());
            const expansion& tx   = expansion_product(t_n, p2.x.rep());
            const expansion& sy   = expansion_product(s_n, p1.y.rep());
            const expansion& ty   = expansion_product(t_n, p2.y.rep());
            const expansion& sz   = expansion_product(s_n, p1.z.rep());
            const expansion& tz   = expansion_product(t_n, p2.z.rep());
            return vec3HE(
                expansion_nt(expansion_nt::SUM, sx,tx),
                expansion_nt(expansion_nt::SUM, sy,ty),
                expansion_nt(expansion_nt::SUM, sz,tz),
                expansion_nt(st_d)
            );
        }
        expansion& st_d_2 = expansion_product(st_d, p2.w.rep());
        st_d_2.optimize();
        const expansion& t_n_2  = expansion_product(t_n, p1.w.rep());
        const expansion& s_n_2  = expansion_product(s_n, p2.w.rep());
        const expansion& sx   = expansion_product(s_n_2, p1.x.rep());
        const expansion& tx   = expansion_product(t_n_2, p2.x.rep());
        const expansion& sy   = expansion_product(s_n_2, p1.y.rep());
        const expansion& ty   = expansion_product(t_n_2, p2.y.rep());
        const expansion& sz   = expansion_product(s_n_2, p1.z.rep());
        const expansion& tz   = expansion_product(t_n_2, p2.z.rep());
        return vec3HE(
            expansion_nt(expansion_nt::SUM, sx,tx),
            expansion_nt(expansion_nt::SUM, sy,ty),
            expansion_nt(expansion_nt::SUM, sz,tz),
            expansion_nt(st_d_2)
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

        bool same_point(const vec2HE& v1, const vec2HE& v2) {
            return (
                ratio_compare(v1.x, v1.w, v2.x, v2.w) == ZERO &&
                ratio_compare(v1.y, v1.w, v2.y, v2.w) == ZERO
            );
        }

        bool same_point(const vec3HE& v1, const vec3HE& v2) {
            return (
                ratio_compare(v1.x, v1.w, v2.x, v2.w) == ZERO &&
                ratio_compare(v1.y, v1.w, v2.y, v2.w) == ZERO &&
                ratio_compare(v1.z, v1.w, v2.z, v2.w) == ZERO 
            );
        }
        
        Sign orient_2d(
            const vec2HE& p0, const vec2HE& p1, const vec2HE& p2
        ) {
            expansion& Delta = expansion_det3x3(
                p0.x.rep(), p0.y.rep(), p0.w.rep(),
                p1.x.rep(), p1.y.rep(), p1.w.rep(),
                p2.x.rep(), p2.y.rep(), p2.w.rep()
            );
            return Sign(
                Delta.sign()*
                p0.w.rep().sign()*
                p1.w.rep().sign()*
                p2.w.rep().sign()
            );
        }
        
        Sign dot_2d(const vec2HE& p0, const vec2HE& p1, const vec2HE& p2) {
            vec2HE U = p1 - p0;
            vec2HE V = p2 - p0;
            const expansion& x1x2 = expansion_product(U.x.rep(), V.x.rep());
            const expansion& y1y2 = expansion_product(U.y.rep(), V.y.rep());
            const expansion& S = expansion_sum(x1x2, y1y2);
            return Sign(S.sign()*U.w.sign()*V.w.sign());
        }

        Sign orient_2dlifted_SOS(
            const vec2HE& p0, const vec2HE& p1,
            const vec2HE& p2, const vec2HE& p3,
            double h0, double h1, double h2, double h3
        ) {
            expansion_nt a13(expansion_nt::DIFF, h0, h1);
            expansion_nt a23(expansion_nt::DIFF, h0, h2);
            expansion_nt a33(expansion_nt::DIFF, h0, h3);                
                
            vec2HE U1 = p1-p0;
            const expansion_nt& w1 = U1.w;
            Sign sw1 = w1.sign();
            
            vec2HE U2 = p2-p0;
            const expansion_nt& w2 = U2.w;
            Sign sw2 = w2.sign();

            vec2HE U3 = p3-p0;
            const expansion_nt& w3 = U3.w;
            Sign sw3 = w3.sign();                

            geo_assert(sw1 != ZERO && sw2 != ZERO && sw3 != ZERO);
            
            expansion_nt w2w3Delta1 = det2x2(U2.x, U2.y, U3.x, U3.y);
            expansion_nt w1w3Delta2 = det2x2(U1.x, U1.y, U3.x, U3.y);
            expansion_nt w1w2Delta3 = det2x2(U1.x, U1.y, U2.x, U2.y);
                
            Sign Delta3_sign = Sign(w1w2Delta3.sign()*sw1*sw2);
            geo_assert(Delta3_sign != ZERO);
            
            expansion_nt w1w2w3R = a13*w1*w2w3Delta1-a23*w2*w1w3Delta2+a33*w3*w1w2Delta3;
            Sign R_sign = Sign(w1w2w3R.sign()*sw1*sw2*sw3);
            
            // Simulation of simplicity
            if(R_sign == ZERO) {
                const vec2HE* p_sort[4] = {
                    &p0, &p1, &p2, &p3
                };
                std::sort(
                    p_sort, p_sort+4,
                    [](const vec2HE* p1, const vec2HE* p2)->bool{
                        vec2HELexicoCompare cmp;
                        return cmp(*p1,*p2);
                    }
                );
                for(index_t i = 0; i < 4; ++i) {
                    if(p_sort[i] == &p0) {
                        expansion_nt w1w2w3Z = w2*w1w3Delta2-w1*w2w3Delta1+w3*w1w2Delta3;
                        Sign Z_sign = Sign(w1w2w3Z.sign()*sw1*sw2*sw3);
                        if(Z_sign != ZERO) {
                            return Sign(Delta3_sign*Z_sign);
                        }
                    } else if(p_sort[i] == &p1) {
                        Sign Delta1_sign = Sign(w2w3Delta1.sign()*sw2*sw3);
                        if(Delta1_sign != ZERO) {
                            return Sign(Delta3_sign * Delta1_sign);
                        }
                    } else if(p_sort[i] == &p2) {
                        Sign Delta2_sign = Sign(w1w3Delta2.sign()*sw1*sw3);
                        if(Delta2_sign != ZERO) {
                            return Sign(-Delta3_sign * Delta2_sign);
                        }
                    } else if(p_sort[i] == &p3) {
                        return NEGATIVE;
                    }
                }
            }
            return Sign(Delta3_sign * R_sign);
        }
    }

    bool get_three_planes_intersection(
        vec3HE& result,
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
        
        result.w = det3x3(
            N1.x, N1.y, N1.z,
            N2.x, N2.y, N2.z,
            N3.x, N3.y, N3.z
        );

        if(result.w.sign() == ZERO) {
            return false;
        }
        
        result.x = det3x3(
            B.x, N1.y, N1.z,
            B.y, N2.y, N2.z,
            B.z, N3.y, N3.z
        );

        result.y = det3x3(
            N1.x, B.x, N1.z,
            N2.x, B.y, N2.z,
            N3.x, B.z, N3.z
        );

        result.z = det3x3(
            N1.x, N1.y, B.x,
            N2.x, N2.y, B.y,
            N3.x, N3.y, B.z
        );

        return true;
    }

}

