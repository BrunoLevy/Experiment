#include <OGF/Experiment/algo/CDT_with_interval.h>
#include <OGF/Experiment/algo/interval_geometry.h>

namespace GEO {

    Sign CDT2d_with_interval::orient2d(
        index_t i, index_t j, index_t k
    ) const {
        ++orient_cnt_total_;
        const vec2& pi = point_[i];
        const vec2& pj = point_[j];
        const vec2& pk = point_[k];

        vec2I piI(pi.x, pi.y);
        vec2I pjI(pj.x, pj.y);
        vec2I pkI(pk.x, pk.y);

        interval_nt Delta = det(pjI-piI,pkI-piI);
        interval_nt::Sign2 Delta_s = Delta.sign();
        Sign result = CDT2d::orient2d(i,j,k);

        if(interval_nt::sign_is_determined(Delta_s)) {
            Sign result_i = interval_nt::convert_sign(Delta_s);
            if(result_i != result) {
                std::cerr << "=============> orient2d interval disagreement: "
                          << int(result) << " " << int(result_i) 
                          << std::endl;
            }
        } else {
            ++orient_cnt_fail_;
        }

        return result;
    }

    Sign CDT2d_with_interval::incircle(
        index_t i,index_t j,index_t k,index_t l
    ) const {
        ++incircle_cnt_total_;
        
        const vec2& pi = point_[i];
        const vec2& pj = point_[j];
        const vec2& pk = point_[k];
        const vec2& pl = point_[l];

        vec2I piI(pi.x, pi.y);
        vec2I pjI(pj.x, pj.y);
        vec2I pkI(pk.x, pk.y);
        vec2I plI(pl.x, pl.y);        

        vec2I U = piI - plI;
        vec2I V = pjI - plI;
        vec2I W = pkI - plI;
        
        interval_nt Delta = det3x3(
            U.x, U.y, dot(U,U),
            V.x, V.y, dot(V,V),
            W.x, W.y, dot(W,W)
        );

        interval_nt::Sign2 Delta_s = Delta.sign();
        Sign result = CDT2d::incircle(i,j,k,l); // calls PCK::in_circle_2d_SOS() -> PCK::side3_2d()
                                                // and side3_2d() has triangle orientation
                                                // --> so we need to multiply with global orientation orient_012_

        if(interval_nt::sign_is_determined(Delta_s) && Delta_s != interval_nt::Sign2(interval_nt::SIGN2_ZERO)) {
            Sign result_i = Sign(interval_nt::convert_sign(Delta_s) * orient_012_);
            if(result_i != result) {
                std::cerr << "=============> incircle2d interval disagreement: "
                          << int(result) << " " << int(result_i)
                          << std::endl;
            } /* else {
                 std::cerr << "=============> incircle2d interval OK "
                           << std::endl;
            } */
        } else {
            ++incircle_cnt_fail_;
        }
        

        return result;
    }
}

