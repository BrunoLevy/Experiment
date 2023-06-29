/*
 *  Copyright (c) 2000-2022 Inria
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *
 *  * Redistributions of source code must retain the above copyright notice,
 *  this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright notice,
 *  this list of conditions and the following disclaimer in the documentation
 *  and/or other materials provided with the distribution.
 *  * Neither the name of the ALICE Project-Team nor the names of its
 *  contributors may be used to endorse or promote products derived from this
 *  software without specific prior written permission.
 * 
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 *  Contact: Bruno Levy
 *
 *     https://www.inria.fr/fr/bruno-levy
 *
 *     Inria,
 *     Domaine de Voluceau,
 *     78150 Le Chesnay - Rocquencourt
 *     FRANCE
 *
 */

#include <iomanip>
#include <limits>

//#define USE_BASIC_INTERVAL 
//#define USE_RH_INTERVAL

#ifdef USE_RH_INTERVAL
// Interval class in "round to nearest" mode, by Richard Harris:
// https://accu.org/journals/overload/19/103/harris_1974/
// Propagates proportional errors at a rate of 1+/-0.5eps
// Handles denormals properly (as a special case)

namespace GEO {
    
    class interval_nt {
    public:
        interval_nt() : lb_(0.0), ub_(0.0), control_(0.0) {
            check();
        }

        interval_nt(const double x) : lb_(x), ub_(x), control_(x) {
            check();            
        }

        interval_nt(const interval_nt& rhs) = default;

        interval_nt(const expansion_nt& rhs) {
            *this = rhs;
            check();
        }
        
        interval_nt& operator=(const interval_nt& rhs) = default;
        
        interval_nt& operator=(double rhs) {
            lb_ = rhs;
            ub_ = rhs;
            control_=expansion_nt(rhs);
            check();
            return *this;
        }

        interval_nt& operator=(const expansion_nt& rhs) {
            *this = rhs.component(0);
            for(index_t i=1; i<rhs.length(); ++i) {
                *this += rhs.component(i);
            }
            control_ = rhs;
            return *this;
        }
        
        double inf() const {
            return lb_;
        }
        
        double sup() const {
            return ub_; 
        }
        
        bool is_nan() const {
            return !(lb_==lb_) || !(ub_==ub_);
        }

	enum Sign2 {
            SIGN2_NEGATIVE=-1,
            SIGN2_ZERO=0,
            SIGN2_POSITIVE=1,
            SIGN2_UNDETERMINED=2
	};

        Sign2 sign() const {
            geo_assert(!is_nan());
            if(lb_ == 0.0 && ub_ == 0.0) {
                return SIGN2_ZERO;
            }
            if(ub_ < 0.0) {
                return SIGN2_NEGATIVE;
            }
            if(lb_ > 0.0) {
                return SIGN2_POSITIVE;
            }
            return SIGN2_UNDETERMINED;
        }

        static bool sign_is_determined(Sign2 s) {
            return
                s == SIGN2_ZERO ||
                s == SIGN2_NEGATIVE ||
                s == SIGN2_POSITIVE ;
        }

        static Sign convert_sign(Sign2 s) {
            geo_assert(sign_is_determined(s));
            if(s == SIGN2_NN) {
                return NEGATIVE;
            }
            if(s == SIGN2_PP) {
                return POSITIVE;
            }
            return ZERO;
        }
        
        interval_nt& negate() {
            lb_ = -lb_;
            ub_ = -ub_;
            std::swap(lb_, ub_);
            check();
            return *this;
        }
        
        interval_nt& operator+=(const interval_nt &x) {
            lb_ += x.lb_;
            ub_ += x.ub_;
            adjust();
            check();
            return *this;
        }
        
        interval_nt& operator-=(const interval_nt &x) {
            lb_ -= x.ub_;
            ub_ -= x.lb_;
            adjust();
            check();
            return *this;
        }
        
        interval_nt& operator*=(const interval_nt &x) {
            if(!is_nan() && !x.is_nan()) {
                double ll = lb_*x.lb_;
                double lu = lb_*x.ub_;
                double ul = ub_*x.lb_;
                double uu = ub_*x.ub_;
                
                if(!(ll==ll)) ll = 0.0;
                if(!(lu==lu)) lu = 0.0;
                if(!(ul==ul)) ul = 0.0;
                if(!(uu==uu)) uu = 0.0;

                if(lu<ll) std::swap(lu, ll);
                if(ul<ll) std::swap(ul, ll);
                if(uu<ll) std::swap(uu, ll);

                if(lu>uu) uu = lu;
                if(ul>uu) uu = ul;
                
                lb_ = ll;
                ub_ = uu;
                
                adjust();
            } else {
                lb_ = std::numeric_limits<double>::quiet_NaN();
                ub_ = std::numeric_limits<double>::quiet_NaN();
            }
            check();
            return *this;            
        }
        
    protected:
        void adjust() {
            add_error(lb_, ub_);
        }
        
        static void add_error(double &lb, double &ub) {
            static const double i = std::numeric_limits<double>::infinity();
            static const double e = std::numeric_limits<double>::epsilon();
            static const double m = std::numeric_limits<double>::min();
            static const double l = 1.0-e;
            static const double u = 1.0+e;

            if(lb==lb && ub==ub && (lb!=ub || (lb!=i && lb!=-i))) {
                if(lb>ub) {
                    std::swap(lb, ub);
                }

                if(lb>m) {
                    lb *= l;
                } else if(lb<-m) {
                    lb *= u;
                } else {
                    lb -= e*m;
                }

                if(ub>m) {
                    ub *= u;
                } else if(ub<-m) {
                    ub *= l;
                } else {
                    ub += e*m;
                }
            } else {
                lb = std::numeric_limits<double>::quiet_NaN();
                ub = std::numeric_limits<double>::quiet_NaN();
            }
        }

        void check() const {
            typedef std::numeric_limits< double > dbl;
            if(inf() > sup()) {
                std::cerr.precision(dbl::max_digits10);
                std::cerr << "inf() > sup() !!" << std::endl;
                std::cerr << "inf()=" << inf() << std::endl;
                std::cerr << "sup()=" << sup() << std::endl;
                geo_assert_not_reached;
            }
            if(control_ < inf() || control_ > sup()) {
                std::cerr.precision(dbl::max_digits10);
                std::cerr << "[" << inf() << "," << sup() << "]"
                          << "   " << control_.estimate() << ":"
                          << control_.rep().length()
                          << std::endl;
                expansion_nt control1 = control_ - inf();
                expansion_nt control2 = sup() - control_;
                std::cerr << control1.estimate() << " "
                          << control2.estimate() << std::endl;
                geo_assert_not_reached;
            }
        }
        
    private:
        double lb_;
        double ub_;
        expansion_nt control_;
    };

    inline interval_nt operator+(const interval_nt& a, const interval_nt& b) {
        interval_nt result = a;
        return result += b;
    }

    inline interval_nt operator-(const interval_nt& a, const interval_nt& b) {
        interval_nt result = a;
        return result -= b;
    }

    inline interval_nt operator*(const interval_nt& a, const interval_nt& b) {
        interval_nt result = a;
        return result *= b;
    }
    
}
#endif




#ifdef USE_BASIC_INTERVAL 
#ifndef GEOGRAM_NUMERICS_INTERVAL_NT
#define GEOGRAM_NUMERICS_INTERVAL_NT
#include <OGF/Experiment/third_party/BasicInterval/interval.h>
namespace GEO {
    typedef numerics::interval interval_nt;
}
#endif
#endif

#ifndef GEOGRAM_NUMERICS_INTERVAL_NT
#define GEOGRAM_NUMERICS_INTERVAL_NT

#include <OGF/Experiment/common/common.h>
#include <geogram/basic/memory.h>
#include <geogram/numerics/expansion_nt.h>

/**
 * \file geogram/numerics/expansion_nt.h
 * \brief Interval number type
 * \details This file provides a "number-type" that implements interval 
 *   arithmetics.
 * Inspiration:
 *   - https://www.codeproject.com/Articles/1040839/Interval-arithmetic-using-round-to-nearest-mode-pa
 * References (thank you Guillaume Moroz):
 *   - https://github.com/goualard-f/GAOL
 *   - Interval Arithmetic: from Principles to Implementation
 *     Timothy J. Hickey, Qun Ju, and Maarten H. van Emden
 *      https://doi.org/10.1145/502102.502106
 *   - Ball arithmetics: https://arblib.org/
 */

// For now, only supported with AVX2 (needs FMA)
#ifdef __AVX2__
#include <immintrin.h> 



namespace GEO {

    class expansion_nt;
    
    class Experiment_API interval_nt {
    public:
        interval_nt() : value_(zero()), control_(0.0) {
            check();
        }
        
        interval_nt(double rhs) : value_(
            (rhs == 0.0) ? zero() : _mm_set_pd(rhs,rhs)
        ), control_(rhs) {
            check();
        }

        interval_nt(const interval_nt& rhs) = default;

        interval_nt(const expansion_nt& rhs) {
            *this = rhs;
            check();
        }
        
        interval_nt& operator=(const interval_nt& rhs) = default;
        
        interval_nt& operator=(double rhs) {
            value_=((rhs==0.0) ? zero() : _mm_set_pd(rhs,rhs));
            control_=expansion_nt(rhs);
            check();
            return *this;
        }

        interval_nt& operator=(const expansion_nt& rhs);        

        static void check(__m128d I) {
            typedef std::numeric_limits< double > dbl;
            double i = inf(I);
            double s = sup(I);
            if(i > s) {
                std::cerr.precision(dbl::max_digits10);
                std::cerr << "i > s !!" << std::endl;
                std::cerr << "i=" << i << std::endl;
                std::cerr << "s=" << s << std::endl;
                geo_assert_not_reached;
            }
        }
        
        void check() const {
            typedef std::numeric_limits< double > dbl;
            if(inf() > sup()) {
                std::cerr.precision(dbl::max_digits10);
                std::cerr << "inf() > sup() !!" << std::endl;
                std::cerr << "inf()=" << inf() << std::endl;
                std::cerr << "sup()=" << sup() << std::endl;
                geo_assert_not_reached;
            }
            if(control_ < inf() || control_ > sup()) {
                std::cerr.precision(dbl::max_digits10);
                std::cerr << "[" << inf() << "," << sup() << "]"
                          << "   " << control_.estimate() << ":" << control_.rep().length()
                          << std::endl;
                expansion_nt control1 = control_ - inf();
                expansion_nt control2 = sup() - control_;
                std::cerr << control1.estimate() << " " << control2.estimate() << std::endl;
                geo_assert_not_reached;
            }
        }
        
	enum Sign2 {
            SIGN2_ERROR = -1,
            SIGN2_ZERO,
            SIGN2_NP,
            SIGN2_PP,
            SIGN2_ZP,
            SIGN2_NN,
            SIGN2_NZ,
            SIGN2_COUNT
	};

        Sign2 sign() const {
            return sign(value_);
        }

        static bool sign_is_determined(Sign2 s) {
            return
                s == SIGN2_ZERO ||
                s == SIGN2_NN ||
                s == SIGN2_PP ;
        }

        static Sign convert_sign(Sign2 s) {
            geo_assert(sign_is_determined(s));
            if(s == SIGN2_NN) {
                return NEGATIVE;
            }
            if(s == SIGN2_PP) {
                return POSITIVE;
            }
            return ZERO;
        }

        double inf() const {
	   return inf(value_);
	}

        double sup() const {
	   return sup(value_);
	}
       
        interval_nt operator-() const {
            // Swap both components
            __m128d swapped = _mm_shuffle_pd(value_, value_, 1);
            // Flip signs
            __m128d sign_mask = _mm_set_pd(-0.0, -0.0);
            interval_nt result = interval_nt(_mm_xor_pd(swapped,sign_mask));
            result.control_ = -control_;
            result.check();
            return result;
        }
        
        interval_nt& operator+=(const interval_nt& rhs) {
            __m128d err;
            value_ = TwoSum(value_, rhs.value_, err);
            check(value_);
            adjust(err);
            control_ += rhs.control_;
            check();
            return *this;
        }
        
        interval_nt& operator-=(const interval_nt& rhs) {
            rhs.check();
            check();

            /*
            __m128d A = value_;
            __m128d B = rhs.value_;
            __m128d E;
            __m128d C = TwoDiff(A,B,E);
            */
            __m128d err;
            value_ = TwoDiff(value_, rhs.value_, err);
            // TODO: check TwoDiff result.
            
            adjust(err);
            control_ -= rhs.control_;
            check();
            return *this;
        }
        
        interval_nt& operator*=(const interval_nt& rhs) {
            __m128d err;
            value_ = TwoProd(value_, rhs.value_, err);
            adjust(err);
            control_ *= rhs.control_;
            check();
            return *this;
        }

        bool operator==(const interval_nt& rhs) const {
            __m128d mask = _mm_cmpeq_pd(value_, rhs.value_);
            int result = _mm_movemask_pd(mask);
            return (result == 3);
        }
        
        bool operator!=(const interval_nt& rhs) const {
            return !operator==(rhs);
        }
        
    protected:
        interval_nt(__m128d rhs) : value_(rhs) {
        }
        
        static Sign2 sign(__m128d value);
        
        static inline __m128d zero() {
            // Note: _mm_set_pd() arguments are inversed,
            // it is _mm_set_pd(e1,e0)
            // (why on earth did they pick this order ?)
            return _mm_set_pd(0.0,-0.0);
        }

        static __m128d TwoSum(__m128d a, __m128d b, __m128d& err) {
            __m128d result = _mm_add_pd(a, b);
            __m128d bb = _mm_sub_pd(result, a);
            err = _mm_add_pd(
                _mm_sub_pd(a, _mm_sub_pd(result, bb)), _mm_sub_pd(b, bb)
            );
            return adjust_zero(result);
        }

        static __m128d TwoDiff(__m128d a, __m128d b, __m128d& err) {
            __m128d result = _mm_sub_pd(a, b);
            __m128d bb = _mm_sub_pd(result, a);
            err = _mm_sub_pd(
                _mm_sub_pd(a, _mm_sub_pd(result, bb)), _mm_add_pd(b, bb)
            );
            return adjust_zero(result);
        }

        static __m128d TwoProd(__m128d a, __m128d b, __m128d& err);

        static __m128d adjust_zero(__m128d x) {
            // Can we do something smarter here ?
            //  (in SSE registers, branchless) ?
            geo_decl_aligned(double unpacked[2]);
            _mm_store_pd(unpacked,x);
            
            if(unpacked[0] == 0.0 && unpacked[1] == 0.0) {
                return zero();
            }
            
            //	[0, b]	- return [+0, b]                
            if(unpacked[0] == 0.0) {
                return _mm_andnot_pd(_mm_set1_pd(-0.0),x);
            }

            //	[a, 0]	- return [a, -0]
            if(unpacked[1] == 0.0) {
                return _mm_or_pd(_mm_set_pd(-0.0,-0.0),x);
            }
            
            return x;
        }

        static double inf(__m128d x) {
            geo_decl_aligned(double retval[2]);
            _mm_store_pd(retval, x);
            return retval[0];
        }

        static double sup(__m128d x) {
            geo_decl_aligned(double retval[2]);
            _mm_store_pd(retval, x);
            return retval[1];
        }
        
        void adjust(__m128d err);
        
    private:
        __m128d value_;
        expansion_nt control_;
    };


    inline interval_nt operator+(const interval_nt& a, const interval_nt& b) {
        interval_nt result = a;
        return result += b;
    }

    inline interval_nt operator-(const interval_nt& a, const interval_nt& b) {
        interval_nt result = a;
        return result -= b;
    }

    inline interval_nt operator*(const interval_nt& a, const interval_nt& b) {
        interval_nt result = a;
        return result *= b;
    }
    
}

#endif
#endif
