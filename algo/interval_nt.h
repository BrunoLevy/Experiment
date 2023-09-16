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

#ifndef GEO_INTERVAL_NT
#define GEO_INTERVAL_NT

#include <geogram/numerics/expansion_nt.h>
#include <iomanip>
#include <limits>

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
            if(s == SIGN2_NEGATIVE) {
                return NEGATIVE;
            }
            if(s == SIGN2_POSITIVE) {
                return POSITIVE;
            }
            return ZERO;
        }
        
        interval_nt& negate() {
            lb_ = -lb_;
            ub_ = -ub_;
            control_.rep().negate();
            std::swap(lb_, ub_);
            check();
            return *this;
        }
        
        interval_nt& operator+=(const interval_nt &x) {
            lb_ += x.lb_;
            ub_ += x.ub_;
            control_ += x.control_;
            adjust();
            check();
            return *this;
        }
        
        interval_nt& operator-=(const interval_nt &x) {
            lb_ -= x.ub_;
            ub_ -= x.lb_;
            control_ -= x.control_;
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
            control_ *= x.control_;
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
            // std::cerr << "Ok" << std::flush;
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
