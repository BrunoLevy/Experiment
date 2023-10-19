/*
 *  Copyright (c) 2000-2023 Inria
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

#include <OGF/Experiment/algo/interval_geometry.h>
#include <geogram/numerics/exact_geometry.h>

namespace GEO {
    
    namespace PCK {
        Sign orient_2d_projected(
            const vec3HI& p0, const vec3HI& p1, const vec3HI& p2,
            coord_index_t axis,
            bool& determined
        ) {
            determined = false;
            
            coord_index_t u = coord_index_t((axis+1)%3);
            coord_index_t v = coord_index_t((axis+2)%3);

            interval_nt::Sign2 sw0 = p0.w.sign();
            interval_nt::Sign2 sw1 = p1.w.sign();
            interval_nt::Sign2 sw2 = p2.w.sign();                        

            if(
                !interval_nt::sign_is_determined(sw0) ||
                !interval_nt::sign_is_determined(sw1) ||
                !interval_nt::sign_is_determined(sw2)
            ) {
                return ZERO;
            }
            
            interval_nt Delta = det3x3(
                p0[u], p0[v], p0.w,
                p1[u], p1[v], p1.w,
                p2[u], p2[v], p2.w
            );

            interval_nt::Sign2 sD = Delta.sign();

            if(!interval_nt::sign_is_determined(sD)) {
                return ZERO;
            }

            determined = true;

            std::cerr << " [" << Delta.inf() << " " << Delta.sup() << "](" << Delta.sup() - Delta.inf() << ")(" << int(Delta.sign()) << ") "
                      << " [" << p0.w.inf() << " " << p0.w.sup() << "](" << p0.w.sup() - p0.w.inf() << ")(" << int(p0.w.sign()) << ") "
                      << " [" << p1.w.inf() << " " << p1.w.sup() << "](" << p1.w.sup() - p1.w.inf() << ")(" << int(p1.w.sign()) << ") "
                      << " [" << p2.w.inf() << " " << p2.w.sup() << "](" << p2.w.sup() - p2.w.inf() << ")(" << int(p2.w.sign()) << ") "
                      << std::endl;

            std::cerr << int(interval_nt::convert_sign(sD)) << " "
                      << int(interval_nt::convert_sign(sw0)) << " "
                      << int(interval_nt::convert_sign(sw1)) << " "
                      << int(interval_nt::convert_sign(sw2)) << std::endl;
            
            return Sign(
                interval_nt::convert_sign(sD)*
                interval_nt::convert_sign(sw0)*
                interval_nt::convert_sign(sw1)*
                interval_nt::convert_sign(sw2)
            );
        }
    }

    
}
