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

#ifndef GEOGRAM_NUMERICS_INTERVAL_GEOMETRY
#define GEOGRAM_NUMERICS_INTERVAL_GEOMETRY

#include <OGF/Experiment/common/common.h>
#include <geogram/basic/geometry.h>
#include <geogram/numerics/interval_nt.h>

#include <geogram/basic/common.h>
#include <geogram/basic/geometry.h>
#include <geogram/numerics/expansion_nt.h>

/**
 * \file geogram/numerics/exact_geometry.h
 * \brief Exact predicates and constructs
 * \details Implements vector types with interval
 *  coordinates (vec2I, vec3I), vector types with
 *  homogeneous expansion coordinates (vec2HI, vec3HI),
 *  and filters for 2d orientation predicate, 
 *  3d-lifted orientation predicate 
 *  (can be used to implement incircle).
 */

namespace GEO {

    struct vec2HE;
    struct vec3HE;
    
    /**
     * \brief vec2 with coordinates as expansions
     * \details Coordinates support +,-,*
     */
    typedef vecng<2,interval_nt> vec2I;

    /**
     * \brief vec3 with coordinates as expansions
     * \details Coordinates support +,-,*
     */
    typedef vecng<3,interval_nt> vec3I;    


    /**
     * \brief 2D vector in homogeneous coordinates
     *  with coordinates as arithmetic expansions
     * \details Coordinates support +,-,* and / by
     *  multiplying w.
     */
    struct GEOGRAM_API vec2HI {

        vec2HI(
            const interval_nt& x_in,
            const interval_nt& y_in,
            const interval_nt& w_in
        ) : x(x_in), y(y_in), w(w_in) {
        }

        vec2HI(const vec2HI& rhs) :
            x(rhs.x), y(rhs.y), w(rhs.w) {
        }

        explicit vec2HI(const vec2& rhs) : 
            x(rhs.x), y(rhs.y), w(1.0) {
        }

        explicit vec2HI(const vec2HE& rhs);
        
        vec2HI& operator=(const vec2HI& rhs) {
            if(&rhs != this) {
                x=rhs.x;
                y=rhs.y;
                w=rhs.w;
            }
            return *this;
        }

        vec2HI& operator=(const vec2HE& rhs);
            
        interval_nt* data() {
            return &x;
        }

        const interval_nt* data() const {
            return &x;
        }
        
        interval_nt& operator[](coord_index_t i) {
            geo_debug_assert(i < 2);
            return data()[i];
        }

        const interval_nt& operator[](coord_index_t i) const {
            geo_debug_assert(i < 2);
            return data()[i];
        }

        interval_nt x;
        interval_nt y;
        interval_nt w;
    };

    /**
     * \brief 3D vector in homogeneous coordinates
     *  with coordinates as arithmetic expansions.
     * \details Coordinates support +,-,* and / by
     *  multiplying w.
     */
    struct GEOGRAM_API vec3HI {

        vec3HI(
            const interval_nt& x_in,
            const interval_nt& y_in,
            const interval_nt& z_in,
            const interval_nt& w_in            
        ) : x(x_in), y(y_in), z(z_in), w(w_in) {
        }

        vec3HI(
            double x_in,
            double y_in,
            double z_in,
            double w_in
        ) : x(x_in), y(y_in), z(z_in), w(w_in) {
        }
        
        vec3HI(const vec3HI& rhs) :
            x(rhs.x), y(rhs.y), z(rhs.z), w(rhs.w) {
        }

        explicit vec3HI(const vec3& rhs) : 
            x(rhs.x), y(rhs.y), z(rhs.z), w(1.0) {
        }

        explicit vec3HI(const vec3HE& rhs);
            
        vec3HI& operator=(const vec3HI& rhs) {
            if(&rhs != this) {
                x=rhs.x;
                y=rhs.y;
                z=rhs.z;                
                w=rhs.w;
            }
            return *this;
        }

        vec3HI& operator=(const vec3HE& rhs);            

        interval_nt* data() {
            return &x;
        }

        const interval_nt* data() const {
            return &x;
        }
        
        interval_nt& operator[](coord_index_t i) {
            geo_debug_assert(i < 3);
            return data()[i];
        }

        const interval_nt& operator[](coord_index_t i) const {
            geo_debug_assert(i < 3);
            return data()[i];
        }

        interval_nt x;
        interval_nt y;
        interval_nt z;        
        interval_nt w;
    };

    namespace PCK {
        
        Sign GEOGRAM_API orient_2d_projected(
            const vec3HI& p0, const vec3HI& p1, const vec3HI& p2,
            coord_index_t axis,
            bool& determined
        );
        
    }
    
}

#endif


