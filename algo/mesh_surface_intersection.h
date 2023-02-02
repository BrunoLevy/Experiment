
/*
 *
 */
 
#ifndef H__OGF_EXPERIMENT_MESH_SURFACE_INTERSECTION__H
#define H__OGF_EXPERIMENT_MESH_SURFACE_INTERSECTION__H

#include <OGF/Experiment/common/common.h>

namespace GEO {

    class Mesh;
    
    struct MeshSurfaceIntersectionParams {
        /**
         * \brief detect and fix duplicated vertices in input 
         */
        bool pre_detect_duplicated_vertices = true;
        
        /**
         * \brief detect and fix duplicated facets in input 
         */
        bool pre_detect_duplicated_facets = true;
        
        /** 
         * \brief detect and compute intersections between facets that share 
         *  a facet or an edge
         */
        bool detect_intersecting_neighbors = true;

        /**
         * \brief connect the facets of the arrangement and split the 
         *  non-manifold edges
         */
        bool post_connect_facets = true;

        /**
         * \brief Do not order the facets in the AABB. Used for debugging,
         *  for keeping facet ids (need to order the facets before, else
         *  the AABB will take ages).
         */
        bool debug_do_not_order_facets = false;

        /**
         * \brief Enable floating point exceptions, to detect overflows and
         *  underflows (shit happens ! [Forrest Gump]).
         */
        bool debug_enable_FPE = true;

        /**
         * \brief Double-check that the constraints do not self intersect,
         *  that no vertex is outside the triangle, and that there is no
         *  duplicated vertex
         */
        bool debug_check_constraints = true;

        /**
         * \brief Use barycentric coordinates in each triangle. The other
         *  option is projection onto the most normal axis.
         */
        bool barycentric = true;
    };

    void Experiment_API mesh_intersect_surface(
        Mesh& M, const MeshSurfaceIntersectionParams& params
    );

    void Experiment_API mesh_classify_intersections(
        Mesh& M, const std::string& expr
    );
}

#endif

