
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
         * \brief if set, assign an operand id to each connected component
         *  of input.
         */
        bool per_component_ids = true;

        /**
         * \brief Display information while computing the intersection
         */
        bool verbose = true;

        /**
         * \brief If set, compute constrained Delaunay triangulation
         *  in the intersected triangles. If there are intersections
         *  in coplanar facets, it guarantees uniqueness of their
         *  triangulation.
         */
        bool delaunay = false;

        /**
         * \brief If set, then Delaunay mode uses approximated incircle
         *  predicate (else it uses exact arithmetics)
         */
        bool approx_incircle = false;
    };

    void Experiment_API mesh_intersect_surface(
        Mesh& M, const MeshSurfaceIntersectionParams& params
    );

    void Experiment_API mesh_classify_intersections(
        Mesh& M, std::function<bool(index_t)> eqn,
        const std::string& attribute="", bool reorder=true
    );

    void Experiment_API mesh_classify_intersections(
        Mesh& M, const std::string& expr,
        const std::string& attribute="", bool reorder=true
    );
    
}

#endif

