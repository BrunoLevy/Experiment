/*
 *  OGF/Graphite: Geometry and Graphics Programming Library + Utilities
 *  Copyright (C) 2000-2015 INRIA - Project ALICE
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact for Graphite: Bruno Levy - Bruno.Levy@inria.fr
 *  Contact for this Plugin: Bruno Levy - Bruno.Levy@inria.fr
 *
 *     Project ALICE
 *     LORIA, INRIA Lorraine,
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX
 *     FRANCE
 *
 *  Note that the GNU General Public License does not permit incorporating
 *  the Software into proprietary programs.
 *
 * As an exception to the GPL, Graphite can be linked with the following
 * (non-GPL) libraries:
 *     Qt, tetgen, SuperLU, WildMagic and CGAL
 */


#ifndef H__OGF_EXPERIMENT_COMMANDS_MESH_GROB_GEOBODIES_COMMANDS__H
#define H__OGF_EXPERIMENT_COMMANDS_MESH_GROB_GEOBODIES_COMMANDS__H

#include <OGF/Experiment/common/common.h>
#include <OGF/mesh/commands/mesh_grob_commands.h>

namespace OGF {

    gom_class Experiment_API MeshGrobGeobodiesCommands :
        public MeshGrobCommands {

    public:
        MeshGrobGeobodiesCommands();
        ~MeshGrobGeobodiesCommands() override;

    gom_slots:

	enum GraphiteGeobodyReconstructionStrategy {
	    GEOBODY_RECON_AUTO = 0,
	    GEOBODY_RECON_TRAINED = 1,
	    GEOBODY_RECON_OBG = 2
	};

        /**
         * \brief Reconstructs a surface from closed contours
	 * \param[in] reconstruction the reconstructed mesh
	 * \param[in] points an optional additional pointset
	 * \param[in] points_weight importance of passing through the points
	 * \param[in] lines an optional additional polygonal line (constraints)
         * \advanced
         * \param[in] resample_l edge length used internally for resampling
	 * \param[in] relative_l if set then l is relative to average edge length
         * \param[in] Poisson_depth Poisson octree depth
         */
        void reconstruct_from_contours(
            const NewMeshGrobName& reconstruction="reconstruction",
	    const NewMeshGrobName& points = "",
	    double points_weight=30.0,
	    const NewMeshGrobName& lines = "",
	    GraphiteGeobodyReconstructionStrategy strat = GEOBODY_RECON_AUTO,
            double resample_l=0.05, bool relative_l = true,
            index_t Nsmooth_iter=1000,
            index_t Poisson_depth=12
        );


	/**
	 * \brief Reconstruct from contours using dynamic programming
	 * \param[in] lines optional mesh with line constraints
	 */
	void reconstruct_from_contours_DP(const NewMeshGrobName& lines = "");

        /**
         * \brief Resamples a polyline made of several closed contours
         * \param[out] resample the name of the resampled mesh
         * \param[in] l the desired edge length
	 * \param[in] relative_l if set then l is relative to average edge length
         */
        void resample_contours(
            const NewMeshGrobName& resample="resample",
	    double l=0.05, bool relative_l=true
        );

        /**
         * \brief Computes the normals of a contour polyline
	 * \param[in] strat strategy for normal estimation
         * \param[in] show if set, creates a polyline with normals display
         * \param[in] scale normals scale used if show is set
         */
        void compute_contours_normals(
	    GraphiteGeobodyReconstructionStrategy strat = GEOBODY_RECON_AUTO,
            index_t Nsmooth_iter=1000, bool show = false, double scale = 30.0
        );


    } ;
}

#endif
