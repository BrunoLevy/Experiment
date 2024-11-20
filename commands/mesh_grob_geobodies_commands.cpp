
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

#include <OGF/Experiment/commands/mesh_grob_geobodies_commands.h>

#ifdef GEOGRAM_WITH_VORPALINE
#include <vorpalib/mesh/mesh_geobodies.h>
#endif


namespace OGF {

    MeshGrobGeobodiesCommands::MeshGrobGeobodiesCommands() {
    }

    MeshGrobGeobodiesCommands::~MeshGrobGeobodiesCommands() {
    }

    void MeshGrobGeobodiesCommands::compute_contours_normals(
        index_t Nsmooth_iter, bool show, double scale
    ) {
#ifdef GEOGRAM_WITH_VORPALINE
        GEO::compute_contours_normals(mesh_grob(), Nsmooth_iter);
        Attribute<double> normal(
            mesh_grob()->vertices.attributes(), "normal"
        );


        mesh_grob()->update();
        // For debugging
        // (use scale = 30 with ore_body)
        if(show) {
            MeshGrob* normals =
                MeshGrob::find_or_create(scene_graph(),"normals");
            normals->clear();
            normals->vertices.set_dimension(3);
            for(index_t v: mesh_grob()->vertices) {
                vec3 p(mesh_grob()->vertices.point_ptr(v));
                vec3 n(normal[3*v],normal[3*v+1],normal[3*v+2]);
                index_t v1 = normals->vertices.create_vertex(p.data());
                index_t v2 = normals->vertices.create_vertex((p+scale*n).data());
                normals->edges.create_edge(v1,v2);
            }
            normals->update();
        }
#else
        geo_argused(Nsmooth_iter);
        geo_argused(show);
        geo_argused(scale);
        Logger::err("Geobodies") << "Needs Tessael's VORPALINE component"
                                 << std::endl;
#endif
    }


    void MeshGrobGeobodiesCommands::resample_contours(
        const NewMeshGrobName& resample_name, double l
    ) {
#ifdef GEOGRAM_WITH_VORPALINE
        if(resample_name == mesh_grob()->name()) {
            Logger::err("Resample") << "resample and input cannot be the same"
                                    << std::endl;
            return;
        }
        MeshGrob* resample = MeshGrob::find_or_create(
            scene_graph(),resample_name
        );
        GEO::resample_contours(mesh_grob(), resample, l);
        resample->update();
#else
        geo_argused(resample_name);
        geo_argused(l);
        Logger::err("Geobodies") << "Needs Tessael's VORPALINE component"
                                 << std::endl;
#endif
    }

    void MeshGrobGeobodiesCommands::reconstruct_from_contours(
        const NewMeshGrobName& recon_name,
        double resample_l, index_t Nsmooth_iter, index_t Poisson_depth
    ) {
#ifdef GEOGRAM_WITH_VORPALINE
        if(recon_name == mesh_grob()->name()) {
            Logger::err("Recon") << "recon and input cannot be the same"
                                    << std::endl;
            return;
        }
        MeshGrob* recon = MeshGrob::find_or_create(scene_graph(),recon_name );
        recon->clear();
        GEO::reconstruct_from_contours(
            mesh_grob(), recon, resample_l, Nsmooth_iter, Poisson_depth
        );
        recon->update();
#else
        geo_argused(recon_name);
        geo_argused(resample_l);
        geo_argused(Nsmooth_iter);
        geo_argused(Poisson_depth);
        Logger::err("Geobodies") << "Needs Tessael's VORPALINE component"
                                 << std::endl;
#endif
    }

    void MeshGrobGeobodiesCommands::reconstruct_from_contours_and_points(
        const NewMeshGrobName& recon_name,
	const MeshGrobName& points_name,
	double points_weight,
        double resample_l, index_t Nsmooth_iter, index_t Poisson_depth
    ) {
#ifdef GEOGRAM_WITH_VORPALINE
        if(recon_name == mesh_grob()->name()) {
            Logger::err("Recon") << "recon and input cannot be the same"
                                    << std::endl;
            return;
        }
        MeshGrob* recon = MeshGrob::find_or_create(scene_graph(),recon_name );
        recon->clear();
	MeshGrob* points = nullptr;
	if(points_name != "") {
	    points = MeshGrob::find(scene_graph(), points_name);
	}
        GEO::reconstruct_from_contours_and_points(
            mesh_grob(), points, recon,
	    points_weight, resample_l, Nsmooth_iter, Poisson_depth
        );
        recon->update();
#else
        geo_argused(recon_name);
	geo_argused(points_name);
	geo_argused(points_weight);
        geo_argused(resample_l);
        geo_argused(Nsmooth_iter);
        geo_argused(Poisson_depth);
        Logger::err("Geobodies") << "Needs Tessael's VORPALINE component"
                                 << std::endl;
#endif
    }

}
