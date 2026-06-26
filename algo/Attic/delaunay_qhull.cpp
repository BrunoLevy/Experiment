/*
 *  OGF/Graphite: Geometry and Graphics Programming Library + Utilities
 *  Copyright (C) 2000 Bruno Levy
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
 *  Contact: Bruno Levy
 *
 *     levy@loria.fr
 *
 *     ISA Project
 *     LORIA, INRIA Lorraine,
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX
 *     FRANCE
 *
 *  Note that the GNU General Public License does not permit incorporating
 *  the Software into proprietary programs.
 */

#include <OGF/Experiment/algo/delaunay_qhull.h>
#include <geogram/basic/process.h>
#include <fstream>

extern "C" {
#include <OGF/Experiment/third_party/qhull/qhull_a.h>
}

namespace OGF {

    //--------------------------------------------------------------------------

    DelaunayNd_QHull::DelaunayNd_QHull(index_t dimension) :
	Delaunay(coord_index_t(dimension)) {
    }

    DelaunayNd_QHull::~DelaunayNd_QHull() {
    }

    void DelaunayNd_QHull::set_vertices(
	index_t nb_vertices, const double* vertices
    ) {
        ogf_assert(!Process::is_running_threads()) ;
	Delaunay::set_vertices(nb_vertices, vertices) ;
	qhull_.compute(
	    dimension(), nb_vertices, vertices, cell_to_v_, cell_to_cell_
	);
	set_arrays(
	    index_t(cell_to_v_.size() / cell_size()),
	    cell_to_v_.data(), cell_to_cell_.data()
	) ;
    }

    //---------------------------------------------------------------------------

    static const char* qhull_command = "qdelaunay QJ" ;
    static const char* hidden_options =
    " d n v H U Qb QB Qc Qf Qg Qi Qm Qr QR Qv Qx TR E V FC Fi Fo Ft Fp FV Q0 Q1 Q2 Q3 Q4 Q5 Q6 Q7 Q8 Q9 ";

    QHull::QHull() {
        ogf_assert(!Process::is_running_threads()) ;
	static boolT firstcall = True ;
	if (firstcall) {
	    qh_meminit(stderr);
	    firstcall= False;
	}
    }

    QHull::~QHull() {
    }

    void QHull::compute(
	index_t dimension, index_t nb_vertices, const double* vertices,
	std::vector<index_t>& cell_to_v, std::vector<index_t>& cell_to_cell
    ) {
	// Step 1: initialize and call QHull

	qh_initqhull_start(NULL, stdout, stderr);
	int exitcode = setjmp(qh errexit);
	if(!exitcode) {
	    qh NOerrexit = False;
	    qh_initflags((char*)qhull_command);
	    qh DELAUNAY = True ;
	    qh PROJECTdelaunay = True ;
	    qh SCALElast= True;    /* 'Qbb' */

	    // It seems that we do not need this one.
            // qh KEEPcoplanar= True; /* 'Qc', to keep coplanars in 'p' */

	    // Behave like command-line version
	    if(dimension < 5) {
		// Option C-0
		qh premerge_centrum = 0 ;
		qh PREmerge = True ;
	    } else {
		// Option Qx
		qh MERGEexact = True ;
	    }

	    qh_checkflags(qh qhull_command, (char*)hidden_options);
	    qh_initflags(qh qhull_command);
	    qh_init_B(
		(double*)vertices, int(nb_vertices), int(dimension), False
	    ) ;
	    qh_qhull() ;

	    // Step 2: get result from QHull

	    // qh_countfacets initializes 'visitid' in each facet.
	    // Note: In QHull parlance, a "facet" is a tet !!
	    //    (a tet is a 3D facet of a 4D convex hull)

	    int numfacets, numsimplicial, totneighbors, numridges,
		numcoplanar, numtricoplanars ;

	    qh_countfacets(
		qh facet_list, NULL, False,
		&numfacets, &numsimplicial,
		&totneighbors, &numridges, &numcoplanar, &numtricoplanars
	    ) ;

	    index_t cell_size = dimension + 1 ;
	    cell_to_v.resize(index_t(numfacets)*cell_size) ;
	    cell_to_cell.resize(index_t(numfacets)*cell_size) ;

	    index_t cur_cell_to_v = 0 ;
	    index_t cur_cell_to_cell = 0 ;

	    facetT* facet ;
	    FORALLfacets {
		if(!qh_skipfacet(facet)) {
		    vertexT* vertex ;
		    vertexT** vertexp ;

		    // Rem: setsize means 'size of the set'
		    // (and not 'set the size')
		    ogf_debug_assert(qh_setsize(facet->vertices) == cell_size) ;
		    FOREACHvertex_(facet->vertices) {
			cell_to_v[cur_cell_to_v] =
			    index_t(qh_pointid(vertex->point)) ;
			cur_cell_to_v++ ;
		    }
		    facetT* neighbor ;
		    facetT** neighborp ;

		    // Rem: setsize means 'size of the set'
		    // (and not 'set the size')
		    ogf_debug_assert(qh_setsize(facet->neighbors) == cell_size) ;
		    FOREACHneighbor_(facet) {
			cell_to_cell[cur_cell_to_cell] = neighbor->visitid - 1 ;
			cur_cell_to_cell++ ;
		    }
		}
	    }

	    // Step 3: Free QHull memory
	    qh_freeqhull(False); // <---- I'm unsure whether the
	                         // argument should be True or False here...
	    int curlong, totlong ;
	    qh_memfreeshort(&curlong, &totlong);
	    if (curlong || totlong) {
		Logger::warn("QHull")  << "did not free " << totlong
				       << " bytes of long memory ("
				       << curlong << " pieces)" << std::endl ;
	    }
	    ogf_assert(cell_to_v.size() == cell_to_cell.size()) ;
	}
	qh NOerrexit = True ;
    }
}
