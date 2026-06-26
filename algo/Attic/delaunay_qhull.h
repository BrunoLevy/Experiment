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


#ifndef OGF_CELLS_TRI_DELAUNAY_QHULL
#define OGF_CELLS_TRI_DELAUNAY_QHULL

#include <OGF/Experiment/common/common.h>
#include <geogram/delaunay/delaunay.h>

namespace OGF {

    //-------------------------------------------------------------------------

    /**
     * Computes convex hull (and Delaunay triangulation) in arbitrary dimension.
     */
    class Experiment_API QHull {
      public:
	QHull() ;
	~QHull() ;
	void compute(
	    index_t dimension,
	    index_t nb_vertices, const double* vertices,
	    std::vector<index_t>& cell_to_v, std::vector<index_t>& cell_to_cell
	) ;
    } ;

    //--------------------------------------------------------------------------

    class Experiment_API DelaunayNd_QHull : public Delaunay {
      public:
	DelaunayNd_QHull(index_t dimension) ;
	virtual ~DelaunayNd_QHull() ;
	virtual void set_vertices(
	    index_t nb_vertices, const double* vertices
	) ;
      protected:
	QHull qhull_ ;
	std::vector<index_t> cell_to_v_ ;
	std::vector<index_t> cell_to_cell_ ;
    } ;

    //--------------------------------------------------------------------------

    class Experiment_API Delaunay2d_QHull : public DelaunayNd_QHull {
      public:
        Delaunay2d_QHull() : DelaunayNd_QHull(2) { }
    } ;

    //---------------------------------------------------------------------------

    class Experiment_API Delaunay3d_QHull : public DelaunayNd_QHull {
      public:
        Delaunay3d_QHull() : DelaunayNd_QHull(3) { }
    } ;

    //--------------------------------------------------------------------------

}

#endif
