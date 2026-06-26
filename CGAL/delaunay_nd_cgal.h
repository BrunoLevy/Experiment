#ifndef H__OGF_CGALND_DELAUNAY_ND_CGAL__H
#define H__OGF_CGALND_DELAUNAY_ND_CGAL__H

#include <OGF/Experiment/common/common.h>
#include <geogram/delaunay/delaunay.h>

namespace OGF {

    // TODO: make it a template of dimension, and register
    // several dimensions.
    class DelaunayNd_CGAL : public Delaunay {
      public:
	DelaunayNd_CGAL(coord_index_t dim);
	~DelaunayNd_CGAL() override;
	void set_vertices(index_t nb_vertices, const double* vertices) override;
      protected:
	std::vector<index_t> cell_to_v_ ;
	std::vector<index_t> cell_to_cell_ ;
    } ;

}

#endif
