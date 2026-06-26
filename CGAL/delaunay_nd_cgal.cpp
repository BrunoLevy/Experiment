/*
  Note: CGAL example for nD Delaunay triangulation:
  https://doc.cgal.org/latest/Triangulation/delaunay_triangulation_8cpp-example.html
*/

#include <OGF/Experiment/CGAL/delaunay_nd_cgal.h>
#define CGAL_EIGEN3_ENABLED
#include <CGAL/Epick_d.h>
#include <CGAL/Delaunay_triangulation.h>
#include <CGAL/Hilbert_sort_d.h>


namespace {
    using namespace GEO;

    template <int DIM> class DelaunayNd_CGAL_Impl {
	// Let us build a nD Delaunay triangulation with an index_t
	// additional field in the vertices and simplices

	// The kernel: Exact predicates inexact constructions dim d
        typedef CGAL::Epick_d< CGAL::Dimension_tag<DIM> > K;

	// The vertices and simplices (full cell = max dim simplex)
	// Additional template argument is for storing indices, that
	// can be set/get through the data() function
	//                                        |
	//                                        v
	typedef CGAL::Triangulation_vertex   <K,index_t> TV;
	typedef CGAL::Triangulation_full_cell<K,index_t> TC;

	// The triangulation data structure stores the combinatorics
	typedef CGAL::Triangulation_data_structure<
		CGAL::Dimension_tag<DIM>, TV, TC
	> TDS;

	// And, finally, here is our nD Delaunay triangulation
	typedef CGAL::Delaunay_triangulation<K,TDS> DT;

	// Some shorthands
	typedef typename K::Point_d Point;
	typedef typename DT::Vertex_handle Vertex_handle;
	typedef typename DT::Full_cell_handle Full_cell_handle;

	// An nD point with an Id. We need to store Ids with the points
	// because we will sort them, and we need to keep relation with
	// initial indices.
	class PointWithId : public Point {
	public:
	    PointWithId(
		const double* begin, const double* end, index_t id_in
	    ) : Point(begin, end), id(id_in) {
	    }
	    PointWithId(const PointWithId& rhs) = default;
	    PointWithId& operator=(const PointWithId& rhs) = default;
	    index_t id;
	};

    public:
	static unsigned int compute(
	    index_t nb_vertices, const double* vertices,
	    std::vector<index_t>& cell_to_v, std::vector<index_t>& cell_to_cell
	) {
	    DT dt(DIM);

	    // --------------------------- Step 1: Hilbert sort

	    std::vector<DelaunayNd_CGAL_Impl::PointWithId> points ;
	    points.reserve(nb_vertices) ;
	    for(unsigned int i=0; i<nb_vertices; i++) {
		points.emplace_back(vertices + i*DIM, vertices + (i+1)*DIM, i);
	    }

	    std::cerr << "Hilbert sort begin" << std::endl ;

	    spatial_sort(
		points.begin(), points.end(), dt.geom_traits(),
		CGAL::Hilbert_sort_median_policy()
	    );

	    std::cerr << "Hilbert sort end" << std::endl ;

	    // --------------------------- Step 2: insert vertices

	    std::cerr << "Delaunay begin" << std::endl ;

	    index_t cnt = 0;
	    Full_cell_handle hint;
	    for(const auto& p: points) {
		std::cerr << "." << std::flush;
		Vertex_handle vx = dt.insert(p,hint);
		vx->data() = p.id;
		hint = vx->full_cell();
		if((cnt % 100) == 0) {
		    std::cerr << cnt;
		}
		++cnt;
	    }

	    std::cerr << std::endl << "Delaunay end" << std::endl ;

	    // --------------------------- Step 3: decode CGAL structure

	    std::cerr << "Decode begin" << std::endl ;

	    index_t cur = 0 ;
	    for(auto it=dt.full_cells_begin(); it!=dt.full_cells_end(); it++) {
		if(dt.is_infinite(it)) {
		    it->data() = NO_INDEX;
		} else {
		    it->data() = cur;
		    cur++ ;
		}
	    }
	    unsigned int nb_cells = cur ;
	    cell_to_v.resize(0);
	    cell_to_v.reserve(nb_cells*(DIM+1)) ;
	    cell_to_cell.resize(0) ;
	    cell_to_cell.reserve(nb_cells*(DIM+1)) ;
	    for(auto it=dt.full_cells_begin(); it!=dt.full_cells_end(); it++) {
		if(!dt.is_infinite(it)) {
		    for(int j=0; j<=DIM; j++) {
			cell_to_v.push_back(it->vertex(j)->data());
		    }
		    for(int j=0; j<=DIM; j++) {
			cell_to_cell.push_back(it->neighbor(j)->data());
		    }
		}
	    }
	    std::cerr << "Decode end" << std::endl ;
	    std::cerr << "Nb cells = " << nb_cells << std::endl;
	    return nb_cells ;
	}
    };
}

namespace OGF {
    DelaunayNd_CGAL::DelaunayNd_CGAL(coord_index_t dim) : Delaunay(dim) {
	ogf_assert(dim == 4 || dim == 5 || dim == 6) ;
    }

    DelaunayNd_CGAL::~DelaunayNd_CGAL() {
    }

    void DelaunayNd_CGAL::set_vertices(
	index_t nb_vertices, const double* vertices
    ) {
	Delaunay::set_vertices(nb_vertices, vertices) ;
	unsigned int nb_cells = 0 ;
	switch(dimension()) {
	    default:
		geo_assert_not_reached;
	    break;
	    case 4:
		nb_cells = DelaunayNd_CGAL_Impl<4>::compute(
		    nb_vertices, vertices, cell_to_v_, cell_to_cell_
		) ;
		break ;
	    case 5:
		nb_cells = DelaunayNd_CGAL_Impl<5>::compute(
		    nb_vertices, vertices, cell_to_v_, cell_to_cell_
		) ;
		break ;
	    case 6:
		nb_cells = DelaunayNd_CGAL_Impl<6>::compute(
		    nb_vertices, vertices, cell_to_v_, cell_to_cell_
		) ;
		break ;
	}
	set_arrays(nb_cells, cell_to_v_.data(), cell_to_cell_.data()) ;
    }
}
