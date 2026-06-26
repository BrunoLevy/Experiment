/*
  New version:
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
        typedef CGAL::Epick_d< CGAL::Dimension_tag<DIM> > K;

	// ------------- store ids ---------------.
	//                                        v
	typedef CGAL::Triangulation_vertex   <K,index_t> TV;
	typedef CGAL::Triangulation_full_cell<K,index_t> TC;

	typedef CGAL::Triangulation_data_structure<
		CGAL::Dimension_tag<DIM>, TV, TC
	> TDS;

	typedef CGAL::Delaunay_triangulation<K,TDS> DT;
	typedef typename K::Point_d Point ;
	typedef typename DT::Vertex_handle Vertex_handle;
	typedef typename DT::Full_cell_handle Full_cell_handle;

	class PointId : public Point{
	public:
	    PointId(
		const double* begin, const double* end, index_t id_in
	    ) : Point(begin, end), id(id_in) {
	    }
	    PointId(const PointId& rhs) = default;
	    PointId& operator=(const PointId& rhs) = default;
	    index_t id;
	};

    public:
	static unsigned int compute(
	    index_t nb_vertices, const double* vertices,
	    std::vector<index_t>& cell_to_v, std::vector<index_t>& cell_to_cell
	) {
	    DT dt(DIM);

	    // --------------------------- Step 1: Hilbert sort

	    std::vector<DelaunayNd_CGAL_Impl::PointId> points ;
	    points.reserve(nb_vertices) ;
	    for(unsigned int i=0; i<nb_vertices; i++) {
		points.push_back(
		     PointId(vertices + i*DIM, vertices + (i+1)*DIM, i)
		) ;
	    }

	    std::cerr << "Hilbert sort begin" << std::endl ;

	    spatial_sort(
		points.begin(), points.end(), dt.geom_traits(),
		CGAL::Hilbert_sort_median_policy()
	    );

	    std::cerr << "Hilbert sort end" << std::endl ;

	    // --------------------------- Step 2: insert vertices

	    std::cerr << "Delaunay begin" << std::endl ;

	    Full_cell_handle hint;
	    for(const auto& p: points) {
		std::cerr << "." << std::flush;
		Vertex_handle vx = dt.insert(p,hint);
		vx->data() = p.id;
		hint = vx->full_cell();
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

/*
#include <CGAL/Pure_complex_data_structure.h>
#include <CGAL/Simple_cartesian_d.h>
#include <CGAL/Filtered_kernel_d.h>
#include <CGAL/Delaunay_complex.h>
#include <CGAL/Hilbert_sort_d.h>



namespace OGF {

    template <int DIM> class DelaunayNd_CGAL_Impl {
	typedef CGAL::Simple_cartesian_d<double,DIM> K ;
	typedef CGAL::Filtered_kernel_d<K> FK ;
	typedef CGAL::Delaunay_complex<FK> PC ;
	typedef typename PC::Point_d Point ;

	class PointId : public Point {
	public:
	    PointId(
		const double* begin, const double* end, int id_in
	    ) : Point(begin, end), id(id_in) {
	    }
	    PointId(const PointId& rhs) : Point((const Point&)rhs), id(rhs.id) {}
	    PointId& operator=(const PointId& rhs) { Point::operator=(rhs) ; id = rhs.id ; }
	    int id ;
	} ;

        public:
	static unsigned int compute(
	    unsigned int nb_vertices, const double* vertices,
	    std::vector<int>& cell_to_v,
	    std::vector<int>& cell_to_cell
	) {

	    // --------------------------- Step 1: Hilbert sort

	    std::cerr << "Hilbert sort begin" << std::endl ;

	    std::vector<DelaunayNd_CGAL_Impl::PointId> points ;
	    points.reserve(nb_vertices) ;
	    for(unsigned int i=0; i<nb_vertices; i++) {
		points.push_back(
		    DelaunayNd_CGAL_Impl::PointId(
			vertices + i*DIM, vertices + (i+1)*DIM, i
		     )
		) ;
	    }

	    CGAL::Hilbert_sort_d<DelaunayNd_CGAL_Impl::FK> hilbert_sort(
		DelaunayNd_CGAL_Impl::FK(),DIM
	    ) ;
	    hilbert_sort(points.begin(), points.end()) ;

	    std::cerr << "Hilbert sort end" << std::endl ;

	    // --------------------------- Step 2: insert vertices

	    std::cerr << "Delaunay begin" << std::endl ;

	    PC pc(DIM) ;
	    typename PC::Simplex_handle s ;
	    for(unsigned int i=0; i<nb_vertices; i++) {
		typename PC::Vertex_handle v = pc.insert(points[i],s) ;
		v->data().id = points[i].id ;
		s = v->simplex() ;
	    }

	    std::cerr << "Delaunay end" << std::endl ;

	    // --------------------------- Step 3: decode CGAL structure

	    std::cerr << "Decode begin" << std::endl ;

	    int cur = 0 ;
	    for(typename PC::Simplex_iterator it = pc.simplices_begin();
		it != pc.simplices_end(); it++) {
		if(pc.is_infinite(it)) {
		    it->data().id = -1 ;
		} else {
		    it->data().id = cur ; cur++ ;
		}
	    }
	    unsigned int nb_cells = cur ;
	    cell_to_v.resize(0);
	    cell_to_v.reserve(nb_cells*(DIM+1)) ;
	    cell_to_cell.resize(0) ;
	    cell_to_cell.reserve(nb_cells*(DIM+1)) ;
	    for(
		typename PC::Simplex_iterator it = pc.simplices_begin();
		it != pc.simplices_end(); it++
            ) {
		if(it->data().id != -1) {
		    for(unsigned int j=0; j<=DIM; j++) {
			cell_to_v.push_back(it->vertex(j)->data().id) ;
		    }
		    for(unsigned int j=0; j<=DIM; j++) {
			cell_to_cell.push_back(it->neighbor(j)->data().id) ;
		    }
		}
	    }
	    std::cerr << "Decode end" << std::endl ;
	    return nb_cells ;
	}
    } ;

    DelaunayNd_CGAL::DelaunayNd_CGAL(int dim) : Delaunay(dim) {
	ogf_assert(dim == 4 || dim == 5 || dim == 6) ;
    }

    DelaunayNd_CGAL::~DelaunayNd_CGAL() {
    }

    void DelaunayNd_CGAL::set_vertices(
	unsigned int nb_vertices, const double* vertices
    ) {
	Delaunay::set_vertices(nb_vertices, vertices) ;
	unsigned int nb_cells = 0 ;
	switch(dimension()) {
	    case 4:
		DelaunayNd_CGAL_Impl<4>::compute(
		    nb_vertices, vertices, cell_to_v_, cell_to_cell_
		) ;
		break ;
	    case 5:
		DelaunayNd_CGAL_Impl<5>::compute(
		    nb_vertices, vertices, cell_to_v_, cell_to_cell_
		) ;
		break ;
	    case 6:
		nb_cells = DelaunayNd_CGAL_Impl<6>::compute(
		    nb_vertices, vertices, cell_to_v_, cell_to_cell_
		) ;
		break ;
	}

	set_arrays(nb_cells, &(cell_to_v_[0]), &(cell_to_cell_[0])) ;
    }
}
*/
