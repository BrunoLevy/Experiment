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

#include <OGF/Experiment/commands/mesh_grob_calabi_yau_commands.h>

#ifdef EXPERIMENT_WITH_CGAL
#include <OGF/Experiment/CGAL/delaunay_nd_cgal.h>
#endif

#include <OGF/Experiment/algo/tiny_cas.h>
#include <geogram/mesh/index.h>
#include <geogram/basic/geometry_nd.h>
#include <complex>

namespace {
    using namespace OGF;

    typedef vecng<6,double> vec6;
    typedef Matrix<6,double> mat6;
    typedef vecng<7,double> vec7;
    typedef std::complex<double> C;
    typedef vecng<2,C> vec2C;

    class BiQuadric {
    public:

	// The complex polynomial
	C P(const C& Z1, const C& Z2) const {
	    return 1.0 + Z2*Z2 + Z1*Z1*Z2;
	}

	// The complex gradient of the complex polynomial
	vec2C gradP(const C& Z1, const C& Z2) const {
	    return vec2C{2.0*Z1*Z2, 2.0*Z2+Z1*Z1};
	}

	// The complex equation, simplest form, rational
	C eval0(const vec6& p) const {
	    double x1 = p[0]; double y1 = p[1]; double z1 = p[2];
	    double x2 = p[3]; double y2 = p[4]; double z2 = p[5];
	    C Z1 = C(x1,y1) / (1.0-z1);
	    C Z2 = C(x2,y2) / (1.0-z2);
	    return P(Z1,Z2);
	}

	// The complex equation, polynomial version
	C eval(const vec6& p) const {
	    double x1 = p[0]; double y1 = p[1]; double z1 = p[2];
	    double x2 = p[3]; double y2 = p[4]; double z2 = p[5];
	    return
		geo_sqr(1 - z1) * geo_sqr(1 - z2)              +
		geo_sqr(1 - z1) * geo_sqr(C(x2,y2)) +
		(1 - z2) * geo_sqr(C(x1,y1)) * C(x2,y2);
	}

	// The implicit equations of the manifold,
	// there are four of them (defines a 2-dim submanifold in IR^6)
	double F(const vec6& p, index_t i) const {
	    switch(i) {
	    case 0:
		return geo_sqr(p[0]) + geo_sqr(p[1]) + geo_sqr(p[2]) - 1.0;
	    case 1:
		return geo_sqr(p[3]) + geo_sqr(p[4]) + geo_sqr(p[5]) - 1.0;
	    case 2:
		return real(eval(p));
	    case 3:
		return imag(eval(p));
	    }
	    geo_assert_not_reached;
	}

	// The gradients of the 4 implicit equations of the manifold
	vec6 gradF(const vec6& p, index_t i) const {
	    double x1 = p[0]; double y1 = p[1]; double z1 = p[2];
	    double x2 = p[3]; double y2 = p[4]; double z2 = p[5];
	    switch(i) {
	    case 0:
		return grad_eval_real(p);
	    case 1:
		return grad_eval_imag(p);
	    case 2: {
		double g[6] = {2.0*x1, 2.0*y1, 2.0*z1, 0.0, 0.0, 0.0};
		return vec6(g);
	    }
	    case 3: {
		double g[6] = {0.0, 0.0, 0.0, 2.0*x2, 2.0*y2, 2.0*z2};
		return vec6(g);
	    }
	    }
	    geo_assert_not_reached;
	}

	// A basis of the 2D tangent space to the manifold
	std::pair<vec6, vec6> tangent_space(const vec6& p) {
	    mat6 M;
	    M.load_zero();
	    for(index_t i=0; i<4; ++i) {
		vec6 g = gradF(p,i);
		for(index_t j=0; j<6; ++j) {
		    M(i,j) = g[j];
		}
	    }

	    // minimum pivot for matrix inversion
	    static constexpr double eps = 1e-6;

	    // Compute U: orthogonal to g1,g2,g3,g4
	    vec6 U;
	    mat6 M1 = M;
	    // Loop until matrix is invertible
	    for(;;) {
		// Generate random coefficients in last two rows
		for(index_t i=4; i<6; ++i) {
		    for(index_t j=0; j<6; ++j) {
			M1(i,j) = 2.0*(Numeric::random_float64() - 0.5);
		    }
		}
		mat6 M1_inv;
		if(M1.compute_inverse(M1_inv, eps)) {
		    U = M1_inv * vec6({0.0, 0.0, 0.0, 0.0, 1.0, 1.0});
		    break;
		}
	    }

	    // Compute V: orthogonal to g1,g2,g3,g4,U
	    vec6 V;
	    mat6 M2 = M;
	    for(index_t j=0; j<6; ++j) {
		M2(4,j) = U[j];
	    }
	    // Loop until matrix is invertible
	    for(;;) {
		// Generate random coefficients in last row
		for(index_t j=0; j<6; ++j) {
		    M2(5,j) = 2.0*(Numeric::random_float64() - 0.5);
		}
		mat6 M2_inv;
		if(M2.compute_inverse(M2_inv, eps)) {
		    V = M2_inv * vec6({0.0, 0.0, 0.0, 0.0, 0.0, 1.0});
		    break;
		}
	    }

	    return std::make_pair(normalize(U),normalize(V));
	}

	/**************************************************************/


	// Written by Claude, does not seem to work ...
	vec2C project_onto_manifold(
	    const C& Z1, const C& Z2,
	    int maxIter = 100, double tol = 1e-12, double mu = 0.0 // 1e-8
	) const {
	    vec2C Z{Z1,Z2};
	    C lambda = 0.0;

	    for (int k = 0; k < maxIter; ++k) {
		C p = P(Z[0], Z[1]);
		if (std::abs(p) < tol) {
		    break;
		}

		vec2C g = gradP(Z[0], Z[1]);

		// ||grad||^2
		// (+ mu for Levenberg-Marquardt damping near singular points)
		double gNorm2 = std::norm(g[0]) + std::norm(g[1]) + mu;

		// grad . (target - Z)
		C gdotDelta = g[0] * (Z1 - Z[0]) + g[1] * (Z2 - Z[1]);

		// Newton update for the Lagrange multiplier
		lambda = -(p + gdotDelta) / gNorm2;

		// Z = target + lambda * conj(grad)
		Z[0] = Z1 + lambda * std::conj(g[0]);
		Z[1] = Z2 + lambda * std::conj(g[1]);
	    }
	    return Z;
	}

	vec6 project_onto_manifold(vec6 p) const {
	    double l1 = length(vec3{p[0],p[1],p[2]});
	    double l2 = length(vec3{p[3],p[4],p[5]});
	    p = vec6{
		p[0]/l1, p[1]/l1, p[2]/l1,
		p[3]/l2, p[4]/l2, p[5]/l2
	    };

	    // project onto the zero set of F (use expr in eval0,
	    // use parameterization in terms of Z1,Z2...)

	    double x1 = p[0]; double y1 = p[1]; double z1 = p[2];
	    double x2 = p[3]; double y2 = p[4]; double z2 = p[5];

	    // Projection onto manifold using formula0 and gradient

	    C Z1 = C(x1,y1) / (1.0-z1);
	    C Z2 = C(x2,y2) / (1.0-z2);

	    vec2C ZZ = project_onto_manifold(Z1,Z2);
	    Z1 = ZZ[0];
	    Z2 = ZZ[1];

	    // Inverse stereographic projection
	    double X1 = real(Z1); double Y1 = imag(Z1);
	    double X2 = real(Z2); double Y2 = imag(Z2);
	    double L1 = (1.0 + X1*X1 + Y1*Y1);
	    double L2 = (1.0 + X2*X2 + Y2*Y2);
	    p = vec6{
		2.0*X1 / L1, 2.0*Y1 / L1, (L1-2.0) / L1,
		2.0*X2 / L2, 2.0*Y2 / L2, (L2-2.0) / L2
	    };

	    // Tiny perturbation to avoid cocyclic points that make
	    // Delaunay triangulation slow
	    l1 = 1.0 + 1e-6*(Numeric::random_float64()-0.5);
	    l2 = 1.0 + 1e-6*(Numeric::random_float64()-0.5);
	    p = vec6{
		p[0]/l1, p[1]/l1, p[2]/l1,
		p[3]/l2, p[4]/l2, p[5]/l2
	    };
	    return p;
	}

    protected:
	// gradient of the real part of the polynomial complex eqn
	// Generated using TinyCAS, see equation_gradients()
	vec6 grad_eval_real(const vec6& p) const {
	    double x1 = p[0]; double y1 = p[1]; double z1 = p[2];
	    double x2 = p[3]; double y2 = p[4]; double z2 = p[5];
	    return vec6({
		+2*x1*x2 -2*x1*x2*z2 -2*y1*y2 +2*y1*y2*z2,
		-2*y1*x2 +2*y1*x2*z2 -2*x1*y2 +2*x1*y2*z2,
		-2 +4*z2 -2*z2*z2 +2*z1 -4*z1*z2 +2*z1*z2*z2
		                   -2*x2*x2 +2*y2*y2 +2*z1*x2*x2 -2*z1*y2*y2,
		+2*x2 -4*z1*x2 +2*z1*z1*x2 +x1*x1 -y1*y1 -x1*x1*z2 +y1*y1*z2,
		-2*y2 +4*z1*y2 -2*z1*z1*y2 -2*x1*y1 +2*x1*y1*z2,
		-2 +2*z2 +4*z1 -4*z1*z2 -2*z1*z1 +2*z1*z1*z2
		                   -x1*x1*x2 +y1*y1*x2 +2*x1*y1*y2
	    });
	}

	// gradient of the imaginary part of the polynomial complex eqn
	// Generated using TinyCAS, see equation_gradients()
	vec6 grad_eval_imag(const vec6& p) const {
	    double x1 = p[0]; double y1 = p[1]; double z1 = p[2];
	    double x2 = p[3]; double y2 = p[4]; double z2 = p[5];
	    return vec6({
		+2*x1*y2 -2*x1*y2*z2 +2*y1*x2 -2*y1*x2*z2,
		-2*y1*y2 +2*y1*y2*z2 +2*x1*x2 -2*x1*x2*z2,
		-4*x2*y2 +4*z1*x2*y2,
		-4*z1*y2 +2*z1*z1*y2 +2*x1*y1 -2*x1*y1*z2,
		+2*x2 -4*z1*x2 +2*z1*z1*x2 +x1*x1 -y1*y1 -x1*x1*z2 +y1*y1*z2,
		-x1*x1*y2 +y1*y1*y2 -2*x1*y1*x2
	    });
	}
    };


    // Chop chop chop !
    template <unsigned int DIM> inline void clip_polygon(
	vector<vecng<DIM,double>>& P,
	vector<index_t>& Pid,
	const vecng<DIM+1,double>& pi,
	index_t piid
    ) {
	typedef vecng<DIM,double> vecd;
	vecd n(pi.data());
	double alpha = pi[DIM];

	vector<vecng<DIM,double>> P2;
	vector<index_t> P2id;

	for(index_t i=0; i<P.size(); ++i) {
	    const vecd& p = P[i];
	    const vecd& q = P[(i+1)%P.size()];
	    index_t pid = Pid[i];

	    Sign sp = geo_sgn(dot(p,n) + alpha);
	    Sign sq = geo_sgn(dot(q,n) + alpha);

	    if(sp >= 0) {
		P2.push_back(p);
		P2id.push_back(pid);
	    }
	    if(sp * sq < 0) {
		double t = - (dot(p,n) + alpha) / dot(q-p,n);
		P2.push_back(p + t*(q-p));
		P2id.push_back(sp > 0 ? piid : pid);
	    }
	}
	std::swap(P,P2);
	std::swap(Pid, P2id);
    }

    // Chop chop chop !
    template <unsigned int DIM> inline vecng<DIM+1,double> bisector(
	const vecng<DIM,double>& p,
	const vecng<DIM,double>& q
    ) {
	vecng<DIM,double> n = p - q;
	vecng<DIM+1,double> pi;
	for(index_t i=0; i<DIM; ++i) {
	    pi[i] = n[i];
	}
	pi[DIM] = -0.5*dot(p+q,n);
	return pi;
    }

    template <unsigned int DIM> inline vecng<DIM,double> centroid(
	const vector<vecng<DIM,double>>& P
    ) {
	double sA = 0.0;
	vecng<DIM,double> result;
	const vecng<DIM,double>& p = P[0];
	for(index_t i=1; i+1<P.size(); ++i) {
	    const vecng<DIM,double>& q = P[i];
	    const vecng<DIM,double>& r = P[i+1];
	    double A = GEO::Geom::triangle_area(p,r,q);
	    result += (A/3.0)*(p+q+r);
	    sA += A;
	}
	return (1.0/sA)*result;
    }

}


namespace OGF {

    MeshGrobCalabiYauCommands::MeshGrobCalabiYauCommands() {
    }

    MeshGrobCalabiYauCommands::~MeshGrobCalabiYauCommands() {
    }

    void MeshGrobCalabiYauCommands::remove_points(index_t nb_points_to_keep) {
	vector<index_t> remove_points(mesh_grob()->vertices.nb());
	for(index_t v: mesh_grob()->vertices) {
	    remove_points[v] = index_t(v >= nb_points_to_keep);
	}
	mesh_grob()->vertices.delete_elements(remove_points);
	mesh_grob()->update();
    }

    void MeshGrobCalabiYauCommands::show_tangent_basis(double size) {
	BiQuadric Q;
	index_t nb_v = mesh_grob()->vertices.nb();
	for(index_t v=0; v<nb_v; ++v) {
	    vec6 p = mesh_grob()->vertices.point<6>(v);
	    auto [U,V] = Q.tangent_space(p);
	    index_t first_v = mesh_grob()->vertices.create_vertices(4);
	    index_t v1 = first_v;
	    index_t v2 = first_v+1;
	    index_t v3 = first_v+2;
	    index_t v4 = first_v+3;
	    mesh_grob()->vertices.point<6>(v1) = p + size * U;
	    mesh_grob()->vertices.point<6>(v2) = p - size * U;
	    mesh_grob()->vertices.point<6>(v3) = p + size * V;
	    mesh_grob()->vertices.point<6>(v4) = p - size * V;
	    mesh_grob()->edges.create_edge(v,v1);
	    mesh_grob()->edges.create_edge(v,v2);
	    mesh_grob()->edges.create_edge(v,v3);
	    mesh_grob()->edges.create_edge(v,v4);
	}
	mesh_grob()->update();
    }

    void MeshGrobCalabiYauCommands::tangent_RVD(
	double r, bool three_only, index_t nb_CVT
    ) {

	if(nb_CVT > 1) {
	    for(index_t i=0; i<nb_CVT; ++i) {
		tangent_RVD(r,three_only,1);
		project_onto_manifold();
	    }
	    return;
	}

#ifdef EXPERIMENT_WITH_CGAL
	MeshGrob* RVD = MeshGrob::find_or_create(scene_graph(), "RVD");
	RVD->clear();
	RVD->vertices.set_dimension(6);

	// Step 1: compute Delaunay
	index_t dim = mesh_grob()->vertices.dimension();
	Delaunay_var del = new DelaunayNd_CGAL(coord_index_t(dim));
	del->set_stores_cicl(true);
	Logger::out("CGAL") << "Computing Delaunay" << std::endl;
	del->set_vertices(
	    mesh_grob()->vertices.nb(), mesh_grob()->vertices.point_ptr(0)
	);
	Logger::out("CGAL") << "Delaunay done" << std::endl;

	vector<trindex> triangles;

	// Step 2: compute clipped tangent planes
	BiQuadric Q;
	index_t nb_v = mesh_grob()->vertices.nb();
	vector<vec6> new_points;
	if(nb_CVT == 1) {
	    new_points.resize(nb_v);
	}

	vector<index_t> infty_count_histo(5,0);
	for(index_t v=0; v<nb_v; ++v) {
	    vec6 p = mesh_grob()->vertices.point<6>(v);
	    auto [U,V] = Q.tangent_space(p);
	    vector<vec6> cell;
	    cell.push_back(p-r*U-r*V);
	    cell.push_back(p-r*U+r*V);
	    cell.push_back(p+r*U+r*V);
	    cell.push_back(p+r*U-r*V);
	    vector<index_t> cell_id(4,NO_INDEX);

	    vector<index_t> Nv;
	    del->get_neighbors(v,Nv);
	    sort_unique(Nv);
	    for(index_t w: Nv) {
		vec6 q = mesh_grob()->vertices.point<6>(w);
		clip_polygon(cell, cell_id, bisector(p,q), w);
	    }

	    index_t infty_count = 0;
	    for(index_t lv=0; lv<cell.size(); ++lv) {
		infty_count += (cell_id[lv] == NO_INDEX);
	    }
	    geo_assert(infty_count <= 4);
	    infty_count_histo[infty_count]++;

	    if(nb_CVT == 1) {
		new_points[v] = (infty_count == 0) ? centroid(cell) : p;
	    }

	    for(index_t lv=0; lv<cell_id.size(); ++lv) {
		index_t i = v;
		index_t j = cell_id[lv];
		index_t k = cell_id[(lv+1)%cell_id.size()];
		if(i != NO_INDEX && j != NO_INDEX && k != NO_INDEX) {
		    triangles.push_back(trindex(i,j,k));
		}
	    }

	    if(infty_count == 0) {
		index_t new_v = RVD->vertices.create_vertices(cell.size());
		index_t new_f = RVD->facets.create_polygon(cell.size());
		for(index_t lv=0; lv<cell.size(); ++lv) {
		    RVD->vertices.point<6>(new_v+lv) = cell[lv];
		    RVD->facets.set_vertex(new_f, lv, new_v+lv);
		}
	    }
	}
	std::cerr << "Infinite cells:" << std::endl;
	for(index_t i=0; i<=4; ++i) {
	    std::cerr << i << ":" << infty_count_histo[i] << std::endl;
	}

	vector<index_t> trgl_count_histo;
	std::sort(triangles.begin(), triangles.end());
	vector<index_t> triangles2;
	triangles2.reserve(triangles.size()*3/3);
	auto b = triangles.begin();
	while(b != triangles.end()) {
	    auto e = b;
	    while(e != triangles.end() && *e == *b) {
		++e;
	    }
	    index_t trgl_count = index_t(e-b);
	    if(trgl_count >= trgl_count_histo.size()) {
		trgl_count_histo.resize(trgl_count+1, 0);
	    }
	    trgl_count_histo[trgl_count]++;
	    if(!three_only || trgl_count == 3) {
		triangles2.push_back(b->indices[0]);
		triangles2.push_back(b->indices[1]);
		triangles2.push_back(b->indices[2]);
	    }
	    b=e;
	}
	mesh_grob()->facets.assign_triangle_mesh(triangles2,true);


	if(nb_CVT == 1) {
	    for(index_t v: mesh_grob()->vertices) {
		mesh_grob()->vertices.point<6>(v) = new_points[v];
	    }
	}

	std::cerr << "Triangles multiplicities (3 is normal):" << std::endl;
	for(index_t i=0; i<trgl_count_histo.size(); ++i) {
	    std::cerr << i << ":" << trgl_count_histo[i] << std::endl;
	}

	mesh_grob()->update();
	RVD->update();

#else
	geo_argused(r);
	geo_argused(three_only);
	geo_argused(nb_CVT);
	Logger::err("CGAL") << "Not compiled with CGAL support" << std::endl;
#endif
    }

    void MeshGrobCalabiYauCommands::compute_Delaunay_highdim(
	index_t nb_pts, bool build_graph
    ) {
#ifdef EXPERIMENT_WITH_CGAL
	index_t dim = mesh_grob()->vertices.dimension();
	Delaunay_var del = new DelaunayNd_CGAL(coord_index_t(dim));
	del->set_stores_cicl(build_graph);
	if(nb_pts == 0) {
	    nb_pts = mesh_grob()->vertices.nb();
	}
	Logger::out("CGAL") << "Computing Delaunay" << std::endl;
	del->set_vertices(
	    nb_pts, mesh_grob()->vertices.point_ptr(0)
	);
	Logger::out("CGAL") << "Delaunay done" << std::endl;
	if(build_graph) {
	    vector<index_t> N;
	    double avg_arity = 0.0;
	    for(index_t v=0; v<del->nb_vertices(); ++v) {
		del->get_neighbors(v,N);
		avg_arity += double(N.size());
	    }
	    avg_arity /= double(del->nb_vertices());
	    Logger::out("CGAL") << "Average arity:" << avg_arity << std::endl;
	}
#else
	geo_argused(nb_pts);
	geo_argused(build_graph);
	Logger::err("CGAL") << "Not compiled with CGAL support" << std::endl;
#endif
    }

    void MeshGrobCalabiYauCommands::equation_check() {
	BiQuadric Q;

	for(index_t v: mesh_grob()->vertices) {
	    const vec6& p = mesh_grob()->vertices.point<6>(v);
	    C f1 = Q.eval0(p);
	    C f2 = Q.eval(p);
	    std::cerr << "f1=" << f1 << " |f1|=" << std::abs(f1)
		      << " f2=" << f2 << std::endl;
	}
    }

    // Z1 = (x1 + iy1) / (1-z1)
    // Z2 = (x2 + iy2) / (1-z2)
    //
    // F1(x1,y1,z1,x2,y2,z2) = || [x1,y1,z1] ||^2 - 1 -> GF3 = 2 [x1,y1,z1,0,0,0]
    // F2(x1,y1,z1,x2,y2,z2) = || [x2,y2,z2] ||^2 - 1 -> GF4 = 2 [0,0,0,x2,y2,z2]
    // P = 1 + Z2^2 + Z1^2*Z2
    // consider Q = (1-z1)^2(1-z2)^2 P (same denominator for the terms of P)
    // Q = (1-z1)^2(1-z2)^2 + (1-z1)^2(x2+iy2)^2 + (1-z2)(x1+iy1)^2(x2+iy2)
    // F3(x1,y1,z1,x2,y2,z2) = Re(Q(Z1,Z2))
    // F4(x1,y1,z1,x2,y2,z2) = Im(Q(Z1,Z2))

    // TODO
    // DONE: Check equation of Q on data
    // DONE: Functions to compute derivatives of F1,F2,F3,F4
    // DONE: Function to compute tangent space basis
    // DONE: Visualize tangent space vectors
    // DONE: Sutherland-Hogdman re-entrant clipping in 6D

    void MeshGrobCalabiYauCommands::equation_gradients() {
	using ::OGF::TinyCAS::Polynom;

	typedef Polynom R;
	typedef std::complex<R> C;

	std::vector<std::string> var_names = {
	    "x1","y1","z1",
	    "x2","y2","z2"
	};

	index_t nb_var = index_t(var_names.size());
	R Rx1(nb_var,0u);
	R Ry1(nb_var,1u);
	R Rz1(nb_var,2u);
	R Rx2(nb_var,3u);
	R Ry2(nb_var,4u);
	R Rz2(nb_var,5u);

	R Rzero(nb_var,0.0);
	R Rone(nb_var,1.0);

	C one(Rone, Rzero);
	C i(Rzero, Rone);
	C x1(Rx1, Rzero);
	C y1(Ry1, Rzero);
	C z1(Rz1, Rzero);
	C x2(Rx2, Rzero);
	C y2(Ry2, Rzero);
	C z2(Rz2, Rzero);

	auto sq = [](const C& Z)->C {
	    return Z*Z;
	};

	C Q =
	    sq(one - z1) * sq(one - z2)              +
	    sq(one - z1) * sq(x2 + i*y2)             +
	      (one - z2) * sq(x1 + i*y1) * (x2+i*y2) ;

	std::cerr
	    << "F3=" << real(Q).to_string(var_names) << std::endl
	    << "F4=" << imag(Q).to_string(var_names) << std::endl;

	for(index_t i=0; i<6; ++i) {
	    std::cerr << "der F3_" << var_names[i] << "="
		      << real(Q).derivative(i).to_string(var_names)
		      << std::endl;
	}

	for(index_t i=0; i<6; ++i) {
	    std::cerr << "der F4_" << var_names[i] << "="
		      << imag(Q).derivative(i).to_string(var_names)
		      << std::endl;
	}

    }

    void MeshGrobCalabiYauCommands::export_tet6(const std::string& filename) {
	std::ofstream out(filename.c_str());
	if(!out) {
	    Logger::err("CY") << "Could not open " << filename << std::endl;
	    return;
	}
	out << mesh_grob()->vertices.nb() << " vertices" << std::endl;
	out << mesh_grob()->facets.nb() << " trgls" << std::endl;
	for(index_t v: mesh_grob()->vertices) {
	    out << mesh_grob()->vertices.point<6>(v) << std::endl;
	}
	for(index_t t: mesh_grob()->facets) {
	    out	<< mesh_grob()->facets.vertex(t,0) << " "
		<< mesh_grob()->facets.vertex(t,1) << " "
		<< mesh_grob()->facets.vertex(t,2)
		<< std::endl;
	}
    }

    void MeshGrobCalabiYauCommands::project_onto_manifold() {
	BiQuadric Q;
	for(index_t v: mesh_grob()->vertices) {
	    vec6 p = mesh_grob()->vertices.point<6>(v);
	    mesh_grob()->vertices.point<6>(v) = Q.project_onto_manifold(p);
	}
	mesh_grob()->update();
    }

}
