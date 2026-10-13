#ifndef __MeshCurvatureQuery_h__
#define __MeshCurvatureQuery_h__

#include <Discregrid/mesh/triangle_mesh.hpp>

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <igl/AABB.h>
#include <igl/boundary_facets.h>
#include <igl/cotmatrix.h>
#include <igl/fast_winding_number.h>
#include <igl/invert_diag.h>
#include <igl/massmatrix.h>
#include <igl/per_vertex_normals.h>

#include <algorithm>
#include <cmath>
#include <vector>

namespace SPH
{
	/** \brief Closest point, signed distance and curvature queries for triangle meshes
	* which are used to generate tube maps (see SimulatorBase::initTubeMap()).
	*
	* Closest points are determined by an AABB tree, the sign of the distance by fast
	* winding numbers, so the mesh is not required to be watertight. Per-vertex curvatures
	* are precomputed (both normalized by the Voronoi area):
	* - mean curvature H: projection of the cotangent Laplacian of the vertex positions onto the
	*   area-weighted vertex normal. The convention of the tube maps is used, i.e., H = -(k1 + k2)
	*   which is negative for convex surfaces (sphere of radius R: H = -2/R).
	* - Gaussian curvature K = k1 * k2: angle defect.
	* The faces of closed meshes may be oriented either way. Open meshes must be oriented such that
	* the face normals point to the fluid side.
	*/
	class MeshCurvatureQuery
	{
	public:
		explicit MeshCurvatureQuery(const Discregrid::TriangleMesh &mesh)
		{
			m_V.resize((int)mesh.nVertices(), 3);
			for (int i = 0; i < (int)mesh.nVertices(); ++i)
				m_V.row(i) = mesh.vertex(i).transpose();

			m_F.resize((int)mesh.nFaces(), 3);
			for (int f = 0; f < (int)mesh.nFaces(); ++f)
			{
				m_F(f, 0) = (int)mesh.faceVertex(f, 0);
				m_F(f, 1) = (int)mesh.faceVertex(f, 1);
				m_F(f, 2) = (int)mesh.faceVertex(f, 2);
			}

			m_tree.init(m_V, m_F);
			igl::fast_winding_number(m_V.cast<float>().eval(), m_F, 2, m_fwnBVH);

			computeVertexCurvatures();
		}

		/** Determine the signed distance sd (negative inside), the closest point cp and the index
		* of the closest triangle. Only triangles closer than maxDistance are considered. If there
		* is no such triangle, |sd| = maxDistance and triangleIndex = -1.
		*/
		void signedDistance(const Eigen::Vector3d &x, const double maxDistance,
			int &triangleIndex, Eigen::Vector3d &cp, double &sd) const
		{
			Eigen::RowVector3d cpRow;
			const Eigen::RowVector3d xRow = x.transpose();
			const double sqrDist = m_tree.squared_distance(m_V, m_F, xRow, maxDistance * maxDistance, triangleIndex, cpRow);
			const double dist = std::sqrt(std::max(0.0, sqrDist));
			cp = cpRow.transpose();

			// winding number: inside ~ 1, outside ~ 0
			const double w = igl::fast_winding_number(m_fwnBVH, 2, xRow.cast<float>().eval());
			sd = (1.0 - 2.0 * std::abs(w)) * dist;
		}

		/** Compute the barycentric coordinates (w0, w1, w2) of the point p (e.g. a closest point) w.r.t. face f. */
		void barycentricCoordinates(const int f, const Eigen::Vector3d &p, double &w0, double &w1, double &w2) const
		{
			const Eigen::Vector3d a = m_V.row(m_F(f, 0));
			const Eigen::Vector3d e0 = Eigen::Vector3d(m_V.row(m_F(f, 1))) - a;
			const Eigen::Vector3d e1 = Eigen::Vector3d(m_V.row(m_F(f, 2))) - a;
			const Eigen::Vector3d e2 = p - a;
			const double d00 = e0.dot(e0);
			const double d01 = e0.dot(e1);
			const double d11 = e1.dot(e1);
			const double d20 = e2.dot(e0);
			const double d21 = e2.dot(e1);
			const double denom = d00 * d11 - d01 * d01;
			if (denom <= 0.0)
			{
				// degenerate triangle
				w0 = w1 = w2 = 1.0 / 3.0;
				return;
			}
			w1 = (d11 * d20 - d01 * d21) / denom;
			w2 = (d00 * d21 - d01 * d20) / denom;
			w0 = 1.0 - w1 - w2;
		}

		/** Interpolate the per-vertex mean curvature on face f using the weights (w0, w1, w2). */
		double meanCurvature(const int f, const double w0, const double w1, const double w2) const
		{
			return w0 * m_H[m_F(f, 0)] + w1 * m_H[m_F(f, 1)] + w2 * m_H[m_F(f, 2)];
		}

		/** Interpolate the per-vertex Gaussian curvature on face f using the weights (w0, w1, w2). */
		double gaussianCurvature(const int f, const double w0, const double w1, const double w2) const
		{
			return w0 * m_K[m_F(f, 0)] + w1 * m_K[m_F(f, 1)] + w2 * m_K[m_F(f, 2)];
		}

	private:
		/** Signed volume of the mesh (positive for closed meshes with outward oriented faces). */
		double signedVolume() const
		{
			double vol = 0.0;
			for (int f = 0; f < m_F.rows(); ++f)
			{
				const Eigen::Vector3d a = m_V.row(m_F(f, 0));
				const Eigen::Vector3d b = m_V.row(m_F(f, 1));
				const Eigen::Vector3d c = m_V.row(m_F(f, 2));
				vol += a.dot(b.cross(c));
			}
			return vol / 6.0;
		}

		Eigen::MatrixXd m_V;
		Eigen::MatrixXi m_F;
		Eigen::VectorXd m_H;
		Eigen::VectorXd m_K;
		igl::AABB<Eigen::MatrixXd, 3> m_tree;
		igl::FastWindingNumberBVH m_fwnBVH;

		void computeVertexCurvatures()
		{
			// mean curvature normal: -M^-1 L V = (k1 + k2) n
			Eigen::SparseMatrix<double> L, M, Minv;
			igl::cotmatrix(m_V, m_F, L);
			igl::massmatrix(m_V, m_F, igl::MASSMATRIX_TYPE_VORONOI, M);
			igl::invert_diag(M, Minv);
			L = -Minv * L;
			const Eigen::MatrixXd HN = L * m_V;

			Eigen::MatrixXi BF;
			igl::boundary_facets(m_F, BF);

			// outward vertex normals (flip them if a closed mesh is oriented inwards)
			Eigen::MatrixXd N;
			igl::per_vertex_normals(m_V, m_F, N);
			if ((BF.rows() == 0) && (signedVolume() < 0.0))
				N = -N;

			// signed mean curvature, negative for convex surfaces
			m_H = -(HN.cwiseProduct(N)).rowwise().sum();

			// Gaussian curvature: angle defect / Voronoi area
			const Eigen::VectorXd A = M.diagonal();
			const double eps = 1e-16;
			std::vector<char> isBoundary((size_t)m_V.rows(), 0);
			for (int r = 0; r < BF.rows(); ++r)
				for (int c = 0; c < BF.cols(); ++c)
					isBoundary[(size_t)BF(r, c)] = 1;

			auto cornerAngle = [&](const int vi, const int vj, const int vk) -> double
			{
				const Eigen::Vector3d pi = m_V.row(vi);
				const Eigen::Vector3d a = Eigen::Vector3d(m_V.row(vj)) - pi;
				const Eigen::Vector3d b = Eigen::Vector3d(m_V.row(vk)) - pi;
				const double na = a.norm();
				const double nb = b.norm();
				if (na < eps || nb < eps)
					return 0.0;
				const double c = std::max(-1.0, std::min(1.0, a.dot(b) / (na * nb)));
				return std::acos(c);
			};

			Eigen::VectorXd angleSum = Eigen::VectorXd::Zero(m_V.rows());
			for (int f = 0; f < m_F.rows(); ++f)
			{
				const int i0 = m_F(f, 0), i1 = m_F(f, 1), i2 = m_F(f, 2);
				angleSum[i0] += cornerAngle(i0, i1, i2);
				angleSum[i1] += cornerAngle(i1, i2, i0);
				angleSum[i2] += cornerAngle(i2, i0, i1);
			}

			m_K.resize(m_V.rows());
			for (int i = 0; i < m_V.rows(); ++i)
			{
				const double target = isBoundary[(size_t)i] ? M_PI : (2.0 * M_PI);
				m_K[i] = (target - angleSum[i]) / std::max(A[i], eps);
			}

			// remove invalid values (e.g. at degenerate triangles)
			for (int i = 0; i < m_V.rows(); ++i)
			{
				if (!std::isfinite(m_H[i]))
					m_H[i] = 0.0;
				if (!std::isfinite(m_K[i]))
					m_K[i] = 0.0;
			}
		}
	};
}

#endif
