#include "RegularSampling2D.h"

#include <vector>

using namespace SPH;

RegularSampling2D::RegularSampling2D()
{
}

void RegularSampling2D::sampleMesh(const Matrix3r& rotation, const Vector3r & translation, const unsigned numVertices, const Vector3r * vertices, const unsigned int numFaces, const unsigned int * faces, const Real maxDistance, std::vector<Vector3r>& samples)
{
	using Vector3ui = Eigen::Matrix<unsigned int, 3, 1>;
	using Vector3b = Eigen::Matrix<bool, 3, 1>;

	// transform vertices
	std::vector<Vector3r> x(numVertices);
	#pragma omp parallel
	{
		#pragma omp for schedule(static)
		for (int i = 0; i < static_cast<int>(numVertices); i++)
		{
			x[i] = rotation * vertices[i] + translation;
		}
		#pragma omp barrier
		#pragma omp for schedule(static)
		for (int i = 0; i < (int) numFaces; i++)
		{
			// get face indices
			const Vector3ui & face = Eigen::Map<const Vector3ui>(faces + 3 * i);
			// find edges that cut the z=0 plane
			Vector3b cutsZ;
			for (unsigned int c = 0; c < 3; c++)
			{
				const Real z1 = x[face[c]].z();
				const Real z2 = x[face[(c + 1) % 3]].z();
				const bool bothGreaterOrEqual = z1 >= 0 && z2 >= 0;
				const bool bothLess = z1 < 0 && z2 < 0;
				cutsZ[c] = !(bothGreaterOrEqual || bothLess);
			}
			// ignore faces that do not cut the z=0 plane
			if (!cutsZ.any())
				continue;

			// compute cuts of edges with the plane
			std::vector<Vector3r> cuts;
			cuts.reserve(2);
			for (unsigned int c = 0; c < 3; c++)
			{
				if (!cutsZ[c])
					continue;
				const Vector3r & a = x[face[c]];
				const Vector3r & b = x[face[(c + 1) % 3]];
				const Vector3r d = b - a;
				const Real alpha = a.z() / d.z();
				cuts.emplace_back(a.x() - alpha * d.x(), a.y() - alpha * d.y(), 0);
			}

			// sample line between cuts
			const Vector3r & v0 = cuts[0];
			const Vector3r dir = cuts[1] - v0;
			const Real l = dir.norm();
			// ceil for oversampling, +1 to get from end to end
			const auto numSamples = static_cast<unsigned int>(std::ceil(l / maxDistance)) + 1u;
			#pragma omp critical
			{
				for (unsigned sample = 0; sample < numSamples; sample++)
					samples.emplace_back(v0 + static_cast<Real>(sample) / (numSamples - 1) * dir);
			}
		}
	}
	
	// remove duplicate points
	const auto less = [](const Vector3r & a, const Vector3r & b) -> bool { return (a[0] < b[0]) || (a[0] == b[0] && a[1] < b[1] || (a[0] == b[0] && a[1] == b[1] && a[2] < b[2])); };
	std::sort(samples.begin(), samples.end(), less);
	
	const Real squaredMinDist = 1e-12f;
	const auto equals = [squaredMinDist](const Vector3r & a, const Vector3r & b) -> bool { return (a - b).squaredNorm() < squaredMinDist; };
	samples.erase(std::unique(samples.begin(), samples.end(), equals), samples.end());
}
