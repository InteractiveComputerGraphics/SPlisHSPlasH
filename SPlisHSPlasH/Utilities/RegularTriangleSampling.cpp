#include "RegularTriangleSampling.h"

#include <set>
#include <vector>

using namespace SPH;

RegularTriangleSampling::RegularTriangleSampling()
{
}

void RegularTriangleSampling::sampleMesh(const unsigned int numVertices, const Vector3r * vertices, const unsigned int numFaces, const unsigned int * faces, const Real maxDistance, std::vector<Vector3r>& samples)
{
	const std::vector<Vector2ui> edges = uniqueEdges(numFaces, faces);
	
	appendVertexSamples(numVertices, vertices, samples);
	appendEdgeSamples(maxDistance, vertices, edges, samples);
	appendFaceSamples(maxDistance, vertices, numFaces, faces, samples);
}

void RegularTriangleSampling::appendVertexSamples(const unsigned int  numVertices, const Vector3r* vertices,
	std::vector<Vector3r>& samples)
{
	samples.insert(samples.end(), vertices, vertices + numVertices);
}

void RegularTriangleSampling::appendEdgeSamples(const Real d, const Vector3r* vertices, const std::vector<Vector2ui> & edges,
	std::vector<Vector3r>& samples, bool skipVertices)
{
	const auto iSkip = static_cast<unsigned int>(skipVertices);
	for (const Vector2ui & edge : edges)
	{
		const Vector3r & v0 = vertices[edge[0]];
		const Vector3r dir = vertices[edge[1]] - v0;
		const Real l = dir.norm();
		// ceil for oversampling, +1 to get from end to end
		const auto numSamples = static_cast<unsigned int>(std::ceil(l / d)) + 1u;
		// check for special case: zero-length edge
		if (numSamples == 1)
		{
			samples.emplace_back(v0);
			continue;
		}
		for (unsigned i = iSkip; i < numSamples - iSkip; i++)
			samples.emplace_back(v0 + static_cast<Real>(i) / (numSamples - 1) * dir);
	}
}

void RegularTriangleSampling::appendFaceSamples(const Real d, const Vector3r* vertices, const unsigned numFaces,
	const unsigned int * faces, std::vector<Vector3r>& samples, bool skipEdges)
{
	const auto iSkip = static_cast<unsigned int>(skipEdges);
	for (unsigned int iFace = 0; iFace < numFaces; iFace++)
	{
		using Matrix3ui = Eigen::Matrix<unsigned int, 3, 1, Eigen::DontAlign>;
		const auto face = Eigen::Map<const Matrix3ui>(faces + 3 * iFace);
		// chose longest edge as base
		Vector3r v0, base, toTop;
		Real l2 = 0;
		for (unsigned int c = 0; c < 3; c++)
		{
			const Vector3r v = vertices[face[(c + 1) % 3]] - vertices[face[c]];
			const Real le2 = v.squaredNorm();
			if (le2 > l2)
			{
				l2 = le2;
				v0 = vertices[face[c]];
				base = v;
				// edge to top of the triangle
				toTop = vertices[face[(c + 2) % 3]] - vertices[face[c]];
			}
		}
		const Real l = std::sqrt(l2);
		const Vector3r bn = base / l;
		// height of the triangle
		const Real h = (toTop - bn.dot(toTop) * bn).norm();
		// ceil for oversampling, +1 to get from end to end
		const auto numLines = static_cast<unsigned int>(std::ceil(h / d)) + 1u;
		for (unsigned int line = iSkip; line < numLines - iSkip; line++)
		{
			Vector3r lineStart = v0;
			Real lineLength = 0;
			if (numLines > 0)
			{
				const Real alpha = static_cast<Real>(line) / (numLines - 1);
				lineStart += alpha * toTop;
				lineLength = (1 - alpha) * l;
			}
			// ceil for oversampling, +1 to get from end to end
			const unsigned int numSamples = static_cast<unsigned int>(std::ceil(lineLength / d)) + 1u;
			for (unsigned i = iSkip; i < numSamples - iSkip; i++)
				samples.emplace_back(lineStart + (i * lineLength) / (numSamples - 1) * bn);
		}
	}
}

std::vector<RegularTriangleSampling::Vector2ui> RegularTriangleSampling::uniqueEdges(unsigned int numFaces, const unsigned int * faces)
{
	std::vector<Vector2ui> edges;
	edges.reserve(3 * numFaces);
	for (unsigned int i = 0; i < numFaces; ++i)
	{
		for (unsigned int c = 0; c < 3; c++)
		{
			unsigned int ia = faces[3 * i + c];
			unsigned int ib = faces[3 * i + (c + 1) % 3];
			// smaller index first
			if (ia > ib)
				std::swap(ia, ib);
			edges.emplace_back(ia, ib);
		}
	}
	const auto less = [](const Vector2ui & a, const Vector2ui & b) -> bool { return (a[0] < b[0]) || (a[0] == b[0] && a[1] < b[1]); };
	std::sort(edges.begin(), edges.end(), less);
	edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
	return edges;
}
