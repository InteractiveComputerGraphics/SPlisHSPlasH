#ifndef __RegularTriangleSampling_H__
#define __RegularTriangleSampling_H__

#include "../Common.h"

#include <vector>

namespace SPH
{
	/** \brief This class implements a per-triangle regular sampling for the surface
	* of 3D models. 
	*/
	class RegularTriangleSampling
	{
	public:
		RegularTriangleSampling();

		/** Performs the poisson sampling with the
		* respective parameters. Compare
		* http://graphics.cs.umass.edu/pubs/sa_2010.pdf
		*
		* @param numVertices number of mesh vertices
		* @param vertices vertex data of sampled data
		* @param numFaces number of faces in the mesh
		* @param faces face data of sampled mesh
		* @param maxDistance maximal distance of sampled vertices
		* @param samples vector to store the samples
		*/
		static void sampleMesh(const unsigned int numVertices, const Vector3r *vertices, const unsigned int numFaces, const unsigned int *faces,
			const Real maxDistance, std::vector<Vector3r> &samples);
		
	private:
		using Vector2ui = Eigen::Matrix<unsigned int, 2, 1, Eigen::DontAlign>;

		static void appendVertexSamples(const unsigned int numVertices, const Vector3r * vertices, std::vector<Vector3r> &samples);
		static void appendEdgeSamples(const Real d, const Vector3r * vertices, const std::vector<Vector2ui> & edges, std::vector<Vector3r> &samples, bool skipVertices = true);
		static void appendFaceSamples(const Real d, const Vector3r * vertices, const unsigned int numFaces, const unsigned int * faces, std::vector<Vector3r> &samples, bool skipEdges = true);
		
		static std::vector<Vector2ui> uniqueEdges(unsigned int numFaces, const unsigned int *faces);
		

	};
}

#endif // __RegularTriangleSampling_H__
