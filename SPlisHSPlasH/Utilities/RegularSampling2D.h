#ifndef __RegularSampling2D_H__
#define __RegularSampling2D_H__

#include "../Common.h"

#include <vector>
#include <Eigen/src/Core/util/ForwardDeclarations.h>

namespace SPH
{
	/** \brief This class implements a per-triangle regular sampling for the surface
	* of 3D models. 
	*/
	class RegularSampling2D
	{
	public:
		RegularSampling2D();

		/** Performs the poisson sampling with the
		* respective parameters. Compare
		* http://graphics.cs.umass.edu/pubs/sa_2010.pdf
		*
		* @param rotation rotation of the mesh
		* @param translation translation of the mesh
		* @param numVertices number of mesh vertices
		* @param vertices vertex data of sampled data
		* @param numFaces number of faces in the mesh
		* @param faces face data of sampled mesh
		* @param maxDistance maximal distance of sampled vertices
		* @param samples vector to store the samples
		*/
		static void sampleMesh(const Matrix3r& rotation, const Vector3r & translation, const unsigned numVertices, const Vector3r *vertices,
		                       const unsigned int numFaces, const unsigned int *faces, const Real maxDistance, std::vector<Vector3r> &samples);
		

	};
}

#endif // __RegularSampling2D_H__
