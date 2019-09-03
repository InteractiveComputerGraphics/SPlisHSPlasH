#ifndef VolumeSampling_H
#define VolumeSampling_H

#include "../Common.h"
#include <vector>


namespace Utilities
{
	/** \brief This class implements a volume sampling of 3D models. 
	*/
	class VolumeSampling
	{
	public:
		/** Performs the volume sampling with the
		* respective parameters. 
		*
		* @param numVertices number of vertices
		* @param vertices vertex data 
		* @param numFaces number of faces
		* @param faces index list of faces
		* @param radius radius of sampled particles
		* @param region defines a subregion of the mesh to be sampled (nullptr if not used)
		* @param resolution resolution of the used SDF
		* @param invert defines if the mesh should be inverted and the outside is sampled
		* @param sampleMode 0=regular, 1=almost dense, 2=dense
		* @param samples sampled vertices that will be returned
		*/
		static void sampleMesh(const unsigned int numVertices, const Vector3r *vertices, 
			const unsigned int numFaces, const unsigned int *faces,
			const Real radius, const AlignedBox3r *region,
			const std::array<unsigned int, 3> &resolution, const bool invert,
			const unsigned int sampleMode,
			std::vector<Vector3r> &samples);		
	};
}

#endif // VolumeSampling_H
