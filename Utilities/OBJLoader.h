#ifndef __OBJLoader_h__
#define __OBJLoader_h__

#include "SPlisHSPlasH/Common.h"
#include <string>
#include "SPlisHSPlasH/TriangleMesh.h"

namespace SPH
{
	/** \brief Struct to store the position and normal indices
	*/
	struct MeshFaceIndices
	{
		int posIndices[3];
		int normalIndices[3];
	};

	/** \brief Read for OBJ files. 
	*/
	class OBJLoader
	{
	public:
		/** This function loads an OBJ file. 
		  * Only triangulated meshes are supported. 
		  */
		static void loadObj(const std::string &filename, TriangleMesh &mesh, const Vector3r scale = Vector3r::Ones());		
	};
}
 
#endif