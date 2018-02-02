#ifndef __PartioReaderWriter_h__
#define __PartioReaderWriter_h__

#include "SPlisHSPlasH/Common.h"
#include <vector>

namespace Utilities
{
	/** \brief Class for reading and writing partio files.
	*/
	class PartioReaderWriter
	{
	public:
		static bool readParticles(const std::string &fileName, const Vector3r &translation, const Matrix3r &rotation, const Real scale,
			std::vector<Vector3r> &pos, std::vector<Vector3r> &vel);

		static bool readParticles(const std::string &fileName, const Vector3r &translation, const Matrix3r &rotation, const Real scale,
			std::vector<Vector3r> &positions, std::vector<Vector3r> &velocities, Real &particleRadius);

		static bool readParticles(const std::string &fileName, const Vector3r &translation, const Matrix3r &rotation, const Real scale,
			std::vector<Vector3r> &pos);

		static void writeParticles(const std::string &fileName, const unsigned int numParticles, const Vector3r *particlePositions,
			const Vector3r *particleVelocities, const Real particleRadius);
	};

}

#endif
