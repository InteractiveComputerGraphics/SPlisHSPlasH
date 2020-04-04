#ifndef __Exporter_h__
#define __Exporter_h__

#include <vector>
#include <array>

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/Utilities/SceneLoader.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "SPlisHSPlasH/TimeManager.h"

// This class handles exporting to external formats.
class Exporter {

public:

	// Flags
	static int PARTIO_EXPORT; static bool enablePartioExport;
	static int VTK_EXPORT; static bool enableVTKExport;
	static int RB_EXPORT; static bool enableRigidBodyExport;
	static int RB_VTK_EXPORT; static bool enableRigidBodyVTKExport;

	// Configs
	static int DATA_EXPORT_FPS; static Real framesPerSecond;
	static int PARTICLE_EXPORT_ATTRIBUTES; static std::string particleAttributes;

private:

	bool m_isFirstFrameVTK;

public:

	std::string m_outputPath;

	Exporter();
	~Exporter();

	void particleExport(std::string exportName = "", std::string temporalIdentifier = "", std::string folder = "", bool partio = enablePartioExport, bool vtk = enableVTKExport);
	void writeParticlesPartio(const std::string& fileName, SPH::FluidModel* model);
	void writeParticlesVTK(const std::string& fileName, SPH::FluidModel* model);

	void rigidBodyExport(std::string scenePath, std::string temporalIdentifier, bool isFirstFrame, Utilities::SceneLoader::Scene);
	void writeRigidBodiesBIN(const std::string& exportPath, std::string scene_path, std::string temporalIdentifier, bool isFirstFrame, Utilities::SceneLoader::Scene);
	void writeRigidBodiesVTK(const std::string& exportPath, std::string temporalIdentifier);

	bool isTrackingParticles() {
		return enablePartioExport || enableVTKExport;
	}

	bool isTrackingRigidBodies() {
		return enableRigidBodyExport || enableRigidBodyVTKExport;
	}

	void saveParticleSnapshot() {
		particleExport("", std::to_string(SPH::TimeManager::getCurrent()->getTime()), "snapshots", true, true);
	}

	void reset() {
		m_isFirstFrameVTK = true;
	}

	// VTK expects big endian
	template<typename T>
	inline void swapByteOrder(T* v)
	{
		constexpr size_t n = sizeof(T);
		uint8_t* bytes = reinterpret_cast<uint8_t*>(v);
		for (unsigned int c = 0u; c < n / 2; c++)
			std::swap(bytes[c], bytes[n - c - 1]);
	}
};

#endif