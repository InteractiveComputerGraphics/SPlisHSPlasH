#ifndef __RigidBodyParticleExporter_VTK_h__
#define __RigidBodyParticleExporter_VTK_h__

#include "ExporterBase.h"
#include "SPlisHSPlasH/FluidModel.h"

namespace SPH {
	/** \brief Rigid body particle exporter for the VTK format.
	*/
	class RigidBodyParticleExporter_VTK : public ExporterBase {
	protected:
		bool m_isFirstFrame;
		std::vector<std::string> m_attributes;
		std::string m_exportPath;

		void writeRigidBodies(const unsigned int frame);

		// VTK expects big endian
		template<typename T>
		inline void swapByteOrder(T* v) {
			constexpr size_t n = sizeof(T);
			uint8_t* bytes = reinterpret_cast<uint8_t*>(v);
			for (unsigned int c = 0u; c < n / 2; c++)
				std::swap(bytes[c], bytes[n - c - 1]);
		}

	public:
		RigidBodyParticleExporter_VTK(SimulatorBase* base);
		RigidBodyParticleExporter_VTK(const RigidBodyParticleExporter_VTK&) = delete;
		RigidBodyParticleExporter_VTK& operator=(const RigidBodyParticleExporter_VTK&) = delete;
		virtual ~RigidBodyParticleExporter_VTK(void);

		virtual void init(const std::string& outputPath);
		virtual void step(const unsigned int frame);
		virtual void reset();
		virtual void setActive(const bool active);
	};
}

#endif
