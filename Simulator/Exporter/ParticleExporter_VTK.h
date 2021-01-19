#ifndef __ParticleExporter_VTK_h__
#define __ParticleExporter_VTK_h__

#include "ExporterBase.h"
#include "SPlisHSPlasH/FluidModel.h"

namespace SPH
{
	/** \brief Particle exporter for the VTK format.
	*/
	class ParticleExporter_VTK : public ExporterBase
	{
	protected: 
		std::string m_exportPath;

		void writeParticles(const std::string& fileName, FluidModel* model);

		// VTK expects big endian
		template<typename T>
		inline void swapByteOrder(T* v)
		{
			constexpr size_t n = sizeof(T);
			uint8_t* bytes = reinterpret_cast<uint8_t*>(v);
			for (unsigned int c = 0u; c < n / 2; c++)
				std::swap(bytes[c], bytes[n - c - 1]);
		}

	public:
		ParticleExporter_VTK(SimulatorBase *base);
		ParticleExporter_VTK(const ParticleExporter_VTK&) = delete;
		ParticleExporter_VTK& operator=(const ParticleExporter_VTK&) = delete;
		virtual ~ParticleExporter_VTK(void);

		virtual void init(const std::string& outputPath);
		virtual void step(const unsigned int frame);
		virtual void reset();
		virtual void setActive(const bool active); 
	};
}

#endif
