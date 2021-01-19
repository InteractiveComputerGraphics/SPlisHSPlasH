#ifndef __ParticleExporter_Partio_h__
#define __ParticleExporter_Partio_h__

#include "ExporterBase.h"
#include "SPlisHSPlasH/FluidModel.h"

namespace SPH
{
	/** \brief Particle exporter for the partio format.
	*/
	class ParticleExporter_Partio : public ExporterBase
	{
	protected: 
		std::string m_exportPath;

		void writeParticlesPartio(const std::string& fileName, FluidModel* model);

	public:
		ParticleExporter_Partio(SimulatorBase *base);
		ParticleExporter_Partio(const ParticleExporter_Partio&) = delete;
        ParticleExporter_Partio& operator=(const ParticleExporter_Partio&) = delete;
		virtual ~ParticleExporter_Partio(void);

		virtual void init(const std::string& outputPath);
		virtual void step(const unsigned int frame);
		virtual void reset();
		virtual void setActive(const bool active); 
	};
}

#endif
