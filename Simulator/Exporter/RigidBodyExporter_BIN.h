#ifndef __RigidBodyExporter_BIN_h__
#define __RigidBodyExporter_BIN_h__

#include "ExporterBase.h"
#include "SPlisHSPlasH/FluidModel.h"

namespace SPH
{
	/** \brief Rigid body exporter for the bin format (own file format).
	*/
	class RigidBodyExporter_BIN : public ExporterBase
	{
	protected: 
		bool m_isFirstFrame;
		std::string m_exportPath;

		void writeRigidBodies(const std::string& fileName);

	public:
		RigidBodyExporter_BIN(SimulatorBase* base);
		RigidBodyExporter_BIN(const RigidBodyExporter_BIN&) = delete;
        RigidBodyExporter_BIN& operator=(const RigidBodyExporter_BIN&) = delete;
		virtual ~RigidBodyExporter_BIN(void);

		virtual void init(const std::string& outputPath);
		virtual void step(const unsigned int frame);
		virtual void reset();
		virtual void setActive(const bool active); 
	};
}

#endif
