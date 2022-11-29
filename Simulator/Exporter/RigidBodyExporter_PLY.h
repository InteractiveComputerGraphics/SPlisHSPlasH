#ifndef __RigidBodyExporter_PLY_h__
#define __RigidBodyExporter_PLY_h__

#include "ExporterBase.h"

namespace SPH
{
	/** \brief Rigid body exporter for the OBJ format.
	*/
	class RigidBodyExporter_PLY : public ExporterBase
	{
	protected: 
		bool m_isFirstFrame;
		std::string m_exportPath;

		void writeRigidBodies(const unsigned int frame);

	public:
		RigidBodyExporter_PLY(SimulatorBase* base);
		RigidBodyExporter_PLY(const RigidBodyExporter_PLY&) = delete;
        RigidBodyExporter_PLY& operator=(const RigidBodyExporter_PLY&) = delete;
		virtual ~RigidBodyExporter_PLY(void);

		virtual void init(const std::string& outputPath);
		virtual void step(const unsigned int frame);
		virtual void reset();
		virtual void setActive(const bool active); 
	};
}

#endif
