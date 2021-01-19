#ifndef __RigidBodyExporter_OBJ_h__
#define __RigidBodyExporter_OBJ_h__

#include "ExporterBase.h"

namespace SPH
{
	/** \brief Rigid body exporter for the OBJ format.
	*/
	class RigidBodyExporter_OBJ : public ExporterBase
	{
	protected: 
		bool m_isFirstFrame;
		std::string m_exportPath;

		void writeRigidBodies(const unsigned int frame);

	public:
		RigidBodyExporter_OBJ(SimulatorBase* base);
		RigidBodyExporter_OBJ(const RigidBodyExporter_OBJ&) = delete;
        RigidBodyExporter_OBJ& operator=(const RigidBodyExporter_OBJ&) = delete;
		virtual ~RigidBodyExporter_OBJ(void);

		virtual void init(const std::string& outputPath);
		virtual void step(const unsigned int frame);
		virtual void reset();
		virtual void setActive(const bool active); 
	};
}

#endif
