#ifndef __RigidBodyExporter_VTK_h__
#define __RigidBodyExporter_VTK_h__

#include "ExporterBase.h"
#include "SPlisHSPlasH/FluidModel.h"

namespace SPH
{
	/** \brief Rigid body exporter for the VTK format.
	*/
	class RigidBodyExporter_VTK : public ExporterBase
	{
	protected: 
		bool m_isFirstFrame;
		std::string m_exportPath;

		void writeRigidBodies(const unsigned int frame);

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
		RigidBodyExporter_VTK(SimulatorBase* base);
		RigidBodyExporter_VTK(const RigidBodyExporter_VTK&) = delete;
        RigidBodyExporter_VTK& operator=(const RigidBodyExporter_VTK&) = delete;
		virtual ~RigidBodyExporter_VTK(void);

		virtual void init(const std::string& outputPath);
		virtual void step(const unsigned int frame);
		virtual void reset();
		virtual void setActive(const bool active); 
	};
}

#endif
