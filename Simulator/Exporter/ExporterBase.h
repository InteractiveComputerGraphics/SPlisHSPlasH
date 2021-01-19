#ifndef __ExporterBaseBase_h__
#define __ExporterBaseBase_h__

#include "SPlisHSPlasH/Common.h"
#include "Simulator/SimulatorBase.h"

namespace SPH
{
	/** \brief Base class for data exporters.
	*/
	class ExporterBase 
	{
	protected:
		SimulatorBase* m_base;
		bool m_active;

	public:
		ExporterBase(SimulatorBase* base) : m_active(false) { m_base = base; };
		ExporterBase(const ExporterBase&) = delete;
        ExporterBase& operator=(const ExporterBase&) = delete;
		virtual ~ExporterBase(void) {};

		virtual void step(const unsigned int frame) = 0;
		virtual void reset() {};

		virtual void init(const std::string &outputPath) {};

		virtual void setActive(const bool active) { m_active = active; }
		virtual bool getActive() const { return m_active; }
	};
}

#endif
