#ifndef __PartioViewer_GUI_Base_h__
#define __PartioViewer_GUI_Base_h__

#include "SPlisHSPlasH/Common.h"
#include "ParameterObject.h"
#include "Simulator/SimulatorBase.h"
#include "extern/cxxopts/cxxopts.hpp"

namespace SPH 
{	
	class PartioViewer_GUI_Base
	{
		protected:
			SimulatorBase *m_simulatorBase;

		public:
			PartioViewer_GUI_Base() {};
			virtual ~PartioViewer_GUI_Base() {};

		public:
			virtual void init() {}
			virtual void initSimulationParameterGUI() {}
			virtual void initParameterGUI() {}
			virtual void render() {}
			virtual void run() {}
			virtual void stop() {}
			virtual void update() {}
			virtual void cleanup() {}
			virtual void addOptions(cxxopts::Options &options) {}
			virtual void parseOptions(cxxopts::ParseResult &result) {}
			virtual unsigned int getWidth() const { return 0; }
			virtual unsigned int getHeight() const { return 0; }

			SPH::SimulatorBase * getSimulatorBase() const { return m_simulatorBase; }
			void setSimulatorBase(SPH::SimulatorBase * val) { m_simulatorBase = val; }
	};
}

#endif