#ifndef __Simulator_GUI_Base_h__
#define __Simulator_GUI_Base_h__

#include "SPlisHSPlasH/Common.h"
#include "ParameterObject.h"
#include "Simulator/SimulatorBase.h"

namespace SPH 
{	
	class Simulator_GUI_Base
	{
		protected:
			SimulatorBase *m_simulatorBase;

		public:
			Simulator_GUI_Base(SimulatorBase *simulatorBase) : m_simulatorBase(simulatorBase) {};
			virtual ~Simulator_GUI_Base() {};

		public:
			virtual void init(int argc, char **argv, const char *name) {}
			virtual void initSimulationParameterGUI() {}
			virtual void initParameterGUI() {}
			virtual void render() {}
			virtual void reset() {}
			virtual void update() {}
			virtual void cleanup() {}
			virtual void run() {}
			virtual void stop() {}
			virtual void addKeyFunc(char k, std::function<void()> const& func) {}

			SPH::SimulatorBase * getSimulatorBase() const { return m_simulatorBase; }
	};
}

#endif