#ifndef __ScriptObject_h__
#define __ScriptObject_h__

#include <string>
#include <ParameterObject.h>

namespace SPH
{
	class SimulatorBase;

	class ScriptObject : public GenParam::ParameterObject
	{
	protected: 
		SimulatorBase* m_base;
		bool m_scriptLoaded;
		std::string m_scriptModule;
		std::string m_scriptFile;

		virtual void initParameters();
		void addFunctionParameters();
		void removeFunctionParameters();

	public:
		static int SCRIPT_FILE;

		ScriptObject(SimulatorBase *base);
		virtual ~ScriptObject(void);

		void init();
		std::string loadScriptFile(const std::string& fileName);

		void execResetFct();
		void execStepFct();

		void updateFunctionParameters();

	};
}

#endif
