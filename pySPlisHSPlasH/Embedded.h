#ifndef __Embedded_h__
#define __Embedded_h__

#include <string>
#include <functional>
#include <iostream>
#include <pybind11/pybind11.h>

namespace SPH
{
	class SimulatorBase;

	class Embedded
	{
	private:
		static Embedded *current;

	protected:
		pybind11::object m_mainScope;
		std::vector<std::string> m_functions;

		void findFunctions(const std::string &moduleName);

	public:
		Embedded ();
		~Embedded ();

		// Singleton
		static Embedded* getCurrent ();
		static void setCurrent (Embedded* tm);
		static bool hasCurrent();


		std::string import_script(const std::string& fileName);
		void exec_script_str(const std::string& str);
		void exec_script_file(const std::string& fileName);
		void exec_fct(const std::string &moduleName, const std::string& fctName);
		void exec_init(const std::string& moduleName, SimulatorBase *base);

		const std::vector<std::string>& getFunctions() const { return m_functions; }

		void reloadModule(const std::string &moduleName);
		void releaseModule(const std::string& moduleName);
	};
}

#endif
