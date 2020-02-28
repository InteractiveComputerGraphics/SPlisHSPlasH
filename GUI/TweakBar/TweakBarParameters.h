#ifndef __TweakBarParameters_h__
#define __TweakBarParameters_h__

#include "SPlisHSPlasH/Common.h"
#include "extern/AntTweakBar/include/AntTweakBar.h"
#include <vector>
#include "ParameterObject.h"

#ifdef USE_DOUBLE
#define TW_TYPE_REAL TW_TYPE_DOUBLE
#define TW_TYPE_DIR3R TW_TYPE_DIR3D
#else
#define TW_TYPE_REAL TW_TYPE_FLOAT
#define TW_TYPE_DIR3R TW_TYPE_DIR3F
#endif

namespace SPH
{
	class TweakBarParameters
	{
	public: 
		typedef std::pair<GenParam::ParameterObject*, unsigned int> ParameterIndex;

		static void createParameterGUI(TwBar *tweakBar);
		static void createParameterObjectGUI(TwBar *tweakBar, GenParam::ParameterObject *paramObj);

		static void TW_CALL setParameterValue(const void *value, void *clientData);
		static void TW_CALL getParameterValue(void *value, void *clientData);

		static void TW_CALL setTimeStepSizeCB(const void *value, void *clientData);
		static void TW_CALL getTimeStepSizeCB(void *value, void *clientData);

		static void cleanup();

	protected:
		static std::vector<std::unique_ptr<ParameterIndex>> m_params;
		static std::vector<std::string> m_objectNames;
	};
}
 
#endif