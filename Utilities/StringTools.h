#ifndef __StringTools_h__
#define __StringTools_h__

#include "SPlisHSPlasH/Common.h"
#include <vector>

namespace SPH
{
	/** \brief Tools to handle std::string objects 
	*/
	class StringTools
	{
	public:
		static void tokenize(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters = " ");
	};
}

#endif
