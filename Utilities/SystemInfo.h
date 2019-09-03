#ifndef __SystemInfo_h__
#define __SystemInfo_h__

#include "SPlisHSPlasH/Common.h"
#if WIN32
#define NOMINMAX
#include "windows.h"
#else
#include <unistd.h>
#include <limits.h>
#endif

namespace Utilities
{
	class SystemInfo
	{
	public:
		static std::string getHostName()
		{
#ifdef WIN32
			const unsigned int bufferSize = 32767;
			TCHAR  *infoBuf = new TCHAR[bufferSize];
			DWORD  bufCharCount = bufferSize;
			if (!GetComputerName(infoBuf, &bufCharCount))
				return "";
			std::string name = infoBuf;
			delete[] infoBuf;
			return name;
#else
			const unsigned int bufferSize = 32767;
			char hostname[bufferSize];
			gethostname(hostname, bufferSize);
			return hostname;
#endif
		}
	};
}

#endif
