#ifndef __FileSystem_h__
#define __FileSystem_h__

#include "SPlisHSPlasH/Common.h"

namespace SPH
{
	/** \brief This class implements different file system functions.
	*/
	class FileSystem
	{
	public:
		static std::string getFilePath(const std::string &path);
		static std::string getFileName(const std::string &path);
		static bool isRelativePath(const std::string &path);
		static int makeDir(const std::string &path);

		/** Make all subdirectories.
		*/
		static int makeDirs(const std::string &path);
		static std::string normalizePath(const std::string &path);
		static bool fileExists(const std::string& fileName);
		static std::string getProgramPath();
	};
}

#endif
