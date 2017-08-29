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
		
		/** Compute the MD5 hash of a file.
		*/
		static std::string getFileMD5(const std::string &filename);

		/** Write the MD5 hash of a file to the md5File.  
		 */
		static bool writeMD5File(const std::string& fileName, const std::string& md5File);

		/** Compare an MD5 hash with the hash stored in an MD5 file.
		*/
		static bool checkMD5(const std::string& md5Hash, const std::string& md5File);
	};
}

#endif
