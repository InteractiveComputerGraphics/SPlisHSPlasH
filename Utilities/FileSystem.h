#ifndef __FileSystem_h__
#define __FileSystem_h__

#include "StringTools.h"
#include "extern/md5/md5.h"
#include <sys/stat.h>
#ifdef WIN32
#include <direct.h>
#define NOMINMAX
#include "windows.h"
#include <commdlg.h>
#else
#include <unistd.h>
#include <dirent.h>
#include <dlfcn.h>
#endif
#ifdef __APPLE__
#include <mach-o/dyld.h>
#endif

#ifndef S_ISDIR
#define S_ISDIR(mode)  (((mode) & S_IFMT) == S_IFDIR)
#endif

#ifndef S_ISREG
#define S_ISREG(mode)  (((mode) & S_IFMT) == S_IFREG)
#endif

namespace Utilities
{
	/** \brief This class implements different file system functions.
	*/
	class FileSystem
	{
	public:

		static std::string getFilePath(const std::string &path)
		{
			std::string npath = normalizePath(path);

			std::string result = npath;
			size_t i = result.rfind('.', result.length());
			if (i != std::string::npos)
			{
				result = result.substr(0, i);
			}
			size_t p1 = result.rfind('\\', result.length());
			size_t p2 = result.rfind('/', result.length());
			if ((p1 != std::string::npos) && (p2 != std::string::npos))
				result = result.substr(0, std::max(p1, p2));
			else if (p1 != std::string::npos)
				result = result.substr(0, p1);
			else if (p2 != std::string::npos)
				result = result.substr(0, p2);
			return result;
		}

		static std::string getFileName(const std::string &path)
		{
			std::string npath = normalizePath(path);

			std::string result = npath;
			size_t i = result.rfind('.', result.length());
			if (i != std::string::npos)
			{
				result = result.substr(0, i);
			}
			size_t p1 = result.rfind('\\', result.length());
			size_t p2 = result.rfind('/', result.length());
			if ((p1 != std::string::npos) && (p2 != std::string::npos))
				result = result.substr(std::max(p1, p2) + 1, result.length());
			else if (p1 != std::string::npos)
				result = result.substr(p1 + 1, result.length());
			else if (p2 != std::string::npos)
				result = result.substr(p2 + 1, result.length());
			return result;
		}

		static std::string getFileNameWithExt(const std::string &path)
		{
			std::string npath = normalizePath(path);

			std::string result = npath;
			size_t p1 = result.rfind('\\', result.length());
			size_t p2 = result.rfind('/', result.length());
			if ((p1 != std::string::npos) && (p2 != std::string::npos))
				result = result.substr(std::max(p1, p2) + 1, result.length());
			else if (p1 != std::string::npos)
				result = result.substr(p1 + 1, result.length());
			else if (p2 != std::string::npos)
				result = result.substr(p2 + 1, result.length());
			return result;
		}

		static std::string getFileExt(const std::string &path)
		{
			std::string npath = normalizePath(path);

			std::string result = npath;
			size_t i = result.rfind('.', result.length());
			if (i != std::string::npos)
			{
				result = result.substr(i + 1, result.length());
			}
			return result;
		}

		static bool isRelativePath(const std::string &path)
		{
			std::string npath = normalizePath(path);

			// Windows
			size_t i = npath.find(":");
			if (i != std::string::npos)
				return false;
			else if (npath[0] == '/')
				return false;
			return true;
		}

		static int makeDir(const std::string &path)
		{
			std::string npath = normalizePath(path);

			struct stat st;
			int status = 0;

			if (stat(path.c_str(), &st) != 0)
			{
#if WIN32
				status = _mkdir(path.c_str());
#else
				status = mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
#endif  
				if (status != 0 && errno != EEXIST)
					status = -1;
			}
			else if (!(S_IFDIR & st.st_mode))
			{
				errno = ENOTDIR;
				status = -1;
			}

			return status;
		}

		/** Make all subdirectories.
		*/
		static int makeDirs(const std::string &path)
		{
			char *pp;
			char *sp;
			int  status;
#ifdef WIN32
			char *copyOfPath = _strdup(path.c_str());
#else
			char *copyOfPath = strdup(path.c_str());
#endif

			status = 0;
			pp = copyOfPath;
			pp = pp + 3;		// Cut away Drive:
			while ((status == 0) && (((sp = strchr(pp, '/')) != 0) || ((sp = strchr(pp, '\\')) != 0)))
			{
				if (sp != pp)
				{
					*sp = '\0';
					status = makeDir(copyOfPath);
					*sp = '/';
				}
				pp = sp + 1;
			}
			if (status == 0)
				status = makeDir(path);
			free(copyOfPath);
			return status;
		}


		static std::string normalizePath(const std::string &path)
		{
			if (path.size() == 0)
				return path;
			std::string result = path;
			std::replace(result.begin(), result.end(), '\\', '/');
			std::vector<std::string> tokens;
			StringTools::tokenize(result, tokens, "/");
			unsigned int index = 0;
			while (index < tokens.size())
			{
				if ((tokens[index] == "..") && (index > 0))
				{
					tokens.erase(tokens.begin() + index - 1, tokens.begin() + index + 1);
					index-=2;
				}
				index++;
			}
			result = "";
			if (path[0] == '/')
				result = "/";
			result = result + tokens[0];
			for (unsigned int i = 1; i < tokens.size(); i++)
				result = result + "/" + tokens[i];

			return result;
		}

		static bool fileExists(const std::string& fileName)
		{
			if (FILE *file = fopen(fileName.c_str(), "r"))
			{
				fclose(file);
				return true;
			}
			else
				return false;
		}

		static std::string getProgramPath()
		{
			std::string executablePath;
			std::string dynamicModulePath;

#ifdef WIN32
			char exeBuffer[MAX_PATH] = { 0 };
			if (GetModuleFileNameA(NULL, exeBuffer, MAX_PATH))
				executablePath = exeBuffer;

			HMODULE moduleHandle = nullptr;
			if (GetModuleHandleExA(GET_MODULE_HANDLE_EX_FLAG_FROM_ADDRESS | GET_MODULE_HANDLE_EX_FLAG_UNCHANGED_REFCOUNT,
				reinterpret_cast<LPCSTR>(&getProgramPath), &moduleHandle) && moduleHandle)
			{
				char moduleBuffer[MAX_PATH] = { 0 };
				if (GetModuleFileNameA(moduleHandle, moduleBuffer, MAX_PATH))
					dynamicModulePath = moduleBuffer;
			}
#elif defined(__APPLE__)
			char exeBuffer[1024] = { 0 };
			uint32_t bufferSize = sizeof(exeBuffer);
			if (_NSGetExecutablePath(exeBuffer, &bufferSize) == 0)
				executablePath = exeBuffer;

			Dl_info moduleInfo{};
			if (dladdr(reinterpret_cast<void*>(&getProgramPath), &moduleInfo) && moduleInfo.dli_fname)
				dynamicModulePath = moduleInfo.dli_fname;
#else
			char exeBuffer[1024] = { 0 };
			char procPath[32];
			sprintf(procPath, "/proc/%d/exe", getpid());
			int exePathLength = std::min((int)readlink(procPath, exeBuffer, sizeof(exeBuffer) - 1), (int)(sizeof(exeBuffer) - 1));
			if (exePathLength > 0)
			{
				exeBuffer[exePathLength] = '\0';
				executablePath = exeBuffer;
			}

			Dl_info moduleInfo{};
			if (dladdr(reinterpret_cast<void*>(&getProgramPath), &moduleInfo) && moduleInfo.dli_fname)
				dynamicModulePath = moduleInfo.dli_fname;
#endif

			auto getParentDirectory = [](const std::string &path) -> std::string
			{
				size_t lastSeparator = path.find_last_of("\\/");
				return (lastSeparator != std::string::npos) ? path.substr(0, lastSeparator) : "";
			};

			// If running as a dynamic module / Python extension (module differs from executable)
			if (!dynamicModulePath.empty() && dynamicModulePath != executablePath)
			{
				std::string searchDirectory = getParentDirectory(dynamicModulePath);

				// Walk up the directory tree (up to 6 levels) to locate resources/data
				for (int level = 0; level < 6 && !searchDirectory.empty(); ++level)
				{
					if (isDirectory(searchDirectory + "/resources") || isDirectory(searchDirectory + "/data"))
						return searchDirectory;

					if (isDirectory(searchDirectory + "/bin/resources") || isDirectory(searchDirectory + "/bin/data"))
						return searchDirectory + "/bin";

					searchDirectory = getParentDirectory(searchDirectory);
				}

				return getParentDirectory(dynamicModulePath);
			}
			else
			{
				// Standalone executable fallback
				return getParentDirectory(executablePath);
			}
		}

		static bool copyFile(const std::string &source, const std::string &dest)
		{
			const size_t bufferSize = 8192;
			char buffer[bufferSize];
			size_t size;

			FILE* sourceFile = fopen(source.c_str(), "rb");
			FILE* destFile = fopen(dest.c_str(), "wb");

			if ((sourceFile == NULL) || (destFile == NULL))
				return false;

			while (size = fread(buffer, 1, bufferSize, sourceFile))
			{
				fwrite(buffer, 1, size, destFile);
			}

			fclose(sourceFile);
			fclose(destFile);

			return true;
		}

		static bool isFile(const std::string &path) 
		{
			struct stat st;
			if (!stat(path.c_str(), &st))
				return S_ISREG(st.st_mode);
			return false;
		}

		static bool isDirectory(const std::string &path)
		{
			struct stat st;
			if (!stat(path.c_str(), &st))
				return S_ISDIR(st.st_mode);
			return false;
		}

		static bool getFilesInDirectory(const std::string& path, std::vector<std::string> &res)
		{
#ifdef WIN32
			std::string p = path + "\\*";
			WIN32_FIND_DATA data;
			HANDLE hFind = FindFirstFile(p.c_str(), &data);
			if (hFind  != INVALID_HANDLE_VALUE) 
			{
				do 
				{
					res.push_back(data.cFileName);
				} 
				while (FindNextFile(hFind, &data) != 0);
				FindClose(hFind);
				return true;
			}
			return false;
#else
			DIR* dir = opendir(path.c_str());
			if (dir != NULL)
			{
				struct dirent *dp;
				while ((dp = readdir(dir)) != NULL) 
					res.push_back(dp->d_name);
				closedir(dir);
				return true;
			}
			return false;
#endif
		}

		
		/** Compute the MD5 hash of a file.
		*/
		static std::string getFileMD5(const std::string &filename)
		{
			std::ifstream file(filename);

			if (!file)
				std::cerr << "Cannot open file: " << filename << std::endl;
			else
			{
				MD5 context(file);
				char *md5hex = context.hex_digest();
				std::string res(md5hex);
				delete[] md5hex;
				return res;
			}
			return "";
		}

		/** Write the MD5 hash of a file to the md5File.  
		 */
		static bool writeMD5File(const std::string& fileName, const std::string& md5File)
		{
			std::ofstream fstream;
			fstream.open(md5File.c_str(), std::ios::out);
			if (fstream.fail())
			{
				std::cerr << "Failed to open file: " << md5File << "\n";
				return false;
			}
			std::string md5 = getFileMD5(fileName);
			if (md5 != "")
				fstream.write(md5.c_str(), md5.size());
			fstream.close();
			return true;
		}

		/** Compare an MD5 hash with the hash stored in an MD5 file.
		*/
		static bool checkMD5(const std::string& md5Hash, const std::string& md5File)
		{
			std::ifstream fstream;
			fstream.open(md5File.c_str(), std::ios::in);
			if (fstream.fail())
			{
				std::cerr << "Failed to open file: " << md5File << "\n";
				return false;
			}
			std::string str((std::istreambuf_iterator<char>(fstream)),
				std::istreambuf_iterator<char>());
			fstream.close();

			return str == md5Hash;
		}

#ifdef WIN32
		/** Open windows file dialog.\n
		* dialogType 0 = Open file dialog\n
		* dialogType 1 = Save file dialog\n
		*/
		static const std::string fileDialog(int dialogType,
			const std::string &initialDir,
			const std::string &filter)
		{
			std::string initDir = normalizePath(initialDir);
			std::replace(initDir.begin(), initDir.end(), '/', '\\');

			OPENFILENAME ofn;       // common dialog box structure
			char fileNameBuffer[512]; 
			fileNameBuffer[0] = '\0';			

			const std::string filterWithEscape = filter + '\0';
			
			ZeroMemory(&ofn, sizeof(ofn));
			ofn.lStructSize = sizeof(ofn);
			ofn.lpstrFile = fileNameBuffer;
			ofn.nMaxFile = sizeof(fileNameBuffer);
			ofn.lpstrFilter = filterWithEscape.c_str();
			ofn.nFilterIndex = 1;
			ofn.lpstrInitialDir = (LPSTR)initDir.c_str();
			ofn.Flags = OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST | OFN_NOCHANGEDIR;

			if (dialogType == 0)
			{
				if (GetOpenFileName(&ofn))
					return std::string(fileNameBuffer);
			}
			else
			{
				if (GetSaveFileName(&ofn))
					return std::string(fileNameBuffer);
			}
			return "";
		}
#endif
	};
}

#endif
