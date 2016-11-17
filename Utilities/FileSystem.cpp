#include "FileSystem.h"
#include "StringTools.h"
#include <vector>
#include <iostream>
#if WIN32
#include <direct.h>
#define NOMINMAX
#include "windows.h"
#else
#include <sys/stat.h>
#endif

using namespace SPH;


std::string FileSystem::getFilePath(const std::string &path)
{
	std::string npath =	normalizePath(path);

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

std::string FileSystem::getFileName(const std::string &path)
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
		result = result.substr(std::max(p1, p2)+1, result.length());
	else if (p1 != std::string::npos)
		result = result.substr(p1+1, result.length());
	else if (p2 != std::string::npos)
		result = result.substr(p2+1, result.length());
	return result;
}

bool FileSystem::isRelativePath(const std::string &path)
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

int FileSystem::makeDir(const std::string &path)
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

int FileSystem::makeDirs(const std::string &path)
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

std::string FileSystem::normalizePath(const std::string &path)
{
	std::string result = path;
	std::replace(result.begin(), result.end(), '\\', '/'); 
	std::vector<std::string> tokens;
	StringTools::tokenize(result, tokens, "/");
	unsigned int index = 0;
	while (index < tokens.size())
	{
		if ((tokens[index] == "..") && (index > 0))
		{
			tokens.erase(tokens.begin() + index-1, tokens.begin() + index + 1);
			index--;
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

bool FileSystem::fileExists(const std::string& fileName) 
{
	if (FILE *file = fopen(fileName.c_str(), "r")) 
	{
		fclose(file);
		return true;
	}
	else
		return false;
}

std::string FileSystem::getProgramPath()
{
	char buffer[1000];
#ifdef WIN32	
	GetModuleFileName(NULL, buffer, 1000);
#else
	char szTmp[32];
	sprintf(szTmp, "/proc/%d/exe", getpid());
	int bytes = std::min((int) readlink(szTmp, buffer, 1000), 999);
	buffer[bytes] = '\0';
#endif
	std::string::size_type pos = std::string(buffer).find_last_of("\\/");
	return std::string(buffer).substr(0, pos);

}
