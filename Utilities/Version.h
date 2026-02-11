#ifndef __Version_h__
#define __Version_h__

#define STRINGIZE_HELPER(x) #x
#define STRINGIZE(x) STRINGIZE_HELPER(x)
#define WARNING(desc) message(__FILE__ "(" STRINGIZE(__LINE__) ") : Warning: " #desc)

#define GIT_SHA1 "9b1e2156bd4dc8265ba8080036fc77c70bc5620b"
#define GIT_REFSPEC "refs/heads/Kee2023"
#define GIT_LOCAL_STATUS "CLEAN"

#define SPLISHSPLASH_VERSION "2.16.0"

#ifdef DL_OUTPUT

#endif

#endif
