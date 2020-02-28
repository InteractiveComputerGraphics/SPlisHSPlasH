add_library(zlib STATIC
	src/adler32.c
	src/compress.c
	src/crc32.c
	src/crc32.h
	src/deflate.c
	src/deflate.h
	src/gzio.c
	src/infback.c
	src/inffast.c
	src/inffast.h
	src/inffixed.h
	src/inflate.c
	src/inflate.h
	src/inftrees.c
	src/inftrees.h
	src/trees.c
	src/trees.h
	src/uncompr.c
	src/zconf.h
	src/zlib.h
	src/zutil.c
	src/zutil.h
)
target_include_directories(zlib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)
set_target_properties(zlib PROPERTIES FOLDER "External Dependencies")
