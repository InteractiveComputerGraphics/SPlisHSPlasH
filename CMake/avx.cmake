macro(set_avx_flags)
    set(AVX_FLAGS)

    include(CheckCXXSourceRuns)
    set(CMAKE_REQUIRED_FLAGS)
	
	# AVX
    if(MSVC AND NOT MSVC_VERSION LESS 1600)
        set(CMAKE_REQUIRED_FLAGS "/arch:AVX")
	elseif(UNIX OR MINGW)
		set(CMAKE_REQUIRED_FLAGS "-mavx")
	endif()
  
	check_cxx_source_runs("
		#include <immintrin.h>
		int main()
		{
			float v[8] = { 1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f };
			float r[8];
			__m256 first = _mm256_loadu_ps(v);
			__m256 second = _mm256_setr_ps(0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f, 0.5f);
			__m256 result = _mm256_add_ps(first, second);
			_mm256_storeu_ps(r, result);

			for( int i = 0; i < 8; i++)
				if ((i+1)*1.0f+0.5f != r[i])
					return -1;
			return 0;
		}"
		FOUND_AVX)

    # AVX2
    if(MSVC AND NOT MSVC_VERSION LESS 1800)
        set(CMAKE_REQUIRED_FLAGS "/arch:AVX2")
	elseif(UNIX OR MINGW)
		set(CMAKE_REQUIRED_FLAGS "-mavx2")
	endif()

	check_cxx_source_runs("
		#include <immintrin.h>
		int main()
		{
			int v[8] = {10, 20, 30, 40, 50, 60, 70, 80};
			int r[8];
			__m256i first = _mm256_loadu_si256((__m256i*)v);
			__m256i second = _mm256_set_epi32(5, 5, 5, 5, 5, 5, 5, 5);
			__m256i result = _mm256_add_epi32(first, second);
			_mm256_storeu_si256((__m256i*)r, result);

			for( int i = 0; i < 8; i++)
				if ((i+1)*10+5 != r[i])
					return -1;
			return 0;
		}"
		FOUND_AVX2)

    # set compiler flags
	if (FOUND_AVX2)
		if(MSVC)
            set(AVX_FLAGS "/arch:AVX2")
		elseif(UNIX OR MINGW)
			set(AVX_FLAGS "-mavx2;-mfma")
		endif()
	elseif (FOUND_AVX)
		if(MSVC)
            set(AVX_FLAGS "/arch:AVX")
		elseif(UNIX OR MINGW)
			set(AVX_FLAGS "-mavx")
		endif()
    endif()
endmacro()