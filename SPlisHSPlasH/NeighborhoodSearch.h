#pragma once

//#define GPU_NEIGHBORHOOD_SEARCH
#if USE_DOUBLE
	#define CUNSEARCH_USE_DOUBLE_PRECISION
#endif

#ifdef GPU_NEIGHBORHOOD_SEARCH
	#include "cuNSearch.h"
	typedef cuNSearch::NeighborhoodSearch NeighborhoodSearch;
#else
	#include "CompactNSearch.h"
	typedef CompactNSearch::NeighborhoodSearch NeighborhoodSearch;
#endif
