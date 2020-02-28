#include "SPlisHSPlasH/Common.h"
#include <iostream>
#include "Utilities/Timing.h"
#include "Utilities/Counting.h"
#include "PartioViewer.h"


// Enable memory leak detection
#ifdef _DEBUG
#ifndef EIGEN_ALIGN
	#define new DEBUG_NEW 
#endif
#endif

INIT_LOGGING
INIT_TIMING
INIT_COUNTING

using namespace SPH;
using namespace std;
using namespace Utilities;


// main 
int main( int argc, char **argv )
{
	REPORT_MEMORY_LEAKS;

	PartioViewer viewer;	
	return viewer.run(argc, argv);
}
