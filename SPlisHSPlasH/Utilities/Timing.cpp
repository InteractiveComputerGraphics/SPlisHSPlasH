#include "Timing.h"

using namespace SPH;

int IDFactory::id = 0;
std::unordered_map<int, AverageTime> Timing::m_averageTimes;
std::stack<TimingHelper> Timing::m_timingStack;
bool Timing::m_dontPrintTimes = false;
unsigned int Timing::m_startCounter = 0;
unsigned int Timing::m_stopCounter = 0;


void SPH::Timing::reset()
{
	while (!m_timingStack.empty())
		m_timingStack.pop();
	m_averageTimes.clear();
	m_startCounter = 0;
	m_stopCounter = 0;
}

