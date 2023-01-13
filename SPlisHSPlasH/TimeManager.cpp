#include "TimeManager.h"

using namespace SPH;
using namespace GenParam;

int TimeManager::TIME_STEP_SIZE = -1;
TimeManager* TimeManager::current = 0;

TimeManager::TimeManager () 
{
	time = 0;
	h = static_cast<Real>(0.001);
}

TimeManager::~TimeManager () 
{
	current = 0;
}

void TimeManager::initParameters()
{
	ParameterObject::initParameters();

	TIME_STEP_SIZE = createNumericParameter("timeStepSize", "Time step size", &h);
	setGroup(TIME_STEP_SIZE, "General|General");
	setDescription(TIME_STEP_SIZE, "The initial time step size used for the time integration. If you use an adaptive time stepping, this size will change during the simulation.");
	static_cast<RealParameter*>(getParameter(TIME_STEP_SIZE))->setMinValue(static_cast<Real>(1e-9));
}

TimeManager* TimeManager::getCurrent ()
{
	if (current == 0)
	{
		current = new TimeManager ();
		current->initParameters();
	}
	return current;
}

void TimeManager::setCurrent (TimeManager* tm)
{
	current = tm;
}

bool TimeManager::hasCurrent()
{
	return (current != 0);
}

Real TimeManager::getTime()
{
	return time;
}

void TimeManager::setTime(Real t)
{
	time = t;
}

Real TimeManager::getTimeStepSize()
{
	return h;
}

void TimeManager::setTimeStepSize(Real tss)
{
	h = tss;
}

void SPH::TimeManager::saveState(BinaryFileWriter &binWriter)
{
	binWriter.write(time);
	binWriter.write(h);
}

void SPH::TimeManager::loadState(BinaryFileReader &binReader)
{
	binReader.read(time);
	binReader.read(h);
}
