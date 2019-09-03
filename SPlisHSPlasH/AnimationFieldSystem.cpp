#include "AnimationFieldSystem.h"
#include "FluidModel.h"
#include "TimeStep.h"
#include "Utilities/Logger.h"
#include "Simulation.h"


using namespace SPH;


AnimationFieldSystem::AnimationFieldSystem() :
	m_fields()
{	
}

AnimationFieldSystem::~AnimationFieldSystem(void)
{
	for (size_t i = 0; i < m_fields.size(); i++)
		delete m_fields[i];
}

void AnimationFieldSystem::step()
{
	for (size_t i = 0; i < m_fields.size(); i++)
	{
		m_fields[i]->step();
	}
}

void AnimationFieldSystem::reset()
{
	for (size_t i = 0; i < m_fields.size(); i++)
	{
		m_fields[i]->reset();
	}
}

void AnimationFieldSystem::addAnimationField(
	const std::string &particleFieldName,
	const Vector3r &pos, const Matrix3r & rotation, const Vector3r &scale,
	const std::string expression[3], const unsigned int type)
{
	m_fields.push_back(new AnimationField(
		particleFieldName,
		pos, rotation, scale,
		expression, type));
}
