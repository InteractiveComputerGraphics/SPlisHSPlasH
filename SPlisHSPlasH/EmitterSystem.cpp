#include "EmitterSystem.h"
#include "FluidModel.h"
#include "TimeStep.h"
#include "FluidModel.h"


using namespace SPH;


EmitterSystem::EmitterSystem()
{	
	m_numReusedParticles = 0;
	m_numberOfEmittedParticles = 0;
	m_reuseParticles = false;
	m_boxMin = Vector3r(-1, -1, -1); 
	m_boxMax = Vector3r(1, 1, 1);

	m_reusedParticles.reserve(m_maxParticlesToReusePerStep);

}

EmitterSystem::~EmitterSystem(void)
{
	for (size_t i = 0; i < m_emitters.size(); i++)
		delete m_emitters[i];

	std::cout << "Sum of emitted particles: " << m_numberOfEmittedParticles << "\n";
	std::cout << "Sum of reused particles: " << m_numReusedParticles << "\n";
}

void EmitterSystem::reuseParticles(FluidModel *model)
{
	if (m_reuseParticles)
	{
		m_reusedParticles.clear();
		for (unsigned int i = 0; i < model->numActiveParticles(); i++)
		{
			Vector3r &x = model->getPosition(0, i);
			if ((x[0] < m_boxMin[0]) || (x[1] < m_boxMin[1]) || (x[2] < m_boxMin[2]) ||
				(x[0] > m_boxMax[0]) || (x[1] > m_boxMax[1]) || (x[2] > m_boxMax[2]))
			{
				m_reusedParticles.push_back(i);
				model->getVelocity(0, i) *= 0.95;	// make particles slow so that they don't influence
													// the CFL condition
													//model->getPosition(0, i) = Vector3r(1000 + i, 1000, 1000);
			}
		}
	}
}

void EmitterSystem::step(TimeStep *timeStep)
{
	if (m_emitters.size() == 0)
		return;

	reuseParticles(timeStep->getModel());
	unsigned int indexReuse = 0;	
	for (size_t i = 0; i < m_emitters.size(); i++)
	{
		unsigned int numEmittedParticles = 0;
		m_emitters[i]->step(timeStep, m_reusedParticles, indexReuse, numEmittedParticles);
		m_numberOfEmittedParticles += numEmittedParticles;
	}
	m_numReusedParticles += indexReuse;
}

void EmitterSystem::reset()
{
	m_reusedParticles.clear();
	m_numReusedParticles = 0;
	m_numberOfEmittedParticles = 0;
	for (size_t i = 0; i < m_emitters.size(); i++)
	{
		m_emitters[i]->reset();
	}
}

void EmitterSystem::addEmitter(const unsigned int width, const unsigned int height,
	const Vector3r &pos, const Vector3r &dir, const Vector3r &initialVel,
	const Real emitsPerSecond, const unsigned int type)
{
	m_emitters.push_back(new Emitter(width, height,
		pos, dir, initialVel, emitsPerSecond, type));
}

void EmitterSystem::enableReuseParticles(const Vector3r &boxMin /*= Vector3r(-1, -1, -1)*/, const Vector3r &boxMax /*= Vector3r(1, 1, 1)*/)
{
	m_reuseParticles = true;
	m_boxMin = boxMin;
	m_boxMax = boxMax;
}

void EmitterSystem::disableReuseParticles()
{
	m_reuseParticles = false;
}
