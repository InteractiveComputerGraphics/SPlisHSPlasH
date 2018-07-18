#include "Emitter.h"
#include "SPHKernels.h"
#include <iostream>
#include "TimeManager.h"
#include "TimeStep.h"
#include "FluidModel.h"
#include "Simulation.h"

using namespace SPH;


Emitter::Emitter(FluidModel *model, const unsigned int width, const unsigned int height,
	const Vector3r &pos, const Vector3r &dir, const Vector3r &initialVel, 
	const Real emitsPerSecond, const unsigned int type)
{	
	m_model = model;
	m_width = width;
	m_height = height;
	m_x = pos;
	m_dir = dir;
	m_v = initialVel;
	m_emitsPerSecond = emitsPerSecond;
	m_type = type;

	m_nextEmitTime = 0.0;
	m_emitCounter = 0;	
}

Emitter::~Emitter(void)
{
}

void Emitter::reset()
{
	m_nextEmitTime = 0.0;
	m_emitCounter = 0;	
}


/** Returns two orthogonal vectors to vec which are also orthogonal to each other.
*/
void Emitter::getOrthogonalVectors(const Vector3r &vec, Vector3r &x, Vector3r &y)
{
	// Get plane vectors x, y
	Vector3r v(1, 0, 0);

	// Check, if v has same direction as vec
	if (fabs(v.dot(vec)) > 0.999)
		v = Vector3r(0, 1, 0);

	x = vec.cross(v);
	y = vec.cross(x);
	x.normalize();
	y.normalize();
}



void Emitter::emitParticles(std::vector <unsigned int> &reusedParticles, unsigned int &indexReuse, unsigned int &numEmittedParticles)
{
	if (TimeManager::getCurrent()->getTime() < m_nextEmitTime)
		return;

	Vector3r axis1, axis2;
	Vector3r d = m_dir;
	getOrthogonalVectors(d, axis1, axis2);

	Simulation *sim = Simulation::getCurrent();

	const Real radius = sim->getValue<Real>(Simulation::PARTICLE_RADIUS);
	const Real diam = static_cast<Real>(2.0)*radius;

	const Real startX = -static_cast<Real>(0.5)*m_width*diam;
	const Real startZ = -static_cast<Real>(0.5)*m_height*diam;

	if ((m_model->numActiveParticles() < m_model->numParticles()) ||
		(reusedParticles.size() > 0))
	{
		unsigned int indexNotReuse = m_model->numActiveParticles();
		for (unsigned int i = 0; i < m_width; i++)
		{
			for (unsigned int j = 0; j < m_height; j++)
			{
				unsigned int index = 0;
				bool reused = false;
				if (indexReuse < reusedParticles.size())
				{
					index = reusedParticles[indexReuse];
					reused = true;
				}
				else
				{
					index = indexNotReuse;
				}

				if (index < m_model->numParticles())
				{
					m_model->getPosition(index) = (i*diam + startX)*axis1 + (j*diam + startZ)*axis2 + m_x;
					m_model->getVelocity(index) = m_v;

					if (reused)
					{
						indexReuse++;
					}
					else
					{
						numEmittedParticles++;
						indexNotReuse++;
					}
					index++;
				}
			}
		}

		if (numEmittedParticles != 0)
		{
			m_model->setNumActiveParticles(m_model->numActiveParticles() + numEmittedParticles);
			Simulation *sim = Simulation::getCurrent();
			sim->emittedParticles(m_model, m_model->numActiveParticles() - numEmittedParticles);
			sim->getNeighborhoodSearch()->resize_point_set(m_model->getPointSetIndex(), &m_model->getPosition(0)[0], m_model->numActiveParticles());
		}
	}
	else
	{
		if (m_model->numActiveParticles() < m_model->numParticles())
		{
			unsigned int index = m_model->numActiveParticles();
			for (unsigned int i = 0; i < m_width; i++)
			{
				for (unsigned int j = 0; j < m_height; j++)
				{
					if (index < m_model->numParticles())
					{
						m_model->getPosition(index) = (i*diam + startX)*axis1 + (j*diam + startZ)*axis2 + m_x;
						m_model->getVelocity(index) = m_v;
						numEmittedParticles++;
					}
					index++;
				}
			}
			m_model->setNumActiveParticles(m_model->numActiveParticles() + numEmittedParticles);
			Simulation *sim = Simulation::getCurrent();
			sim->emittedParticles(m_model, m_model->numActiveParticles() - numEmittedParticles);
			sim->getNeighborhoodSearch()->resize_point_set(m_model->getPointSetIndex(), &m_model->getPosition(0)[0], m_model->numActiveParticles());
		}
	}

	m_nextEmitTime += static_cast<Real>(1.0) / m_emitsPerSecond;
	m_emitCounter++;
}



void Emitter::emitParticlesCircle(std::vector <unsigned int> &reusedParticles, unsigned int &indexReuse, unsigned int &numEmittedParticles)
{
	if (TimeManager::getCurrent()->getTime() < m_nextEmitTime)
		return;

	Vector3r axis1, axis2;
	Vector3r d = m_dir;
	getOrthogonalVectors(d, axis1, axis2);

	Simulation *sim = Simulation::getCurrent();
	const Real r = sim->getValue<Real>(Simulation::PARTICLE_RADIUS);
	const Real diam = static_cast<Real>(2.0)*r;

	const Real radius = (static_cast<Real>(0.5) * (Real)m_width * diam);
	const Real radius2 = radius*radius;

	const Real startX = -static_cast<Real>(0.5)*(m_width - 1)*diam;
	const Real startZ = -static_cast<Real>(0.5)*(m_width - 1)*diam;

	if ((m_model->numActiveParticles() < m_model->numParticles()) ||
		(reusedParticles.size() > 0))
	{
		unsigned int indexNotReuse = m_model->numActiveParticles();
		for (unsigned int i = 0; i < m_width; i++)
		{
			for (unsigned int j = 0; j < m_width; j++)
			{
				const Real x = (i*diam + startX);
				const Real y = (j*diam + startZ);

				unsigned int index = 0;
				bool reused = false;
				if (indexReuse < reusedParticles.size())
				{
					index = reusedParticles[indexReuse];
					reused = true;
				}
				else
				{
					index = indexNotReuse;
				}

				if ((index < m_model->numParticles()) && (x*x + y*y <= radius2))
				{
					m_model->getPosition(index) = x*axis1 + y*axis2 + m_x;
					m_model->getVelocity(index) = m_v;

					if (reused)
					{
						indexReuse++;
					}
					else
					{
						numEmittedParticles++;
						indexNotReuse++;
					}
					index++;
				}
			}
		}

		if (numEmittedParticles != 0)
		{
			m_model->setNumActiveParticles(m_model->numActiveParticles() + numEmittedParticles);
			Simulation *sim = Simulation::getCurrent();
			sim->emittedParticles(m_model, m_model->numActiveParticles() - numEmittedParticles);
			sim->getNeighborhoodSearch()->resize_point_set(m_model->getPointSetIndex(), &m_model->getPosition(0)[0], m_model->numActiveParticles());
		}
	}
	else
	{
		if (m_model->numActiveParticles() < m_model->numParticles())
		{
			unsigned int index = m_model->numActiveParticles();
			for (unsigned int i = 0; i < m_width; i++)
			{
				for (unsigned int j = 0; j < m_width; j++)
				{
					const Real x = (i*diam + startX);
					const Real y = (j*diam + startZ);
					if ((index < m_model->numParticles()) && (x*x + y*y <= radius2))
					{
						m_model->getPosition(index) = x*axis1 + y*axis2 + m_x;
						m_model->getVelocity(index) = m_v;
						numEmittedParticles++;
						index++;
					}
				}
			}
			m_model->setNumActiveParticles(m_model->numActiveParticles() + numEmittedParticles);
			Simulation *sim = Simulation::getCurrent();
			sim->emittedParticles(m_model, m_model->numActiveParticles() - numEmittedParticles);
			sim->getNeighborhoodSearch()->resize_point_set(m_model->getPointSetIndex(), &m_model->getPosition(0)[0], m_model->numActiveParticles());
		}
	}

	m_nextEmitTime += static_cast<Real>(1.0) / m_emitsPerSecond;
	m_emitCounter++;
}

void Emitter::step(std::vector <unsigned int> &reusedParticles, unsigned int &indexReuse, unsigned int &numEmittedParticles)
{
	if (m_type == 1)
		emitParticlesCircle(reusedParticles, indexReuse, numEmittedParticles);
	else
		emitParticles(reusedParticles, indexReuse, numEmittedParticles);
}

