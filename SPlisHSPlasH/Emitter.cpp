#include "Emitter.h"
#include "SPHKernels.h"
#include <iostream>
#include "TimeManager.h"
#include "TimeStep.h"

using namespace SPH;


Emitter::Emitter(const unsigned int width, const unsigned int height,
	const Vector3r &pos, const Vector3r &dir, const Vector3r &initialVel, 
	const Real emitsPerSecond, const unsigned int type)
{	
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



void Emitter::emitParticles(TimeStep *timeStep, std::vector <unsigned int> &reusedParticles, unsigned int &indexReuse, unsigned int &numEmittedParticles)
{
	if (TimeManager::getCurrent()->getTime() < m_nextEmitTime)
		return;

	Vector3r axis1, axis2;
	Vector3r d = m_dir;
	getOrthogonalVectors(d, axis1, axis2);

	FluidModel *model = timeStep->getModel();

	const Real diam = 2.0*model->getParticleRadius();

	const Real startX = -0.5*m_width*diam;
	const Real startZ = -0.5*m_height*diam;

	if ((model->numActiveParticles() < model->numParticles()) ||
		(reusedParticles.size() > 0))
	{
		unsigned int indexNotReuse = model->numActiveParticles();
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

				if (index < model->numParticles())
				{
					model->getPosition(0, index) = (i*diam + startX)*axis1 + (j*diam + startZ)*axis2 + m_x;
					model->getVelocity(0, index) = m_v;

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
			model->setNumActiveParticles(model->numActiveParticles() + numEmittedParticles);
			timeStep->emittedParticles(model->numActiveParticles() - numEmittedParticles);
			model->getNeighborhoodSearch()->resize_point_set(0, &model->getPosition(0, 0)[0], model->numActiveParticles());
		}
	}
	else
	{
		if (model->numActiveParticles() < model->numParticles())
		{
			unsigned int index = model->numActiveParticles();
			for (unsigned int i = 0; i < m_width; i++)
			{
				for (unsigned int j = 0; j < m_height; j++)
				{
					if (index < model->numParticles())
					{
						model->getPosition(0, index) = (i*diam + startX)*axis1 + (j*diam + startZ)*axis2 + m_x;
						model->getVelocity(0, index) = m_v;
						numEmittedParticles++;
					}
					index++;
				}
			}
			model->setNumActiveParticles(model->numActiveParticles() + numEmittedParticles);
			timeStep->emittedParticles(model->numActiveParticles() - numEmittedParticles);
			model->getNeighborhoodSearch()->resize_point_set(0, &model->getPosition(0, 0)[0], model->numActiveParticles());
		}
	}

	m_nextEmitTime += 1.0 / m_emitsPerSecond;
	m_emitCounter++;
}



void Emitter::emitParticlesCircle(TimeStep *timeStep, std::vector <unsigned int> &reusedParticles, unsigned int &indexReuse, unsigned int &numEmittedParticles)
{
	if (TimeManager::getCurrent()->getTime() < m_nextEmitTime)
		return;

	Vector3r axis1, axis2;
	Vector3r d = m_dir;
	getOrthogonalVectors(d, axis1, axis2);

	FluidModel *model = timeStep->getModel();
	const Real diam = 2.0*model->getParticleRadius();

	const Real radius = (0.5 * (Real)m_width * diam);
	const Real radius2 = radius*radius;

	const Real startX = -0.5*(m_width - 1)*diam;
	const Real startZ = -0.5*(m_width - 1)*diam;

	if ((model->numActiveParticles() < model->numParticles()) ||
		(reusedParticles.size() > 0))
	{
		unsigned int indexNotReuse = model->numActiveParticles();
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

				if ((index < model->numParticles()) && (x*x + y*y <= radius2))
				{
					model->getPosition(0, index) = x*axis1 + y*axis2 + m_x;
					model->getVelocity(0, index) = m_v;

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
			model->setNumActiveParticles(model->numActiveParticles() + numEmittedParticles);
			timeStep->emittedParticles(model->numActiveParticles() - numEmittedParticles);
			model->getNeighborhoodSearch()->resize_point_set(0, &model->getPosition(0, 0)[0], model->numActiveParticles());
		}
	}
	else
	{
		if (model->numActiveParticles() < model->numParticles())
		{
			unsigned int index = model->numActiveParticles();
			for (unsigned int i = 0; i < m_width; i++)
			{
				for (unsigned int j = 0; j < m_width; j++)
				{
					const Real x = (i*diam + startX);
					const Real y = (j*diam + startZ);
					if ((index < model->numParticles()) && (x*x + y*y <= radius2))
					{
						model->getPosition(0, index) = x*axis1 + y*axis2 + m_x;
						model->getVelocity(0, index) = m_v;
						numEmittedParticles++;
						index++;
					}
				}
			}
			model->setNumActiveParticles(model->numActiveParticles() + numEmittedParticles);
			timeStep->emittedParticles(model->numActiveParticles() - numEmittedParticles);
			model->getNeighborhoodSearch()->resize_point_set(0, &model->getPosition(0, 0)[0], model->numActiveParticles());
		}
	}

	m_nextEmitTime += 1.0 / m_emitsPerSecond;
	m_emitCounter++;
}

void Emitter::step(TimeStep *timeStep, std::vector <unsigned int> &reusedParticles, unsigned int &indexReuse, unsigned int &numEmittedParticles)
{
	if (m_type == 1)
		emitParticlesCircle(timeStep, reusedParticles, indexReuse, numEmittedParticles);
	else
		emitParticles(timeStep, reusedParticles, indexReuse, numEmittedParticles);
}

