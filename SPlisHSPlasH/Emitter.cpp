#include "Emitter.h"
#include "SPHKernels.h"
#include <iostream>
#include "TimeManager.h"
#include "TimeStep.h"
#include "FluidModel.h"
#include "Simulation.h"

using namespace SPH;

Emitter::Emitter(FluidModel *model,
	const unsigned int width, const unsigned int height,
	const Vector3r &pos, const Matrix3r & rotation,
	const Real velocity,
	const unsigned int type)
	: m_model(model)
	, m_width(width)
	, m_height(height)
	, m_x(pos)
	, m_rotation(rotation)
	, m_velocity(velocity)
	, m_type(type)
	, m_nextEmitTime(0)
	, m_emitStartTime(0)
	, m_emitEndTime(std::numeric_limits<Real>::max())
	, m_emitCounter(0)
{
}

Emitter::~Emitter(void)
{
}

void Emitter::reset()
{
	m_nextEmitTime = m_emitStartTime;
	m_emitCounter = 0;	
}

Vector3r Emitter::getSize(const Real width, const Real height, const int type)
{
	Simulation *sim = Simulation::getCurrent();
	const Real radius = sim->getParticleRadius();
	const Real diam = static_cast<Real>(2.0)*radius;

	const Real supportRadius = sim->getSupportRadius();
	const Real animationMarginAround = diam;
	Vector3r size;
	if (type == 0)
	{
		size = {
			2 * supportRadius,
			height * diam + 2 * animationMarginAround,
			width * diam + 2 * animationMarginAround
		};
	}
	else
	{
		// height and radius of cylinder
		const Real h = 2 * supportRadius;
		const Real r = 0.5f * width * diam + animationMarginAround;
		size = { h, 2*r, 2*r };
	}

	return size;
}

void Emitter::emitParticles(std::vector <unsigned int> &reusedParticles, unsigned int &indexReuse, unsigned int &numEmittedParticles)
{
	TimeManager *tm = TimeManager::getCurrent();
	const Real t = tm->getTime();
	const Real timeStepSize = tm->getTimeStepSize();
	const Vector3r & emitDir = m_rotation.col(0);
	Vector3r emitVel = m_velocity * emitDir;
	Simulation *sim = Simulation::getCurrent();
	const Real radius = sim->getParticleRadius();
	const Real diam = static_cast<Real>(2.0)*radius;

	// shortly before the emitter starts, cleanup the emitter from particles
	if (t < m_emitStartTime || t > m_emitEndTime)
 		emitVel = emitDir * radius * 10 / 0.25;

	if (t >= m_emitStartTime - 0.25 && t <= m_emitEndTime)
	{
		// animate emitted particles
		const Vector3r & x0 = m_x;

		const Real animationMarginAhead = sim->getSupportRadius();
		const Vector3r size = getSize(m_width, m_height, m_type);
		const Vector3r halfSize = 0.5 * size;
		const Vector3r pos = x0 + 0.5f * animationMarginAhead * emitDir;

		const unsigned int nModels = sim->numberOfFluidModels();
		for (unsigned int m = 0; m < nModels; m++)
		{
			FluidModel *fm = sim->getFluidModel(m);
			const unsigned int numParticles = fm->numActiveParticles();
			#pragma omp parallel for schedule(static) default(shared)
			for (int i = 0; i < (int)numParticles; i++)
			{
				Vector3r &xi = fm->getPosition(i);
				if (inBox(xi, pos, m_rotation, halfSize))
				{
					fm->getVelocity(i) = emitVel;
					fm->getPosition(i) += timeStepSize * emitVel;
					fm->setParticleState(i, ParticleState::AnimatedByEmitter);
				}
			}
		}
	}
	if (t < m_nextEmitTime || t > m_emitEndTime)
	{
		return;
	}

	const Vector3r axisHeight = m_rotation.col(1);
	const Vector3r axisWidth = m_rotation.col(2);
	
	const Real startX = -static_cast<Real>(0.5)*(m_width  - 1)*diam;
	const Real startZ = -static_cast<Real>(0.5)*(m_height - 1)*diam;

	const Real dt = t - m_nextEmitTime;

	const Vector3r velocityOffset = dt * emitVel;
	const Vector3r offset = m_x + velocityOffset;

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
					m_model->getPosition(index) = (i*diam + startX)*axisWidth + (j*diam + startZ)*axisHeight + offset;
					m_model->getVelocity(index) = emitVel;
					m_model->setParticleState(index, ParticleState::AnimatedByEmitter);

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
						m_model->getPosition(index) = (i*diam + startX)*axisWidth + (j*diam + startZ)*axisHeight + offset;
						m_model->getVelocity(index) = emitVel;
						m_model->setParticleState(index, ParticleState::AnimatedByEmitter);
						numEmittedParticles++;
					}
					index++;
				}
			}
			m_model->setNumActiveParticles(m_model->numActiveParticles() + numEmittedParticles);
			sim->emittedParticles(m_model, m_model->numActiveParticles() - numEmittedParticles);
			sim->getNeighborhoodSearch()->resize_point_set(m_model->getPointSetIndex(), &m_model->getPosition(0)[0], m_model->numActiveParticles());
		}
	}

	m_nextEmitTime += diam / m_velocity;
	m_emitCounter++;
}



void Emitter::emitParticlesCircle(std::vector <unsigned int> &reusedParticles, unsigned int &indexReuse, unsigned int &numEmittedParticles)
{
	TimeManager *tm = TimeManager::getCurrent();
	const Real t = tm->getTime();
	const Real timeStepSize = tm->getTimeStepSize();
	const Vector3r & emitDir = m_rotation.col(0);
	Simulation *sim = Simulation::getCurrent();
	const Real particleRadius = sim->getParticleRadius();
	const Real diam = static_cast<Real>(2.0)*particleRadius;
	Vector3r emitVel = m_velocity * emitDir;

	// shortly before the emitter starts, cleanup the emitter from particles
	if (t < m_emitStartTime || t > m_emitEndTime)
		emitVel = emitDir * particleRadius * 10 / 0.25;

	if (t >= m_emitStartTime - 0.25 && t <= m_emitEndTime)
	{
		// animate emitted particles		
		const Vector3r & x0 = m_x;

		const Real animationMarginAhead = sim->getSupportRadius();
		const Vector3r size = getSize(m_width, m_height, m_type);
		const Real h = size[0];
 		const Real r = 0.5f * size[1];
		const Real r2 = r * r;
		const Vector3r pos = x0 + 0.5f * animationMarginAhead * emitDir;


		const unsigned int nModels = sim->numberOfFluidModels();
		for (unsigned int m = 0; m < nModels; m++)
		{
			FluidModel *fm = sim->getFluidModel(m);
			const unsigned int numParticles = fm->numActiveParticles();
 			#pragma omp parallel for schedule(static) default(shared)
 			for (int i = 0; i < (int)numParticles; i++)
 			{
 				Vector3r &xi = fm->getPosition(i);
 				if (inCylinder(xi, pos, m_rotation, h, r2))
 				{
 					fm->getVelocity(i) = emitVel;
 					fm->getPosition(i) += timeStepSize * emitVel;
 					fm->setParticleState(i, ParticleState::AnimatedByEmitter);
 				}
 			}
		}
 	}
	if (t < m_nextEmitTime || t > m_emitEndTime)
	{
		return;
	}

	Vector3r axisHeight = m_rotation.col(1);
	Vector3r axisWidth = m_rotation.col(2);

	const Real radius = (static_cast<Real>(0.5) * (Real)m_width * diam);
	const Real radius2 = radius*radius;

	const Real startX = -static_cast<Real>(0.5)*(m_width - 1)*diam;
	const Real startZ = -static_cast<Real>(0.5)*(m_width - 1)*diam;

	const Real dt = t - m_nextEmitTime;
	const Vector3r velocity = emitVel;
	const Vector3r velocityOffset = dt * velocity;
	const Vector3r offset = m_x + velocityOffset;

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
					m_model->getPosition(index) = x*axisWidth + y*axisHeight + offset;
					m_model->getVelocity(index) = velocity;
					m_model->setParticleState(index, ParticleState::AnimatedByEmitter);

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
						m_model->getPosition(index) = x*axisWidth + y*axisHeight + offset;
						m_model->getVelocity(index) = velocity;
						m_model->setParticleState(index, ParticleState::AnimatedByEmitter);
						numEmittedParticles++;
						index++;
					}
				}
			}
			m_model->setNumActiveParticles(m_model->numActiveParticles() + numEmittedParticles);
			sim->emittedParticles(m_model, m_model->numActiveParticles() - numEmittedParticles);
			sim->getNeighborhoodSearch()->resize_point_set(m_model->getPointSetIndex(), &m_model->getPosition(0)[0], m_model->numActiveParticles());
		}
	}

	m_nextEmitTime += diam / m_velocity;
	m_emitCounter++;
}

void Emitter::step(std::vector <unsigned int> &reusedParticles, unsigned int &indexReuse, unsigned int &numEmittedParticles)
{
	if (m_type == 1)
		emitParticlesCircle(reusedParticles, indexReuse, numEmittedParticles);
	else
		emitParticles(reusedParticles, indexReuse, numEmittedParticles);
}

void SPH::Emitter::saveState(BinaryFileWriter &binWriter)
{
	binWriter.write(m_nextEmitTime);
	binWriter.write(m_emitCounter);
}

void SPH::Emitter::loadState(BinaryFileReader &binReader)
{
	binReader.read(m_nextEmitTime);
	binReader.read(m_emitCounter);
}

