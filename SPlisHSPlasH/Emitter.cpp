#include "Emitter.h"
#include "FluidModel.h"
#include "SPHKernels.h"
#include "Simulation.h"
#include "TimeManager.h"
#include "TimeStep.h"
#include "Utilities/Logger.h"
#include <iostream>

using namespace SPH;

Emitter::Emitter(FluidModel* model, const unsigned int width, const unsigned int height, const Vector3r& pos,
                 const Matrix3r& rotation, const Real velocity, const unsigned int type, const bool useBoundary,
                 const unsigned int velocityProfile)
    : m_model(model), m_width(width), m_x(pos), m_rotation(rotation), m_velocity(velocity), m_type(type),
      m_useBoundary(useBoundary), m_velocityProfile(velocityProfile)
{
    Simulation* sim = Simulation::getCurrent();
    m_depth = getDepth();

    if (m_type == 1) {
        // for cylindrical emitters, the height must not be smaller than the width, to spawn the initial particles.
        m_height = width;
    }
    else {
        m_height = height;
    }

    m_size = getSize(static_cast<Real>(m_width), static_cast<Real>(m_height), m_type);
}

Emitter::~Emitter(void) {}

void Emitter::reset() { m_emitCounter = 0; }

int Emitter::getDepth()
{
    // This is its own function to not repeat the calculation in the constructor but still have it available for the
    // static getSize().
    Simulation* sim = Simulation::getCurrent();
    return static_cast<int>(std::ceil(sim->getSupportRadius() / sim->getParticleRadius()));
}

Vector3r Emitter::getSize(const Real width, const Real height, const int type)
{

    const Simulation* sim = Simulation::getCurrent();
    const Real particleDiameter = 2 * sim->getParticleRadius();
    const Real depth = static_cast<Real>(getDepth());
    Vector3r size;

    // The emitter must be larger in emit direction, else particles may spawn just outside.
    if (type == 0) {
        // box emitter
        size = {depth + 1, height, width};
    }
    else if (type == 1) {
        // cylindrical emitter
        size = {depth + 1, width, width};
    }

    return size * particleDiameter;
}

Vector3r Emitter::getSizeExtraMargin(const Real width, const Real height, const int type)
{
    // The box or cylinder around the emitter needs some padding.
    const Simulation* sim = Simulation::getCurrent();
    // const BoundaryHandlingMethods* boundaryMethod = sim->getBoundaryHandlingMethod();
    Real extraMargin;
    if (type == 0) {
        // box emitter
        if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012) {
            extraMargin = static_cast<Real>(2.0);
        }
        else {
            extraMargin = static_cast<Real>(2.5);
        }
    }
    else if (type == 1) {
        // cylindrical emitter
        if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012) {
            extraMargin = static_cast<Real>(2.0);
        }
        else {
            extraMargin = static_cast<Real>(2.5);
        }
    }
    else {
        extraMargin = 0.0;
    }

    return getSize(width + extraMargin, height + extraMargin, type);
}

void Emitter::emitParticles(std::vector<unsigned int>& reusedParticles, unsigned int& indexReuse,
                            unsigned int& numEmittedParticles)
{
    TimeManager* tm = TimeManager::getCurrent();
    const Real t = tm->getTime();
    const Real timeStepSize = tm->getTimeStepSize();
    const Vector3r axisDepth = m_rotation.col(0);
    const Vector3r axisHeight = m_rotation.col(1);
    const Vector3r axisWidth = m_rotation.col(2);
    Vector3r bulkEmitVel = m_velocity * axisDepth;
    Simulation* sim = Simulation::getCurrent();
    const Real particleRadius = sim->getParticleRadius();
    const Real particleDiameter = static_cast<Real>(2.0) * particleRadius;

    unsigned int indexNextNewParticle = m_model->numActiveParticles();

    // fill the emitter at the start of the simulation
    // The emitter is filled starting from the plane the particles leave the emitter.
    // TODO: maybe this should be done in the constructor or somewhere else?
    const Vector3r startPos = m_x + m_depth * particleRadius * axisDepth - (m_width - 1) * particleRadius * axisWidth
        - (m_height - 1) * particleRadius * axisHeight;
    if (t == 0.0) {
        for (unsigned int i = 0; i < m_depth; i++) {
            for (unsigned int j = 0; j < m_width; j++) {
                for (unsigned int k = 0; k < m_height; k++) {
                    const Vector3r spawnPos = startPos - i * axisDepth * particleDiameter
                        + j * axisWidth * particleDiameter + k * axisHeight * particleDiameter;

                    // Spawn particles in a grid.
                    // Skip the particles outside when a cylindrical emitter is used.
                    if (m_type == 0
                        || (m_type == 1
                            && inCylinder(spawnPos, m_x, m_rotation, m_size[0], (m_size[1] * m_size[1]) / 4)))
                    {
                        m_model->setPosition(indexNextNewParticle, spawnPos);
                        m_model->setParticleState(indexNextNewParticle, ParticleState::AnimatedByEmitter);
                        m_model->setObjectId(indexNextNewParticle, m_objectId);

                        indexNextNewParticle++;
                        numEmittedParticles++;
                    }
                }
            }
        }
    }

    if (t >= m_emitStartTime && t < m_emitEndTime) {
        // emitter is active
        for (unsigned int i = 0; i < indexNextNewParticle; i++) {
            const Vector3r tempPos = m_model->getPosition(i);
            // Check if the particle is inside the emitter.
            const bool insideEmitter = (m_type == 0 && inBox(tempPos, m_x, m_rotation, 0.5 * m_size))
                || (m_type == 1 && inCylinder(tempPos, m_x, m_rotation, m_size[0], (m_size[1] * m_size[1]) / 4));

            if (insideEmitter) {
                // Advect ALL particles inside the emitter.
                //
                // TODO: Doesn't get all particles within the rigid body. Perhaps use getSizeExtraMargin()?
                Vector3r localEmitVel;
                if (m_type == 1 && m_velocityProfile != 0) {
                    // different velocity profiles for circular emitters
                    const Vector3r localPos = tempPos - m_x;
                    const Real r = sqrt(pow(axisWidth.dot(localPos), 2.0) + (pow(axisHeight.dot(localPos), 2.0)));
                    const Vector3r maxEmitVel =
                        (1.0 / (1.0 - (2.0 / (static_cast<Real>(m_velocityProfile) + 2.0)))) * bulkEmitVel;
                    const Real relativeRadius = r / (m_size[1] / 2.0);
                    localEmitVel = maxEmitVel * (1.0 - pow(relativeRadius, static_cast<Real>(m_velocityProfile)));
                }
                else {// box shaped emitter or uniform velocity profile
                    localEmitVel = bulkEmitVel;
                }

                m_model->setPosition(i, tempPos + timeStepSize * localEmitVel);
                m_model->setVelocity(i, localEmitVel);
            }
            if (!insideEmitter && m_model->getParticleState(i) == ParticleState::AnimatedByEmitter
                && m_model->getObjectId(i) == m_objectId)
            {
                // particle has left the emitter during last step
                m_model->setParticleState(i, ParticleState::Active);
                m_model->setObjectId(i, 0);
                // reuse or spawn a new particle upstream
                if (indexReuse < reusedParticles.size()) {
                    // reuse a particle
                    m_model->setPosition(reusedParticles[indexReuse], tempPos - axisDepth * particleDiameter * m_depth);
                    m_model->setParticleState(reusedParticles[indexReuse], ParticleState::AnimatedByEmitter);
                    m_model->setVelocity(reusedParticles[indexReuse], bulkEmitVel);
                    m_model->setObjectId(reusedParticles[indexReuse], m_objectId);
                    indexReuse++;
                }
                else if (m_model->numActiveParticles() < m_model->numParticles()) {
                    // spawn a new particle
                    m_model->setPosition(indexNextNewParticle, tempPos - axisDepth * particleDiameter * m_depth);
                    m_model->setParticleState(indexNextNewParticle, ParticleState::AnimatedByEmitter);
                    m_model->setObjectId(indexNextNewParticle, m_objectId);

                    indexNextNewParticle++;
                    numEmittedParticles++;
                }
                else {
                    LOG_INFO << "No particles left for the emitter to reuse or activate!";
                }
            }
        }
    }

    m_model->setNumActiveParticles(m_model->numActiveParticles() + numEmittedParticles);
    sim->emittedParticles(m_model, m_model->numActiveParticles() - numEmittedParticles);
    sim->getNeighborhoodSearch()->resize_point_set(m_model->getPointSetIndex(), &m_model->getPosition(0)[0],
                                                   m_model->numActiveParticles());
}

void Emitter::step(std::vector<unsigned int>& reusedParticles, unsigned int& indexReuse,
                   unsigned int& numEmittedParticles)
{
    emitParticles(reusedParticles, indexReuse, numEmittedParticles);
}

void SPH::Emitter::saveState(BinaryFileWriter& binWriter) { binWriter.write(m_emitCounter); }

void SPH::Emitter::loadState(BinaryFileReader& binReader) { binReader.read(m_emitCounter); }
