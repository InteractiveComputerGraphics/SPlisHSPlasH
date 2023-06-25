#include "SPlisHSPlasH/BoundaryModel_Akinci2012.h"
#include "StrongCouplingBoundarySolver.h"
#include "TimeManager.h"
#include "DynamicRigidBody.h"
#include "Simulator/DynamicBoundarySimulator.h"
#include "SPlisHSPlasH/SPHKernels.h"

using namespace SPH;

StrongCouplingBoundarySolver* StrongCouplingBoundarySolver::current = nullptr;

StrongCouplingBoundarySolver::StrongCouplingBoundarySolver() :
	m_density(),
	m_pressureGrad(),
	m_v_s(),
	m_s(),
	m_pressure(),
	m_v_rr(),
	m_minus_rho_div_v_s(),
	m_minus_rho_div_v_rr(),
	m_diagonalElement(),
	m_artificialVolume(),
	m_lastPressure(),
	m_predictedVelocity(),
	m_predictedPosition(),
	m_v_rr_body(),
	m_omega_rr_body(),
	m_bodyContacts()
{
	m_maxIterations = 100;
	m_minIterations = 2;
	m_maxDensityDeviation = 0.001;

	m_kernelFct = CubicKernel::W;
	m_gradKernelFct = CubicKernel::gradW;
	m_W_zero = CubicKernel::W_zero();

	Simulation* sim = Simulation::getCurrent();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();

	m_restDensity = 1;
	m_density.resize(nBoundaries);
	m_v_s.resize(nBoundaries);
	m_s.resize(nBoundaries);
	m_pressure.resize(nBoundaries);
	m_v_rr.resize(nBoundaries);
	m_minus_rho_div_v_s.resize(nBoundaries);
	m_minus_rho_div_v_rr.resize(nBoundaries);
	m_diagonalElement.resize(nBoundaries);
	m_pressureGrad.resize(nBoundaries);
	m_artificialVolume.resize(nBoundaries);
	m_lastPressure.resize(nBoundaries);
	m_predictedVelocity.resize(nBoundaries);
	m_predictedPosition.resize(nBoundaries);
	m_v_rr_body.resize(nBoundaries, Vector3r::Zero());
	m_omega_rr_body.resize(nBoundaries, Vector3r::Zero());
	m_bodyContacts.resize(nBoundaries, 0);
	m_contactsAllBodies = 0;

	for (unsigned int i = 0; i < nBoundaries; i++) {
		BoundaryModel_Akinci2012* bm = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModel(i));
		m_density[i].resize(bm->numberOfParticles(), m_restDensity);
		m_v_s[i].resize(bm->numberOfParticles(), Vector3r::Zero());
		m_s[i].resize(bm->numberOfParticles(), 0.0);
		m_pressure[i].resize(bm->numberOfParticles(), 0.0);
		m_v_rr[i].resize(bm->numberOfParticles(), Vector3r::Zero());
		m_minus_rho_div_v_s[i].resize(bm->numberOfParticles(), 0.0);
		m_minus_rho_div_v_rr[i].resize(bm->numberOfParticles(), 0.0);
		m_diagonalElement[i].resize(bm->numberOfParticles(), 0.0);
		m_pressureGrad[i].resize(bm->numberOfParticles(), Vector3r::Zero());
		m_artificialVolume[i].resize(bm->numberOfParticles(), 0.0);
		m_lastPressure[i].resize(bm->numberOfParticles(), 0.0);
		m_predictedVelocity[i].resize(bm->numberOfParticles(), Vector3r::Zero());
		m_predictedPosition[i].resize(bm->numberOfParticles(), Vector3r::Zero());
		
		bm->addField({ "density", FieldType::Scalar, [this, i](const unsigned int j) -> Real* { return &getDensity(i, j); } });
		bm->addField({ "sourceTerm", FieldType::Scalar, [this, i](const unsigned int j)->Real* {return &getSourceTerm(i,j); } });
		bm->addField({ "v_s", FieldType::Vector3, [this,i](const unsigned int j)->Vector3r* {return &getV_s(i,j); } });
		bm->addField({ "v_rr", FieldType::Vector3, [this,i](const unsigned int j)->Vector3r* {return &getV_rr(i,j); } });
		bm->addField({ "pressure", FieldType::Scalar,[this,i](const unsigned int j)->Real* {return &getPressure(i,j); } });
		bm->addField({ "rigidbodyVelocity", FieldType::Vector3,[bm](const unsigned int j)->Vector3r* {return &(static_cast<DynamicRigidBody*>(bm->getRigidBodyObject())->getVelocity()); } });
		bm->addField({ "v_rr_body", FieldType::Vector3,[this,i](const unsigned int j)->Vector3r* {return &getV_rr_body(i); } });
	}

}


void SPH::StrongCouplingBoundarySolver::resize(unsigned int size) {
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();

	m_restDensity = 1000;
	m_density.resize(nBoundaries);
	m_v_s.resize(nBoundaries);
	m_s.resize(nBoundaries);
	m_pressure.resize(nBoundaries);
	m_v_rr.resize(nBoundaries);
	m_minus_rho_div_v_s.resize(nBoundaries);
	m_minus_rho_div_v_rr.resize(nBoundaries);
	m_diagonalElement.resize(nBoundaries);
	m_pressureGrad.resize(nBoundaries);
	m_artificialVolume.resize(nBoundaries);
	m_lastPressure.resize(nBoundaries);
	m_predictedVelocity.resize(nBoundaries);
	m_predictedPosition.resize(nBoundaries);
	m_v_rr_body.resize(nBoundaries, Vector3r::Zero());
	m_omega_rr_body.resize(nBoundaries, Vector3r::Zero());
	m_bodyContacts.resize(nBoundaries, 0);

	for (unsigned int i = 0; i < nBoundaries; i++) {
		BoundaryModel_Akinci2012* bm = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModel(i));
		m_density[i].resize(bm->numberOfParticles(), m_restDensity);
		m_v_s[i].resize(bm->numberOfParticles(), Vector3r::Zero());
		m_s[i].resize(bm->numberOfParticles(), 0.0);
		m_pressure[i].resize(bm->numberOfParticles(), 0.0);
		m_v_rr[i].resize(bm->numberOfParticles(), Vector3r::Zero());
		m_minus_rho_div_v_s[i].resize(bm->numberOfParticles(), 0.0);
		m_minus_rho_div_v_rr[i].resize(bm->numberOfParticles(), 0.0);
		m_diagonalElement[i].resize(bm->numberOfParticles(), 0.0);
		m_pressureGrad[i].resize(bm->numberOfParticles(), Vector3r::Zero());
		m_artificialVolume[i].resize(bm->numberOfParticles(), 0.0);
		m_lastPressure[i].resize(bm->numberOfParticles(), 0.0);
		m_predictedVelocity[i].resize(bm->numberOfParticles(), Vector3r::Zero());
		m_predictedPosition[i].resize(bm->numberOfParticles(), Vector3r::Zero());
	}
}


void SPH::StrongCouplingBoundarySolver::reset() {
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();

	for (unsigned int i = 0; i < nBoundaries; i++) {
		BoundaryModel_Akinci2012* bm = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModel(i));
		for (int j = 0; j < bm->numberOfParticles(); j++) {
			m_v_s[i][j] = Vector3r::Zero();
			m_s[i][j] = 0.0;
			m_pressure[i][j] = 0.0;
			m_v_rr[i][j] = Vector3r::Zero();
			m_minus_rho_div_v_s[i][j] = 0.0;
			m_minus_rho_div_v_rr[i][j] = 0.0;
			m_diagonalElement[i][j] = 0.0;
			m_pressureGrad[i][j] = Vector3r::Zero();
			m_artificialVolume[i][j] = 0.0;
			m_lastPressure[i][j] = 0.0;
			m_predictedVelocity[i][j] = Vector3r::Zero();
			m_predictedPosition[i][j] = Vector3r::Zero();
		}
		m_v_rr_body[i] = Vector3r::Zero();
		m_omega_rr_body[i] = Vector3r::Zero();
		m_bodyContacts[i] = 0;
	}
}

StrongCouplingBoundarySolver::~StrongCouplingBoundarySolver() {
	current = nullptr;
	for (unsigned int i = 0; i < Simulation::getCurrent()->numberOfBoundaryModels(); i++) {
		m_density[i].clear();
		m_v_s[i].clear();
		m_s[i].clear();
		m_pressure[i].clear();
		m_v_rr[i].clear();
		m_minus_rho_div_v_s[i].clear();
		m_minus_rho_div_v_rr[i].clear();
		m_diagonalElement[i].clear();
		m_pressureGrad[i].clear();
		m_artificialVolume[i].clear();
		m_lastPressure[i].clear();
		m_predictedVelocity[i].clear();
		m_predictedPosition[i].clear();
	}
	m_density.clear();
	m_v_s.clear();
	m_s.clear();
	m_pressure.clear();
	m_v_rr.clear();
	m_minus_rho_div_v_s.clear();
	m_minus_rho_div_v_rr.clear();
	m_diagonalElement.clear();
	m_pressureGrad.clear();
	m_artificialVolume.clear();
	m_lastPressure.clear();
	m_predictedVelocity.clear();
	m_predictedPosition.clear();
	m_v_rr_body.clear();
	m_omega_rr_body.clear();
	m_bodyContacts.clear();
}

StrongCouplingBoundarySolver* StrongCouplingBoundarySolver::getCurrent() {
	if (current == nullptr) {
		current = new StrongCouplingBoundarySolver();
	} else {
		if (Simulation::getCurrent()->numberOfBoundaryModels() != current->m_density.size()) {
			current->resize(Simulation::getCurrent()->numberOfBoundaryModels());
		}
	}
	return current;
}

void StrongCouplingBoundarySolver::performNeighborhoodSearchSort() {
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfBoundaryModels();

	for (unsigned int i = 0; i < nModels; i++) {
		BoundaryModel_Akinci2012* bm = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModel(i));

		const unsigned int numPart = bm->numberOfParticles();
		if (numPart != 0) {
			auto const& d = sim->getNeighborhoodSearch()->point_set(bm->getPointSetIndex());
			d.sort_field(&m_density[i][0]);
			d.sort_field(&m_v_s[i][0]);
			d.sort_field(&m_s[i][0]);
			d.sort_field(&m_pressure[i][0]);
			d.sort_field(&m_v_rr[i][0]);
			d.sort_field(&m_minus_rho_div_v_rr[i][0]);
			d.sort_field(&m_diagonalElement[i][0]);
			d.sort_field(&m_pressureGrad[i][0]);
			d.sort_field(&m_artificialVolume[i][0]);
			d.sort_field(&m_lastPressure[i][0]);
			d.sort_field(&m_predictedVelocity[i][0]);
			d.sort_field(&m_predictedPosition[i][0]);
		}
	}
}

void StrongCouplingBoundarySolver::computeContacts() {
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nFluids = sim->numberOfFluidModels();
	m_contactsAllBodies = 0;
	for (unsigned int boundaryPointSetIndex = nFluids; boundaryPointSetIndex < sim->numberOfPointSets(); boundaryPointSetIndex++) {
		unsigned int contacts = 0;
		BoundaryModel_Akinci2012* bm = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModelFromPointSet(boundaryPointSetIndex));
		int bmIndex = boundaryPointSetIndex - nFluids;
		if (bm->getRigidBodyObject()->isDynamic()) {
            #pragma omp parallel default(shared)
			{
                #pragma omp for schedule(static) reduction(+:contacts)
				for (int r = 0; r < bm->numberOfParticles(); r++) {
					for (unsigned int pid = nFluids; pid < sim->numberOfPointSets(); pid++) {
						if (pid != boundaryPointSetIndex) {
							if (sim->numberOfNeighbors(boundaryPointSetIndex, pid, r) > 0) {
								contacts++;
								break;
							}
						}
					}
				}
			}
		}
		m_bodyContacts[bmIndex] = contacts;
		m_contactsAllBodies += contacts;
	}
}


void StrongCouplingBoundarySolver::computeDensityAndVolume() {
	Simulation* sim = Simulation::getCurrent();
	TimeManager* tm = TimeManager::getCurrent();
	const Real dt = tm->getTimeStepSize();
	const unsigned int nFluids = sim->numberOfFluidModels();

	// Compute density and artificial volume for all particles
	for (unsigned int boundaryPointSetIndex = nFluids; boundaryPointSetIndex < sim->numberOfPointSets(); boundaryPointSetIndex++) {
		BoundaryModel_Akinci2012* bm = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModelFromPointSet(boundaryPointSetIndex));
		int bmIndex = boundaryPointSetIndex - nFluids;
		RigidBodyObject* rb = bm->getRigidBodyObject();
        #pragma omp parallel default(shared)
		{
            #pragma omp for schedule(static)  
			for (int r = 0; r < bm->numberOfParticles(); r++) {
				if (rb->isDynamic()) {
					// compute density for particle r
					Real particleDensity = getRestDensity() * bm->getVolume(r) * W_zero();
					// iterate over all rigid bodies
					for (unsigned int pid = nFluids; pid < sim->numberOfPointSets(); pid++) {
						BoundaryModel_Akinci2012* bm_neighbor = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModelFromPointSet(pid));
						for (unsigned int j = 0; j < sim->numberOfNeighbors(boundaryPointSetIndex, pid, r); j++) {
							const unsigned int k = sim->getNeighbor(boundaryPointSetIndex, pid, r, j);
							particleDensity += getRestDensity() * bm->getVolume(r) * W(bm->getPosition(r) - bm_neighbor->getPosition(k));
						}
					}
					setDensity(bmIndex, r, std::max(particleDensity, getRestDensity()));
					//gamma = 0.7, see the paper
					const Real gamma = static_cast<Real>(0.7); 
					if (getDensity(bmIndex, r) > getRestDensity()) {						
						setArtificialVolume(bmIndex, r, gamma * getRestDensity() * bm->getVolume(r) / getDensity(bmIndex, r));
					} else {
						setArtificialVolume(bmIndex, r, gamma * bm->getVolume(r));
					}
				} else {
					setArtificialVolume(bmIndex, r, bm->getVolume(r));
				}
			}
		}
	}
}

void StrongCouplingBoundarySolver::computeV_s() {
	Simulation* sim = Simulation::getCurrent();
	TimeManager* tm = TimeManager::getCurrent();
	const Real dt = tm->getTimeStepSize();
	const unsigned int nFluids = sim->numberOfFluidModels();

	// Compute v_s for all particles
	for (unsigned int boundaryPointSetIndex = nFluids; boundaryPointSetIndex < sim->numberOfPointSets(); boundaryPointSetIndex++) {
		BoundaryModel_Akinci2012* bm = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModelFromPointSet(boundaryPointSetIndex));
		DynamicRigidBody* rb = dynamic_cast<DynamicRigidBody*>(bm->getRigidBodyObject());
		if (rb != nullptr && rb->isDynamic()) {
            #pragma omp parallel default(shared)
			{
                #pragma omp for schedule(static)  
				for (int r = 0; r < bm->numberOfParticles(); r++) {
					Vector3r gravForce = rb->getMass() * Vector3r(sim->getVecValue<Real>(Simulation::GRAVITATION));
					// the rb->getForce() contains only fluid-rigid forces
					Vector3r F_R = rb->getForce() + gravForce;
					Vector3r pos = bm->getPosition(r) - rb->getPosition();
					// compute v_s
					setV_s(boundaryPointSetIndex - nFluids, r, rb->getVelocity() + dt * rb->getInvMass() * F_R +
						(rb->getAngularVelocity() + rb->getInertiaTensorInverseW() * (dt * rb->getTorque() + dt * (rb->getInertiaTensorW() * rb->getAngularVelocity()).cross(rb->getAngularVelocity()))).cross(pos));
					// compute the predicted velocity first here, since for cases without rigid-rigid contacts, the computeSourceTermRHS will not be called and V_rr is 0
					// but we still need the predicted velocity for fluid-rigid strong coupling
					setPredictedVelocity(boundaryPointSetIndex - nFluids, r, getV_s(boundaryPointSetIndex - nFluids, r));
					setPredictedPosition(boundaryPointSetIndex - nFluids, r, bm->getPosition(r) + (bm->getVelocity(r) + getV_s(boundaryPointSetIndex - nFluids, r)) / 2 * dt);
				}
			}
		}
	}
}

void StrongCouplingBoundarySolver::computeSourceTerm() {
	Simulation* sim = Simulation::getCurrent();
	TimeManager* tm = TimeManager::getCurrent();
	const Real dt = tm->getTimeStepSize();
	const unsigned int nFluids = sim->numberOfFluidModels();

	computeV_s();

	// compute source term s for all particles
	for (unsigned int boundaryPointSetIndex = nFluids; boundaryPointSetIndex < sim->numberOfPointSets(); boundaryPointSetIndex++) {
		BoundaryModel_Akinci2012* bm = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModelFromPointSet(boundaryPointSetIndex));
		int bmIndex = boundaryPointSetIndex - nFluids;
		if (bm->getRigidBodyObject()->isDynamic()) {
            #pragma omp parallel default(shared)
			{
                #pragma omp for schedule(static)  
				for (int r = 0; r < bm->numberOfParticles(); r++) {
					Real rho_div_v_s = 0;
					// iterate over all rigid bodies except bm, since the divergence of particles in the same rigid body should be 0.
					for (unsigned int pid = nFluids; pid < sim->numberOfPointSets(); pid++) {
						if (boundaryPointSetIndex != pid) {
							BoundaryModel_Akinci2012* bm_neighbor = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModelFromPointSet(pid));
							int neighborIndex = pid - nFluids;
							for (unsigned int j = 0; j < sim->numberOfNeighbors(boundaryPointSetIndex, pid, r); j++) {
								const unsigned int k = sim->getNeighbor(boundaryPointSetIndex, pid, r, j);
								// sum up divergence of v_s
								rho_div_v_s += getArtificialVolume(neighborIndex, k) * getDensity(neighborIndex, k) * (getV_s(neighborIndex, k) - getV_s(bmIndex, r)).dot(gradW(bm->getPosition(r) - bm_neighbor->getPosition(k)));
							}
						}
					}
					Real s = (getRestDensity() - getDensity(bmIndex, r)) / dt + rho_div_v_s;
					setSourceTerm(bmIndex, r, s);
					// used to computed the density deviation
					setMinus_rho_div_v_s(bmIndex, r, -rho_div_v_s);
				}
			}
		}
	}
}

void StrongCouplingBoundarySolver::computeDiagonalElement() {
	Simulation* sim = Simulation::getCurrent();
	TimeManager* tm = TimeManager::getCurrent();
	const Real dt = tm->getTimeStepSize();
	const unsigned int nFluids = sim->numberOfFluidModels();
	for (unsigned int pointSetIndex_r = nFluids; pointSetIndex_r < sim->numberOfPointSets(); pointSetIndex_r++) {
		BoundaryModel_Akinci2012* bm_r = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModelFromPointSet(pointSetIndex_r));
		int index_r = pointSetIndex_r - nFluids;
		DynamicRigidBody* rb = dynamic_cast<DynamicRigidBody*>(bm_r->getRigidBodyObject());
		if (rb != nullptr && rb->isDynamic()) {
            #pragma omp parallel default(shared)
			{
                #pragma omp for schedule(static)      
				for (int r = 0; r < bm_r->numberOfParticles(); r++) {
					Vector3r pos_r = bm_r->getPosition(r) - rb->getPosition();
					Vector3r grad_p_b = Vector3r::Zero();
					const Real density_r2 = getDensity(index_r, r) * getDensity(index_r, r);
					// rk are neighboring rigid particles of r of other rigid bodies k
					for (unsigned int pointSetIndex_rk = nFluids; pointSetIndex_rk < sim->numberOfPointSets(); pointSetIndex_rk++) {
						if (pointSetIndex_rk != pointSetIndex_r) {
							BoundaryModel_Akinci2012* bm_rk = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModelFromPointSet(pointSetIndex_rk));
							int index_rk = pointSetIndex_rk - nFluids;
							DynamicRigidBody* rb_rk = static_cast<DynamicRigidBody*>(bm_rk->getRigidBodyObject());
							for (unsigned int n = 0; n < sim->numberOfNeighbors(pointSetIndex_r, pointSetIndex_rk, r); n++) {
								const unsigned int rk = sim->getNeighbor(pointSetIndex_r, pointSetIndex_rk, r, n);
								grad_p_b += getArtificialVolume(index_rk, rk) * getDensity(index_rk, rk) / density_r2 * gradW(bm_r->getPosition(r) - bm_rk->getPosition(rk));
							}
						}
					}
					grad_p_b *= getDensity(index_r, r);
					Vector3r v_b = -dt * rb->getInvMass() * getArtificialVolume(index_r, r) * grad_p_b;
					Vector3r omega_b = -dt * rb->getInertiaTensorInverseW() * getArtificialVolume(index_r, r) * pos_r.cross(grad_p_b);
					Vector3r v_b_r = v_b + omega_b.cross(pos_r);
					Real b_r = 0;
					// for all neighboring rigid bodies
					for (unsigned int pointSetIndex_rk = nFluids; pointSetIndex_rk < sim->numberOfPointSets(); pointSetIndex_rk++) {
						if (pointSetIndex_rk != pointSetIndex_r) {
							Vector3r v_b_k_body = Vector3r::Zero();
							Vector3r omega_b_k_body = Vector3r::Zero();
							BoundaryModel_Akinci2012* bm_rk = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModelFromPointSet(pointSetIndex_rk));
							int index_rk = pointSetIndex_rk - nFluids;
							DynamicRigidBody* rb_rk = static_cast<DynamicRigidBody*>(bm_rk->getRigidBodyObject());
							for (unsigned int n = 0; n < sim->numberOfNeighbors(pointSetIndex_r, pointSetIndex_rk, r); n++) {
								const unsigned int rk = sim->getNeighbor(pointSetIndex_r, pointSetIndex_rk, r, n);
								Vector3r pos_rk = bm_rk->getPosition(rk) - rb_rk->getPosition();
								Vector3r grad_p_b_rkr = getDensity(index_rk, rk) * getArtificialVolume(index_r, r) * getDensity(index_r, r) / density_r2 * gradW(bm_rk->getPosition(rk) - bm_r->getPosition(r));
								v_b_k_body += -dt * rb_rk->getInvMass() * getArtificialVolume(index_rk, rk) * grad_p_b_rkr;
								omega_b_k_body += -dt * rb_rk->getInertiaTensorInverseW() * getArtificialVolume(index_rk, rk) * pos_rk.cross(grad_p_b_rkr);
							}
							// divergence
							Real sum_rk = 0;
							for (unsigned int n = 0; n < sim->numberOfNeighbors(pointSetIndex_r, pointSetIndex_rk, r); n++) {
								const unsigned int rk = sim->getNeighbor(pointSetIndex_r, pointSetIndex_rk, r, n);
								Vector3r pos_rk = bm_rk->getPosition(rk) - rb_rk->getPosition();
								Vector3r v_b_rk = v_b_k_body + omega_b_k_body.cross(pos_rk);
								sum_rk += getArtificialVolume(index_rk, rk) * getDensity(index_rk, rk) * (v_b_rk - v_b_r).dot(gradW(bm_r->getPosition(r) - bm_rk->getPosition(rk)));
							}
							b_r += sum_rk;
						}
					}
					b_r = -b_r;
					//setDiagonalElement(index_r, r, b_r < 0 ? b_r : 0);
					setDiagonalElement(index_r, r, b_r);
				}
			}
		}
	}
}

void StrongCouplingBoundarySolver::computeSourceTermRHS() {
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nFluids = sim->numberOfFluidModels();
	for (unsigned int boundaryPointSetIndex = nFluids; boundaryPointSetIndex < sim->numberOfPointSets(); boundaryPointSetIndex++) {
		computeSourceTermRHSForBody(boundaryPointSetIndex);
	}
}

void SPH::StrongCouplingBoundarySolver::computeSourceTermRHSForBody(const unsigned int& boundaryPointSetIndex) {
	Simulation* sim = Simulation::getCurrent();
	TimeManager* tm = TimeManager::getCurrent();
	const Real dt = tm->getTimeStepSize();
	const unsigned int nFluids = sim->numberOfFluidModels();
	BoundaryModel_Akinci2012* bm = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModelFromPointSet(boundaryPointSetIndex));

	// compute pressure gradient, v_rr and omega_rr for rigid bodies
	int bmIndex = boundaryPointSetIndex - nFluids;
	DynamicRigidBody* rb = dynamic_cast<DynamicRigidBody*>(bm->getRigidBodyObject());
	Real vx = 0;
	Real vy = 0;
	Real vz = 0;
	Real omegax = 0;
	Real omegay = 0;
	Real omegaz = 0;
	if (rb != nullptr && rb->isDynamic()) {
        #pragma omp parallel default(shared)
		{
            #pragma omp for schedule(static) reduction(+:vx) reduction(+:vy) reduction(+:vz) reduction(+:omegax) reduction(+:omegay) reduction(+:omegaz)                                                    
			for (int r = 0; r < bm->numberOfParticles(); r++) {
				Vector3r pressureGrad_r = Vector3r::Zero();
				const Real density_r = getDensity(bmIndex, r);
				const Real pressure_r = getPressure(bmIndex, r);
				for (unsigned int pid = nFluids; pid < sim->numberOfPointSets(); pid++) {
					if (pid != boundaryPointSetIndex) {
						BoundaryModel_Akinci2012* bm_neighbor = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModelFromPointSet(pid));
						int neighborIndex = pid - nFluids;
						for (unsigned int j = 0; j < sim->numberOfNeighbors(boundaryPointSetIndex, pid, r); j++) {
							const unsigned int k = sim->getNeighbor(boundaryPointSetIndex, pid, r, j);
							const Real density_k = getDensity(neighborIndex, k);
							const Real volume_k = getArtificialVolume(neighborIndex, k);
							const Real pressure_k = getPressure(neighborIndex, k);
							pressureGrad_r += volume_k * density_k * (pressure_r / (density_r * density_r) + pressure_k / (density_k * density_k)) * gradW(bm->getPosition(r) - bm_neighbor->getPosition(k));
						}
					}
				}
				pressureGrad_r *= density_r;
				setPressureGrad(bmIndex, r, pressureGrad_r);

				Vector3r v_rr_body_tmp = -dt * rb->getInvMass() * getArtificialVolume(bmIndex, r) * getPressureGrad(bmIndex, r);
				vx += v_rr_body_tmp.x();
				vy += v_rr_body_tmp.y();
				vz += v_rr_body_tmp.z();
				Vector3r pos = bm->getPosition(r) - rb->getPosition();
				Vector3r omega_rr_body_tmp = -dt * rb->getInertiaTensorInverseW() * getArtificialVolume(bmIndex, r) * pos.cross(getPressureGrad(bmIndex, r));
				omegax += omega_rr_body_tmp.x();
				omegay += omega_rr_body_tmp.y();
				omegaz += omega_rr_body_tmp.z();
			}
		}
	}
	Vector3r v_rr_body = Vector3r(vx, vy, vz);
	Vector3r omega_rr_body = Vector3r(omegax, omegay, omegaz);
	setV_rr_body(bmIndex, v_rr_body);
	setOmega_rr_body(bmIndex, omega_rr_body);
	
	// compute v_rr and predicted velocity (v_rr+v_s) for all particles
	if (rb->isDynamic()) {
        #pragma omp parallel default(shared)
		{
            #pragma omp for schedule(static)  
			for (int r = 0; r < bm->numberOfParticles(); r++) {
				Vector3r pos = bm->getPosition(r) - rb->getPosition();
				setV_rr(bmIndex, r, getV_rr_body(bmIndex) + getOmega_rr_body(bmIndex).cross(pos));
				setPredictedVelocity(bmIndex, r, getV_rr(bmIndex, r) + getV_s(bmIndex, r));
				setPredictedPosition(bmIndex, r, bm->getPosition(r) + (bm->getVelocity(r) + getPredictedVelocity(bmIndex, r)) / 2 * dt);

				Real minus_rho_div_v_rr = 0.0;
				const Vector3r v_rr_r = getV_rr(bmIndex, r);
				for (unsigned int pid = nFluids; pid < sim->numberOfPointSets(); pid++) {
					// divergence of particles in the same body should be 0;
					if (boundaryPointSetIndex != pid) {
						BoundaryModel_Akinci2012* bm_neighbor = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModelFromPointSet(pid));
						int neighborIndex = pid - nFluids;
						for (unsigned int j = 0; j < sim->numberOfNeighbors(boundaryPointSetIndex, pid, r); j++) {
							const unsigned int k = sim->getNeighbor(boundaryPointSetIndex, pid, r, j);
							const Real density_k = getDensity(neighborIndex, k);
							const Real volume_k = getArtificialVolume(neighborIndex, k);
							const Vector3r v_rr_k = getV_rr(pid - nFluids, k);
							minus_rho_div_v_rr += volume_k * density_k * (v_rr_k - v_rr_r).dot(gradW(bm->getPosition(r) - bm_neighbor->getPosition(k)));
						}
					}
				}
				minus_rho_div_v_rr = -minus_rho_div_v_rr;
				setSourceTermRHS(bmIndex, r, minus_rho_div_v_rr);
			}
		}
	}
}


void StrongCouplingBoundarySolver::pressureSolveIteration(Real& avgDensityDeviation) {
	Simulation* sim = Simulation::getCurrent();
	const Real dt = TimeManager::getCurrent()->getTimeStepSize();
	const unsigned int nFluids = sim->numberOfFluidModels();
	Real avgDensityDev = 0.0;
	for (unsigned int boundaryPointSetIndex = nFluids; boundaryPointSetIndex < sim->numberOfPointSets(); boundaryPointSetIndex++) {
		BoundaryModel_Akinci2012* bm = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModelFromPointSet(boundaryPointSetIndex));
		int bmIndex = boundaryPointSetIndex - nFluids;
		if (bm->getRigidBodyObject()->isDynamic()) {
			int numContactsBody = getBodyContacts(bmIndex);
			computeSourceTermRHSForBody(boundaryPointSetIndex);
			// beta_r_RJ
			Real relaxation = 0.5 / numContactsBody;
            #pragma omp parallel default(shared)
			{
                #pragma omp for schedule(static) reduction(+:avgDensityDev)
				for (int r = 0; r < bm->numberOfParticles(); r++) {
					if (std::abs(getDiagonalElement(bmIndex, r)) > 1e-10) {
						Real pressureNextIter = getPressure(bmIndex, r) + relaxation / getDiagonalElement(bmIndex, r) * (getSourceTerm(bmIndex, r) - getSourceTermRHS(bmIndex, r));
						setPressure(bmIndex, r, std::max(pressureNextIter, static_cast<Real>(0.0)));
						avgDensityDev -= (getSourceTerm(bmIndex, r) - getMinus_rho_div_v_s(bmIndex, r) - getSourceTermRHS(bmIndex, r)) * dt;
					} else {
						setPressure(bmIndex, r, static_cast<Real>(0.0));
					}
				}
			}
		}
	}
	avgDensityDev /= m_contactsAllBodies;
	avgDensityDev /= getRestDensity();
	avgDensityDeviation = std::abs(avgDensityDev);
}

void StrongCouplingBoundarySolver::applyForce() {
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nFluids = sim->numberOfFluidModels();
	for (unsigned int boundaryPointSetIndex = nFluids; boundaryPointSetIndex < sim->numberOfPointSets(); boundaryPointSetIndex++) {
		BoundaryModel_Akinci2012* bm = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModelFromPointSet(boundaryPointSetIndex));
		int bmIndex = boundaryPointSetIndex - nFluids;
		DynamicRigidBody* rb = dynamic_cast<DynamicRigidBody*>(bm->getRigidBodyObject());
		if (rb != nullptr && rb->isDynamic()) {
            #pragma omp parallel default(shared)
			{
                #pragma omp for schedule(static)
				for (int r = 0; r < bm->numberOfParticles(); r++) {
					const Vector3r f = -getArtificialVolume(bmIndex, r) * getPressureGrad(bmIndex, r);				
					bm->addForce(bm->getPosition(r), f);
					if (f.norm() > 1e-10) {
						const Vector3r fric = copmuteParticleFriction(boundaryPointSetIndex, r, f);
						bm->addForce(bm->getPosition(r), fric);
					}
					setPressure(bmIndex, r, static_cast<Real>(0.0));
					setPressureGrad(bmIndex, r, Vector3r::Zero());
				}
			}
		}
	}
}

Vector3r StrongCouplingBoundarySolver::copmuteParticleFriction(const unsigned int& boundaryPointSetIndex, const unsigned int& index, const Vector3r& force) {
	Simulation* sim = Simulation::getCurrent();
	DynamicBoundarySimulator* boundarySimulator = sim->getDynamicBoundarySimulator();
	TimeManager* tm = TimeManager::getCurrent();
	const Real dt = tm->getTimeStepSize();
	const unsigned int nFluids = sim->numberOfFluidModels();
	BoundaryModel_Akinci2012* bm = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModelFromPointSet(boundaryPointSetIndex));
	DynamicRigidBody* rb = dynamic_cast<DynamicRigidBody*>(bm->getRigidBodyObject());
	if (rb == nullptr) {
		return Vector3r::Zero();
	}
	int bmIndex = boundaryPointSetIndex - nFluids;

	// compute normal
	Vector3r numerator = bm->getPosition(index) * W_zero();
	Real denom = W_zero();
	for (unsigned int n = 0; n < sim->numberOfNeighbors(boundaryPointSetIndex, boundaryPointSetIndex, index); n++) {
		const unsigned int k = sim->getNeighbor(boundaryPointSetIndex, boundaryPointSetIndex, index, n);
		numerator += bm->getPosition(k) * W(bm->getPosition(index) - bm->getPosition(k));
		denom += W(bm->getPosition(index) - bm->getPosition(k));
	}
	Vector3r normal = bm->getPosition(index) - numerator / denom;
	normal.normalize();

	// compute averaged relative predicted velocity
	Vector3r numeratorV = Vector3r::Zero();
	Real denomV = 0;
	for (unsigned int pid = nFluids; pid < sim->numberOfPointSets(); pid++) {
		BoundaryModel_Akinci2012* bm_k = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModelFromPointSet(pid));
		if (pid != boundaryPointSetIndex) {
			for (unsigned int n = 0; n < sim->numberOfNeighbors(boundaryPointSetIndex, pid, index); n++) {
				const unsigned int k = sim->getNeighbor(boundaryPointSetIndex, pid, index, n);
				numeratorV += getPredictedVelocity(pid - nFluids, k) * W(bm->getPosition(index) - bm_k->getPosition(k));
				denomV += W(bm->getPosition(index) - bm_k->getPosition(k));
			}
		}
	}
	Vector3r neighborV = numeratorV / denomV;
	Vector3r relativeV = getPredictedVelocity(bmIndex, index) - neighborV;
	Vector3r normalVPred = relativeV.dot(normal) * normal;
	Vector3r tangentialVPred = relativeV - normalVPred;

	Vector3r tangent = tangentialVPred;
	tangent.normalize();
	Vector3r normalF = -force.dot(normal) * normal; 
	Vector3r tangentF = -rb->getFrictionCoeff() * normalF.norm() * tangent;
	Vector3r pos = bm->getPosition(index);
	Matrix3r productMat;
	productMat << 0, -pos.z(), pos.y(),
		          pos.z(), 0, -pos.x(),
		          -pos.y(), pos.x(), 0;
	Matrix3r K = rb->getInvMass() * Matrix3r::Identity() - productMat * rb->getInertiaTensorInverseW() * productMat;
	Vector3r max = -tangentialVPred / (dt * getBodyContacts(bmIndex) * tangent.transpose() * K * tangent);
	if (max.norm() < tangentF.norm()) {
		tangentF = max;
	}
	return tangentF;
}

Real StrongCouplingBoundarySolver::W(const Vector3r& r) {
	return m_kernelFct(r);
}

Vector3r StrongCouplingBoundarySolver::gradW(const Vector3r& r) {
	return m_gradKernelFct(r);
}

Real StrongCouplingBoundarySolver::W_zero() {
	return m_W_zero;
}
