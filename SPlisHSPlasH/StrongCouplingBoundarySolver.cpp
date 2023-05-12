#include "SPlisHSPlasH/BoundaryModel_Akinci2012.h"
#include "StrongCouplingBoundarySolver.h"
#include "TimeManager.h"
#include "DynamicRigidBody.h"

using namespace SPH;

StrongCouplingBoundarySolver* StrongCouplingBoundarySolver::current = nullptr;

StrongCouplingBoundarySolver::StrongCouplingBoundarySolver() :
	m_density(),
	m_pressureGrad(),
	m_v_s(),
	m_s(),
	m_pressure(),
	m_v_rr(),
	m_minus_rho_div_v_rr(),
	m_diagonalElement(),
	m_artificialVolume(),
	m_lastPressure(),
	m_predictVelocity()
{
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();

	m_density.resize(nBoundaries);
	m_v_s.resize(nBoundaries);
	m_s.resize(nBoundaries);
	m_pressure.resize(nBoundaries);
	m_v_rr.resize(nBoundaries);
	m_minus_rho_div_v_rr.resize(nBoundaries);
	m_diagonalElement.resize(nBoundaries);
	m_pressureGrad.resize(nBoundaries);
	m_artificialVolume.resize(nBoundaries);
	m_lastPressure.resize(nBoundaries);
	m_predictVelocity.resize(nBoundaries);
	m_restDensity.resize(nBoundaries, 1);
	m_v_rr_body.resize(nBoundaries, Vector3r::Zero());
	m_omega_rr_body.resize(nBoundaries, Vector3r::Zero());

	for (unsigned int i = 0; i < nBoundaries; i++) {
		BoundaryModel_Akinci2012* bm = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModel(i));
		
		m_density[i].resize(bm->numberOfParticles(), m_restDensity[i]);
		m_v_s[i].resize(bm->numberOfParticles(), Vector3r::Zero());
		m_s[i].resize(bm->numberOfParticles(), 0);
		m_pressure[i].resize(bm->numberOfParticles(), 0);
		m_v_rr[i].resize(bm->numberOfParticles(), Vector3r::Zero());
		m_minus_rho_div_v_rr[i].resize(bm->numberOfParticles(), 0);
		m_diagonalElement[i].resize(bm->numberOfParticles(), 0);
		m_pressureGrad[i].resize(bm->numberOfParticles(), Vector3r::Zero());
		m_artificialVolume[i].resize(bm->numberOfParticles(), 0);
		m_lastPressure[i].resize(bm->numberOfParticles(), 0);
		m_predictVelocity[i].resize(bm->numberOfParticles(), Vector3r::Zero());
		
		bm->addField({ "density", FieldType::Scalar, [this, i](const unsigned int j) -> Real* { return &getDensity(i, j); } });
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
		m_minus_rho_div_v_rr[i].clear();
		m_diagonalElement[i].clear();
		m_pressureGrad[i].clear();
		m_artificialVolume[i].clear();
		m_lastPressure[i].clear();
		m_predictVelocity[i].clear();
	}
	m_density.clear();
	m_v_s.clear();
	m_s.clear();
	m_pressure.clear();
	m_v_rr.clear();
	m_minus_rho_div_v_rr.clear();
	m_diagonalElement.clear();
	m_pressureGrad.clear();
	m_artificialVolume.clear();
	m_lastPressure.clear();
	m_predictVelocity.clear();
}

void StrongCouplingBoundarySolver::computeV_s() {
	// Use dynamic boundary simulator and Akinci2012
	Simulation* sim = Simulation::getCurrent();
	DynamicBoundarySimulator* boundarySimulator = sim->getDynamicBoundarySimulator();
	TimeManager* tm = TimeManager::getCurrent();
	const Real dt = tm->getTimeStepSize();
	const unsigned int nFluids = sim->numberOfFluidModels();

	// Compute v_s for all particles
	for (unsigned int boundaryPointSetIndex = nFluids; boundaryPointSetIndex < sim->numberOfPointSets(); boundaryPointSetIndex++) {
		BoundaryModel_Akinci2012* bm = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModelFromPointSet(boundaryPointSetIndex));
		DynamicRigidBody* rb = static_cast<DynamicRigidBody*>(bm->getRigidBodyObject());
		if (rb->isDynamic()) {
            #pragma omp parallel default(shared)
			{
                #pragma omp for schedule(static)  
				for (int r = 0; r < bm->numberOfParticles(); r++) {
					// compute v_s
					Vector3r gravForce = rb->getMass() * Vector3r(sim->getVecValue<Real>(Simulation::GRAVITATION));
					// the rb->getForce() contains only fluid-rigid forces
					Vector3r F_R = rb->getForce() + gravForce;
					setV_s(boundaryPointSetIndex - nFluids, r, rb->getVelocity() + dt * rb->getInvMass() * F_R +
						(rb->getAngularVelocity() + rb->getInertiaTensorInverseW() * (dt * rb->getTorque() + dt * (rb->getInertiaTensorW() * rb->getAngularVelocity()).cross(rb->getAngularVelocity()))).cross(bm->getPosition(r)));
				}
			}
		}
	}
}

void StrongCouplingBoundarySolver::computeSourceTerm() {
	Simulation* sim = Simulation::getCurrent();
	DynamicBoundarySimulator* boundarySimulator = sim->getDynamicBoundarySimulator();
	TimeManager* tm = TimeManager::getCurrent();
	const Real dt = tm->getTimeStepSize();
	const unsigned int nFluids = sim->numberOfFluidModels();

	// compute source term s for all particles
	for (unsigned int boundaryPointSetIndex = nFluids; boundaryPointSetIndex < sim->numberOfPointSets(); boundaryPointSetIndex++) {
		BoundaryModel_Akinci2012* bm = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModelFromPointSet(boundaryPointSetIndex));
		int bmIndex = boundaryPointSetIndex - nFluids;
		if (bm->getRigidBodyObject()->isDynamic()) {
            #pragma omp parallel default(shared)
			{
                #pragma omp for schedule(static)  
				for (int r = 0; r < bm->numberOfParticles(); r++) {
					Real s = 0;
					// iterate over all rigid bodies except bm, since the divergence of particles in the same rigid body should be 0.
					for (unsigned int pid = nFluids; pid < sim->numberOfPointSets(); pid++) {
						if (boundaryPointSetIndex != pid) {
							BoundaryModel_Akinci2012* bm_neighbor = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModelFromPointSet(pid));
							int neighborIndex = pid - nFluids;
							for (unsigned int j = 0; j < sim->numberOfNeighbors(boundaryPointSetIndex, pid, r); j++) {
								const unsigned int k = sim->getNeighbor(boundaryPointSetIndex, pid, r, j);
								// sum up divergence of v_s
								s += getArtificialVolume(neighborIndex, k) * getDensity(neighborIndex, k) * (getV_s(neighborIndex, k) - getV_s(bmIndex, r)).dot(sim->gradW(bm->getPosition(r) - bm_neighbor->getPosition(k)));
							}
						}
					}
					s += (getRestDensity(bmIndex) - getDensity(bmIndex, r)) / dt;
					setSourceTerm(bmIndex, r, s);
				}
			}
		}
	}
}

void StrongCouplingBoundarySolver::computeDiagonalElement() {
	Simulation* sim = Simulation::getCurrent();
	DynamicBoundarySimulator* boundarySimulator = sim->getDynamicBoundarySimulator();
	TimeManager* tm = TimeManager::getCurrent();
	const Real dt = tm->getTimeStepSize();
	const unsigned int nFluids = sim->numberOfFluidModels();
	for (unsigned int pointSetIndex_r = nFluids; pointSetIndex_r < sim->numberOfPointSets(); pointSetIndex_r++) {
		BoundaryModel_Akinci2012* bm_r = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModelFromPointSet(pointSetIndex_r));
		int indexR = pointSetIndex_r - nFluids;
		DynamicRigidBody* rb = static_cast<DynamicRigidBody*>(bm_r->getRigidBodyObject());
		if (rb->isDynamic()) {
            #pragma omp parallel default(shared)
			{
                #pragma omp for schedule(static)      
				for (int r = 0; r < bm_r->numberOfParticles(); r++) {
					Matrix3r productMat_r = Matrix3r().setZero();
					Vector3r pos_r = bm_r->getPosition(r);
					productMat_r << 0, -pos_r.z(), pos_r.y(),
						            pos_r.z(), 0, -pos_r.x(),
						            -pos_r.y(), pos_r.x(), 0;
					Real diagonal = 0;
					for (unsigned int pointSetIndex_i = nFluids; pointSetIndex_i < sim->numberOfPointSets(); pointSetIndex_i++) {
						BoundaryModel_Akinci2012* bm_i = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModelFromPointSet(pointSetIndex_i));
						int indexI = pointSetIndex_i - nFluids;
						DynamicRigidBody* rb_i = static_cast<DynamicRigidBody*>(bm_i->getRigidBodyObject());

						// neighbors of r in other rigid bodies, denoted as i
						if (pointSetIndex_i != pointSetIndex_r) {
							for (unsigned int m = 0; m < sim->numberOfNeighbors(pointSetIndex_r, pointSetIndex_i, r); m++) {
								const unsigned int i = sim->getNeighbor(pointSetIndex_r, pointSetIndex_i, r, m);
								Vector3r coeffSum_ik = Vector3r().setZero();

								Vector3r gradWri = sim->gradW(bm_r->getPosition(r) - bm_i->getPosition(i));
								Matrix3r productMat_i = Matrix3r().setZero();
								Vector3r pos_i = bm_i->getPosition(i);
								productMat_i << 0, -pos_i.z(), pos_i.y(),
									            pos_i.z(), 0, -pos_i.x(),
									            -pos_i.y(), pos_i.x(), 0;

								// neighbors of i in the same rigid body, denoted as k (k cannot be r)
								for (unsigned int n = 0; n < sim->numberOfNeighbors(pointSetIndex_i, pointSetIndex_i, i); n++) {
									const unsigned int k = sim->getNeighbor(pointSetIndex_i, pointSetIndex_i, i, n);
									Vector3r coeffSum_kj = Vector3r().setZero();

									Vector3r pos_k = bm_i->getPosition(k);
									Matrix3r productMat_k = Matrix3r().setZero();
									productMat_k << 0, -pos_k.z(), pos_k.y(),
										            pos_k.z(), 0, -pos_k.x(),
										            -pos_k.y(), pos_k.x(), 0;
									Real invMass = 0;
									if (rb_i->isDynamic()) {
										invMass = rb_i->getInvMass();
									} else {
										invMass = 1;
									}
									Matrix3r K_ik = invMass * Matrix3r::Identity() - productMat_i * rb_i->getInertiaTensorInverseW() * productMat_k;

									// neighbors of k in all bodies, denoted as j (loop can be optimized?)
									//for (unsigned int pointSetIndex_j = nFluids; pointSetIndex_j < sim->numberOfPointSets(); pointSetIndex_j++) {
									//	BoundaryModel_Akinci2012* bm_j = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModelFromPointSet(pointSetIndex_j));
									//	for (unsigned int o = 0; o < sim->numberOfNeighbors(pointSetIndex_i, pointSetIndex_j, k); o++) {
									//		const unsigned int j = sim->getNeighbor(pointSetIndex_i, pointSetIndex_j, k, o);
									//		int indexJ = pointSetIndex_j - nFluids;
									//		Vector3r gradWkj = sim->gradW(bm_i->getPosition(k) - bm_j->getPosition(j));
									//		//if (pointSetIndex_i == pointSetIndex_r && bm_i->getParticleID(k) == bm->getParticleID(r)) {
									//		//	coeffSum_kj += bm_j->getArtificialVolume(j) * bm_j->getDensity(j) / (bm_i->getDensity(k) * bm_i->getDensity(k)) * gradWkj;
									//		//}
									//		if (pointSetIndex_j == pointSetIndex_r && j == r) {
									//			coeffSum_kj += getArtificialVolume(indexJ, j) * getDensity(indexJ, j) / (getDensity(indexJ, j) * getDensity(indexJ, j)) * gradWkj;
									//		}
									//	}
									//}

									//if (coeffSum_kj != Vector3r().setZero()) {
									//	coeffSum_kj = getArtificialVolume(indexI, k) * K_ik * getDensity(indexI, k) * coeffSum_kj;
									//	coeffSum_ik += coeffSum_kj;
									//}
									Vector3r gradWkr = sim->gradW(bm_i->getPosition(k) - bm_r->getPosition(r));
									coeffSum_ik += getArtificialVolume(indexI, k) * K_ik * getDensity(indexI, k) * getArtificialVolume(indexR, r) * getDensity(indexR, r) / (getDensity(indexR, r) * getDensity(indexR, r)) * gradWkr;
								}
								coeffSum_ik *= dt;

								Vector3r coeffSum_rk = Vector3r().setZero();
								// neighbors of r in the same rigid body, denoted as k, k != r
								for (unsigned int n = 0; n < sim->numberOfNeighbors(pointSetIndex_r, pointSetIndex_r, r); n++) {
									const unsigned int k = sim->getNeighbor(pointSetIndex_r, pointSetIndex_r, r, n);

									Vector3r pos_k = bm_r->getPosition(k);
									Matrix3r productMat_k = Matrix3r().setZero();
									productMat_k << 0, -pos_k.z(), pos_k.y(),
										            pos_k.z(), 0, -pos_k.x(),
										            -pos_k.y(), pos_k.x(), 0;
									Matrix3r K_rk = rb->getInvMass() * Matrix3r::Identity() - productMat_r * rb->getInertiaTensorInverseW() * productMat_k;

									// neighbors of k in all bodies, denoted as j
									//for (unsigned int pointSetIndex_j = nFluids; pointSetIndex_j < sim->numberOfPointSets(); pointSetIndex_j++) {
									//	BoundaryModel_Akinci2012* bm_j = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModelFromPointSet(pointSetIndex_j));
									//	int indexJ = pointSetIndex_j - nFluids;
									//	for (unsigned int o = 0; o < sim->numberOfNeighbors(pointSetIndex_r, pointSetIndex_j, k); o++) {
									//		const unsigned int j = sim->getNeighbor(pointSetIndex_r, pointSetIndex_j, k, o);
									//		Vector3r gradWkj = sim->gradW(bm_r->getPosition(k) - bm_j->getPosition(j));
									//		if (k == r) {
									//			coeffSum_kj += getArtificialVolume(indexJ, j) * getDensity(indexJ, j) / (getDensity(indexR, k) * getDensity(indexR, k)) * gradWkj;
									//		}

									//		if (pointSetIndex_j == pointSetIndex_r && j == r) {
									//			coeffSum_kj += getArtificialVolume(indexJ, j) * getDensity(indexJ, j) / (getDensity(indexJ, j) * getDensity(indexJ, j)) * gradWkj;
									//		}
									//	}
									//}
									Vector3r gradWkr = sim->gradW(bm_r->getPosition(k) - bm_r->getPosition(r));
										
									coeffSum_rk += getArtificialVolume(indexR, k) * K_rk * getDensity(indexR, k) * getArtificialVolume(indexR, r) * getDensity(indexR, r) / (getDensity(indexR, r) * getDensity(indexR, r)) * gradWkr;

									//if (coeffSum_kj != Vector3r().setZero()) {
									//	coeffSum_kj = getArtificialVolume(indexR, k) * K_rk * getDensity(indexR, k) * coeffSum_kj;
									//	coeffSum_rk += coeffSum_kj;
									//}
								}

								Vector3r coeffSum_rj = Vector3r().setZero();
								// k = r, for all neighbors of r in all bodies, denoted as j
								for (unsigned int pointSetIndex_j = nFluids; pointSetIndex_j < sim->numberOfPointSets(); pointSetIndex_j++) {
									BoundaryModel_Akinci2012* bm_j = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModelFromPointSet(pointSetIndex_j));
									int indexJ = pointSetIndex_j - nFluids;
									for (unsigned int o = 0; o < sim->numberOfNeighbors(pointSetIndex_r, pointSetIndex_j, r); o++) {
										const unsigned int j = sim->getNeighbor(pointSetIndex_r, pointSetIndex_j, r, o);
										Vector3r gradWrj = sim->gradW(bm_r->getPosition(r) - bm_j->getPosition(j));
										coeffSum_rj += getArtificialVolume(indexJ, j) * getDensity(indexJ, j) / (getDensity(indexR, r) * getDensity(indexR, r)) * gradWrj;
									}
								}


								Matrix3r K_rr = rb->getInvMass() * Matrix3r::Identity() - productMat_r * rb->getInertiaTensorInverseW() * productMat_r;

								coeffSum_rj = getArtificialVolume(indexR, r) * K_rr * getDensity(indexR, r) * coeffSum_rj;
								coeffSum_rk += coeffSum_rj;
								coeffSum_rk *= dt;
								diagonal += getArtificialVolume(indexI, i) * getDensity(indexI, i) * (coeffSum_ik - coeffSum_rk).dot(gradWri);

							}
						}					
					}
					setDiagonalElement(indexR, r, diagonal);
					//std::cout << getDiagonalElement(indexR, r) << std::endl;
				}
			}
		}
	}
}

void StrongCouplingBoundarySolver::computeDensityAndVolume() {
	// Use dynamic boundary simulator and Akinci2012
	Simulation* sim = Simulation::getCurrent();
	DynamicBoundarySimulator* boundarySimulator = sim->getDynamicBoundarySimulator();
	TimeManager* tm = TimeManager::getCurrent();
	const Real dt = tm->getTimeStepSize();
	const unsigned int nFluids = sim->numberOfFluidModels();

	// Compute density and artificial volume for all particles
	for (unsigned int boundaryPointSetIndex = nFluids; boundaryPointSetIndex < sim->numberOfPointSets(); boundaryPointSetIndex++) {
		BoundaryModel_Akinci2012* bm = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModelFromPointSet(boundaryPointSetIndex));
		int bmIndex = boundaryPointSetIndex - nFluids;
		DynamicRigidBody* rb = static_cast<DynamicRigidBody*>(bm->getRigidBodyObject());
        #pragma omp parallel default(shared)
		{
            #pragma omp for schedule(static)  
			for (int r = 0; r < bm->numberOfParticles(); r++) {
				if (rb->isDynamic()) {
					// compute density for particle r
					Real particleDensity = getRestDensity(bmIndex) * bm->getVolume(r) * sim->W_zero();
					// iterate over all rigid bodies
					for (unsigned int pid = nFluids; pid < sim->numberOfPointSets(); pid++) {
						BoundaryModel_Akinci2012* bm_neighbor = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModelFromPointSet(pid));
						for (unsigned int j = 0; j < sim->numberOfNeighbors(boundaryPointSetIndex, pid, r); j++) {
							const unsigned int k = sim->getNeighbor(boundaryPointSetIndex, pid, r, j);
							particleDensity += getRestDensity(bmIndex) * bm->getVolume(r) * sim->W(bm->getPosition(r) - bm_neighbor->getPosition(k));
						}
					}
					setDensity(bmIndex, r, particleDensity);
					if (getDensity(bmIndex, r) > getRestDensity(bmIndex)) {
						setArtificialVolume(bmIndex, r, getRestDensity(bmIndex) * bm->getVolume(r) / getDensity(bmIndex, r));
					}else{
						setArtificialVolume(bmIndex, r, bm->getVolume(r));
					}
				} else {
					setArtificialVolume(bmIndex, r, bm->getVolume(r));
				}
			}
		}
		
	}
}

void StrongCouplingBoundarySolver::computePressureGrad() {
	Simulation* sim = Simulation::getCurrent();
	DynamicBoundarySimulator* boundarySimulator = sim->getDynamicBoundarySimulator();
	TimeManager* tm = TimeManager::getCurrent();
	const Real dt = tm->getTimeStepSize();
	const unsigned int nFluids = sim->numberOfFluidModels();

	for (unsigned int boundaryPointSetIndex = nFluids; boundaryPointSetIndex < sim->numberOfPointSets(); boundaryPointSetIndex++) {
		BoundaryModel_Akinci2012* bm = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModelFromPointSet(boundaryPointSetIndex));
		int bmIndex = boundaryPointSetIndex - nFluids;
		DynamicRigidBody* rb = static_cast<DynamicRigidBody*>(bm->getRigidBodyObject());
		if (rb->isDynamic()) {
            #pragma omp parallel default(shared)
			{
                #pragma omp for schedule(static)  
				for (int r = 0; r < bm->numberOfParticles(); r++) {
					// pressure gradient for particle r
					Vector3r pressureGrad_r = Vector3r().setZero();
					const Real density_r = getDensity(bmIndex, r);
					const Real pressure_r = getPressure(bmIndex, r);
					for (unsigned int pid = nFluids; pid < sim->numberOfPointSets(); pid++) {
						BoundaryModel_Akinci2012* bm_neighbor = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModelFromPointSet(pid));
						int neighborIndex = pid - nFluids;
						for (unsigned int j = 0; j < sim->numberOfNeighbors(boundaryPointSetIndex, pid, r); j++) {
							const unsigned int k = sim->getNeighbor(boundaryPointSetIndex, pid, r, j);
							const Real density_k = getDensity(neighborIndex, k);
							const Real volume_k = getArtificialVolume(neighborIndex, k);
							const Real pressure_k = getPressure(neighborIndex, k);
							pressureGrad_r += volume_k * density_k * (pressure_r / (density_r * density_r) + pressure_k / (density_k * density_k)) * sim->gradW(bm->getPosition(r) - bm_neighbor->getPosition(k));
						}
					}
					pressureGrad_r *= density_r;
					setPressureGrad(bmIndex, r, pressureGrad_r);
				}
			}
		}
	}
}

void StrongCouplingBoundarySolver::computeSourceTermRHS() {
	Simulation* sim = Simulation::getCurrent();
	DynamicBoundarySimulator* boundarySimulator = sim->getDynamicBoundarySimulator();
	TimeManager* tm = TimeManager::getCurrent();
	const Real dt = tm->getTimeStepSize();
	const unsigned int nFluids = sim->numberOfFluidModels();

	// compute v_rr and omega_rr for rigid body using the pressure gradient
	for (unsigned int boundaryPointSetIndex = nFluids; boundaryPointSetIndex < sim->numberOfPointSets(); boundaryPointSetIndex++) {
		BoundaryModel_Akinci2012* bm = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModelFromPointSet(boundaryPointSetIndex));
		int bmIndex = boundaryPointSetIndex - nFluids;
		DynamicRigidBody* rb = static_cast<DynamicRigidBody*>(bm->getRigidBodyObject());
		Vector3r v_rr_body = Vector3r().setZero();
		Vector3r omega_rr_body = Vector3r().setZero();
		if (rb->isDynamic()) {
            #pragma omp parallel default(shared)
			{
                #pragma omp for schedule(static)  
				for (int r = 0; r < bm->numberOfParticles(); r++) {
					v_rr_body += -dt * rb->getInvMass() * getArtificialVolume(bmIndex, r) * getPressureGrad(bmIndex, r);
					omega_rr_body += -dt * rb->getInertiaTensorInverseW() * getArtificialVolume(bmIndex, r) * bm->getPosition(r).cross(getPressureGrad(bmIndex, r));
				}
			}
		}
		setV_rr_body(bmIndex, v_rr_body);
		setOmega_rr_body(bmIndex, omega_rr_body);
	}
	// compute v_rr and predicted velocity (v_rr+v_s) for all particles
	for (unsigned int boundaryPointSetIndex = nFluids; boundaryPointSetIndex < sim->numberOfPointSets(); boundaryPointSetIndex++) {
		BoundaryModel_Akinci2012* bm = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModelFromPointSet(boundaryPointSetIndex));
		int bmIndex = boundaryPointSetIndex - nFluids;
		DynamicRigidBody* rb = static_cast<DynamicRigidBody*>(bm->getRigidBodyObject());
		if (rb->isDynamic()) {
            #pragma omp parallel default(shared)
			{
                #pragma omp for schedule(static)  
				for (int r = 0; r < bm->numberOfParticles(); r++) {
					setV_rr(bmIndex, r, getV_rr_body(bmIndex) + getOmega_rr_body(bmIndex).cross(bm->getPosition(r)));
					setPredictVelocity(bmIndex, r, getV_rr(bmIndex, r) + getV_s(bmIndex, r));
				}
			}
		}
	}
	// compute the -rho * (div v_rr), which is the RHS to the source term
	for (unsigned int boundaryPointSetIndex = nFluids; boundaryPointSetIndex < sim->numberOfPointSets(); boundaryPointSetIndex++) {
		BoundaryModel_Akinci2012* bm = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModelFromPointSet(boundaryPointSetIndex));
		int bmIndex = boundaryPointSetIndex - nFluids;
        #pragma omp parallel default(shared)
		{
            #pragma omp for schedule(static)  
			for (int r = 0; r < bm->numberOfParticles(); r++) {
				Real minus_rho_div_v_rr = 0;
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
							const Vector3r v_rr_k = getV_rr(neighborIndex, k);
							minus_rho_div_v_rr += volume_k * density_k * (v_rr_k - v_rr_r).dot(sim->gradW(bm->getPosition(r) - bm_neighbor->getPosition(k)));
						}
					}
				}
				minus_rho_div_v_rr = -minus_rho_div_v_rr;
				setSourceTermRHS(bmIndex, r, minus_rho_div_v_rr);
			}
		}
	}
}

void SPH::StrongCouplingBoundarySolver::performNeighborhoodSearchSort() {
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
			d.sort_field(&m_predictVelocity[i][0]);
		}
	}
}