#include "SPlisHSPlasH/BoundaryModel_Akinci2012.h"
#include "StrongCouplingBoundarySolver.h"
#include "TimeManager.h"
#include "DynamicRigidBody.h"

using namespace SPH;

StrongCouplingBoundarySolver* StrongCouplingBoundarySolver::current = nullptr;

void StrongCouplingBoundarySolver::computeSourceTerm() {
	Simulation* sim = Simulation::getCurrent();
	DynamicBoundarySimulator* boundarySimulator = sim->getDynamicBoundarySimulator();
	TimeManager* tm = TimeManager::getCurrent();
	const Real dt = tm->getTimeStepSize();
	const unsigned int nFluids = sim->numberOfFluidModels();

	// compute source term s for all particles
	for (unsigned int boundaryPointSetIndex = nFluids; boundaryPointSetIndex < sim->numberOfPointSets(); boundaryPointSetIndex++) {
		BoundaryModel_Akinci2012* bm = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModelFromPointSet(boundaryPointSetIndex));
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
							for (unsigned int j = 0; j < sim->numberOfNeighbors(boundaryPointSetIndex, pid, r); j++) {
								const unsigned int k = sim->getNeighbor(boundaryPointSetIndex, pid, r, j);
								// sum up divergence of v_s
								s += bm_neighbor->getArtificialVolume(k) * bm_neighbor->getDensity(k) * (bm_neighbor->getV_s(k) - bm->getV_s(r)).dot(sim->gradW(bm->getPosition(r) - bm_neighbor->getPosition(k)));
							}
						}
					}
					s += (bm->getRestDensity() - bm->getDensity(r)) / dt;
					bm->setSourceTerm(r, s);
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
		BoundaryModel_Akinci2012* bm = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModelFromPointSet(pointSetIndex_r));
		DynamicRigidBody* rb = static_cast<DynamicRigidBody*>(bm->getRigidBodyObject());
		if (rb->isDynamic()) {
            #pragma omp parallel default(shared)
			{
                #pragma omp for schedule(static)      
				for (int r = 0; r < bm->numberOfParticles(); r++) {
					Matrix3r productMat_r;
					Vector3r pos_r = bm->getPosition(r);
					productMat_r << 0, -pos_r.z(), pos_r.y(),
						pos_r.z(), 0, -pos_r.x(),
						-pos_r.y(), pos_r.x(), 0;
					Real diagonal = 0;
					for (unsigned int pointSetIndex_i = nFluids; pointSetIndex_i < sim->numberOfPointSets(); pointSetIndex_i++) {
						BoundaryModel_Akinci2012* bm_i = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModelFromPointSet(pointSetIndex_i));
						DynamicRigidBody* rb_i = static_cast<DynamicRigidBody*>(bm_i->getRigidBodyObject());

						// neighbors of r
						for (unsigned int m = 0; m < sim->numberOfNeighbors(pointSetIndex_r, pointSetIndex_i, r); m++) {

							const unsigned int i = sim->getNeighbor(pointSetIndex_r, pointSetIndex_i, r, m);
							Vector3r coeffSum_ik;

							Vector3r gradWri = sim->gradW(bm->getPosition(r) - bm_i->getPosition(i));
							Matrix3r productMat_i;
							Vector3r pos_i = bm_i->getPosition(i);
							productMat_i << 0, -pos_i.z(), pos_i.y(),
								pos_i.z(), 0, -pos_i.x(),
								-pos_i.y(), pos_i.x(), 0;

							// neighbors of i in the same rigid body, denoted as k
							for (unsigned int n = 0; n < sim->numberOfNeighbors(pointSetIndex_i, pointSetIndex_i, i); n++) {
								const unsigned int k = sim->getNeighbor(pointSetIndex_i, pointSetIndex_i, i, n);
								Vector3r coeffSum_kj;

								Vector3r pos_k = bm_i->getPosition(k);
								Matrix3r productMat_k;
								productMat_k << 0, -pos_k.z(), pos_k.y(),
									pos_k.z(), 0, -pos_k.x(),
									-pos_k.y(), pos_k.x(), 0;
								Matrix3r K_ik = rb_i->getInvMass() * Matrix3r::Identity() - productMat_i * rb_i->getInertiaTensorInverseW() * productMat_k;

								// neighbors of k in all bodies, denoted as j
								for (unsigned int pointSetIndex_j = nFluids; pointSetIndex_j < sim->numberOfPointSets(); pointSetIndex_j++) {
									BoundaryModel_Akinci2012* bm_j = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModelFromPointSet(pointSetIndex_j));
									for (unsigned int o = 0; o < sim->numberOfNeighbors(pointSetIndex_i, pointSetIndex_j, k); o++) {
										const unsigned int j = sim->getNeighbor(pointSetIndex_i, pointSetIndex_j, k, o);
										Vector3r gradWkj = sim->gradW(bm_i->getPosition(k) - bm_j->getPosition(j));
										if (pointSetIndex_i == pointSetIndex_r && bm_i->getParticleID(k) == bm->getParticleID(r)) {
											coeffSum_kj += bm_j->getArtificialVolume(j) * bm_j->getDensity(j) / (bm_i->getDensity(k) * bm_i->getDensity(k)) * gradWkj;
										}
										if (pointSetIndex_j == pointSetIndex_r && bm_j->getParticleID(j) == bm->getParticleID(r)) {
											coeffSum_kj += bm_j->getArtificialVolume(j) * bm_j->getDensity(j) / (bm_j->getDensity(j) * bm_j->getDensity(j)) * gradWkj;
										}
									}
								}
								coeffSum_kj = bm_i->getArtificialVolume(k) * K_ik * bm_i->getDensity(k) * coeffSum_kj;
								coeffSum_ik += coeffSum_kj;
							}

							coeffSum_ik *= dt;

							Vector3r coeffSum_rk;
							// neighbors of r in the same rigid body, denoted as k
							for (unsigned int n = 0; n < sim->numberOfNeighbors(pointSetIndex_r, pointSetIndex_r, r); n++) {
								const unsigned int k = sim->getNeighbor(pointSetIndex_r, pointSetIndex_r, r, n);
								Vector3r coeffSum_kj;

								Vector3r pos_k = bm->getPosition(k);
								Matrix3r productMat_k;
								productMat_k << 0, -pos_k.z(), pos_k.y(),
									pos_k.z(), 0, -pos_k.x(),
									-pos_k.y(), pos_k.x(), 0;
								Matrix3r K_rk = rb->getInvMass() * Matrix3r::Identity() - productMat_r * rb->getInertiaTensorInverseW() * productMat_k;

								// neighbors of k in all bodies, denoted as j
								for (unsigned int pointSetIndex_j = nFluids; pointSetIndex_j < sim->numberOfPointSets(); pointSetIndex_j++) {
									BoundaryModel_Akinci2012* bm_j = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModelFromPointSet(pointSetIndex_j));
									for (unsigned int o = 0; o < sim->numberOfNeighbors(pointSetIndex_r, pointSetIndex_j, k); o++) {
										const unsigned int j = sim->getNeighbor(pointSetIndex_r, pointSetIndex_j, k, o);
										Vector3r gradWkj = sim->gradW(bm->getPosition(k) - bm_j->getPosition(j));
										if (bm->getParticleID(k) == bm->getParticleID(r)) {
											coeffSum_kj += bm_j->getArtificialVolume(j) * bm_j->getDensity(j) / (bm->getDensity(k) * bm->getDensity(k)) * gradWkj;
										}
										if (pointSetIndex_j == pointSetIndex_r && bm_j->getParticleID(j) == bm->getParticleID(r)) {
											coeffSum_kj += bm_j->getArtificialVolume(j) * bm_j->getDensity(j) / (bm_j->getDensity(j) * bm_j->getDensity(j)) * gradWkj;
										}
									}
								}
								coeffSum_kj = bm->getArtificialVolume(k) * K_rk * bm->getDensity(k) * coeffSum_kj;
								coeffSum_rk += coeffSum_kj;
							}

							coeffSum_rk *= dt;

							diagonal += bm_i->getArtificialVolume(i) * bm_i->getDensity(i) * (coeffSum_ik - coeffSum_rk).dot(gradWri);
						}
					}
					bm->setDiagonalElement(r, diagonal);
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
		DynamicRigidBody* rb = static_cast<DynamicRigidBody*>(bm->getRigidBodyObject());
		if (rb->isDynamic()) {
            #pragma omp parallel default(shared)
			{
                #pragma omp for schedule(static)  
				for (int r = 0; r < bm->numberOfParticles(); r++) {
					// compute density
					Real particleDensity = 0;
					// iterate over all rigid bodies
					for (unsigned int pid = nFluids; pid < sim->numberOfPointSets(); pid++) {
						BoundaryModel_Akinci2012* bm_neighbor = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModelFromPointSet(pid));
						for (unsigned int j = 0; j < sim->numberOfNeighbors(boundaryPointSetIndex, pid, r); j++) {
							const unsigned int k = sim->getNeighbor(boundaryPointSetIndex, pid, r, j);
							particleDensity += bm->getRestDensity() * bm->getVolume(r) * sim->W(bm->getPosition(r) - bm_neighbor->getPosition(k));
						}
					}
					// 0.8 for compensation
					bm->setDensity(r, particleDensity / 0.8);
					bm->setArtificialVolume(r, bm->getRestDensity() * bm->getVolume(r) / bm->getDensity(r));
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
		DynamicRigidBody* rb = static_cast<DynamicRigidBody*>(bm->getRigidBodyObject());
		if (rb->isDynamic()) {
            #pragma omp parallel default(shared)
			{
                #pragma omp for schedule(static)  
				for (int r = 0; r < bm->numberOfParticles(); r++) {
					Vector3r pressureGrad_r = Vector3r().setZero();
					const Real density_r = bm->getDensity(r);
					const Real pressure_r = bm->getPressure(r);
					for (unsigned int pid = nFluids; pid < sim->numberOfPointSets(); pid++) {
						BoundaryModel_Akinci2012* bm_neighbor = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModelFromPointSet(pid));
						for (unsigned int j = 0; j < sim->numberOfNeighbors(boundaryPointSetIndex, pid, r); j++) {
							const unsigned int k = sim->getNeighbor(boundaryPointSetIndex, pid, r, j);
							const Real density_k = bm_neighbor->getDensity(k);
							const Real volume_k = bm_neighbor->getRestDensity() / density_k * bm_neighbor->getVolume(k);
							const Real pressure_k = bm_neighbor->getPressure(k);
							pressureGrad_r += volume_k * density_k * (pressure_r / (density_r * density_r) + pressure_k / (density_k * density_k)) * sim->gradW(bm->getPosition(r) - bm_neighbor->getPosition(k));
						}
					}
					pressureGrad_r *= density_r;
					bm->setPressureGrad(r, pressureGrad_r);
				}
			}
		}
	}

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
					bm->setV_s(r, rb->getVelocity() + dt * rb->getInvMass() * F_R +
						(rb->getAngularVelocity() + rb->getInertiaTensorInverseW() * (dt * rb->getTorque() + dt * (rb->getInertiaTensorW() * rb->getAngularVelocity()).cross(rb->getAngularVelocity()))).cross(bm->getPosition(r)));
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
		DynamicRigidBody* rb = static_cast<DynamicRigidBody*>(bm->getRigidBodyObject());
		Vector3r v_rr_body = Vector3r().setZero();
		Vector3r omega_rr_body = Vector3r().setZero();
		if (rb->isDynamic()) {
            #pragma omp parallel default(shared)
			{
                #pragma omp for schedule(static)  
				for (int r = 0; r < bm->numberOfParticles(); r++) {
					v_rr_body += -dt * rb->getInvMass() * bm->getArtificialVolume(r) * bm->getPressureGrad(r);
					omega_rr_body += -dt * rb->getInertiaTensorInverseW() * bm->getArtificialVolume(r) * bm->getPosition(r).cross(bm->getPressureGrad(r));
				}
			}
		}
		bm->setV_rr_body(v_rr_body);
		bm->setOmega_rr_body(omega_rr_body);
	}
	// compute v_rr for all particles (predicted velocity)
	for (unsigned int boundaryPointSetIndex = nFluids; boundaryPointSetIndex < sim->numberOfPointSets(); boundaryPointSetIndex++) {
		BoundaryModel_Akinci2012* bm = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModelFromPointSet(boundaryPointSetIndex));
		DynamicRigidBody* rb = static_cast<DynamicRigidBody*>(bm->getRigidBodyObject());
		if (rb->isDynamic()) {
            #pragma omp parallel default(shared)
			{
                #pragma omp for schedule(static)  
				for (int r = 0; r < bm->numberOfParticles(); r++) {
					bm->setV_rr(r, bm->getV_rr_body() + bm->getOmega_rr_body().cross(bm->getPosition(r)));
				}
			}
		}
	}
	// compute the -rho * (div v_rr), which is the RHS to the source term
	for (unsigned int boundaryPointSetIndex = nFluids; boundaryPointSetIndex < sim->numberOfPointSets(); boundaryPointSetIndex++) {
		BoundaryModel_Akinci2012* bm = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModelFromPointSet(boundaryPointSetIndex));
        #pragma omp parallel default(shared)
		{
            #pragma omp for schedule(static)  
			for (int r = 0; r < bm->numberOfParticles(); r++) {
				Real minus_rho_div_v_rr = 0;
				const Vector3r v_rr_r = bm->getV_rr(r);
				for (unsigned int pid = nFluids; pid < sim->numberOfPointSets(); pid++) {
					BoundaryModel_Akinci2012* bm_neighbor = static_cast<BoundaryModel_Akinci2012*>(sim->getBoundaryModelFromPointSet(pid));
					for (unsigned int j = 0; j < sim->numberOfNeighbors(boundaryPointSetIndex, pid, r); j++) {
						const unsigned int k = sim->getNeighbor(boundaryPointSetIndex, pid, r, j);
						const Real density_k = bm_neighbor->getDensity(k);
						const Real volume_k = bm_neighbor->getArtificialVolume(k);
						const Vector3r v_rr_k = bm_neighbor->getV_rr(k);
						minus_rho_div_v_rr += volume_k * density_k * (v_rr_k - v_rr_r).dot(sim->gradW(bm->getPosition(r) - bm_neighbor->getPosition(k)));
					}
				}
				minus_rho_div_v_rr = -minus_rho_div_v_rr;
				bm->setSourceTermRHS(r, minus_rho_div_v_rr);
			}
		}
	}
}