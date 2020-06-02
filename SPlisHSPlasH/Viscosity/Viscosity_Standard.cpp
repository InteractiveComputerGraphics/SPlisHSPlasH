#include "Viscosity_Standard.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "../Simulation.h"
#include "SPlisHSPlasH/BoundaryModel_Akinci2012.h"
#include "SPlisHSPlasH/BoundaryModel_Koschier2017.h"
#include "SPlisHSPlasH/BoundaryModel_Bender2019.h"
#include "SPlisHSPlasH/Utilities/MathFunctions.h"


using namespace SPH;
using namespace GenParam;

int Viscosity_Standard::VISCOSITY_COEFFICIENT_BOUNDARY = -1;


Viscosity_Standard::Viscosity_Standard(FluidModel *model) :
	ViscosityBase(model)
{
	m_boundaryViscosity = 0.0;
}

Viscosity_Standard::~Viscosity_Standard(void)
{
}

void Viscosity_Standard::initParameters()
{
	ViscosityBase::initParameters();

	VISCOSITY_COEFFICIENT_BOUNDARY = createNumericParameter("viscosityBoundary", "Viscosity coefficient (Boundary)", &m_boundaryViscosity);
	setGroup(VISCOSITY_COEFFICIENT_BOUNDARY, "Viscosity");
	setDescription(VISCOSITY_COEFFICIENT_BOUNDARY, "Coefficient for the viscosity force computation at the boundary.");
	RealParameter* rparam = static_cast<RealParameter*>(getParameter(VISCOSITY_COEFFICIENT_BOUNDARY));
	rparam->setMinValue(0.0);
}

#ifdef USE_AVX

void Viscosity_Standard::step()
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int numParticles = m_model->numActiveParticles();
	const Real h = sim->getSupportRadius();
	const Real h2 = h*h;
	const Real sphereVolume = static_cast<Real>(4.0 / 3.0 * M_PI) * h2*h;
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();
	const unsigned int fluidModelIndex = m_model->getPointSetIndex();
	const Real density0 = m_model->getDensity0();
	FluidModel *model = m_model;
	Real d = 10.0;
	if (sim->is2DSimulation())
		d = 8.0;

	const Scalarf8 dvisc(d * m_viscosity * density0);
	const Scalarf8 dbvisc(d * m_boundaryViscosity);
	const Scalarf8 h2_avx(0.01f*h2);
	const Scalarf8 density0_avx(density0);

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			const Vector3r &xi = m_model->getPosition(i);
			const Vector3r &vi = m_model->getVelocity(i);
			Vector3r &ai = m_model->getAcceleration(i);
			const Real density_i = m_model->getDensity(i);

			const Vector3f8 xi_avx(xi);
			const Vector3f8 vi_avx(vi);
			const Scalarf8 mi_avx(m_model->getMass(i));
			const Scalarf8 density_i_avx(density_i);

			Vector3f8 delta_ai;
			delta_ai.setZero();

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors_in_same_phase_avx(
				compute_Vj(model);
				compute_Vj_gradW_samephase();
 				const Vector3f8 vj_avx = convertVec_zero(&sim->getNeighborList(fluidModelIndex, fluidModelIndex, i)[j], &model->getVelocity(0), count);
 
 				// Viscosity
 				const Scalarf8 density_j_avx = convert_one(&sim->getNeighborList(fluidModelIndex, fluidModelIndex, i)[j], &model->getDensity(0), count);
 				const Vector3f8 xixj = xi_avx - xj_avx;
				delta_ai += V_gradW * ((dvisc / density_j_avx) * (vi_avx - vj_avx).dot(xixj) / (xixj.squaredNorm() + h2_avx));
 			);

			//////////////////////////////////////////////////////////////////////////
			// Boundary
			//////////////////////////////////////////////////////////////////////////
 			if (m_boundaryViscosity != 0.0)
 			{
 				if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
 				{
 					forall_boundary_neighbors_avx(
						const Vector3f8 vj_avx = convertVec_zero(&sim->getNeighborList(fluidModelIndex, pid, i)[j], &bm_neighbor->getVelocity(0), count);
						const Scalarf8 Vj_avx = convert_zero(&sim->getNeighborList(fluidModelIndex, pid, i)[j], &bm_neighbor->getVolume(0), count);
						const Vector3f8 xixj = xi_avx - xj_avx;
 						const Vector3f8 delta_a = CubicKernel_AVX::gradW(xixj) * (dbvisc * (density0_avx * Vj_avx / density_i_avx) * (vi_avx - vj_avx).dot(xixj) / (xixj.squaredNorm() + h2_avx));
						delta_ai += delta_a;
 						bm_neighbor->addForce(xj_avx, -delta_a * mi_avx, count);
 					);
 				}
 				else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
 				{
 					forall_density_maps(
						const Vector3r xixj = xi - xj;
						Vector3r normal = -xixj;
						const Real nl = normal.norm();
						if (nl > static_cast<Real>(0.0001))
						{
							normal /= nl;

							Vector3r t1;
							Vector3r t2;
							MathFunctions::getOrthogonalVectors(normal, t1, t2);
							
							const Real dist = (static_cast<Real>(1.0) - nl / h) * h;
							//const Real dist = 0.5*h;
							const Vector3r x1 = xj - t1*dist;
							const Vector3r x2 = xj + t1*dist;
							const Vector3r x3 = xj - t2*dist;
							const Vector3r x4 = xj + t2*dist;

							const Vector3r xix1 = xi - x1;
							const Vector3r xix2 = xi - x2;
							const Vector3r xix3 = xi - x3;
							const Vector3r xix4 = xi - x4;

							const Vector3r gradW1 = sim->gradW(xix1);
							const Vector3r gradW2 = sim->gradW(xix2);
							const Vector3r gradW3 = sim->gradW(xix3);
							const Vector3r gradW4 = sim->gradW(xix4);

							// each sample point represents the quarter of the volume inside of the boundary
							const Real vol = static_cast<Real>(0.25) * rho * sphereVolume;

							Vector3r v1;
							Vector3r v2;
							Vector3r v3;
							Vector3r v4;
							bm_neighbor->getPointVelocity(x1, v1);
							bm_neighbor->getPointVelocity(x2, v2);
							bm_neighbor->getPointVelocity(x3, v3);
							bm_neighbor->getPointVelocity(x4, v4);

							// compute forces for both sample point
							const Vector3r a1 = (d * m_boundaryViscosity * vol * (vi - v1).dot(xix1) / (xix1.squaredNorm() + static_cast<Real>(0.01)*h2)) * gradW1;
							const Vector3r a2 = (d * m_boundaryViscosity * vol * (vi - v2).dot(xix2) / (xix2.squaredNorm() + static_cast<Real>(0.01)*h2)) * gradW2;
							const Vector3r a3 = (d * m_boundaryViscosity * vol * (vi - v3).dot(xix3) / (xix3.squaredNorm() + static_cast<Real>(0.01)*h2)) * gradW3;
							const Vector3r a4 = (d * m_boundaryViscosity * vol * (vi - v4).dot(xix4) / (xix4.squaredNorm() + static_cast<Real>(0.01)*h2)) * gradW4;
							ai += a1 + a2 + a3 + a4;

							bm_neighbor->addForce(x1, -m_model->getMass(i) * a1);
							bm_neighbor->addForce(x2, -m_model->getMass(i) * a2);
							bm_neighbor->addForce(x3, -m_model->getMass(i) * a3);
							bm_neighbor->addForce(x4, -m_model->getMass(i) * a4);
						}
 					);
 				}
 				else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
 				{
 					forall_volume_maps(
						const Vector3r xixj = xi - xj;
						Vector3r normal = -xixj;
						const Real nl = normal.norm();
						if (nl > static_cast<Real>(0.0001))
						{
							normal /= nl;

							Vector3r t1;
							Vector3r t2;
							MathFunctions::getOrthogonalVectors(normal, t1, t2);

							const Real dist = (static_cast<Real>(1.0) - nl / h) * h;
							//const Real dist = 0.5*h;
							const Vector3r x1 = xj - t1*dist;
							const Vector3r x2 = xj + t1*dist;
							const Vector3r x3 = xj - t2*dist;
							const Vector3r x4 = xj + t2*dist;

							const Vector3r xix1 = xi - x1;
							const Vector3r xix2 = xi - x2;
							const Vector3r xix3 = xi - x3;
							const Vector3r xix4 = xi - x4;

							const Vector3r gradW1 = sim->gradW(xix1);
							const Vector3r gradW2 = sim->gradW(xix2);
							const Vector3r gradW3 = sim->gradW(xix3);
							const Vector3r gradW4 = sim->gradW(xix4);

							// each sample point represents the quarter of the volume inside of the boundary
							const Real vol = static_cast<Real>(0.25) * Vj;

							Vector3r v1;
							Vector3r v2;
							Vector3r v3;
							Vector3r v4;
							bm_neighbor->getPointVelocity(x1, v1);
							bm_neighbor->getPointVelocity(x2, v2);
							bm_neighbor->getPointVelocity(x3, v3);
							bm_neighbor->getPointVelocity(x4, v4);

							// compute forces for both sample point
							const Vector3r a1 = (d * m_boundaryViscosity * vol * (vi - v1).dot(xix1) / (xix1.squaredNorm() + static_cast<Real>(0.01)*h2)) * gradW1;
							const Vector3r a2 = (d * m_boundaryViscosity * vol * (vi - v2).dot(xix2) / (xix2.squaredNorm() + static_cast<Real>(0.01)*h2)) * gradW2;
							const Vector3r a3 = (d * m_boundaryViscosity * vol * (vi - v3).dot(xix3) / (xix3.squaredNorm() + static_cast<Real>(0.01)*h2)) * gradW3;
							const Vector3r a4 = (d * m_boundaryViscosity * vol * (vi - v4).dot(xix4) / (xix4.squaredNorm() + static_cast<Real>(0.01)*h2)) * gradW4;
							ai += a1 + a2 + a3 + a4;

							bm_neighbor->addForce(x1, -m_model->getMass(i) * a1);
							bm_neighbor->addForce(x2, -m_model->getMass(i) * a2);
							bm_neighbor->addForce(x3, -m_model->getMass(i) * a3);
							bm_neighbor->addForce(x4, -m_model->getMass(i) * a4);
						}
 					);
 				}
 			}

			ai[0] += delta_ai.x().reduce();
			ai[1] += delta_ai.y().reduce();
			ai[2] += delta_ai.z().reduce();
		}
	}
}

#else

void Viscosity_Standard::step()
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int numParticles = m_model->numActiveParticles();
	const Real h = sim->getSupportRadius();
	const Real h2 = h*h;
	const Real sphereVolume = static_cast<Real>(4.0 / 3.0 * M_PI) * h2*h;
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();
	const unsigned int fluidModelIndex = m_model->getPointSetIndex();
	const Real density0 = m_model->getDensity0();
	FluidModel *model = m_model;
	Real d = 10.0;
	if (sim->is2DSimulation())
		d = 8.0;

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			const Vector3r &xi = m_model->getPosition(i);
			const Vector3r &vi = m_model->getVelocity(i);
			Vector3r &ai = m_model->getAcceleration(i);
			const Real density_i = m_model->getDensity(i);

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors_in_same_phase(
				const Vector3r &vj = model->getVelocity(neighborIndex);

				// Viscosity
				const Real density_j = model->getDensity(neighborIndex);
				const Vector3r xixj = xi - xj;
				ai += d * m_viscosity * (model->getMass(neighborIndex) / density_j) * (vi - vj).dot(xixj) / (xixj.squaredNorm() + 0.01*h2) * sim->gradW(xi - xj);
			)

			//////////////////////////////////////////////////////////////////////////
			// Boundary
			//////////////////////////////////////////////////////////////////////////
 			if (m_boundaryViscosity != 0.0)
 			{
 				if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
 				{
 					forall_boundary_neighbors(
 						const Vector3r &vj = bm_neighbor->getVelocity(neighborIndex);
 						const Vector3r xixj = xi - xj;
 						const Vector3r a = d * m_boundaryViscosity * (density0 * bm_neighbor->getVolume(neighborIndex) / density_i) * (vi - vj).dot(xixj) / (xixj.squaredNorm() + 0.01*h2) * sim->gradW(xi - xj);
 						ai += a;
 						bm_neighbor->addForce(xj, -m_model->getMass(i) * a);
 					);
 				}
 				else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
 				{
 					forall_density_maps(
						const Vector3r xixj = xi - xj;
						Vector3r normal = -xixj;
						const Real nl = normal.norm();
						if (nl > static_cast<Real>(0.0001))
						{
							normal /= nl;

							Vector3r t1;
							Vector3r t2;
							MathFunctions::getOrthogonalVectors(normal, t1, t2);

							const Real dist = (static_cast<Real>(1.0) - nl / h) * h;
							//const Real dist = 0.25*h;
							const Vector3r x1 = xj - t1*dist;
							const Vector3r x2 = xj + t1*dist;
							const Vector3r x3 = xj - t2*dist;
							const Vector3r x4 = xj + t2*dist;

							const Vector3r xix1 = xi - x1;
							const Vector3r xix2 = xi - x2;
							const Vector3r xix3 = xi - x3;
							const Vector3r xix4 = xi - x4;

							const Vector3r gradW1 = sim->gradW(xix1);
							const Vector3r gradW2 = sim->gradW(xix2);
							const Vector3r gradW3 = sim->gradW(xix3);
							const Vector3r gradW4 = sim->gradW(xix4);

							// each sample point represents the quarter of the volume inside of the boundary
							const Real vol = static_cast<Real>(0.25) * rho * sphereVolume;

							Vector3r v1;
							Vector3r v2;
							Vector3r v3;
							Vector3r v4;
							bm_neighbor->getPointVelocity(x1, v1);
							bm_neighbor->getPointVelocity(x2, v2);
							bm_neighbor->getPointVelocity(x3, v3);
							bm_neighbor->getPointVelocity(x4, v4);

							// compute forces for both sample point
							const Vector3r a1 = d * m_boundaryViscosity * vol * (vi - v1).dot(xix1) / (xix1.squaredNorm() + 0.01*h2) * gradW1;
							const Vector3r a2 = d * m_boundaryViscosity * vol * (vi - v2).dot(xix2) / (xix2.squaredNorm() + 0.01*h2) * gradW2;
							const Vector3r a3 = d * m_boundaryViscosity * vol * (vi - v3).dot(xix3) / (xix3.squaredNorm() + 0.01*h2) * gradW3;
							const Vector3r a4 = d * m_boundaryViscosity * vol * (vi - v4).dot(xix4) / (xix4.squaredNorm() + 0.01*h2) * gradW4;
							ai += a1 + a2 + a3 + a4;

							bm_neighbor->addForce(x1, -m_model->getMass(i) * a1);
							bm_neighbor->addForce(x2, -m_model->getMass(i) * a2);
							bm_neighbor->addForce(x3, -m_model->getMass(i) * a3);
							bm_neighbor->addForce(x4, -m_model->getMass(i) * a4);
						}
 					);
 				}
 				else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
 				{
 					forall_volume_maps(
						const Vector3r xixj = xi - xj;
						Vector3r normal = -xixj;
						const Real nl = normal.norm();
						if (nl > static_cast<Real>(0.0001))
						{
							normal /= nl;

							Vector3r t1;
							Vector3r t2;
							MathFunctions::getOrthogonalVectors(normal, t1, t2);

							const Real dist = (static_cast<Real>(1.0) - nl / h) * h;
							//const Real dist = 0.25*h;
							const Vector3r x1 = xj - t1*dist;
							const Vector3r x2 = xj + t1*dist;
							const Vector3r x3 = xj - t2*dist;
							const Vector3r x4 = xj + t2*dist;

							const Vector3r xix1 = xi - x1;
							const Vector3r xix2 = xi - x2;
							const Vector3r xix3 = xi - x3;
							const Vector3r xix4 = xi - x4;

							const Vector3r gradW1 = sim->gradW(xix1);
							const Vector3r gradW2 = sim->gradW(xix2);
							const Vector3r gradW3 = sim->gradW(xix3);
							const Vector3r gradW4 = sim->gradW(xix4);

							// each sample point represents the quarter of the volume inside of the boundary
							const Real vol = static_cast<Real>(0.25) * Vj;

							Vector3r v1;
							Vector3r v2;
							Vector3r v3;
							Vector3r v4;
							bm_neighbor->getPointVelocity(x1, v1);
							bm_neighbor->getPointVelocity(x2, v2);
							bm_neighbor->getPointVelocity(x3, v3);
							bm_neighbor->getPointVelocity(x4, v4);

							// compute forces for both sample point
							const Vector3r a1 = d * m_boundaryViscosity * vol * (vi - v1).dot(xix1) / (xix1.squaredNorm() + 0.01*h2) * gradW1;
							const Vector3r a2 = d * m_boundaryViscosity * vol * (vi - v2).dot(xix2) / (xix2.squaredNorm() + 0.01*h2) * gradW2;
							const Vector3r a3 = d * m_boundaryViscosity * vol * (vi - v3).dot(xix3) / (xix3.squaredNorm() + 0.01*h2) * gradW3;
							const Vector3r a4 = d * m_boundaryViscosity * vol * (vi - v4).dot(xix4) / (xix4.squaredNorm() + 0.01*h2) * gradW4;
							ai += a1 + a2 + a3 + a4;

							bm_neighbor->addForce(x1, -m_model->getMass(i) * a1);
							bm_neighbor->addForce(x2, -m_model->getMass(i) * a2);
							bm_neighbor->addForce(x3, -m_model->getMass(i) * a3);
							bm_neighbor->addForce(x4, -m_model->getMass(i) * a4);
						}
 					);
 				}
 			}
		}
	}
}

#endif

void Viscosity_Standard::reset()
{
}

