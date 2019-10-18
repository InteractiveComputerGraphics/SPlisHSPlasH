#include "Viscosity_Takahashi2015.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "Utilities/Timing.h"
#include "Utilities/Counting.h"
#include "../Simulation.h"
#include "SPlisHSPlasH/BoundaryModel_Akinci2012.h"
#include "SPlisHSPlasH/BoundaryModel_Koschier2017.h"
#include "SPlisHSPlasH/BoundaryModel_Bender2019.h"

using namespace SPH;
using namespace GenParam;

int Viscosity_Takahashi2015::ITERATIONS = -1;
int Viscosity_Takahashi2015::MAX_ITERATIONS = -1;
int Viscosity_Takahashi2015::MAX_ERROR = -1;

Viscosity_Takahashi2015::Viscosity_Takahashi2015(FluidModel *model) :
	ViscosityBase(model)
{
	m_viscousStress.resize(model->numParticles(), Matrix3r::Zero());
	m_accel.resize(model->numParticles(), Vector3r::Zero());

	m_maxIter = 100;
	m_maxError = 0.01;
	m_iterations = 0;

	model->addField({ "viscous stress", FieldType::Matrix3, [&](const unsigned int i) -> Real* { return &m_viscousStress[i](0,0); } });
	model->addField({ "accel (visco)", FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_accel[i][0]; } });
}

Viscosity_Takahashi2015::~Viscosity_Takahashi2015(void)
{
	m_model->removeFieldByName("viscous stress");
	m_model->removeFieldByName("accel (visco)");

	m_viscousStress.clear();
	m_accel.clear();
}

void Viscosity_Takahashi2015::initParameters()
{
	ViscosityBase::initParameters();

	ITERATIONS = createNumericParameter("viscoIterations", "Iterations", &m_iterations);
	setGroup(ITERATIONS, "Viscosity");
	setDescription(ITERATIONS, "Iterations required by the viscosity solver.");
	getParameter(ITERATIONS)->setReadOnly(true);

	MAX_ITERATIONS = createNumericParameter("viscoMaxIter", "Max. iterations (visco)", &m_maxIter);
	setGroup(MAX_ITERATIONS, "Viscosity");
	setDescription(MAX_ITERATIONS, "Max. iterations of the viscosity solver.");
	static_cast<NumericParameter<unsigned int>*>(getParameter(MAX_ITERATIONS))->setMinValue(1);

	MAX_ERROR = createNumericParameter("viscoMaxError", "Max. visco error", &m_maxError);
	setGroup(MAX_ERROR, "Viscosity");
	setDescription(MAX_ERROR, "Max. error of the viscosity solver.");
	RealParameter* rparam = static_cast<RealParameter*>(getParameter(MAX_ERROR));
	rparam->setMinValue(1e-6);
}

void Viscosity_Takahashi2015::matrixVecProd(const Real* vec, Real *result, void *userData)
{
	Viscosity_Takahashi2015 *visco = (Viscosity_Takahashi2015*)userData;
	FluidModel *model = visco->getModel();
	const unsigned int numParticles = model->numActiveParticles();
	const Real h = TimeManager::getCurrent()->getTimeStepSize();

	computeViscosityAcceleration(visco, vec);

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static) 
		for (int i = 0; i < (int)numParticles; i++)
		{
			const Vector3r &ai = visco->getAccel(i);
			result[3*i] = vec[3*i] - h*ai[0];
			result[3*i+1] = vec[3*i+1] - h*ai[1];
			result[3*i+2] = vec[3*i+2] - h*ai[2];
		}
	}
}

void Viscosity_Takahashi2015::computeViscosityAcceleration(Viscosity_Takahashi2015 *visco, const Real* v)
{
	Simulation *sim = Simulation::getCurrent();
	FluidModel *model = visco->getModel();
	const unsigned int numParticles = model->numActiveParticles();
	const unsigned int fluidModelIndex = model->getPointSetIndex();
	
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static) 
		for (int i = 0; i < (int)numParticles; i++)
		{
			const Vector3r &xi = model->getPosition(i);
			const Vector3r &vi = Eigen::Map<const Vector3r>(&v[3*i]);
			const Real density_i = model->getDensity(i);

			Matrix3r nablaV;
			nablaV.setZero();

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors_in_same_phase(
				const Vector3r &vj = Eigen::Map<const Vector3r>(&v[3 * neighborIndex]);
				const Vector3r gradW = sim->gradW(xi - xj);

				const Matrix3r dyad = (vj - vi) * gradW.transpose();

				nablaV += (model->getMass(neighborIndex) / model->getDensity(neighborIndex)) * dyad;
			)
			visco->getViscousStress(i) = visco->m_viscosity * (nablaV + nablaV.transpose());
		}

		#pragma omp for schedule(static) 
		for (int i = 0; i < (int)numParticles; i++)
		{
			const Vector3r &xi = model->getPosition(i);
			Vector3r &ai = visco->getAccel(i);
			ai.setZero();
			const Real density_i = model->getDensity(i);
			const Real density_i_2 = density_i*density_i;

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors_in_same_phase(
				const Real density_j = model->getDensity(neighborIndex);
				const Real density_j_2 = density_j*density_j;
				const Vector3r gradW = sim->gradW(xi - xj);

				ai += model->getMass(neighborIndex) * (visco->getViscousStress(i) / density_i_2 + visco->getViscousStress(neighborIndex) / density_j_2) * gradW;
			)
		}
	}
}


void Viscosity_Takahashi2015::step()
{
	Simulation *sim = Simulation::getCurrent();
	const int numParticles = (int) m_model->numActiveParticles();
	if (numParticles == 0)
		return;

	const Real density0 = m_model->getDensity0();
	const Real h = TimeManager::getCurrent()->getTimeStepSize();
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();
	const unsigned int fluidModelIndex = m_model->getPointSetIndex();

	//////////////////////////////////////////////////////////////////////////
	// Init linear system solver and preconditioner
	//////////////////////////////////////////////////////////////////////////
	MatrixReplacement A(3*m_model->numActiveParticles(), matrixVecProd, (void*) this);

	m_solver.setTolerance(m_maxError);
	m_solver.setMaxIterations(m_maxIter);
	m_solver.compute(A);

	VectorXr b(3*numParticles);
	VectorXr x(3*numParticles);
	x.setZero();

	//////////////////////////////////////////////////////////////////////////
	// Compute RHS
	//////////////////////////////////////////////////////////////////////////
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static) nowait 
		for (int i = 0; i < (int)numParticles; i++)
		{
			const Vector3r &vi = m_model->getVelocity(i);
			b[3*i] = vi[0];
			b[3*i+1] = vi[1];
			b[3*i+2] = vi[2];
		}
	}

	//////////////////////////////////////////////////////////////////////////
	// Solve linear system 
	//////////////////////////////////////////////////////////////////////////
	START_TIMING("CG solve");
	x = m_solver.solve(b);
	m_iterations = (int)m_solver.iterations();
	STOP_TIMING_AVG;
	INCREASE_COUNTER("Visco iterations", static_cast<Real>(m_iterations));

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			Vector3r &ai = m_model->getAcceleration(i);
			const Vector3r newV(x[3 * i], x[3 * i + 1], x[3 * i + 2]);
			ai += (1.0 / h) * (newV - m_model->getVelocity(i));
		}
	}

	// Compute viscosity forces (XSPH) with boundary to simulate simple friction
	const Real invH = (static_cast<Real>(1.0) / h);
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
			// Boundary
			//////////////////////////////////////////////////////////////////////////
			if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012)
			{
				forall_boundary_neighbors(
					const Vector3r &vj = bm_neighbor->getVelocity(neighborIndex);
					const Vector3r a = -invH * 0.1 * (density0 * bm_neighbor->getVolume(neighborIndex) / density_i) * (vi - vj)* sim->W(xi - xj);
					ai += a;
					bm_neighbor->addForce(xj, -m_model->getMass(i) * a);
				);
			}
			else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
			{
				forall_density_maps(
					Vector3r vj;
					bm_neighbor->getPointVelocity(xi, vj);
					const Vector3r a = -invH * 0.1 * (density0 / density_i) * (vi-vj)* rho;
					ai += a;
					bm_neighbor->addForce(xj, -m_model->getMass(i) * a);
				);
			}
			else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Bender2019)
			{
				forall_volume_maps(
					Vector3r vj;
					bm_neighbor->getPointVelocity(xj, vj);
					const Vector3r a = -invH * 0.1 * (density0 * Vj / density_i) * (vi-vj)* sim->W(xi - xj);
					ai += a;
					bm_neighbor->addForce(xj, -m_model->getMass(i) * a);
				);
			}
		}
	}
}


void Viscosity_Takahashi2015::reset()
{
}

void Viscosity_Takahashi2015::performNeighborhoodSearchSort()
{
}
