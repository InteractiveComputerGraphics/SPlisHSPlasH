#include "Viscosity_Peer2016.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "Utilities/Timing.h"
#include "Utilities/Counting.h"
#include "../Simulation.h"

using namespace SPH;
using namespace GenParam;

int Viscosity_Peer2016::ITERATIONS_V = -1;
int Viscosity_Peer2016::ITERATIONS_OMEGA = -1;
int Viscosity_Peer2016::MAX_ITERATIONS_V = -1;
int Viscosity_Peer2016::MAX_ERROR_V = -1;
int Viscosity_Peer2016::MAX_ITERATIONS_OMEGA = -1;
int Viscosity_Peer2016::MAX_ERROR_OMEGA = -1;


Viscosity_Peer2016::Viscosity_Peer2016(FluidModel *model) :
	ViscosityBase(model)
{
	m_targetNablaV.resize(model->numParticles(), Matrix3r::Zero());
	m_omega.resize(model->numParticles(), Vector3r::Zero());

	m_iterationsV = 0;
	m_iterationsOmega = 0;
	m_maxIterV = 50;
	m_maxErrorV = 0.01;
	m_maxIterOmega = 50;
	m_maxErrorOmega = 0.01;

	model->addField({ "target nablaV", FieldType::Matrix3, [&](const unsigned int i) -> Real* { return &m_targetNablaV[i](0,0); } });
	model->addField({ "omega (visco)", FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_omega[i][0]; } });
}

Viscosity_Peer2016::~Viscosity_Peer2016(void)
{
	m_model->removeFieldByName("target nablaV");
	m_model->removeFieldByName("omega (visco)");

	m_targetNablaV.clear();
	m_omega.clear();
}

void Viscosity_Peer2016::initParameters()
{
	ViscosityBase::initParameters();

	ITERATIONS_V = createNumericParameter("viscoIterationsV", "Iterations (velocity field)", &m_iterationsV);
	setGroup(ITERATIONS_V, "Viscosity");
	setDescription(ITERATIONS_V, "Iterations required by the viscosity solver.");
	getParameter(ITERATIONS_V)->setReadOnly(true);

	ITERATIONS_OMEGA = createNumericParameter("viscoIterationsOmega", "Iterations (vorticity diffusion)", &m_iterationsOmega);
	setGroup(ITERATIONS_OMEGA, "Viscosity");
	setDescription(ITERATIONS_OMEGA, "Iterations required by the viscosity solver.");
	getParameter(ITERATIONS_OMEGA)->setReadOnly(true);

	MAX_ITERATIONS_V = createNumericParameter("viscoMaxIter", "Max. iterations", &m_maxIterV);
	setGroup(MAX_ITERATIONS_V, "Viscosity");
	setDescription(MAX_ITERATIONS_V, "Max. iterations of the viscosity solver.");
	static_cast<NumericParameter<unsigned int>*>(getParameter(MAX_ITERATIONS_V))->setMinValue(1);

	MAX_ERROR_V = createNumericParameter("viscoMaxError", "Max. error", &m_maxErrorV);
	setGroup(MAX_ERROR_V, "Viscosity");
	setDescription(MAX_ERROR_V, "Max. error of the viscosity solver.");
	RealParameter* rparam = static_cast<RealParameter*>(getParameter(MAX_ERROR_V));
	rparam->setMinValue(1e-6);

	MAX_ITERATIONS_OMEGA = createNumericParameter("viscoMaxIterOmega", "Max. iterations (vorticity diffusion)", &m_maxIterOmega);
	setGroup(MAX_ITERATIONS_OMEGA, "Viscosity");
	setDescription(MAX_ITERATIONS_OMEGA, "Max. iterations of the vorticity diffusion solver.");
	static_cast<NumericParameter<unsigned int>*>(getParameter(MAX_ITERATIONS_OMEGA))->setMinValue(1);

	MAX_ERROR_OMEGA = createNumericParameter("viscoMaxErrorOmega", "Max. vorticity diffusion error", &m_maxErrorOmega);
	setGroup(MAX_ERROR_OMEGA, "Viscosity");
	setDescription(MAX_ERROR_OMEGA, "Max. error of the vorticity diffusion solver.");
	rparam = static_cast<RealParameter*>(getParameter(MAX_ERROR_OMEGA));
	rparam->setMinValue(1e-6);
}

void Viscosity_Peer2016::matrixVecProdV(const Real* vec, Real *result, void *userData)
{
	Simulation *sim = Simulation::getCurrent();
	FluidModel *model = (FluidModel*)userData;
	const unsigned int numParticles = model->numActiveParticles();
	if (numParticles == 0)
		return;

	const unsigned int fluidModelIndex = model->getPointSetIndex();
	const unsigned int nFluids = sim->numberOfFluidModels();

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static) 
		for (int i = 0; i < (int)numParticles; i++)
		{
			// Diagonal element
			const Vector3r &xi = model->getPosition(i);
			result[i] = (model->getDensity(i) - model->getMass(i) * sim->W_zero()) * vec[i];

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors_in_same_phase(
				result[i] -= model->getMass(neighborIndex) * sim->W(xi - xj) * vec[neighborIndex];
			)
		}
	}
}

void Viscosity_Peer2016::diagonalMatrixElementV(const unsigned int i, Real &result, void *userData)
{
	// Diagonal element
	FluidModel *model = (FluidModel*)userData;
	Simulation *sim = Simulation::getCurrent();
	result = model->getDensity(i) - model->getMass(i) * sim->W_zero();
}

void Viscosity_Peer2016::matrixVecProdOmega(const Real* vec, Real *result, void *userData)
{
	Simulation *sim = Simulation::getCurrent();
	FluidModel *model = (FluidModel*)userData;
	const unsigned int numParticles = model->numActiveParticles();
	if (numParticles == 0)
		return;

	const unsigned int fluidModelIndex = model->getPointSetIndex();
	const unsigned int nFluids = sim->numberOfFluidModels();

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static) 
		for (int i = 0; i < (int)numParticles; i++)
		{
			// Diagonal element
			const Vector3r &xi = model->getPosition(i);


			// Compute current fluid density for particle i
			Real density_i = model->getMass(i) * sim->W_zero();

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors_in_same_phase(
				density_i += model->getMass(neighborIndex) * sim->W(xi - xj);
			)


			result[i] = (density_i - model->getMass(i) * sim->W_zero()) * vec[i];

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors_in_same_phase(
				result[i] -= model->getMass(neighborIndex) * sim->W(xi - xj) * vec[neighborIndex];
			)
		}
	}
}

void Viscosity_Peer2016::diagonalMatrixElementOmega(const unsigned int i, Real &result, void *userData)
{
	// Diagonal element
	Simulation *sim = Simulation::getCurrent();
	FluidModel *model = (FluidModel*)userData;
	const unsigned int fluidModelIndex = model->getPointSetIndex();
	const unsigned int nFluids = sim->numberOfFluidModels();

	const Vector3r &xi = model->getPosition(i);
	// Compute current fluid density for particle i
	Real density_i = model->getMass(i) * sim->W_zero();

	//////////////////////////////////////////////////////////////////////////
	// Fluid
	//////////////////////////////////////////////////////////////////////////
	forall_fluid_neighbors_in_same_phase(
		density_i += model->getMass(neighborIndex) * sim->W(xi - xj);
	)

	result = density_i - model->getMass(i) * sim->W_zero();
}

void Viscosity_Peer2016::step()
{
	Simulation *sim = Simulation::getCurrent();
	const int numParticles = (int) m_model->numActiveParticles();
	if (numParticles == 0)
		return;

	const Real viscosity = static_cast<Real>(1.0) - m_viscosity;
	const Real density0 = m_model->getDensity0();
	const unsigned int fluidModelIndex = m_model->getPointSetIndex();
	const unsigned int nFluids = sim->numberOfFluidModels();
	FluidModel *model = m_model;

	const Real h = TimeManager::getCurrent()->getTimeStepSize();

	//////////////////////////////////////////////////////////////////////////
	// Init linear system solver and preconditioner
	//////////////////////////////////////////////////////////////////////////
	MatrixReplacement A(m_model->numActiveParticles(), matrixVecProdV, (void*)m_model);
	m_solverV.preconditioner().init(m_model->numActiveParticles(), diagonalMatrixElementV, (void*)m_model);
	m_solverV.setTolerance(m_maxErrorV);
	m_solverV.setMaxIterations(m_maxIterV);
	m_solverV.compute(A);

	MatrixReplacement A2(m_model->numActiveParticles(), matrixVecProdOmega, (void*)m_model);
	m_solverOmega.preconditioner().init(m_model->numActiveParticles(), diagonalMatrixElementOmega, (void*)m_model);
	m_solverOmega.setTolerance(m_maxErrorOmega);
	m_solverOmega.setMaxIterations(m_maxIterOmega);
	m_solverOmega.compute(A2);

	VectorXr b0(numParticles);
	VectorXr b1(numParticles);
	VectorXr b2(numParticles);
	VectorXr x0(numParticles);
	VectorXr x1(numParticles);
	VectorXr x2(numParticles);
	VectorXr g0(numParticles);
	VectorXr g1(numParticles);
	VectorXr g2(numParticles);


	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static) nowait 
		for (int i = 0; i < numParticles; i++)
		{
			const Vector3r &xi = m_model->getPosition(i);
			const Vector3r &vi = m_model->getVelocity(i);
			const Real density_i = m_model->getDensity(i);

			//////////////////////////////////////////////////////////////////////////
			// compute nabla v
			//////////////////////////////////////////////////////////////////////////
			Matrix3r nablaV;
			nablaV.setZero();

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors_in_same_phase(
				const Vector3r &vj = m_model->getVelocity(neighborIndex);
				const Vector3r gradW = sim->gradW(xi - xj);

				Matrix3r dyad = (vj - vi) * gradW.transpose();

				nablaV += (1.0 / density_i) * m_model->getMass(neighborIndex) * dyad;
			)

			//////////////////////////////////////////////////////////////////////////
			// decomposition of velocity gradient
			//////////////////////////////////////////////////////////////////////////
			Matrix3r &target = getTargetNablaV(i);
			Matrix3r R = 0.5 * (nablaV - nablaV.transpose());
			const Real divergence = nablaV(0, 0) + nablaV(1, 1) + nablaV(2, 2);
			const Matrix3r V = (1.0 / 3.0) * divergence * Matrix3r::Identity();
			const Matrix3r S = 0.5 * (nablaV + nablaV.transpose()) - V;

			//////////////////////////////////////////////////////////////////////////
			// extract omega
			//////////////////////////////////////////////////////////////////////////
			Vector3r &omega = getOmega(i);
			omega[0] = static_cast<Real>(2.0) * R(2, 1);
			omega[1] = static_cast<Real>(2.0) * R(0, 2);
			omega[2] = static_cast<Real>(2.0) * R(1, 0);

			//////////////////////////////////////////////////////////////////////////
			// compute target nabla v without spin tensor
			//////////////////////////////////////////////////////////////////////////
			if (density_i >= density0)
			{
				target = V + viscosity * S;
			}
			else
			{
				if (-divergence < 0.0)
					target = V + viscosity * S;
				else
					target = viscosity * S;
			}
		}
	}

	//////////////////////////////////////////////////////////////////////////
	// Compute RHS of vorticity diffusion system
	//////////////////////////////////////////////////////////////////////////
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static) nowait 
		for (int i = 0; i < (int)numParticles; i++)
		{
			const Vector3r &xi = m_model->getPosition(i);
			Vector3r rhs;
			rhs.setZero();

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors_in_same_phase(
				const Real m = m_model->getMass(neighborIndex);
				const Vector3r xij = xi - xj;
				const Real W = sim->W(xij);

				rhs += m *(m_omega[i] - m_omega[neighborIndex]) * W;
			)

			const Vector3r &omegai = getOmega(i);
 			g0[i] = omegai[0];
 			g1[i] = omegai[1];
 			g2[i] = omegai[2];
			b0[i] = viscosity * rhs[0];
			b1[i] = viscosity * rhs[1];
			b2[i] = viscosity * rhs[2];
		}
	}

	//////////////////////////////////////////////////////////////////////////
	// Solve linear system 
	//////////////////////////////////////////////////////////////////////////
	START_TIMING("CG solve omega");
	m_iterationsOmega = 0;
	x0 = m_solverOmega.solveWithGuess(b0, g0);
	//x0 = m_solverOmega.solve(b0);
	if (m_solverOmega.iterations() == 0)
		x0 = g0;
	m_iterationsOmega += (int)m_solverOmega.iterations();

	x1 = m_solverOmega.solveWithGuess(b1, g1);
	//x1 = m_solverOmega.solve(b1);
	if (m_solverOmega.iterations() == 0)
		x1 = g1;
	m_iterationsOmega += (int)m_solverOmega.iterations();

	x2 = m_solverOmega.solveWithGuess(b2, g2);
	//x2 = m_solverOmega.solve(b2);
	if (m_solverOmega.iterations() == 0)
		x2 = g2;
	m_iterationsOmega += (int)m_solverOmega.iterations();
	STOP_TIMING_AVG;
	INCREASE_COUNTER("Visco iterations - Omega", static_cast<Real>(m_iterationsOmega));

	//////////////////////////////////////////////////////////////////////////
	// Determine new spin tensor and add it to target nabla v
	//////////////////////////////////////////////////////////////////////////
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static) nowait 
		for (int i = 0; i < (int)numParticles; i++)
		{
			Matrix3r R;
			R << static_cast<Real>(0.0), -static_cast<Real>(0.5)*x2[i], static_cast<Real>(0.5)*x1[i],
				static_cast<Real>(0.5)*x2[i], 0.0, -static_cast<Real>(0.5)*x0[i],
				-static_cast<Real>(0.5)*x1[i], static_cast<Real>(0.5)*x0[i], static_cast<Real>(0.0);

			Matrix3r &target = getTargetNablaV(i);
			target += R;
		}
	}


	//////////////////////////////////////////////////////////////////////////
	// Compute RHS
	//////////////////////////////////////////////////////////////////////////
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static) nowait 
		for (int i = 0; i < (int)numParticles; i++)
		{
			const Vector3r &xi = m_model->getPosition(i);
			Vector3r rhs;
			rhs.setZero();

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors_in_same_phase(
				const Real m = m_model->getMass(neighborIndex);
				const Vector3r xij = xi - xj;
				const Real W = sim->W(xij);

				rhs += m * 0.5 * (getTargetNablaV(i) + getTargetNablaV(neighborIndex)) * xij * W;
			)

			const Vector3r &vi = m_model->getVelocity(i);
			g0[i] = vi[0];
			g1[i] = vi[1];
			g2[i] = vi[2];
			b0[i] = rhs[0];
			b1[i] = rhs[1];
			b2[i] = rhs[2];
		}
	}

	//////////////////////////////////////////////////////////////////////////
	// Solve linear system 
	//////////////////////////////////////////////////////////////////////////
	START_TIMING("CG solve");
	m_iterationsV = 0;
	x0 = m_solverV.solveWithGuess(b0, g0);
	if (m_solverV.iterations() == 0)
		x0 = g0;
	m_iterationsV += (int)m_solverV.iterations();
	x1 = m_solverV.solveWithGuess(b1, g1);
	if (m_solverV.iterations() == 0)
		x1 = g1;
	m_iterationsV += (int)m_solverV.iterations();
	x2 = m_solverV.solveWithGuess(b2, g2);
	if (m_solverV.iterations() == 0)
		x2 = g2;
	m_iterationsV += (int)m_solverV.iterations();
	STOP_TIMING_AVG;
	INCREASE_COUNTER("Visco iterations - V", static_cast<Real>(m_iterationsV));

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static) nowait 
		for (int i = 0; i < (int)numParticles; i++)
		{
			Vector3r &vi = m_model->getVelocity(i);
			vi[0] = x0[i];
			vi[1] = x1[i];
			vi[2] = x2[i];
		}
	}
}


void Viscosity_Peer2016::reset()
{
}

void Viscosity_Peer2016::performNeighborhoodSearchSort()
{
}
