#include "Viscosity_Peer2015.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "Utilities/Timing.h"
#include "Utilities/Counting.h"
#include "../Simulation.h"
#include "SPlisHSPlasH/BoundaryModel_Akinci2012.h"
#include "SPlisHSPlasH/BoundaryModel_Koschier2017.h"
#include "SPlisHSPlasH/BoundaryModel_Bender2019.h"


using namespace SPH;
using namespace GenParam;

int Viscosity_Peer2015::ITERATIONS = -1;
int Viscosity_Peer2015::MAX_ITERATIONS = -1;
int Viscosity_Peer2015::MAX_ERROR = -1;


Viscosity_Peer2015::Viscosity_Peer2015(FluidModel *model) :
	ViscosityBase(model)
{
	m_density.resize(model->numParticles(), 0.0);
	m_targetNablaV.resize(model->numParticles(), Matrix3r::Zero());

	m_iterations = 0;
	m_maxIter = 50;
	m_maxError = 0.01;

	model->addField({ "target nablaV", FieldType::Matrix3, [&](const unsigned int i) -> Real* { return &m_targetNablaV[i](0,0); } });
}

Viscosity_Peer2015::~Viscosity_Peer2015(void)
{
	m_model->removeFieldByName("target nablaV");

	m_density.clear();
	m_targetNablaV.clear();
}

void Viscosity_Peer2015::initParameters()
{
	ViscosityBase::initParameters();

	ITERATIONS = createNumericParameter("viscoIterations", "Iterations", &m_iterations);
	setGroup(ITERATIONS, "Fluid Model|Viscosity");
	setDescription(ITERATIONS, "Iterations required by the viscosity solver.");
	getParameter(ITERATIONS)->setReadOnly(true);

	MAX_ITERATIONS = createNumericParameter("viscoMaxIter", "Max. iterations (visco)", &m_maxIter);
	setGroup(MAX_ITERATIONS, "Fluid Model|Viscosity");
	setDescription(MAX_ITERATIONS, "Max. iterations of the viscosity solver.");
	static_cast<NumericParameter<unsigned int>*>(getParameter(MAX_ITERATIONS))->setMinValue(1);

	MAX_ERROR = createNumericParameter("viscoMaxError", "Max. visco error", &m_maxError);
	setGroup(MAX_ERROR, "Fluid Model|Viscosity");
	setDescription(MAX_ERROR, "Max. error of the viscosity solver.");
	RealParameter* rparam = static_cast<RealParameter*>(getParameter(MAX_ERROR));
	rparam->setMinValue(1e-6);
}

void Viscosity_Peer2015::computeDensities()
{
	Simulation* sim = Simulation::getCurrent();
	FluidModel* model = m_model;
	const unsigned int fluidModelIndex = model->getPointSetIndex();
	const Real density0 = model->getDensity0();
	const unsigned int numParticles = model->numActiveParticles();
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int nBoundaries = sim->numberOfBoundaryModels();

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			Real& density = m_density[i];

			// Compute current density for particle i
			density = model->getVolume(i) * sim->W_zero();
			const Vector3r& xi = model->getPosition(i);

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors(
				density += fm_neighbor->getVolume(neighborIndex) * sim->W(xi - xj);
			);

			//////////////////////////////////////////////////////////////////////////
			// Boundary
			//////////////////////////////////////////////////////////////////////////
			if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Akinci2012 || sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Gissler2019)
			{
				forall_boundary_neighbors(
					// Boundary: Akinci2012
					density += bm_neighbor->getVolume(neighborIndex) * sim->W(xi - xj);
				);
			}
			else if (sim->getBoundaryHandlingMethod() == BoundaryHandlingMethods::Koschier2017)
			{
				forall_density_maps(
					density += rho;
				);
			}
			else   // Bender2019
			{
				forall_volume_maps(
					density += Vj * sim->W(xi - xj);
				);
			}

			density *= density0;
		}
	}
}

void Viscosity_Peer2015::matrixVecProd(const Real* vec, Real *result, void *userData)
{
	Simulation *sim = Simulation::getCurrent();
	Viscosity_Peer2015* visco = (Viscosity_Peer2015*)userData;
	FluidModel* model = visco->getModel();
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
			result[i] = (visco->m_density[i] - model->getMass(i) * sim->W_zero()) * vec[i];

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors_in_same_phase(
				result[i] -= model->getMass(neighborIndex) * sim->W(xi - xj) * vec[neighborIndex];
			)
		}
	}
}

void Viscosity_Peer2015::diagonalMatrixElement(const unsigned int row, Real &result, void *userData)
{
	// Diagonal element
	Viscosity_Peer2015* visco = (Viscosity_Peer2015*)userData;
	FluidModel* model = visco->getModel();
	Simulation *sim = Simulation::getCurrent();
	result = visco->m_density[row] - model->getDensity0() * model->getVolume(row) * sim->W_zero();
}

void Viscosity_Peer2015::step()
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

	// this method computes its own density values since
	// it is very sensitive to small deviations (like they
	// appear when using AVX)
	computeDensities();

	// Compute target
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static) nowait 
		for (int i = 0; i < numParticles; i++)
		{
			const Vector3r &xi = m_model->getPosition(i);
			const Vector3r &vi = m_model->getVelocity(i);
			const Real density_i = m_model->getDensity(i);

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

			Matrix3r &target = getTargetNablaV(i);
			const Matrix3r R = 0.5 * (nablaV - nablaV.transpose());
			const Real divergence = nablaV(0, 0) + nablaV(1, 1) + nablaV(2, 2);
			const Matrix3r V = (1.0 / 3.0) * divergence * Matrix3r::Identity();
			const Matrix3r S = 0.5 * (nablaV + nablaV.transpose()) - V;
			if (density_i >= density0)
			{
				target = R + V + viscosity * S;
			}
			else
			{
				if (-divergence < 0.0)
					target = R + V + viscosity * S;
				else
					target = R + viscosity * S;
			}
		}
	}

	//////////////////////////////////////////////////////////////////////////
	// Init linear system solver and preconditioner
	//////////////////////////////////////////////////////////////////////////
	MatrixReplacement A(m_model->numActiveParticles(), matrixVecProd, (void*) this);
	m_solver.preconditioner().init(m_model->numActiveParticles(), diagonalMatrixElement, (void*)this);

	m_solver.setTolerance(m_maxError);
	m_solver.setMaxIterations(m_maxIter);
	m_solver.compute(A);

	VectorXr b0(numParticles);
	VectorXr b1(numParticles);
	VectorXr b2(numParticles);
	VectorXr x0(numParticles);
	VectorXr x1(numParticles);
	VectorXr x2(numParticles);
	VectorXr g0(numParticles);
	VectorXr g1(numParticles);
	VectorXr g2(numParticles);

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
	m_iterations = 0;
	x0 = m_solver.solveWithGuess(b0, g0);
	if (m_solver.iterations() == 0)
		x0 = g0;
	m_iterations += (int)m_solver.iterations();
	x1 = m_solver.solveWithGuess(b1, g1);
	if (m_solver.iterations() == 0)
		x1 = g1;
	m_iterations += (int)m_solver.iterations();
	x2 = m_solver.solveWithGuess(b2, g2);
	if (m_solver.iterations() == 0)
		x2 = g2;
	m_iterations += (int)m_solver.iterations();
	STOP_TIMING_AVG;
	INCREASE_COUNTER("Visco iterations", static_cast<Real>(m_iterations));

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


void Viscosity_Peer2015::reset()
{
}

void Viscosity_Peer2015::performNeighborhoodSearchSort()
{
}
