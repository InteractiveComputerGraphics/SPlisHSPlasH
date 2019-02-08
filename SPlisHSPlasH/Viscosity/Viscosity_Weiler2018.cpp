#include "Viscosity_Weiler2018.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "Utilities/Timing.h"
#include "Utilities/Counting.h"
#include "../Simulation.h"

using namespace SPH;
using namespace GenParam;

int Viscosity_Weiler2018::ITERATIONS = -1;
int Viscosity_Weiler2018::MAX_ITERATIONS = -1;
int Viscosity_Weiler2018::MAX_ERROR = -1;
int Viscosity_Weiler2018::VISCOSITY_COEFFICIENT_BOUNDARY = -1;

Viscosity_Weiler2018::Viscosity_Weiler2018(FluidModel *model) :
	ViscosityBase(model), m_vDiff()
{
	m_maxIter = 100;
	m_maxError = 0.01;
	m_iterations = 0;
	m_boundaryViscosity = 0.0;

	m_vDiff.resize(model->numParticles(), Vector3r::Zero());

	model->addField({ "velocity difference", FieldType::Vector3, [&](const unsigned int i) -> Real* { return &m_vDiff[i][0]; } });
}

Viscosity_Weiler2018::~Viscosity_Weiler2018(void)
{
	m_model->removeFieldByName("velocity difference");

	m_vDiff.clear();
}

void Viscosity_Weiler2018::initParameters()
{
	ViscosityBase::initParameters();

	VISCOSITY_COEFFICIENT_BOUNDARY = createNumericParameter("viscosityBoundary", "Viscosity coefficient (Boundary)", &m_boundaryViscosity);
	setGroup(VISCOSITY_COEFFICIENT_BOUNDARY, "Viscosity");
	setDescription(VISCOSITY_COEFFICIENT_BOUNDARY, "Coefficient for the viscosity force computation at the boundary.");
	RealParameter* rparam = static_cast<RealParameter*>(getParameter(VISCOSITY_COEFFICIENT_BOUNDARY));
	rparam->setMinValue(0.0);

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
	rparam = static_cast<RealParameter*>(getParameter(MAX_ERROR));
	rparam->setMinValue(1e-6);
}

void Viscosity_Weiler2018::matrixVecProd(const Real* vec, Real *result, void *userData)
{
	Simulation *sim = Simulation::getCurrent();
	Viscosity_Weiler2018 *visco = (Viscosity_Weiler2018*)userData;
	FluidModel *model = visco->getModel();
	const unsigned int numParticles = model->numActiveParticles();
	const unsigned int fluidModelIndex = model->getPointSetIndex();
	const unsigned int nFluids = sim->numberOfFluidModels();

	const Real h = sim->getSupportRadius();
	const Real h2 = h*h;
	const Real dt = TimeManager::getCurrent()->getTimeStepSize();
	const Real mu = visco->m_viscosity;
	const Real mub = visco->m_boundaryViscosity;
	const Real density0 = model->getDensity0();

	Real d = 10.0;
	if (sim->is2DSimulation())
		d = 6.0;

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static) 
		for (int i = 0; i < (int)numParticles; i++)
		{
			const Vector3r &xi = model->getPosition(i);
			Vector3r ai;
			ai.setZero();
			const Real density_i = model->getDensity(i);
			const Vector3r &vi = Eigen::Map<const Vector3r>(&vec[3 * i]);

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			forall_fluid_neighbors_in_same_phase(
				const Real density_j = model->getDensity(neighborIndex);
				const Vector3r gradW = sim->gradW(xi - xj);

				const Vector3r &vj = Eigen::Map<const Vector3r>(&vec[3 * neighborIndex]);
				const Vector3r xixj = xi - xj;

				ai += d * mu * (model->getMass(neighborIndex) / density_j) * (vi - vj).dot(xixj) / (xixj.squaredNorm() + 0.01*h2) * gradW;
			)

			//////////////////////////////////////////////////////////////////////////
			// Boundary
			//////////////////////////////////////////////////////////////////////////
			forall_boundary_neighbors(
				const Vector3r &vj = bm_neighbor->getVelocity(neighborIndex);
				const Vector3r gradW = sim->gradW(xi - xj);

				const Vector3r xixj = xi - xj;
				ai += d * mub * (density0 * bm_neighbor->getVolume(neighborIndex) / density_i) * (vi - vj).dot(xixj) / (xixj.squaredNorm() + 0.01*h2) * gradW;
			)

			result[3 * i] = vec[3 * i] - dt / density_i*ai[0];
			result[3 * i + 1] = vec[3 * i + 1] - dt / density_i*ai[1];
			result[3 * i + 2] = vec[3 * i + 2] - dt / density_i*ai[2];
		}
	}
}

#ifdef USE_BLOCKDIAGONAL_PRECONDITIONER
void Viscosity_Weiler2018::diagonalMatrixElement(const unsigned int i, Matrix3r &result, void *userData)
{
	// Diagonal element
	Simulation *sim = Simulation::getCurrent();
	Viscosity_Weiler2018 *visco = (Viscosity_Weiler2018*)userData;
	FluidModel *model = visco->getModel();
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int fluidModelIndex = model->getPointSetIndex();
	const Real density0 = model->getDensity0();

	Real d = 10.0;
	if (sim->is2DSimulation())
		d = 6.0;

	const Real h = sim->getSupportRadius();
	const Real h2 = h*h;
	const Real dt = TimeManager::getCurrent()->getTimeStepSize();
	const Real mu = visco->m_viscosity;
	const Real mub = visco->m_boundaryViscosity;

	const Real density_i = model->getDensity(i);

	result.setZero();
	
	const Vector3r &xi = model->getPosition(i);

	//////////////////////////////////////////////////////////////////////////
	// Fluid
	//////////////////////////////////////////////////////////////////////////
	forall_fluid_neighbors_in_same_phase(
		const Real density_j = model->getDensity(neighborIndex);
		const Vector3r gradW = sim->gradW(xi - xj);
		const Vector3r xixj = xi - xj;
		result += d * mu * (model->getMass(neighborIndex) / density_j) / (xixj.squaredNorm() + 0.01*h2) * (gradW * xixj.transpose());
	)

	//////////////////////////////////////////////////////////////////////////
	// Boundary
	//////////////////////////////////////////////////////////////////////////
	forall_boundary_neighbors(
		const Vector3r &vj = bm_neighbor->getVelocity(neighborIndex);
		const Vector3r gradW = sim->gradW(xi - xj);

		const Vector3r xixj = xi - xj;
		result += d * mub * (density0 * bm_neighbor->getVolume(neighborIndex) / density_i) / (xixj.squaredNorm() + 0.01*h2) * (gradW * xixj.transpose());
	)
	result = Matrix3r::Identity() - (dt / density_i) * result;
}

#else

void Viscosity_Weiler2018::diagonalMatrixElement(const unsigned int i, Vector3r &result, void *userData)
{
	// Diagonal element
	Viscosity_Weiler2018 *visco = (Viscosity_Weiler2018*)userData;
	FluidModel *model = visco->getModel();

	const Real h = model->getSupportRadius();
	const Real h2 = h*h;
	const Real dt = TimeManager::getCurrent()->getTimeStepSize();
	const Real mu = visco->m_viscosity;
	const Real mub = visco->m_boundaryViscosity;
	const unsigned int nFluids = sim->numberOfFluidModels();
	const unsigned int fluidModelIndex = m_model->getPointSetIndex();
	const Real density0 = model->getDensity0();

	Real d = 10.0;
	if (sim->is2DSimulation())
		d = 6.0;

	const Real density_i = model->getDensity(i);

	result.setZero();

	const Vector3r &xi = model->getPosition(0, i);

	//////////////////////////////////////////////////////////////////////////
	// Fluid
	//////////////////////////////////////////////////////////////////////////
	forall_fluid_neighbors_in_same_phase(
		const Real density_j = model->getDensity(neighborIndex);
		const Vector3r gradW = model->gradW(xi - xj);
		const Vector3r xixj = xi - xj;
		Matrix3r r = d * mu * (model->getMass(neighborIndex) / density_j) / (xixj.squaredNorm() + 0.01*h2) * (gradW * xixj.transpose());
		result += r.diagonal();
	)

	//////////////////////////////////////////////////////////////////////////
	// Boundary
	//////////////////////////////////////////////////////////////////////////
	forall_boundary_neighbors(
		const Vector3r &vj = model->getVelocity(pid, neighborIndex);
		const Vector3r gradW = model->gradW(xi - xj);

		const Vector3r xixj = xi - xj;
		Matrix3r r = d * mub * (density0 * bm_neighbor->getVolume(neighborIndex) / density_i) / (xixj.squaredNorm() + 0.01*h2) * (gradW * xixj.transpose());
		result += r.diagonal();
	)
	result = Vector3r::Ones() - (dt / density_i) * result;
}

#endif


void Viscosity_Weiler2018::step()
{
	const int numParticles = (int) m_model->numActiveParticles();
	// prevent solver from running with a zero-length vector
	if (numParticles == 0)
		return;
	const Real density0 = m_model->getDensity0();
	const Real h = TimeManager::getCurrent()->getTimeStepSize();

	//////////////////////////////////////////////////////////////////////////
	// Init linear system solver and preconditioner
	//////////////////////////////////////////////////////////////////////////
	MatrixReplacement A(3*m_model->numActiveParticles(), matrixVecProd, (void*) this);
	m_solver.preconditioner().init(m_model->numActiveParticles(), diagonalMatrixElement, (void*)this);

	m_solver.setTolerance(m_maxError);
	m_solver.setMaxIterations(m_maxIter);
	m_solver.compute(A);

	VectorXr b(3*numParticles);
	VectorXr x(3*numParticles);
	VectorXr g(3*numParticles);

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

			// Warmstart
			g[3 * i] = vi[0]+m_vDiff[i][0];
			g[3 * i + 1] = vi[1]+m_vDiff[i][1];
			g[3 * i + 2] = vi[2]+m_vDiff[i][2];
		}
	}

	//////////////////////////////////////////////////////////////////////////
	// Solve linear system 
	//////////////////////////////////////////////////////////////////////////
	START_TIMING("CG solve");
	x = m_solver.solveWithGuess(b, g);
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
			m_vDiff[i] = (newV - m_model->getVelocity(i));
		}
	}
}


void Viscosity_Weiler2018::reset()
{
}

void Viscosity_Weiler2018::performNeighborhoodSearchSort()
{
}
