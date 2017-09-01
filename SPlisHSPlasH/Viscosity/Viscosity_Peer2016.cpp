#include "Viscosity_Peer2016.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "SPlisHSPlasH/Utilities/Timing.h"

using namespace SPH;

Viscosity_Peer2016::Viscosity_Peer2016(FluidModel *model) :
	ViscosityBase(model)
{
	m_targetNablaV.resize(model->numParticles(), Matrix3r::Zero());
	m_omega.resize(model->numParticles(), Vector3r::Zero());

	m_maxIter = 100;
	m_maxError = 0.01;
}

Viscosity_Peer2016::~Viscosity_Peer2016(void)
{
	m_targetNablaV.clear();
	m_omega.clear();
}


void Viscosity_Peer2016::matrixVecProdV(const Real* vec, Real *result, void *userData)
{
	FluidModel *model = (FluidModel*)userData;
	const unsigned int numParticles = model->numActiveParticles();

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static) 
		for (int i = 0; i < (int)numParticles; i++)
		{
			// Diagonal element
			const Vector3r &xi = model->getPosition(0, i);
			result[i] = (model->getDensity(i) - model->getMass(i) * model->W_zero()) * vec[i];

			for (unsigned int j = 0; j < model->numberOfNeighbors(0, i); j++)
			{
				const unsigned int neighborIndex = model->getNeighbor(0, i, j);
				const Vector3r &xj = model->getPosition(0, neighborIndex);
				result[i] -= model->getMass(neighborIndex) * model->W(xi - xj) * vec[neighborIndex];
			}
		}
	}
}

void Viscosity_Peer2016::diagonalMatrixElementV(const unsigned int i, Real &result, void *userData)
{
	// Diagonal element
	FluidModel *model = (FluidModel*)userData;
	result = model->getDensity(i) - model->getMass(i) * model->W_zero();
}

void Viscosity_Peer2016::matrixVecProdOmega(const Real* vec, Real *result, void *userData)
{
	FluidModel *model = (FluidModel*)userData;
	const unsigned int numParticles = model->numActiveParticles();

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static) 
		for (int i = 0; i < (int)numParticles; i++)
		{
			// Diagonal element
			const Vector3r &xi = model->getPosition(0, i);


			// Compute current fluid density for particle i
			Real density_i = model->getMass(i) * model->W_zero();

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			for (unsigned int j = 0; j < model->numberOfNeighbors(0, i); j++)
			{
				const unsigned int neighborIndex = model->getNeighbor(0, i, j);
				const Vector3r &xj = model->getPosition(0, neighborIndex);
				density_i += model->getMass(neighborIndex) * model->W(xi - xj);
			}


			result[i] = (density_i - model->getMass(i) * model->W_zero()) * vec[i];

			for (unsigned int j = 0; j < model->numberOfNeighbors(0, i); j++)
			{
				const unsigned int neighborIndex = model->getNeighbor(0, i, j);
				const Vector3r &xj = model->getPosition(0, neighborIndex);
				result[i] -= model->getMass(neighborIndex) * model->W(xi - xj) * vec[neighborIndex];
			}
		}
	}
}

void Viscosity_Peer2016::diagonalMatrixElementOmega(const unsigned int i, Real &result, void *userData)
{
	// Diagonal element
	FluidModel *model = (FluidModel*)userData;

	const Vector3r &xi = model->getPosition(0, i);
	// Compute current fluid density for particle i
	Real density_i = model->getMass(i) * model->W_zero();

	//////////////////////////////////////////////////////////////////////////
	// Fluid
	//////////////////////////////////////////////////////////////////////////
	for (unsigned int j = 0; j < model->numberOfNeighbors(0, i); j++)
	{
		const unsigned int neighborIndex = model->getNeighbor(0, i, j);
		const Vector3r &xj = model->getPosition(0, neighborIndex);
		density_i += model->getMass(neighborIndex) * model->W(xi - xj);
	}

	result = density_i - model->getMass(i) * model->W_zero();
}

void Viscosity_Peer2016::step()
{
	const int numParticles = (int) m_model->numActiveParticles();
	const Real viscosity = 1.0 - m_viscosity;
	const Real density0 = m_model->getDensity0();

	const Real h = TimeManager::getCurrent()->getTimeStepSize();

	//////////////////////////////////////////////////////////////////////////
	// Init linear system solver and preconditioner
	//////////////////////////////////////////////////////////////////////////
	MatrixReplacement A(m_model->numActiveParticles(), matrixVecProdV, (void*)m_model);
	m_solverV.preconditioner().init(m_model->numActiveParticles(), diagonalMatrixElementV, (void*)m_model);
	m_solverV.setTolerance(m_maxError);
	m_solverV.setMaxIterations(m_maxIter);
	m_solverV.compute(A);

	MatrixReplacement A2(m_model->numActiveParticles(), matrixVecProdOmega, (void*)m_model);
	m_solverOmega.preconditioner().init(m_model->numActiveParticles(), diagonalMatrixElementOmega, (void*)m_model);
	m_solverOmega.setTolerance(m_maxError);
	m_solverOmega.setMaxIterations(m_maxIter);
	m_solverOmega.compute(A2);

	Eigen::VectorXd b0(numParticles);
	Eigen::VectorXd b1(numParticles);
	Eigen::VectorXd b2(numParticles);
	Eigen::VectorXd x0(numParticles);
	Eigen::VectorXd x1(numParticles);
	Eigen::VectorXd x2(numParticles);
	Eigen::VectorXd g0(numParticles);
	Eigen::VectorXd g1(numParticles);
	Eigen::VectorXd g2(numParticles);


	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static) nowait 
		for (int i = 0; i < numParticles; i++)
		{
			const Vector3r &xi = m_model->getPosition(0, i);
			const Vector3r &vi = m_model->getVelocity(0, i);
			const Real density_i = m_model->getDensity(i);

			//////////////////////////////////////////////////////////////////////////
			// compute nabla v
			//////////////////////////////////////////////////////////////////////////
			Matrix3r nablaV;
			nablaV.setZero();
			for (unsigned int j = 0; j < m_model->numberOfNeighbors(0, i); j++)
			{
				const unsigned int neighborIndex = m_model->getNeighbor(0, i, j);

				const Vector3r &xj = m_model->getPosition(0, neighborIndex);
				const Vector3r &vj = m_model->getVelocity(0, neighborIndex);
				const Vector3r gradW = m_model->gradW(xi - xj);

				Matrix3r dyad = (vj - vi) * gradW.transpose();

				nablaV += (1.0 / density_i) * m_model->getMass(neighborIndex) * dyad;
			}

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
			omega[0] = 2.0 * R(2, 1);
			omega[1] = 2.0 * R(0, 2);
			omega[2] = 2.0 * R(1, 0);

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
			const Vector3r &xi = m_model->getPosition(0, i);
			Vector3r rhs;
			rhs.setZero();

			for (unsigned int j = 0; j < m_model->numberOfNeighbors(0, i); j++)
			{
				const unsigned int neighborIndex = m_model->getNeighbor(0, i, j);

				const Real m = m_model->getMass(neighborIndex);
				const Vector3r &xj = m_model->getPosition(0, neighborIndex);
				const Vector3r xij = xi - xj;
				const Real W = m_model->W(xij);

				rhs += m *(m_omega[i] - m_omega[neighborIndex]) * W;
			}

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
	x0 = m_solverOmega.solveWithGuess(b0, g0);
	//x0 = m_solver.solve(b0);
	if (m_solverOmega.iterations() == 0)
		x0 = g0;
	x1 = m_solverOmega.solveWithGuess(b1, g1);
	//x1 = m_solver.solve(b1);
	if (m_solverOmega.iterations() == 0)
		x1 = g1;
	x2 = m_solverOmega.solveWithGuess(b2, g2);
	//x2 = m_solver.solve(b2);
	if (m_solverOmega.iterations() == 0)
		x2 = g2;
	STOP_TIMING_AVG;

	//////////////////////////////////////////////////////////////////////////
	// Determine new spin tensor and add it to target nabla v
	//////////////////////////////////////////////////////////////////////////
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static) nowait 
		for (int i = 0; i < (int)numParticles; i++)
		{
			Matrix3r R;
			R << 0.0, -0.5*x2[i], 0.5*x1[i],
				0.5*x2[i], 0.0, -0.5*x0[i],
				-0.5*x1[i], 0.5*x0[i], 0.0;

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
			const Vector3r &xi = m_model->getPosition(0, i);
			Vector3r rhs;
			rhs.setZero();

			for (unsigned int j = 0; j < m_model->numberOfNeighbors(0, i); j++)
			{
				const unsigned int neighborIndex = m_model->getNeighbor(0, i, j);

				const Real m = m_model->getMass(neighborIndex);
				const Vector3r &xj = m_model->getPosition(0, neighborIndex);
				const Vector3r xij = xi - xj;
				const Real W = m_model->W(xij);

				rhs += m * 0.5 * (getTargetNablaV(i) + getTargetNablaV(neighborIndex)) * xij * W;
			}

			const Vector3r &vi = m_model->getVelocity(0, i);
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
	int iter = 0;
	x0 = m_solverV.solveWithGuess(b0, g0);
	if (m_solverV.iterations() == 0)
		x0 = g0;
	iter += (int)m_solverV.iterations();
	x1 = m_solverV.solveWithGuess(b1, g1);
	if (m_solverV.iterations() == 0)
		x1 = g1;
	iter += (int)m_solverV.iterations();
	x2 = m_solverV.solveWithGuess(b2, g2);
	if (m_solverV.iterations() == 0)
		x2 = g2;
	iter += (int)m_solverV.iterations();
	STOP_TIMING_AVG;

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static) nowait 
		for (int i = 0; i < (int)numParticles; i++)
		{
			Vector3r &vi = m_model->getVelocity(0, i);
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
