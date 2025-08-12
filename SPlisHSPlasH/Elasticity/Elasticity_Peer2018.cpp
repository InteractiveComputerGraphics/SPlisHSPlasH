#include "Elasticity_Peer2018.h"
#include "SPlisHSPlasH/Simulation.h"
#include "SPlisHSPlasH/Utilities/MathFunctions.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "Utilities/Timing.h"
#include "Utilities/Counting.h"

using namespace SPH;
using namespace GenParam;

int Elasticity_Peer2018::ITERATIONS = -1;
int Elasticity_Peer2018::MAX_ITERATIONS = -1;
int Elasticity_Peer2018::MAX_ERROR = -1;
int Elasticity_Peer2018::ALPHA = -1;
int Elasticity_Peer2018::MAX_NEIGHBORS = -1;


Elasticity_Peer2018::Elasticity_Peer2018(FluidModel *model) :
	ElasticityBase(model)
{
	const unsigned int numParticles = model->numActiveParticles();
	m_restVolumes.resize(numParticles);
	m_current_to_initial_index.resize(numParticles);
	m_initial_to_current_index.resize(numParticles);
	m_initialNeighbors.resize(numParticles);
	m_rotations.resize(numParticles, Matrix3r::Identity());
	m_stress.resize(numParticles);
	m_L.resize(numParticles);
	m_RL.resize(numParticles);
	m_F.resize(numParticles);

	m_iterations = 0;
	m_maxIter = 100;
	m_maxError = static_cast<Real>(1.0e-4);
	m_alpha = 0.0; 
	m_maxNeighbors = -1;

	model->addField({ "rest volume", FieldType::Scalar, [&](const unsigned int i) -> Real* { return &m_restVolumes[i]; }, true });
	model->addField({ "rotation", FieldType::Matrix3, [&](const unsigned int i) -> Real* { return &m_rotations[i](0,0); } });
	model->addField({ "stress", FieldType::Matrix3, [&](const unsigned int i) -> Real* { return &m_stress[i](0,0); } });
	model->addField({ "deformation gradient", FieldType::Matrix3, [&](const unsigned int i) -> Real* { return &m_F[i](0,0); } });
	model->addField({ "correction matrix", FieldType::Matrix3, [&](const unsigned int i) -> Real* { return &m_L[i](0,0); } });
}

Elasticity_Peer2018::~Elasticity_Peer2018(void)
{
	m_model->removeFieldByName("rest volume");
	m_model->removeFieldByName("rotation");
	m_model->removeFieldByName("stress");
	m_model->removeFieldByName("deformation gradient");
	m_model->removeFieldByName("correction matrix");
}

void Elasticity_Peer2018::deferredInit()
{
	initValues();
}

void Elasticity_Peer2018::initParameters()
{
	ElasticityBase::initParameters();

	ITERATIONS = createNumericParameter("elasticityIterations", "Iterations", &m_iterations);
	setGroup(ITERATIONS, "Fluid Model|Elasticity");
	setDescription(ITERATIONS, "Iterations required by the elasticity solver.");
	getParameter(ITERATIONS)->setReadOnly(true);

	MAX_ITERATIONS = createNumericParameter("elasticityMaxIter", "Max. iterations (elasticity)", &m_maxIter);
	setGroup(MAX_ITERATIONS, "Fluid Model|Elasticity");
	setDescription(MAX_ITERATIONS, "Coefficient for the elasticity force computation");
	static_cast<NumericParameter<unsigned int>*>(getParameter(MAX_ITERATIONS))->setMinValue(1);

	MAX_ERROR = createNumericParameter("elasticityMaxError", "Max. elasticity error", &m_maxError);
	setGroup(MAX_ERROR, "Fluid Model|Elasticity");
	setDescription(MAX_ERROR, "Coefficient for the elasticity force computation");
	RealParameter * rparam = static_cast<RealParameter*>(getParameter(MAX_ERROR));
	rparam->setMinValue(1e-7);

	ALPHA = createNumericParameter("alpha", "Zero-energy modes suppression", &m_alpha);
	setGroup(ALPHA, "Fluid Model|Elasticity");
	setDescription(ALPHA, "Coefficent for zero-energy modes suppression method");
	rparam = static_cast<RealParameter*>(getParameter(ALPHA));
	rparam->setMinValue(0.0);

	MAX_NEIGHBORS = createNumericParameter("maxNeighbors", "Max. neighbors", &m_maxNeighbors);
	setGroup(MAX_NEIGHBORS, "Fluid Model|Elasticity");
	setDescription(MAX_NEIGHBORS, "Maximum number of neighbors that are considered.");
}


void Elasticity_Peer2018::initValues()
{
	Simulation *sim = Simulation::getCurrent();
	sim->getNeighborhoodSearch()->find_neighbors();

	FluidModel *model = m_model;
	const unsigned int numParticles = model->numActiveParticles();
	const unsigned int fluidModelIndex = model->getPointSetIndex();

	// Store the neighbors in the reference configurations and
	// compute the volume of each particle in rest state
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static) 
		for (int i = 0; i < (int)numParticles; i++)
		{
			m_current_to_initial_index[i] = i;
			m_initial_to_current_index[i] = i;

			// reset particle state
			if (model->getParticleState(i) == ParticleState::Fixed)
				model->setParticleState(i, ParticleState::Active);

			// only neighbors in same phase will influence elasticity
			unsigned int numNeighbors = sim->numberOfNeighbors(fluidModelIndex, fluidModelIndex, i);
			m_initialNeighbors[i].resize(numNeighbors);
			for (unsigned int j = 0; j < numNeighbors; j++)
				m_initialNeighbors[i][j] = sim->getNeighbor(fluidModelIndex, fluidModelIndex, i, j);

			// if maxNeighbors is set, then sort all neighbors wrt. to their distance to xi 
			// and only take the maxNeighbors next ones.
			if (m_maxNeighbors > 0)
			{
				struct Comparator {
					Comparator(const Vector3r& xi, Vector3r* x) : m_xi(xi), m_x(x) {};
					bool operator()(unsigned int a, unsigned int b)
					{
						return (m_x[a] - m_xi).squaredNorm() < (m_x[b] - m_xi).squaredNorm();
					}

					Vector3r m_xi;
					Vector3r* m_x;
				};

				std::sort(m_initialNeighbors[i].begin(), m_initialNeighbors[i].end(), Comparator(model->getPosition0(i), &model->getPosition0(0)));
				if (m_initialNeighbors[i].size() > m_maxNeighbors)
				{
					numNeighbors = m_maxNeighbors;
					m_initialNeighbors[i].resize(m_maxNeighbors);
				}
			}

			// compute volume
			Real density = model->getMass(i) * sim->W_zero();
			const Vector3r &xi = model->getPosition(i);
			for (size_t j = 0; j < m_initialNeighbors[i].size(); j++)
			{
				const unsigned int neighborIndex = m_initialNeighbors[i][j];
				const Vector3r& xj = model->getPosition(neighborIndex); 
				density += model->getMass(neighborIndex) * sim->W(xi - xj);
			}
			m_restVolumes[i] = model->getMass(i) / density;
			m_rotations[i].setIdentity();
		}
	}

	// mark all particles in the bounding box as fixed
	determineFixedParticles();

	computeMatrixL();
}


void Elasticity_Peer2018::step()
{
	START_TIMING("Elasticity_Peer2018")
	const unsigned int numParticles = m_model->numActiveParticles();
	if (numParticles == 0)
		return;
	const Real dt = TimeManager::getCurrent()->getTimeStepSize();

	//////////////////////////////////////////////////////////////////////////
	// Init linear system solver and preconditioner
	//////////////////////////////////////////////////////////////////////////
	MatrixReplacement A(3 * m_model->numActiveParticles(), matrixVecProd, (void*)this);
	//m_solver.preconditioner().init(m_model->numActiveParticles(), diagonalMatrixElement, (void*)this);

	m_solver.setTolerance(m_maxError);
	m_solver.setMaxIterations(m_maxIter);
	m_solver.compute(A);

	VectorXr b(3 * numParticles);
	VectorXr x(3 * numParticles);
	VectorXr g(3 * numParticles);

	computeRotations();
	computeRHS(b);

	// warmstart
	#pragma omp parallel for schedule(static) 
	for (int i = 0; i < (int)numParticles; i++)
	{
		if (m_model->getParticleState(i) == ParticleState::Active)
			g.segment<3>(3 * i) = m_model->getVelocity(i) + dt * m_model->getAcceleration(i);
		else
			g.segment<3>(3 * i) = m_model->getVelocity(i);
	}

	//////////////////////////////////////////////////////////////////////////
	// Solve linear system 
	//////////////////////////////////////////////////////////////////////////
	START_TIMING("Elasticity - CG solve");
	x = m_solver.solveWithGuess(b, g);
	m_iterations = (int)m_solver.iterations();
	STOP_TIMING_AVG;
	INCREASE_COUNTER("Elasticity - CG iterations", static_cast<Real>(m_iterations));

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			if (m_model->getParticleState(i) == ParticleState::Active)
			{
				Vector3r& ai = m_model->getAcceleration(i);
				ai += (1.0 / dt) * (x.segment<3>(3 * i) - m_model->getVelocity(i));
			}
		}
	}
	STOP_TIMING_AVG
}


void Elasticity_Peer2018::reset()
{
	initValues();
}

void Elasticity_Peer2018::performNeighborhoodSearchSort()
{
	const unsigned int numPart = m_model->numActiveParticles();
	if (numPart == 0)
		return;

	Simulation *sim = Simulation::getCurrent();
	auto const& d = sim->getNeighborhoodSearch()->point_set(m_model->getPointSetIndex());
	d.sort_field(&m_restVolumes[0]);
	d.sort_field(&m_rotations[0]);
	d.sort_field(&m_current_to_initial_index[0]);
	d.sort_field(&m_L[0]);

	for (unsigned int i = 0; i < numPart; i++)
		m_initial_to_current_index[m_current_to_initial_index[i]] = i;
}

void Elasticity_Peer2018::computeRotations()
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int numParticles = m_model->numActiveParticles();
	const unsigned int fluidModelIndex = m_model->getPointSetIndex();
	FluidModel *model = m_model;

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			const unsigned int i0 = m_current_to_initial_index[i];
			const Vector3r &xi = m_model->getPosition(i);
			const Vector3r &xi0 = m_model->getPosition0(i0);
			Matrix3r F;
			F.setZero();

			const size_t numNeighbors = m_initialNeighbors[i0].size();

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			for (unsigned int j = 0; j < numNeighbors; j++)
			{
				const unsigned int neighborIndex = m_initial_to_current_index[m_initialNeighbors[i0][j]];
				// get initial neighbor index considering the current particle order 
				const unsigned int neighborIndex0 = m_initialNeighbors[i0][j];

				const Vector3r &xj = model->getPosition(neighborIndex);
				const Vector3r &xj0 = m_model->getPosition0(neighborIndex0);
				const Vector3r xj_xi = xj - xi;
				const Vector3r xi_xj_0 = xi0 - xj0;
				const Vector3r correctedKernel = m_L[i] * sim->gradW(xi_xj_0);
				F += m_restVolumes[neighborIndex] * xj_xi * correctedKernel.transpose();
			}

			if (sim->is2DSimulation())
				F(2, 2) = 1.0;

//  			Vector3r sigma; 
//  			Matrix3r U, VT;
//  			MathFunctions::svdWithInversionHandling(F, sigma, U, VT);
//  			m_rotations[i] = U * VT;
 			Quaternionr q(m_rotations[i]);
 			MathFunctions::extractRotation(F, q, 10);
 			m_rotations[i] = q.matrix();

			m_RL[i] = m_rotations[i] * m_L[i];
		}
	}
}

void Elasticity_Peer2018::computeMatrixL()
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int numParticles = m_model->numActiveParticles();
	const unsigned int fluidModelIndex = m_model->getPointSetIndex();

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			const unsigned int i0 = m_current_to_initial_index[i];
			const Vector3r &xi0 = m_model->getPosition0(i0);
			Matrix3r L;
			L.setZero();

			const size_t numNeighbors = m_initialNeighbors[i0].size();

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			for (unsigned int j = 0; j < numNeighbors; j++)
			{
				const unsigned int neighborIndex = m_initial_to_current_index[m_initialNeighbors[i0][j]];
				// get initial neighbor index considering the current particle order 
				const unsigned int neighborIndex0 = m_initialNeighbors[i0][j];

				const Vector3r &xj0 = m_model->getPosition0(neighborIndex0);
				const Vector3r xj_xi_0 = xj0 - xi0;
				const Vector3r gradW = sim->gradW(xj_xi_0);

				// minus because gradW(xij0) == -gradW(xji0)
				L -= m_restVolumes[neighborIndex] * gradW * xj_xi_0.transpose();
			}

			// add 1 to z-component. otherwise we get a singular matrix in 2D
			if (sim->is2DSimulation())
				L(2, 2) = 1.0;

			bool invertible = false;
			L.computeInverseWithCheck(m_L[i], invertible, 1e-9);
			if (!invertible)
			{
				//MathFunctions::pseudoInverse(L, m_L[i]);
				m_L[i] = Matrix3r::Identity();
			}
		}
	}
}

#ifdef USE_AVX

void Elasticity_Peer2018::computeRHS(VectorXr & rhs)
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int numParticles = m_model->numActiveParticles();
	const unsigned int fluidModelIndex = m_model->getPointSetIndex();
	FluidModel *model = m_model;
	const Real dt = TimeManager::getCurrent()->getTimeStepSize();

	Real mu = m_youngsModulus / (static_cast<Real>(2.0) * (static_cast<Real>(1.0) + m_poissonRatio));
	Real lambda = m_youngsModulus * m_poissonRatio / ((static_cast<Real>(1.0) + m_poissonRatio) * (static_cast<Real>(1.0) - static_cast<Real>(2.0) * m_poissonRatio));

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			if (model->getParticleState(i) != ParticleState::Active)
				continue;

			const unsigned int i0 = m_current_to_initial_index[i];
			const Vector3r &xi = m_model->getPosition(i);
			const Vector3r &xi0 = m_model->getPosition0(i0);
			const unsigned int numNeighbors = m_initialNeighbors[i0].size();

 			//////////////////////////////////////////////////////////////////////////
 			// compute corotated deformation gradient (Eq. 18)
 			//////////////////////////////////////////////////////////////////////////
			Matrix3f8 F_avx;
			F_avx.setZero();
			const Vector3f8 xi_avx(xi);
			const Vector3f8 xi0_avx(xi0);
			const Matrix3f8 Ri(m_rotations[i]);
			const Matrix3f8 RLi(m_RL[i]);
 
  			//////////////////////////////////////////////////////////////////////////
 			// Fluid
 			//////////////////////////////////////////////////////////////////////////
			for (unsigned int j = 0; j < numNeighbors; j += 8)
 			{
				const unsigned int count = std::min(numNeighbors - j, 8u);
				std::array<unsigned int, 8> indices;
				generateIndices(m_initial_to_current_index.data(), &m_initialNeighbors[i0][j], indices, count);
 
				const Vector3f8 xj0_avx = convertVec_zero(&m_initialNeighbors[i0][j], &model->getPosition0(0), count);
				const Vector3f8 xj_avx = convertVec_zero(&indices[0], &model->getPosition(0), count);
				const Scalarf8 Vj0_avx = convert_zero(&indices[0], &m_restVolumes[0], count);

				const Vector3f8 xj_xi = xj_avx - xi_avx;
				const Vector3f8 xi_xj_0 = xi0_avx - xj0_avx;
				const Vector3f8 correctedRotatedKernel = RLi * CubicKernel_AVX::gradW(xi_xj_0);

				Matrix3f8 dyad;
				dyadicProduct((xj_xi - Ri * (xj0_avx - xi0_avx)), correctedRotatedKernel, dyad);
				F_avx += dyad * Vj0_avx;
 			}

			m_F[i] = F_avx.reduce();
			m_F[i] += Matrix3r::Identity();

			if (sim->is2DSimulation())
				m_F[i](2, 2) = 1.0;

 			//////////////////////////////////////////////////////////////////////////
 			// compute Cauchy strain: epsilon = 0.5 (F + F^T) - I
 			//////////////////////////////////////////////////////////////////////////
 			Vector6r strain;
 			strain[0] = m_F[i](0, 0) - static_cast<Real>(1.0);						// \epsilon_{00}
 			strain[1] = m_F[i](1, 1) - static_cast<Real>(1.0);						// \epsilon_{11}
 			strain[2] = m_F[i](2, 2) - static_cast<Real>(1.0);						// \epsilon_{22}
 			strain[3] = static_cast<Real>(0.5) * (m_F[i](0, 1) + m_F[i](1, 0));			// \epsilon_{01}
 			strain[4] = static_cast<Real>(0.5) * (m_F[i](0, 2) + m_F[i](2, 0));			// \epsilon_{02}
 			strain[5] = static_cast<Real>(0.5) * (m_F[i](1, 2) + m_F[i](2, 1));			// \epsilon_{12}

			//////////////////////////////////////////////////////////////////////////
			// First Piola Kirchhoff stress = 2 mu epsilon + lambda trace(epsilon) I
			//////////////////////////////////////////////////////////////////////////
			const Real trace = strain[0] + strain[1] + strain[2];
			const Real ltrace = lambda*trace;
			Matrix3r& stress = m_stress[i];
			stress(0, 0) = static_cast<Real>(2.0) * mu * strain[0] + ltrace;
			stress(1, 1) = static_cast<Real>(2.0) * mu * strain[1] + ltrace;
			stress(2, 2) = static_cast<Real>(2.0) * mu * strain[2] + ltrace;

			stress(0, 1) = static_cast<Real>(2.0) * mu * strain[3];
			stress(1, 0) = static_cast<Real>(2.0) * mu * strain[3];
			stress(0, 2) = static_cast<Real>(2.0) * mu * strain[4];
			stress(2, 0) = static_cast<Real>(2.0) * mu * strain[4];
			stress(1, 2) = static_cast<Real>(2.0) * mu * strain[5];
			stress(2, 1) = static_cast<Real>(2.0) * mu * strain[5];
		}
	}

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			if (model->getParticleState(i) != ParticleState::Active)
				continue;

			const unsigned int i0 = m_current_to_initial_index[i];
			const Vector3r &xi0 = m_model->getPosition0(i0);
			const unsigned int numNeighbors = m_initialNeighbors[i0].size();

			//////////////////////////////////////////////////////////////////////////
			// Compute elastic force
			//////////////////////////////////////////////////////////////////////////
			Vector3f8 force_avx;
			force_avx.setZero();
			const Scalarf8 Vi0_avx(m_restVolumes[i]);
			const Vector3f8 xi0_avx(xi0);
			const Matrix3f8 RLi(m_RL[i]);
			const Matrix3f8 stress_i(m_stress[i]);

			for (unsigned int j = 0; j < numNeighbors; j += 8)
			{
				const unsigned int count = std::min(numNeighbors - j, 8u);
				std::array<unsigned int, 8> indices;
				generateIndices(m_initial_to_current_index.data(), &m_initialNeighbors[i0][j], indices, count);

				const Matrix3f8& RLj = convertMat_zero(&indices[0], &m_RL[0], count);
				const Vector3f8 xj0_avx = convertVec_zero(&m_initialNeighbors[i0][j], &model->getPosition0(0), count);
				const Scalarf8 Vj0_avx = convert_zero(&indices[0], &m_restVolumes[0], count);
				const Vector3f8 xi_xj_0 = xi0_avx - xj0_avx;
				const Vector3f8 gradW = CubicKernel_AVX::gradW(xi_xj_0);
				const Vector3f8 correctedRotatedKernel_i = RLi * gradW;
				const Vector3f8 correctedRotatedKernel_j = RLj * gradW;

				const Matrix3f8& stress_j = convertMat_zero(&indices[0], &m_stress[0], count);
				Vector3f8 PWi = stress_i * correctedRotatedKernel_i;
				Vector3f8 PWj = stress_j * correctedRotatedKernel_j;
				force_avx += (PWi + PWj) * Vi0_avx * Vj0_avx;
			}
			Vector3r force = force_avx.reduce();

			if (m_alpha != 0.0)
			{
				//////////////////////////////////////////////////////////////////////////
				// Ganzenmüller, G.C. 2015. An hourglass control algorithm for Lagrangian
				// Smooth Particle Hydrodynamics. Computer Methods in Applied Mechanics and 
				// Engineering 286, 87.106.
				//////////////////////////////////////////////////////////////////////////
				Vector3r fi_hg;
				fi_hg.setZero();
				const Vector3r &xi = m_model->getPosition(i);
				for (unsigned int j = 0; j < numNeighbors; j++)
				{
					const unsigned int neighborIndex = m_initial_to_current_index[m_initialNeighbors[i0][j]];
					// get initial neighbor index considering the current particle order 
					const unsigned int neighborIndex0 = m_initialNeighbors[i0][j];

					const Vector3r &xj = model->getPosition(neighborIndex);
					const Vector3r &xj0 = m_model->getPosition0(neighborIndex0);

					// Note: Ganzenmueller defines xij = xj-xi
					const Vector3r xi_xj = -(xi - xj);
					const Real xixj_l = xi_xj.norm();
					if (xixj_l > 1.0e-6)
					{
						// Note: Ganzenmueller defines xij = xj-xi
						const Vector3r xi_xj_0 = -(xi0 - xj0);
						const Real xixj0_l2 = xi_xj_0.squaredNorm();
						const Real W0 = sim->W(xi_xj_0);

						const Vector3r xij_i = m_F[i] * m_rotations[i] * xi_xj_0;
						const Vector3r xji_j = -m_F[neighborIndex] * m_rotations[neighborIndex] * xi_xj_0;
						const Vector3r epsilon_ij_i = xij_i - xi_xj;
						const Vector3r epsilon_ji_j = xji_j + xi_xj;

						const Real delta_ij_i = epsilon_ij_i.dot(xi_xj) / xixj_l;
						const Real delta_ji_j = -epsilon_ji_j.dot(xi_xj) / xixj_l;

						fi_hg -= m_restVolumes[neighborIndex] * W0 / xixj0_l2 * (delta_ij_i + delta_ji_j) * xi_xj / xixj_l;
					}
				}
				fi_hg *= m_alpha * m_youngsModulus * m_restVolumes[i];
				model->getAcceleration(i) += fi_hg / model->getMass(i);
			}
			if (model->getParticleState(i) == ParticleState::Active)
				rhs.segment<3>(3 * i) = model->getVelocity(i) + dt * (model->getAcceleration(i) + 1.0 / model->getMass(i) * force);
			else
				rhs.segment<3>(3 * i) = model->getVelocity(i);
		}
	}
}

void Elasticity_Peer2018::matrixVecProd(const Real* vec, Real *result, void *userData)
{
	Simulation *sim = Simulation::getCurrent();
	Elasticity_Peer2018 * elasticity = static_cast<Elasticity_Peer2018*>(userData);
	FluidModel *model = elasticity->getModel();
	const unsigned int numParticles = model->numActiveParticles();
	const Real dt = TimeManager::getCurrent()->getTimeStepSize();

	const auto &current_to_initial_index = elasticity->m_current_to_initial_index;
	const auto &initial_to_current_index = elasticity->m_initial_to_current_index;
	const auto &initialNeighbors = elasticity->m_initialNeighbors;
	const auto &restVolumes = elasticity->m_restVolumes;
	const auto &rotations = elasticity->m_rotations;
	const auto &L = elasticity->m_L;
	const auto &RL = elasticity->m_RL;
	auto &stress = elasticity->m_stress;
	const Real youngsModulus = elasticity->m_youngsModulus;
	const Real poissonRatio = elasticity->m_poissonRatio;

	Real mu = youngsModulus / (static_cast<Real>(2.0) * (static_cast<Real>(1.0) + poissonRatio));
	Real lambda = youngsModulus * poissonRatio / ((static_cast<Real>(1.0) + poissonRatio) * (static_cast<Real>(1.0) - static_cast<Real>(2.0) * poissonRatio));

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			const unsigned int i0 = current_to_initial_index[i];
			const Vector3r &pi = Eigen::Map<const Vector3r>(&vec[3 * i], 3);
			const Vector3r &xi0 = model->getPosition0(i0);
			const unsigned int numNeighbors = initialNeighbors[i0].size();

 			//////////////////////////////////////////////////////////////////////////
 			// compute corotated deformation gradient (Eq. 18)
 			//////////////////////////////////////////////////////////////////////////
			Matrix3f8 nablaU_avx;
			nablaU_avx.setZero();
			const Vector3f8 pi_avx(pi);
			const Vector3f8 xi0_avx(xi0);
			const Matrix3f8 RLi(RL[i]);

  			//////////////////////////////////////////////////////////////////////////
 			// Fluid
 			//////////////////////////////////////////////////////////////////////////
			for (unsigned int j = 0; j < numNeighbors; j += 8)
 			{
				const unsigned int count = std::min(numNeighbors - j, 8u);
				std::array<unsigned int, 8> indices;
				elasticity->generateIndices(initial_to_current_index.data(), &initialNeighbors[i0][j], indices, count);

				const Vector3f8 xj0_avx = convertVec_zero(&initialNeighbors[i0][j], &model->getPosition0(0), count);
				const Vector3f8 pj_avx = convertVec_zero(&indices[0], &vec[0], count);
				const Scalarf8 Vj0_avx = convert_zero(&indices[0], &restVolumes[0], count);
				const Vector3f8 pj_pi = pj_avx - pi_avx;
				const Vector3f8 xi_xj_0 = xi0_avx - xj0_avx;
				const Vector3f8 correctedRotatedKernel = RLi * CubicKernel_AVX::gradW(xi_xj_0);

				Matrix3f8 dyad;
				dyadicProduct(pj_pi, correctedRotatedKernel, dyad);
				nablaU_avx += dyad * Vj0_avx;
 			}
			Matrix3r nablaU = nablaU_avx.reduce();
			nablaU *= dt;
 			//////////////////////////////////////////////////////////////////////////
 			// compute Cauchy strain: epsilon = 0.5 (nablaU + nablaU^T)
 			//////////////////////////////////////////////////////////////////////////
 			Vector6r strain;
			strain[0] = nablaU(0, 0);									// \epsilon_{00}
			strain[1] = nablaU(1, 1);									// \epsilon_{11}
			strain[2] = nablaU(2, 2);									// \epsilon_{22}
 			strain[3] = static_cast<Real>(0.5) * (nablaU(0, 1) + nablaU(1, 0));			// \epsilon_{01}
 			strain[4] = static_cast<Real>(0.5) * (nablaU(0, 2) + nablaU(2, 0));			// \epsilon_{02}
 			strain[5] = static_cast<Real>(0.5) * (nablaU(1, 2) + nablaU(2, 1));			// \epsilon_{12}

			//////////////////////////////////////////////////////////////////////////
			// First Piola Kirchhoff stress = 2 mu epsilon + lambda trace(epsilon) I
			//////////////////////////////////////////////////////////////////////////
			const Real trace = strain[0] + strain[1] + strain[2];
			const Real ltrace = lambda*trace;
			stress[i](0, 0) = static_cast<Real>(2.0) * mu * strain[0] + ltrace;
			stress[i](1, 1) = static_cast<Real>(2.0) * mu * strain[1] + ltrace;
			stress[i](2, 2) = static_cast<Real>(2.0) * mu * strain[2] + ltrace;

			stress[i](0, 1) = static_cast<Real>(2.0) * mu * strain[3];
			stress[i](1, 0) = static_cast<Real>(2.0) * mu * strain[3];
			stress[i](0, 2) = static_cast<Real>(2.0) * mu * strain[4];
			stress[i](2, 0) = static_cast<Real>(2.0) * mu * strain[4];
			stress[i](1, 2) = static_cast<Real>(2.0) * mu * strain[5];
			stress[i](2, 1) = static_cast<Real>(2.0) * mu * strain[5];
		}
	}

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			if (model->getParticleState(i) == ParticleState::Active)
			{
				const unsigned int i0 = current_to_initial_index[i];
				const Vector3r& xi0 = model->getPosition0(i0);
				const unsigned int numNeighbors = initialNeighbors[i0].size();

				//////////////////////////////////////////////////////////////////////////
				// Compute elastic force
				//////////////////////////////////////////////////////////////////////////
				Vector3f8 force_avx;
				force_avx.setZero();
				const Scalarf8 Vi0_avx(restVolumes[i]);
				const Vector3f8 xi0_avx(xi0);
				const Matrix3f8 RLi(RL[i]);
				const Matrix3f8 stress_i(stress[i]);
				for (unsigned int j = 0; j < numNeighbors; j += 8)
				{
					const unsigned int count = std::min(numNeighbors - j, 8u);
					std::array<unsigned int, 8> indices;
					elasticity->generateIndices(initial_to_current_index.data(), &initialNeighbors[i0][j], indices, count);

					const Matrix3f8& RLj = convertMat_zero(&indices[0], &RL[0], count);
					const Vector3f8 xj0_avx = convertVec_zero(&initialNeighbors[i0][j], &model->getPosition0(0), count);
					const Scalarf8 Vj0_avx = convert_zero(&indices[0], &restVolumes[0], count);
					const Vector3f8 xi_xj_0 = xi0_avx - xj0_avx;
					const Vector3f8 gradW = CubicKernel_AVX::gradW(xi_xj_0);
					const Vector3f8 correctedRotatedKernel_i = RLi * gradW;
					const Vector3f8 correctedRotatedKernel_j = RLj * gradW;
					
					const Matrix3f8& stress_j = convertMat_zero(&indices[0], &stress[0], count);
					Vector3f8 PWi = stress_i * correctedRotatedKernel_i;
					Vector3f8 PWj = stress_j * correctedRotatedKernel_j;
					force_avx += (PWi + PWj) * Vi0_avx * Vj0_avx;
				}

				const Vector3r force = force_avx.reduce();
				const Real factor = dt / model->getMass(i);
				result[3 * i] = vec[3 * i] - factor * force[0];
				result[3 * i + 1] = vec[3 * i + 1] - factor * force[1];
				result[3 * i + 2] = vec[3 * i + 2] - factor * force[2];
			}
			else
			{
				result[3 * i] = vec[3 * i];
				result[3 * i + 1] = vec[3 * i + 1];
				result[3 * i + 2] = vec[3 * i + 2];
			}
		}
	}
}


#else


void Elasticity_Peer2018::computeRHS(VectorXr & rhs)
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int numParticles = m_model->numActiveParticles();
	const unsigned int fluidModelIndex = m_model->getPointSetIndex();
	FluidModel *model = m_model;
	const Real dt = TimeManager::getCurrent()->getTimeStepSize();

	Real mu = m_youngsModulus / (static_cast<Real>(2.0) * (static_cast<Real>(1.0) + m_poissonRatio));
	Real lambda = m_youngsModulus * m_poissonRatio / ((static_cast<Real>(1.0) + m_poissonRatio) * (static_cast<Real>(1.0) - static_cast<Real>(2.0) * m_poissonRatio));

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			const unsigned int i0 = m_current_to_initial_index[i];
			const Vector3r &xi = m_model->getPosition(i);
			const Vector3r &xi0 = m_model->getPosition0(i0);
			const size_t numNeighbors = m_initialNeighbors[i0].size();

 			//////////////////////////////////////////////////////////////////////////
 			// compute corotated deformation gradient (Eq. 18)
 			//////////////////////////////////////////////////////////////////////////
			m_F[i].setZero();
 
  			//////////////////////////////////////////////////////////////////////////
 			// Fluid
 			//////////////////////////////////////////////////////////////////////////
 			for (unsigned int j = 0; j < numNeighbors; j++)
 			{
 				const unsigned int neighborIndex = m_initial_to_current_index[m_initialNeighbors[i0][j]];
 				// get initial neighbor index considering the current particle order 
 				const unsigned int neighborIndex0 = m_initialNeighbors[i0][j];
 
 				const Vector3r &xj = model->getPosition(neighborIndex);
 				const Vector3r &xj0 = m_model->getPosition0(neighborIndex0);
 				const Vector3r xj_xi = xj - xi;
 				const Vector3r xi_xj_0 = xi0 - xj0;
 				const Vector3r correctedRotatedKernel = m_RL[i] * sim->gradW(xi_xj_0);
 				m_F[i] += m_restVolumes[neighborIndex] * (xj_xi - m_rotations[i]*(xj0-xi0)) * correctedRotatedKernel.transpose();
 			}

			m_F[i] += Matrix3r::Identity();

			if (sim->is2DSimulation())
				m_F[i](2, 2) = 1.0;

 			//////////////////////////////////////////////////////////////////////////
 			// compute Cauchy strain: epsilon = 0.5 (F + F^T) - I
 			//////////////////////////////////////////////////////////////////////////
 			Vector6r strain;
 			strain[0] = m_F[i](0, 0) - static_cast<Real>(1.0);						// \epsilon_{00}
 			strain[1] = m_F[i](1, 1) - static_cast<Real>(1.0);						// \epsilon_{11}
 			strain[2] = m_F[i](2, 2) - static_cast<Real>(1.0);						// \epsilon_{22}
 			strain[3] = static_cast<Real>(0.5) * (m_F[i](0, 1) + m_F[i](1, 0));			// \epsilon_{01}
 			strain[4] = static_cast<Real>(0.5) * (m_F[i](0, 2) + m_F[i](2, 0));			// \epsilon_{02}
 			strain[5] = static_cast<Real>(0.5) * (m_F[i](1, 2) + m_F[i](2, 1));			// \epsilon_{12}

			//////////////////////////////////////////////////////////////////////////
			// First Piola Kirchhoff stress = 2 mu epsilon + lambda trace(epsilon) I
			//////////////////////////////////////////////////////////////////////////
			const Real trace = strain[0] + strain[1] + strain[2];
			const Real ltrace = lambda*trace;
			Matrix3r& stress = m_stress[i];
			stress(0, 0) = static_cast<Real>(2.0) * mu * strain[0] + ltrace;
			stress(1, 1) = static_cast<Real>(2.0) * mu * strain[1] + ltrace;
			stress(2, 2) = static_cast<Real>(2.0) * mu * strain[2] + ltrace;

			stress(0, 1) = static_cast<Real>(2.0) * mu * strain[3];
			stress(1, 0) = static_cast<Real>(2.0) * mu * strain[3];
			stress(0, 2) = static_cast<Real>(2.0) * mu * strain[4];
			stress(2, 0) = static_cast<Real>(2.0) * mu * strain[4];
			stress(1, 2) = static_cast<Real>(2.0) * mu * strain[5];
			stress(2, 1) = static_cast<Real>(2.0) * mu * strain[5];
		}
	}

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			const unsigned int i0 = m_current_to_initial_index[i];
			const Vector3r &xi0 = m_model->getPosition0(i0);
			const size_t numNeighbors = m_initialNeighbors[i0].size();

			//////////////////////////////////////////////////////////////////////////
			// Compute elastic force
			//////////////////////////////////////////////////////////////////////////
			Vector3r force;
			force.setZero();
			for (unsigned int j = 0; j < numNeighbors; j++)
			{
				const unsigned int neighborIndex = m_initial_to_current_index[m_initialNeighbors[i0][j]];
				// get initial neighbor index considering the current particle order 
				const unsigned int neighborIndex0 = m_initialNeighbors[i0][j];

				const Vector3r &xj0 = m_model->getPosition0(neighborIndex0);
				const Vector3r xi_xj_0 = xi0 - xj0;
				const Vector3r correctedRotatedKernel_i = m_RL[i] * sim->gradW(xi_xj_0);
				const Vector3r correctedRotatedKernel_j = -m_RL[neighborIndex] * sim->gradW(xi_xj_0);
				const Vector3r PWi = m_stress[i] * correctedRotatedKernel_i;
				const Vector3r PWj = m_stress[neighborIndex] * correctedRotatedKernel_j;
				force += m_restVolumes[i] * m_restVolumes[neighborIndex] * (PWi - PWj);
			}

			if (m_alpha != 0.0)
			{
				//////////////////////////////////////////////////////////////////////////
				// Ganzenmüller, G.C. 2015. An hourglass control algorithm for Lagrangian
				// Smooth Particle Hydrodynamics. Computer Methods in Applied Mechanics and 
				// Engineering 286, 87.106.
				//////////////////////////////////////////////////////////////////////////
				Vector3r fi_hg;
				fi_hg.setZero();
				const Vector3r &xi = m_model->getPosition(i);
				for (unsigned int j = 0; j < numNeighbors; j++)
				{
					const unsigned int neighborIndex = m_initial_to_current_index[m_initialNeighbors[i0][j]];
					// get initial neighbor index considering the current particle order 
					const unsigned int neighborIndex0 = m_initialNeighbors[i0][j];

					const Vector3r &xj = model->getPosition(neighborIndex);
					const Vector3r &xj0 = m_model->getPosition0(neighborIndex0);

					// Note: Ganzenmueller defines xij = xj-xi
					const Vector3r xi_xj = -(xi - xj);
					const Real xixj_l = xi_xj.norm();
					if (xixj_l > 1.0e-6)
					{
						// Note: Ganzenmueller defines xij = xj-xi
						const Vector3r xi_xj_0 = -(xi0 - xj0);
						const Real xixj0_l2 = xi_xj_0.squaredNorm();
						const Real W0 = sim->W(xi_xj_0);

						const Vector3r xij_i = m_F[i] * m_rotations[i] * xi_xj_0;
						const Vector3r xji_j = -m_F[neighborIndex] * m_rotations[neighborIndex] * xi_xj_0;
						const Vector3r epsilon_ij_i = xij_i - xi_xj;
						const Vector3r epsilon_ji_j = xji_j + xi_xj;

						const Real delta_ij_i = epsilon_ij_i.dot(xi_xj) / xixj_l;
						const Real delta_ji_j = -epsilon_ji_j.dot(xi_xj) / xixj_l;

						fi_hg -= m_restVolumes[neighborIndex] * W0 / xixj0_l2 * (delta_ij_i + delta_ji_j) * xi_xj / xixj_l;
					}
				}
				fi_hg *= m_alpha * m_youngsModulus * m_restVolumes[i];
				model->getAcceleration(i) += fi_hg / model->getMass(i);
			}
			if (model->getParticleState(i) == ParticleState::Active)
				rhs.segment<3>(3 * i) = model->getVelocity(i) + dt * (model->getAcceleration(i) + 1.0 / model->getMass(i) * force);
			else
				rhs.segment<3>(3 * i) = model->getVelocity(i);
		}
	}
}

void Elasticity_Peer2018::matrixVecProd(const Real* vec, Real *result, void *userData)
{
	Simulation *sim = Simulation::getCurrent();
	Elasticity_Peer2018 * elasticity = static_cast<Elasticity_Peer2018*>(userData);
	FluidModel *model = elasticity->getModel();
	const unsigned int numParticles = model->numActiveParticles();
	const Real dt = TimeManager::getCurrent()->getTimeStepSize();

	const auto &current_to_initial_index = elasticity->m_current_to_initial_index;
	const auto &initial_to_current_index = elasticity->m_initial_to_current_index;
	const auto &initialNeighbors = elasticity->m_initialNeighbors;
	const auto &restVolumes = elasticity->m_restVolumes;
	const auto &rotations = elasticity->m_rotations;
	const auto &L = elasticity->m_L;
	const auto &RL = elasticity->m_RL;
	auto &stress = elasticity->m_stress;
	const Real youngsModulus = elasticity->m_youngsModulus;
	const Real poissonRatio = elasticity->m_poissonRatio;

	Real mu = youngsModulus / (static_cast<Real>(2.0) * (static_cast<Real>(1.0) + poissonRatio));
	Real lambda = youngsModulus * poissonRatio / ((static_cast<Real>(1.0) + poissonRatio) * (static_cast<Real>(1.0) - static_cast<Real>(2.0) * poissonRatio));

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			const unsigned int i0 = current_to_initial_index[i];
			const Vector3r &pi = Eigen::Map<const Vector3r>(&vec[3 * i], 3);
			const Vector3r &xi0 = model->getPosition0(i0);
			const size_t numNeighbors = initialNeighbors[i0].size();

 			//////////////////////////////////////////////////////////////////////////
 			// compute corotated deformation gradient (Eq. 18)
 			//////////////////////////////////////////////////////////////////////////
 			Matrix3r nablaU;
			nablaU.setZero();

  			//////////////////////////////////////////////////////////////////////////
 			// Fluid
 			//////////////////////////////////////////////////////////////////////////
 			for (unsigned int j = 0; j < numNeighbors; j++)
 			{
 				const unsigned int neighborIndex = initial_to_current_index[initialNeighbors[i0][j]];
 				// get initial neighbor index considering the current particle order 
 				const unsigned int neighborIndex0 = initialNeighbors[i0][j];
 
 				const Vector3r &pj = Eigen::Map<const Vector3r>(&vec[3 * neighborIndex], 3);
 				const Vector3r &xj0 = model->getPosition0(neighborIndex0);
 				const Vector3r pj_pi = pj - pi;
 				const Vector3r xi_xj_0 = xi0 - xj0;
 				const Vector3r correctedRotatedKernel = RL[i] * sim->gradW(xi_xj_0);
				nablaU += restVolumes[neighborIndex] * pj_pi * correctedRotatedKernel.transpose();
 			}
			nablaU *= dt;
 			//////////////////////////////////////////////////////////////////////////
 			// compute Cauchy strain: epsilon = 0.5 (nablaU + nablaU^T)
 			//////////////////////////////////////////////////////////////////////////
 			Vector6r strain;
			strain[0] = nablaU(0, 0);									// \epsilon_{00}
			strain[1] = nablaU(1, 1);									// \epsilon_{11}
			strain[2] = nablaU(2, 2);									// \epsilon_{22}
 			strain[3] = static_cast<Real>(0.5) * (nablaU(0, 1) + nablaU(1, 0));			// \epsilon_{01}
 			strain[4] = static_cast<Real>(0.5) * (nablaU(0, 2) + nablaU(2, 0));			// \epsilon_{02}
 			strain[5] = static_cast<Real>(0.5) * (nablaU(1, 2) + nablaU(2, 1));			// \epsilon_{12}

			//////////////////////////////////////////////////////////////////////////
			// First Piola Kirchhoff stress = 2 mu epsilon + lambda trace(epsilon) I
			//////////////////////////////////////////////////////////////////////////
			const Real trace = strain[0] + strain[1] + strain[2];
			const Real ltrace = lambda*trace;
			stress[i](0, 0) = static_cast<Real>(2.0) * mu * strain[0] + ltrace;
			stress[i](1, 1) = static_cast<Real>(2.0) * mu * strain[1] + ltrace;
			stress[i](2, 2) = static_cast<Real>(2.0) * mu * strain[2] + ltrace;

			stress[i](0, 1) = static_cast<Real>(2.0) * mu * strain[3];
			stress[i](1, 0) = static_cast<Real>(2.0) * mu * strain[3];
			stress[i](0, 2) = static_cast<Real>(2.0) * mu * strain[4];
			stress[i](2, 0) = static_cast<Real>(2.0) * mu * strain[4];
			stress[i](1, 2) = static_cast<Real>(2.0) * mu * strain[5];
			stress[i](2, 1) = static_cast<Real>(2.0) * mu * strain[5];
		}
	}

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			if (model->getParticleState(i) == ParticleState::Active)
			{
				const unsigned int i0 = current_to_initial_index[i];
				const Vector3r& xi0 = model->getPosition0(i0);
				const size_t numNeighbors = initialNeighbors[i0].size();

				//////////////////////////////////////////////////////////////////////////
				// Compute elastic force
				//////////////////////////////////////////////////////////////////////////
				Vector3r force;
				force.setZero();
				for (unsigned int j = 0; j < numNeighbors; j++)
				{
					const unsigned int neighborIndex = initial_to_current_index[initialNeighbors[i0][j]];
					// get initial neighbor index considering the current particle order 
					const unsigned int neighborIndex0 = initialNeighbors[i0][j];

					const Vector3r& xj0 = model->getPosition0(neighborIndex0);
					const Vector3r xi_xj_0 = xi0 - xj0;
					const Vector3r gradW = sim->gradW(xi_xj_0);
					const Vector3r correctedRotatedKernel_i = RL[i] * gradW;
					const Vector3r correctedRotatedKernel_j = -RL[neighborIndex] * gradW;
					const Vector3r PWi = stress[i] * correctedRotatedKernel_i;
					const Vector3r PWj = stress[neighborIndex] * correctedRotatedKernel_j;
					force += restVolumes[i] * restVolumes[neighborIndex] * (PWi - PWj);
				}

				const Real factor = dt / model->getMass(i);
				result[3 * i] = vec[3 * i] - factor * force[0];
				result[3 * i + 1] = vec[3 * i + 1] - factor * force[1];
				result[3 * i + 2] = vec[3 * i + 2] - factor * force[2];
			}
			else
			{
				result[3 * i] = vec[3 * i];
				result[3 * i + 1] = vec[3 * i + 1];
				result[3 * i + 2] = vec[3 * i + 2];
			}
		}
	}
}

#endif





void SPH::Elasticity_Peer2018::saveState(BinaryFileWriter &binWriter)
{
	binWriter.writeBuffer((char*)m_current_to_initial_index.data(), m_current_to_initial_index.size() * sizeof(unsigned int));
	binWriter.writeBuffer((char*)m_initial_to_current_index.data(), m_initial_to_current_index.size() * sizeof(unsigned int));
	binWriter.writeBuffer((char*)m_L.data(), m_L.size() * sizeof(Matrix3r));
}

void SPH::Elasticity_Peer2018::loadState(BinaryFileReader &binReader)
{
	binReader.readBuffer((char*)m_current_to_initial_index.data(), m_current_to_initial_index.size() * sizeof(unsigned int));
	binReader.readBuffer((char*)m_initial_to_current_index.data(), m_initial_to_current_index.size() * sizeof(unsigned int));
	binReader.readBuffer((char*)m_L.data(), m_L.size() * sizeof(Matrix3r));
}
