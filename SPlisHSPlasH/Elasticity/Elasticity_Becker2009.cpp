#include "Elasticity_Becker2009.h"
#include "SPlisHSPlasH/Simulation.h"
#include "SPlisHSPlasH/Utilities/MathFunctions.h"

using namespace SPH;
using namespace GenParam;

std::string Elasticity_Becker2009::METHOD_NAME = "Becker et al. 2009";
int Elasticity_Becker2009::YOUNGS_MODULUS = -1;
int Elasticity_Becker2009::POISSON_RATIO = -1;
int Elasticity_Becker2009::FIXED_BOX_MIN = -1;
int Elasticity_Becker2009::FIXED_BOX_MAX = -1;
int Elasticity_Becker2009::ALPHA = -1;


Elasticity_Becker2009::Elasticity_Becker2009(FluidModel *model) :
	NonPressureForceBase(model)
{
	const unsigned int numParticles = model->numActiveParticles();
	m_restVolumes.resize(numParticles);
	m_current_to_initial_index.resize(numParticles);
	m_initial_to_current_index.resize(numParticles);
	m_initialNeighbors.resize(numParticles);
	m_rotations.resize(numParticles, Matrix3r::Identity());
	m_stress.resize(numParticles);
	m_F.resize(numParticles);

	m_youngsModulus = static_cast<Real>(100000.0);
	m_poissonRatio = static_cast<Real>(0.3);
	m_alpha = 0.0;
	m_fixedBoxMin.setZero();
	m_fixedBoxMax.setZero();

	model->addField({ "rest volume", METHOD_NAME, FieldType::Scalar, [&](const unsigned int i) -> Real* { return &m_restVolumes[i]; }, true });
	model->addField({ "rotation", METHOD_NAME, FieldType::Matrix3, [&](const unsigned int i) -> Real* { return &m_rotations[i](0,0); } });
	model->addField({ "stress", METHOD_NAME, FieldType::Vector6, [&](const unsigned int i) -> Real* { return &m_stress[i][0]; } });
	model->addField({ "deformation gradient", METHOD_NAME, FieldType::Matrix3, [&](const unsigned int i) -> Real* { return &m_F[i](0,0); } });
}

Elasticity_Becker2009::~Elasticity_Becker2009(void)
{
	m_model->removeFieldByName("rest volume");
	m_model->removeFieldByName("rotation");
	m_model->removeFieldByName("stress");
	m_model->removeFieldByName("deformation gradient");
}

void Elasticity_Becker2009::deferredInit()
{
	initValues();
}

void Elasticity_Becker2009::initParameters()
{
	NonPressureForceBase::initParameters();

	YOUNGS_MODULUS = createNumericParameter("youngsModulus", "Young`s modulus", &m_youngsModulus);
	setGroup(YOUNGS_MODULUS, "Fluid Model|Elasticity");
	setDescription(YOUNGS_MODULUS, "Stiffness of the elastic material");
	RealParameter* rparam = static_cast<RealParameter*>(getParameter(YOUNGS_MODULUS));
	rparam->setMinValue(0.0);

	POISSON_RATIO = createNumericParameter("poissonsRatio", "Poisson`s ratio", &m_poissonRatio);
	setGroup(POISSON_RATIO, "Fluid Model|Elasticity");
	setDescription(POISSON_RATIO, "Ratio of transversal expansion and axial compression");
	rparam = static_cast<RealParameter*>(getParameter(POISSON_RATIO));
	rparam->setMinValue(static_cast<Real>(-1.0 + 1e-4));
	rparam->setMaxValue(static_cast<Real>(0.5 - 1e-4));

	ParameterBase::GetVecFunc<Real> getFct = [&]()-> Real* { return m_fixedBoxMin.data(); };
	ParameterBase::SetVecFunc<Real> setFct = [&](Real* val)
	{
		m_fixedBoxMin = Vector3r(val[0], val[1], val[2]);
		determineFixedParticles();
	};
	FIXED_BOX_MIN = createVectorParameter("fixedBoxMin", "Fixed box min", 3u, getFct, setFct);
	setGroup(FIXED_BOX_MIN, "Fluid Model|Elasticity");
	setDescription(FIXED_BOX_MIN, "Minimum point of box of which contains the fixed particles.");
	getParameter(FIXED_BOX_MIN)->setReadOnly(true);


	ParameterBase::GetVecFunc<Real> getFct2 = [&]()-> Real* { return m_fixedBoxMax.data(); };
	ParameterBase::SetVecFunc<Real> setFct2 = [&](Real* val)
	{
		m_fixedBoxMax = Vector3r(val[0], val[1], val[2]);
		determineFixedParticles();
	};
	FIXED_BOX_MAX = createVectorParameter("fixedBoxMax", "Fixed box max", 3u, getFct2, setFct2);
	setGroup(FIXED_BOX_MAX, "Fluid Model|Elasticity");
	setDescription(FIXED_BOX_MAX, "Maximum point of box of which contains the fixed particles.");
	getParameter(FIXED_BOX_MAX)->setReadOnly(true);

	ALPHA = createNumericParameter("alpha", "Zero-energy modes suppression", &m_alpha);
	setGroup(ALPHA, "Fluid Model|Elasticity");
	setDescription(ALPHA, "Coefficent for zero-energy modes suppression method");
	rparam = static_cast<RealParameter*>(getParameter(ALPHA));
	rparam->setMinValue(0.0);
}

/** Mark all particles in the bounding box as fixed.
*/
void Elasticity_Becker2009::determineFixedParticles()
{
	const unsigned int numParticles = m_model->numActiveParticles();

	if (!m_fixedBoxMin.isZero() || !m_fixedBoxMax.isZero())
	{
		for (int i = 0; i < (int)numParticles; i++)
		{
			const Vector3r& x = m_model->getPosition0(i);
			if ((x[0] > m_fixedBoxMin[0]) && (x[1] > m_fixedBoxMin[1]) && (x[2] > m_fixedBoxMin[2]) &&
				(x[0] < m_fixedBoxMax[0]) && (x[1] < m_fixedBoxMax[1]) && (x[2] < m_fixedBoxMax[2]))
			{
				m_model->setParticleState(i, ParticleState::Fixed);
			}
		}
	}
}

void Elasticity_Becker2009::initValues()
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
			const unsigned int numNeighbors = sim->numberOfNeighbors(fluidModelIndex, fluidModelIndex, i);
			m_initialNeighbors[i].resize(numNeighbors);
			for (unsigned int j = 0; j < numNeighbors; j++)
				m_initialNeighbors[i][j] = sim->getNeighbor(fluidModelIndex, fluidModelIndex, i, j);

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
}

void Elasticity_Becker2009::step()
{
	computeRotations();
	computeStress();
	computeForces();
}


void Elasticity_Becker2009::reset()
{
	initValues();
}

void Elasticity_Becker2009::performNeighborhoodSearchSort()
{
	const unsigned int numPart = m_model->numActiveParticles();
	if (numPart == 0)
		return;

	Simulation *sim = Simulation::getCurrent();
	auto const& d = sim->getNeighborhoodSearch()->point_set(m_model->getPointSetIndex());
	d.sort_field(&m_restVolumes[0]);
	d.sort_field(&m_current_to_initial_index[0]);

	for (unsigned int i = 0; i < numPart; i++)
		m_initial_to_current_index[m_current_to_initial_index[i]] = i;
}

void Elasticity_Becker2009::computeRotations()
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
			Matrix3r Apq;
			Apq.setZero();

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
				const Vector3r xj_xi_0 = xj0 - xi0;
				Apq += m_model->getMass(neighborIndex) * sim->W(xj_xi_0) * (xj_xi * xj_xi_0.transpose());
			}

// 			Vector3r sigma;
// 			Matrix3r U, VT;
// 			MathFunctions::svdWithInversionHandling(Apq, sigma, U, VT);
// 			m_rotations[i] = U * VT;
			Quaternionr q(m_rotations[i]);
			MathFunctions::extractRotation(Apq, q, 10);
			m_rotations[i] = q.matrix();
		}
	}
}

void Elasticity_Becker2009::computeStress()
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int numParticles = m_model->numActiveParticles();
	const unsigned int fluidModelIndex = m_model->getPointSetIndex();
	FluidModel *model = m_model;

	Real mu = m_youngsModulus / (static_cast<Real>(2.0) * (static_cast<Real>(1.0) + m_poissonRatio));
	Real lambda = m_youngsModulus * m_poissonRatio / ((static_cast<Real>(1.0) + m_poissonRatio) * (static_cast<Real>(1.0) - static_cast<Real>(2.0) * m_poissonRatio));

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			const unsigned int i0 = m_current_to_initial_index[i];
			const Vector3r& xi = m_model->getPosition(i);
			const Vector3r& xi0 = m_model->getPosition0(i0);

			Matrix3r nablaU;
			nablaU.setZero();
			const size_t numNeighbors = m_initialNeighbors[i0].size();

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			for (unsigned int j = 0; j < numNeighbors; j++)
			{
				const unsigned int neighborIndex = m_initial_to_current_index[m_initialNeighbors[i0][j]];
				// get initial neighbor index considering the current particle order 
				const unsigned int neighborIndex0 = m_initialNeighbors[i0][j];

				const Vector3r& xj = model->getPosition(neighborIndex);
				const Vector3r& xj0 = m_model->getPosition0(neighborIndex0);

				const Vector3r xj_xi = xj - xi;
				const Vector3r xj_xi_0 = xj0 - xi0;

				const Vector3r uji = m_rotations[i].transpose() * xj_xi - xj_xi_0;
				// subtract because kernel gradient is taken in direction of xji0 instead of xij0
				nablaU -= (m_restVolumes[neighborIndex] * uji) * sim->gradW(xj_xi_0).transpose();
			}
			m_F[i] = nablaU + Matrix3r::Identity();

			// compute Cauchy strain: epsilon = 0.5 (nabla u + nabla u^T)
			Vector6r strain;
			strain[0] = nablaU(0, 0);						// \epsilon_{00}
			strain[1] = nablaU(1, 1);						// \epsilon_{11}
			strain[2] = nablaU(2, 2);						// \epsilon_{22}
			strain[3] = static_cast<Real>(0.5) * (nablaU(0, 1) + nablaU(1, 0)); // \epsilon_{01}
			strain[4] = static_cast<Real>(0.5) * (nablaU(0, 2) + nablaU(2, 0)); // \epsilon_{02}
			strain[5] = static_cast<Real>(0.5) * (nablaU(1, 2) + nablaU(2, 1)); // \epsilon_{12}

			//////////////////////////////////////////////////////////////////////////
			// First Piola Kirchhoff stress = 2 mu epsilon + lambda trace(epsilon) I
			//////////////////////////////////////////////////////////////////////////
			const Real trace = strain[0] + strain[1] + strain[2];
			const Real ltrace = lambda * trace;
			Vector6r& stress = m_stress[i];
			stress[0] = static_cast<Real>(2.0) * mu * strain[0] + ltrace;
			stress[1] = static_cast<Real>(2.0) * mu * strain[1] + ltrace;
			stress[2] = static_cast<Real>(2.0) * mu * strain[2] + ltrace;
			stress[3] = static_cast<Real>(2.0) * mu * strain[3];
			stress[4] = static_cast<Real>(2.0) * mu * strain[4];
			stress[5] = static_cast<Real>(2.0) * mu * strain[5];
		}
	}
}

void Elasticity_Becker2009::computeForces()
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
			if (model->getParticleState(i) == ParticleState::Active)
			{
				const unsigned int i0 = m_current_to_initial_index[i];
				const Vector3r& xi0 = m_model->getPosition0(i0);

				const size_t numNeighbors = m_initialNeighbors[i0].size();
				Vector3r fi;
				fi.setZero();

				//////////////////////////////////////////////////////////////////////////
				// Fluid
				//////////////////////////////////////////////////////////////////////////
				for (unsigned int j = 0; j < numNeighbors; j++)
				{
					const unsigned int neighborIndex = m_initial_to_current_index[m_initialNeighbors[i0][j]];
					// get initial neighbor index considering the current particle order 
					const unsigned int neighborIndex0 = m_initialNeighbors[i0][j];

					const Vector3r& xj0 = m_model->getPosition0(neighborIndex0);

					const Vector3r xj_xi_0 = xj0 - xi0;
					const Vector3r gradW0 = sim->gradW(xj_xi_0);

					const Vector3r dji = m_restVolumes[i] * gradW0;
					const Vector3r dij = -m_restVolumes[neighborIndex] * gradW0;

					Vector3r sdji, sdij;
					symMatTimesVec(m_stress[neighborIndex], dji, sdji);
					symMatTimesVec(m_stress[i], dij, sdij);

					const Vector3r fij = -m_restVolumes[neighborIndex] * sdji;
					const Vector3r fji = -m_restVolumes[i] * sdij;

					fi += m_rotations[neighborIndex] * fij - m_rotations[i] * fji;
				}
				fi = 0.5*fi;

				if (m_alpha != 0.0)
				{
					//////////////////////////////////////////////////////////////////////////
					// Ganzenmüller, G.C. 2015. An hourglass control algorithm for Lagrangian
					// Smooth Particle Hydrodynamics. Computer Methods in Applied Mechanics and 
					// Engineering 286, 87.106.
					//////////////////////////////////////////////////////////////////////////
					Vector3r fi_hg;
					fi_hg.setZero();
					const Vector3r& xi = m_model->getPosition(i);
					for (unsigned int j = 0; j < numNeighbors; j++)
					{
						const unsigned int neighborIndex = m_initial_to_current_index[m_initialNeighbors[i0][j]];
						// get initial neighbor index considering the current particle order 
						const unsigned int neighborIndex0 = m_initialNeighbors[i0][j];

						const Vector3r& xj = model->getPosition(neighborIndex);
						const Vector3r& xj0 = m_model->getPosition0(neighborIndex0);

						// Note: Ganzenm�ller defines xij = xj-xi
						const Vector3r xi_xj = -(xi - xj);
						const Real xixj_l = xi_xj.norm();
						if (xixj_l > 1.0e-6)
						{
							// Note: Ganzenm�ller defines xij = xj-xi
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

				// elastic acceleration
				Vector3r& ai = model->getAcceleration(i);
				ai += fi / model->getMass(i);
			}
		}
	}
}

void SPH::Elasticity_Becker2009::saveState(BinaryFileWriter &binWriter)
{
	binWriter.writeBuffer((char*)m_current_to_initial_index.data(), m_current_to_initial_index.size() * sizeof(unsigned int));
	binWriter.writeBuffer((char*)m_initial_to_current_index.data(), m_initial_to_current_index.size() * sizeof(unsigned int));
}

void SPH::Elasticity_Becker2009::loadState(BinaryFileReader &binReader)
{
	binReader.readBuffer((char*)m_current_to_initial_index.data(), m_current_to_initial_index.size() * sizeof(unsigned int));
	binReader.readBuffer((char*)m_initial_to_current_index.data(), m_initial_to_current_index.size() * sizeof(unsigned int));
}

