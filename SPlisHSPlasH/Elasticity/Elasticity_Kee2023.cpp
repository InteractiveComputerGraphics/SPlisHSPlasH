#include "Elasticity_Kee2023.h"
#include "SPlisHSPlasH/Simulation.h"
#include "SPlisHSPlasH/Utilities/MathFunctions.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "Utilities/Timing.h"
#include "Utilities/Counting.h"
#include <extern/md5/md5.h>
#include "Utilities/BinaryFileReaderWriter.h"
#include "Utilities/StringTools.h"
#include <Utilities/FileSystem.h>
#include <array>
#include <chrono>

using namespace SPH;
using namespace GenParam;

std::string Elasticity_Kee2023::METHOD_NAME = "Kee et al. 2023";
int Elasticity_Kee2023::YOUNGS_MODULUS = -1;
int Elasticity_Kee2023::POISSON_RATIO = -1;
int Elasticity_Kee2023::FIXED_BOX_MIN = -1;
int Elasticity_Kee2023::FIXED_BOX_MAX = -1;
int Elasticity_Kee2023::FIXED_BOX2_MIN = -1;
int Elasticity_Kee2023::FIXED_BOX2_MAX = -1;
int Elasticity_Kee2023::ALPHA = -1;
int Elasticity_Kee2023::MAX_NEIGHBORS = -1;
int Elasticity_Kee2023::SOLVER_TYPE = -1;
int Elasticity_Kee2023::LBFGS_WINDOW_SIZE = -1;
int Elasticity_Kee2023::MATERIAL_TYPE = -1;
int Elasticity_Kee2023::MAX_ITER = -1;
int Elasticity_Kee2023::MAX_ERROR = -1;
int Elasticity_Kee2023::MAX_ITER_CG = -1;
int Elasticity_Kee2023::TOL_CG = -1;
int Elasticity_Kee2023::MAX_LS_ITER = -1;
int Elasticity_Kee2023::LS_ARMIJO_PARAM = -1;
int Elasticity_Kee2023::LS_BETA = -1;
int Elasticity_Kee2023::USE_LINE_SEARCH = -1;
int Elasticity_Kee2023::ENUM_SOLVER_NEWTON = -1;
int Elasticity_Kee2023::ENUM_SOLVER_LBFGS = -1;
int Elasticity_Kee2023::ENUM_MATERIAL_STABLE_NEOHOOKEAN = -1;
int Elasticity_Kee2023::ENUM_MATERIAL_COROTATED = -1;


Elasticity_Kee2023::Elasticity_Kee2023(FluidModel *model) :
	NonPressureForceBase(model)
{
	const unsigned int numParticles = model->numActiveParticles();
	m_restVolumes.resize(numParticles);
	m_current_to_initial_index.resize(numParticles);
	m_initial_to_current_index.resize(numParticles);
	m_initialNeighbors.resize(numParticles);
	m_rotations.resize(numParticles, Matrix3r::Identity());
	m_stress.resize(numParticles);
	m_fixedGroupId.resize(numParticles, 0);
	m_L.resize(numParticles);								// kernel gradient correction matrix L
	m_F.resize(numParticles);								// deformation gradient
	m_PL.resize(numParticles);								// stores the rotation matrix times the matrix L

	m_youngsModulus = static_cast<Real>(5000000.0);
	m_poissonRatio = static_cast<Real>(0.45);
	
	m_alpha = 0.0;
	m_maxNeighbors = -1;
	m_fixedBoxMin.setZero();
	m_fixedBoxMax.setZero();
	m_fixedBox2Min.setZero();
	m_fixedBox2Max.setZero();
	m_solverType = 1;
	m_lbfgsWindowSize = 5;
	m_materialType = 1;  // Co-rotated
	m_maxIter = 100;
	m_maxError = static_cast<Real>(1e-6);
	m_maxIterCG = 100;
	m_tolCG = static_cast<Real>(1e-4);
	m_maxLSIter = 20;
	m_lsArmijoParam = static_cast<Real>(1e-4);
	m_lsBeta = static_cast<Real>(0.5);
	m_useLineSearch = true;

	model->addField({ "rest volume", METHOD_NAME, FieldType::Scalar, [&](const unsigned int i) -> Real* { return &m_restVolumes[i]; }, true });
	model->addField({ "rotation", METHOD_NAME, FieldType::Matrix3, [&](const unsigned int i) -> Real* { return &m_rotations[i](0,0); } });
	model->addField({ "stress", METHOD_NAME, FieldType::Scalar, [&](const unsigned int i) -> Real* { return &m_stress[i]; } });
	model->addField({ "deformation gradient", METHOD_NAME, FieldType::Matrix3, [&](const unsigned int i) -> Real* { return &m_F[i](0,0); } });
	model->addField({ "correction matrix", METHOD_NAME, FieldType::Matrix3, [&](const unsigned int i) -> Real* { return &m_L[i](0,0); } });
}

Elasticity_Kee2023::~Elasticity_Kee2023(void)
{
	m_model->removeFieldByName("rest volume");
	m_model->removeFieldByName("rotation");
	m_model->removeFieldByName("stress");
	m_model->removeFieldByName("deformation gradient");
	m_model->removeFieldByName("correction matrix");

	for (auto objIndex = 0; objIndex < m_objects.size(); objIndex++)
	{
		delete m_objects[objIndex];
	}
	m_objects.clear();
}

void Elasticity_Kee2023::deferredInit()
{
	initValues();
}

void Elasticity_Kee2023::initParameters()
{
	NonPressureForceBase::initParameters();

	ParameterBase::GetFunc<Real> getFctYM = [&]()-> Real { return m_youngsModulus; };
	ParameterBase::SetFunc<Real> setFctYM = [&](Real val) 
	{ 
		m_youngsModulus = val; 
		m_mu = m_youngsModulus / (static_cast<Real>(2.0) * (static_cast<Real>(1.0) + m_poissonRatio));
		m_lambda = m_youngsModulus * m_poissonRatio / ((static_cast<Real>(1.0) + m_poissonRatio) * (static_cast<Real>(1.0) - static_cast<Real>(2.0) * m_poissonRatio));

		if (Simulation::getCurrent()->isSimulationInitialized())		// if Young's modulus has changed, recompute the factorization
			Simulation::getCurrent()->reset();
	};
	YOUNGS_MODULUS = createNumericParameter("youngsModulus", "Young`s modulus", getFctYM, setFctYM);
	setGroup(YOUNGS_MODULUS, "Fluid Model|Elasticity");
	setDescription(YOUNGS_MODULUS, "Stiffness of the elastic material");
	RealParameter* rparam = static_cast<RealParameter*>(getParameter(YOUNGS_MODULUS));
	rparam->setMinValue(0.0);

	ParameterBase::GetFunc<Real> getFctPR = [&]()-> Real { return m_poissonRatio; };
	ParameterBase::SetFunc<Real> setFctPR = [&](Real val)
	{
		m_poissonRatio = val;
		m_mu = m_youngsModulus / (static_cast<Real>(2.0) * (static_cast<Real>(1.0) + m_poissonRatio));
		m_lambda = m_youngsModulus * m_poissonRatio / ((static_cast<Real>(1.0) + m_poissonRatio) * (static_cast<Real>(1.0) - static_cast<Real>(2.0) * m_poissonRatio));

		if (Simulation::getCurrent()->isSimulationInitialized())		// if Poisson ration has changed, recompute the factorization
			Simulation::getCurrent()->reset();
	};
	POISSON_RATIO = createNumericParameter("poissonsRatio", "Poisson`s ratio", getFctPR, setFctPR);
	setGroup(POISSON_RATIO, "Fluid Model|Elasticity");
	setDescription(POISSON_RATIO, "Ratio of transversal expansion and axial compression");
	rparam = static_cast<RealParameter*>(getParameter(POISSON_RATIO));
	rparam->setMinValue(static_cast<Real>(-1.0 + 1e-4));
	rparam->setMaxValue(static_cast<Real>(0.5 - 1e-4));

	ParameterBase::GetFunc<Real> getFctAlpha = [&]()-> Real { return m_alpha; };
	ParameterBase::SetFunc<Real> setFctAlpha = [&](Real val)
	{
		m_alpha = val;
		if (Simulation::getCurrent()->isSimulationInitialized())		// if value has changed, recompute the factorization
			Simulation::getCurrent()->reset();
	};
	ALPHA = createNumericParameter("alpha", "Zero-energy modes suppression", getFctAlpha, setFctAlpha);
	setGroup(ALPHA, "Fluid Model|Elasticity");
	setDescription(ALPHA, "Coefficent for zero-energy modes suppression method");
	rparam = static_cast<RealParameter*>(getParameter(ALPHA));
	rparam->setMinValue(0.0);

	ParameterBase::GetFunc<int> getFct5 = [&]()-> int { return m_maxNeighbors; };
	ParameterBase::SetFunc<int> setFct5 = [&](int val)
	{
		m_maxNeighbors = val;
		if (Simulation::getCurrent()->isSimulationInitialized())		// if value has changed, recompute the factorization
			Simulation::getCurrent()->reset();
	};
	MAX_NEIGHBORS = createNumericParameter("maxNeighbors", "Max. neighbors", getFct5, setFct5);
	setGroup(MAX_NEIGHBORS, "Fluid Model|Elasticity");
	setDescription(MAX_NEIGHBORS, "Maximum number of neighbors that are considered.");

	ParameterBase::GetVecFunc<Real> getFctFMin = [&]()-> Real* { return m_fixedBoxMin.data(); };
	ParameterBase::SetVecFunc<Real> setFctFMin = [&](Real* val)
	{
		m_fixedBoxMin = Vector3r(val[0], val[1], val[2]);
		determineFixedParticles();
	};
	FIXED_BOX_MIN = createVectorParameter("fixedBoxMin", "Fixed box min", 3u, getFctFMin, setFctFMin);
	setGroup(FIXED_BOX_MIN, "Fluid Model|Elasticity");
	setDescription(FIXED_BOX_MIN, "Minimum point of box of which contains the fixed particles.");
	getParameter(FIXED_BOX_MIN)->setReadOnly(true);


	ParameterBase::GetVecFunc<Real> getFctFMax = [&]()-> Real* { return m_fixedBoxMax.data(); };
	ParameterBase::SetVecFunc<Real> setFctFMax = [&](Real* val)
	{
		m_fixedBoxMax = Vector3r(val[0], val[1], val[2]);
		determineFixedParticles();
	};
	FIXED_BOX_MAX = createVectorParameter("fixedBoxMax", "Fixed box max", 3u, getFctFMax, setFctFMax);
	setGroup(FIXED_BOX_MAX, "Fluid Model|Elasticity");
	setDescription(FIXED_BOX_MAX, "Maximum point of box of which contains the fixed particles.");
	getParameter(FIXED_BOX_MAX)->setReadOnly(true);

	ParameterBase::GetVecFunc<Real> getFctF2Min = [&]()-> Real* { return m_fixedBox2Min.data(); };
	ParameterBase::SetVecFunc<Real> setFctF2Min = [&](Real* val)
	{
		m_fixedBox2Min = Vector3r(val[0], val[1], val[2]);
		determineFixedParticles();
	};
	FIXED_BOX2_MIN = createVectorParameter("fixedBox2Min", "Fixed box 2 min", 3u, getFctF2Min, setFctF2Min);
	setGroup(FIXED_BOX2_MIN, "Fluid Model|Elasticity");
	setDescription(FIXED_BOX2_MIN, "Minimum point of second box containing fixed particles.");
	getParameter(FIXED_BOX2_MIN)->setReadOnly(true);

	ParameterBase::GetVecFunc<Real> getFctF2Max = [&]()-> Real* { return m_fixedBox2Max.data(); };
	ParameterBase::SetVecFunc<Real> setFctF2Max = [&](Real* val)
	{
		m_fixedBox2Max = Vector3r(val[0], val[1], val[2]);
		determineFixedParticles();
	};
	FIXED_BOX2_MAX = createVectorParameter("fixedBox2Max", "Fixed box 2 max", 3u, getFctF2Max, setFctF2Max);
	setGroup(FIXED_BOX2_MAX, "Fluid Model|Elasticity");
	setDescription(FIXED_BOX2_MAX, "Maximum point of second box containing fixed particles.");
	getParameter(FIXED_BOX2_MAX)->setReadOnly(true);

	SOLVER_TYPE = createEnumParameter("solverType", "Solver type", &m_solverType);
	setGroup(SOLVER_TYPE, "Fluid Model|Elasticity");
	setDescription(SOLVER_TYPE, "Solver used for the elasticity computation.");
	EnumParameter* enumParam = static_cast<EnumParameter*>(getParameter(SOLVER_TYPE));
	enumParam->addEnumValue("Newton", ENUM_SOLVER_NEWTON);
	enumParam->addEnumValue("L-BFGS", ENUM_SOLVER_LBFGS);

	LBFGS_WINDOW_SIZE = createNumericParameter("lbfgsWindowSize", "L-BFGS window size", &m_lbfgsWindowSize);
	setGroup(LBFGS_WINDOW_SIZE, "Fluid Model|Elasticity");
	setDescription(LBFGS_WINDOW_SIZE, "Number of past iterations stored for L-BFGS approximation (only for L-BFGS).");
	static_cast<NumericParameter<int>*>(getParameter(LBFGS_WINDOW_SIZE))->setMinValue(0);

	MAX_ITER = createNumericParameter("maxIterations", "Max. iterations", &m_maxIter);
	setGroup(MAX_ITER, "Fluid Model|Elasticity");
	setDescription(MAX_ITER, "Maximum number of solver iterations per time step.");
	static_cast<NumericParameter<int>*>(getParameter(MAX_ITER))->setMinValue(1);

	MAX_ERROR = createNumericParameter("maxError", "Max. error", &m_maxError);
	setGroup(MAX_ERROR, "Fluid Model|Elasticity");
	setDescription(MAX_ERROR, "Convergence threshold on infinity norm of position update.");
	rparam = static_cast<RealParameter*>(getParameter(MAX_ERROR));
	rparam->setMinValue(static_cast<Real>(1e-15));

	MAX_ITER_CG = createNumericParameter("maxIterationsCG", "Max. CG iterations", &m_maxIterCG);
	setGroup(MAX_ITER_CG, "Fluid Model|Elasticity");
	setDescription(MAX_ITER_CG, "Maximum number of CG iterations for Newton linear solve.");
	static_cast<NumericParameter<int>*>(getParameter(MAX_ITER_CG))->setMinValue(1);

	TOL_CG = createNumericParameter("tolCG", "CG tolerance", &m_tolCG);
	setGroup(TOL_CG, "Fluid Model|Elasticity");
	setDescription(TOL_CG, "Convergence tolerance for CG solver.");
	rparam = static_cast<RealParameter*>(getParameter(TOL_CG));
	rparam->setMinValue(static_cast<Real>(1e-15));

	MAX_LS_ITER = createNumericParameter("maxLSIterations", "Max. LS iterations", &m_maxLSIter);
	setGroup(MAX_LS_ITER, "Fluid Model|Elasticity");
	setDescription(MAX_LS_ITER, "Maximum number of backtracking line search iterations.");
	static_cast<NumericParameter<int>*>(getParameter(MAX_LS_ITER))->setMinValue(1);

	LS_ARMIJO_PARAM = createNumericParameter("lsArmijoParam", "LS Armijo c1", &m_lsArmijoParam);
	setGroup(LS_ARMIJO_PARAM, "Fluid Model|Elasticity");
	setDescription(LS_ARMIJO_PARAM, "Armijo sufficient decrease parameter for line search.");
	rparam = static_cast<RealParameter*>(getParameter(LS_ARMIJO_PARAM));
	rparam->setMinValue(static_cast<Real>(1e-10));
	rparam->setMaxValue(static_cast<Real>(0.5));

	LS_BETA = createNumericParameter("lsBeta", "LS backtracking factor", &m_lsBeta);
	setGroup(LS_BETA, "Fluid Model|Elasticity");
	setDescription(LS_BETA, "Step size reduction factor per line search iteration.");
	rparam = static_cast<RealParameter*>(getParameter(LS_BETA));
	rparam->setMinValue(static_cast<Real>(0.01));
	rparam->setMaxValue(static_cast<Real>(0.99));

	USE_LINE_SEARCH = createBoolParameter("useLineSearch", "Use line search", &m_useLineSearch);
	setGroup(USE_LINE_SEARCH, "Fluid Model|Elasticity");
	setDescription(USE_LINE_SEARCH, "Enable backtracking Armijo line search.");

	MATERIAL_TYPE = createEnumParameter("materialType", "Material type", &m_materialType);
	setGroup(MATERIAL_TYPE, "Fluid Model|Elasticity");
	setDescription(MATERIAL_TYPE, "Constitutive model for the elastic material.");
	enumParam = static_cast<EnumParameter*>(getParameter(MATERIAL_TYPE));
	enumParam->addEnumValue("Stable Neo-Hookean", ENUM_MATERIAL_STABLE_NEOHOOKEAN);
	enumParam->addEnumValue("Co-rotated", ENUM_MATERIAL_COROTATED);
}

/** Mark all particles in the bounding box as fixed.
*/
void Elasticity_Kee2023::determineFixedParticles()
{
	const unsigned int numParticles = m_model->numActiveParticles();

	auto inBox = [](const Vector3r& x, const Vector3r& bmin, const Vector3r& bmax) -> bool
	{
		return (x[0] > bmin[0]) && (x[1] > bmin[1]) && (x[2] > bmin[2]) &&
			   (x[0] < bmax[0]) && (x[1] < bmax[1]) && (x[2] < bmax[2]);
	};

	const bool hasBox1 = !m_fixedBoxMin.isZero() || !m_fixedBoxMax.isZero();
	const bool hasBox2 = !m_fixedBox2Min.isZero() || !m_fixedBox2Max.isZero();

	if (hasBox1 || hasBox2)
	{
		for (int i = 0; i < (int)numParticles; i++)
		{
			const Vector3r& x = m_model->getPosition0(i);
			if (hasBox1 && inBox(x, m_fixedBoxMin, m_fixedBoxMax))
			{
				m_model->setParticleState(i, ParticleState::Fixed);
				m_fixedGroupId[i] = 1;
			}
			else if (hasBox2 && inBox(x, m_fixedBox2Min, m_fixedBox2Max))
			{
				m_model->setParticleState(i, ParticleState::Fixed);
				m_fixedGroupId[i] = 2;
			}
		}
	}
}

/** Compute an MD4 check sum using the neighborhood structure in order to 
* recognize known particle models (cache).
*/
std::string Elasticity_Kee2023::computeMD5(const unsigned int objIndex)
{
	ElasticObject* obj = m_objects[objIndex];

	auto& group = obj->m_particleIndices;
	auto numParticles = group.size();
	std::vector<unsigned int> tempN;
	tempN.resize(obj->m_particleIndices.size() * 2);
	for (size_t i = 0; i < numParticles; i++)
	{
		const unsigned int particleIndex = group[i];
		const size_t numNeighbors = m_initialNeighbors[particleIndex].size();
		tempN[2 * i] = static_cast<unsigned int>(numNeighbors);
		tempN[2 * i + 1] = 0;
		for (auto j = 0; j < numNeighbors; j++)
			tempN[2 * i + 1] += m_initialNeighbors[particleIndex][j] - group[0];
	}

	// compute MD5 checksum for all particle positions, this is used for the cache file
	MD5 context((unsigned char*)&tempN[0], static_cast<unsigned int>(2 * numParticles * sizeof(unsigned int)));
	char* md5hex = context.hex_digest();

	tempN.clear();
	return std::string(md5hex);
}

/** Initialize the particle neighborhoods in the reference configuration.
* Fix particles which lie in the user-defined bounding box. 
* Find out if there are multiple separate objects in the phase. 
* Finally, compute kernel gradient correction matrices and factorization.
*/
void Elasticity_Kee2023::initValues()
{
	Simulation *sim = Simulation::getCurrent();
	sim->getNeighborhoodSearch()->find_neighbors();

	FluidModel *model = m_model;
	const unsigned int numParticles = model->numActiveParticles();
	const unsigned int fluidModelIndex = model->getPointSetIndex();
	m_totalNeighbors = 0u;

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
			if (m_model->getParticleState(i) == ParticleState::Fixed)
				m_model->setParticleState(i, ParticleState::Active);

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

				// sort the neighbors according to their distance
				std::sort(m_initialNeighbors[i].begin(), m_initialNeighbors[i].end(), Comparator(model->getPosition0(i), &model->getPosition0(0)));

				// take only the next maxNeighbors
				if (m_initialNeighbors[i].size() > m_maxNeighbors)
				{
					numNeighbors = m_maxNeighbors;
					m_initialNeighbors[i].resize(m_maxNeighbors);
				}
			}
			m_totalNeighbors += numNeighbors;

			m_rotations[i].setIdentity();
			m_stress[i] = 0.0;
			m_F[i].setIdentity();
			m_PL[i].setIdentity();
		}
	}

	// Symmetrize neighbor lists: if j∈N(i) but i∉N(j), add i to N(j)
	if (m_maxNeighbors > 0)
	{
		for (unsigned int i = 0; i < numParticles; i++)
		{
			for (size_t jn = 0; jn < m_initialNeighbors[i].size(); jn++)
			{
				unsigned int j = m_initialNeighbors[i][jn];
				bool found = false;
				for (size_t kn = 0; kn < m_initialNeighbors[j].size(); kn++)
				{
					if (m_initialNeighbors[j][kn] == i) { found = true; break; }
				}
				if (!found)
					m_initialNeighbors[j].push_back(i);
			}
		}
	}

	// Check neighbor list symmetry
	{
		unsigned int maxActualNeighbors = 0;
		int asymmetricPairs = 0;
		for (unsigned int i = 0; i < numParticles; i++)
		{
			if (m_initialNeighbors[i].size() > maxActualNeighbors)
				maxActualNeighbors = (unsigned int)m_initialNeighbors[i].size();
			for (size_t jn = 0; jn < m_initialNeighbors[i].size(); jn++)
			{
				unsigned int j = m_initialNeighbors[i][jn];
				bool found = false;
				for (size_t kn = 0; kn < m_initialNeighbors[j].size(); kn++)
				{
					if (m_initialNeighbors[j][kn] == i) { found = true; break; }
				}
				if (!found) asymmetricPairs++;
			}
		}
		LOG_INFO << "Neighbor stats: maxNeighbors=" << maxActualNeighbors
			<< ", asymmetricPairs=" << asymmetricPairs;
	}

	// Compute rest volumes using final (symmetrized) neighbor lists
	for (unsigned int i = 0; i < numParticles; i++)
	{
		Real density = model->getMass(i) * sim->W_zero();
		const Vector3r &xi0 = model->getPosition0(i);
		for (size_t j = 0; j < m_initialNeighbors[i].size(); j++)
		{
			const unsigned int neighborIndex0 = m_initialNeighbors[i][j];
			const Vector3r& xj0 = model->getPosition0(neighborIndex0);
			density += model->getMass(neighborIndex0) * sim->W(xi0 - xj0);
		}
		m_restVolumes[i] = model->getMass(i) / density;
	}

	// mark all particles in the bounding box as fixed
	determineFixedParticles();

	// find separate objects
	START_TIMING("findObjects")
	findObjects();
	STOP_TIMING_AVG;

	// if we find the same object, copy the neighborhood info
	size_t numObjects = m_objects.size();
	auto& fluidInfos = sim->getFluidInfos();
	for (auto objIndex = 1; objIndex < numObjects; objIndex++)
	{
		bool foundSameObj = false;
		int objIndex2;
		for (objIndex2 = 0; objIndex2 < objIndex; objIndex2++)
		{
			if (fluidInfos[objIndex].hasSameParticleSampling(fluidInfos[objIndex2]))
			{
				foundSameObj = true;
				break;
			}
		}
		if (foundSameObj)
		{
			ElasticObject* obj = m_objects[objIndex];
			ElasticObject* obj0 = m_objects[objIndex2];
			const std::vector<unsigned int>& group = obj->m_particleIndices;
			const std::vector<unsigned int>& group0 = obj0->m_particleIndices;
			int numParticles = (int)group.size();
			int offset = group[0];
			int offset0 = group0[0];

			for (int i = 0; i < (int)numParticles; i++)
			{
				int particleIndex = group[i];
				int particleIndex0 = group0[i];

				const unsigned int i0 = m_current_to_initial_index[particleIndex];
				const unsigned int i00 = m_current_to_initial_index[particleIndex0];
				const size_t numNeighbors = m_initialNeighbors[i00].size();
				m_initialNeighbors[i0].resize(numNeighbors);
				for (int j = 0; j < numNeighbors; j++)
				{
					m_initialNeighbors[i0][j] = m_initialNeighbors[i00][j] - offset0 + offset;
				}
			}
		}
	}

	// compute kernel gradient correction matrix
	START_TIMING("computeMatrixL")
	computeMatrixL();
	STOP_TIMING_AVG;

	// init factorization
	START_TIMING("initSystem")
	initSystem();
	STOP_TIMING_AVG;
}

/** Find separate objects by object id. */
void Elasticity_Kee2023::findObjects()
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int numParticles = m_model->numActiveParticles();
	const unsigned int fluidModelIndex = m_model->getPointSetIndex();
	FluidModel *model = m_model;
	m_objects.clear();
	std::map<unsigned int, unsigned int> obj2group;

	for (unsigned int i = 0; i < numParticles; i++)
	{
		const unsigned int objId = m_model->getObjectId(i);

		// search for object id in map
		if (obj2group.find(objId) == obj2group.end())
		{
			// if id was not found, then we have a new object
			m_objects.push_back(new ElasticObject());
			const unsigned int groupIndex = (unsigned int)(m_objects.size()-1);
			obj2group[objId] = groupIndex;

			const unsigned int i0 = m_current_to_initial_index[i];
			m_objects[groupIndex]->m_particleIndices.push_back(i0);
		}
		else
		{
			// object already exists
			const unsigned int groupIndex = obj2group[objId];
			const unsigned int i0 = m_current_to_initial_index[i];
			m_objects[groupIndex]->m_particleIndices.push_back(i0);
		}
	}

	// For each object sort the particles so that all fixed particles are at the end of the list.
	// This is needed for the factorization to exclude fixed particles.
	for (size_t groupIndex = 0; groupIndex < m_objects.size(); groupIndex++)
	{
		struct Comparator {
			Comparator(Elasticity_Kee2023* _this) : m_this(_this) {};
			bool operator()(unsigned int a, unsigned int b) 
			{
				if ((m_this->m_model->getParticleState(a) != ParticleState::Active) && (m_this->m_model->getParticleState(b) == ParticleState::Active))
					return false;
				else if ((m_this->m_model->getParticleState(a) == ParticleState::Active) && (m_this->m_model->getParticleState(b) != ParticleState::Active))
					return true;
				else if ((m_this->m_model->getParticleState(a) != ParticleState::Active) && (m_this->m_model->getParticleState(b) != ParticleState::Active))
					return a < b;
				else
					return a < b;
			}

			Elasticity_Kee2023* m_this;
		};

		std::sort(m_objects[groupIndex]->m_particleIndices.begin(), m_objects[groupIndex]->m_particleIndices.end(), Comparator(this));	

		m_objects[groupIndex]->m_nFixed = 0;
		for (size_t i = 0; i < m_objects[groupIndex]->m_particleIndices.size(); i++)
		{
			if (m_model->getParticleState(m_objects[groupIndex]->m_particleIndices[i]) != ParticleState::Active)
				m_objects[groupIndex]->m_nFixed++;
		}
		LOG_INFO << "Object " << groupIndex << " - fixed particles: " << m_objects[groupIndex]->m_nFixed;
	}
}

/** Initialize the solver for the linear system by either computing a factorization
 * or loading a factorization from the cache.
 */
void Elasticity_Kee2023::initSystem()
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int fluidModelIndex = m_model->getPointSetIndex();
	const Real dt = TimeManager::getCurrent()->getTimeStepSize();

	// Compute Lamé parameters 
	m_mu = m_youngsModulus / (static_cast<Real>(2.0) * (static_cast<Real>(1.0) + m_poissonRatio));
	m_lambda = m_youngsModulus * m_poissonRatio / ((static_cast<Real>(1.0) + m_poissonRatio) * (static_cast<Real>(1.0) - static_cast<Real>(2.0) * m_poissonRatio));
	FluidModel *model = m_model;

	size_t numObjects = m_objects.size();

	auto& fluidInfos = sim->getFluidInfos();
	for (auto objIndex = 0; objIndex < numObjects; objIndex++)
	{
		ElasticObject* obj = m_objects[objIndex];

		// compute MD5 check sum
		std::string md5 = computeMD5(objIndex);

		// check if object with same md5 already exists
		bool foundFactorization = false;
		for (size_t i = 0; i < objIndex; i++)
		{
			// reuse the factorization for all objects with the same particle sampling to reduce the memory consumption
			if (fluidInfos[objIndex].hasSameParticleSampling(fluidInfos[i]))
			{
				m_objects[objIndex]->m_factorization = m_objects[i]->m_factorization;
				foundFactorization = true;
				LOG_INFO << "Object " << objIndex << " is using the factorization of object " << i;
				break;
			}
		}
		// if no factorization was found, create a new one
		if (obj->m_factorization == nullptr)
			obj->m_factorization = std::make_shared<Factorization>();

		// generate file name for cache file
		std::string baseName = Utilities::FileSystem::getFileName(fluidInfos[objIndex].samplesFile);
		std::string ext = Utilities::FileSystem::getFileExt(fluidInfos[objIndex].samplesFile);
		std::string cacheFileName = sim->getCachePath() + "/" + baseName + "_" + ext +  
			"_" + std::to_string(fluidInfos[objIndex].mode) +
			"_" + md5 + "_" + Utilities::StringTools::real2String(dt) +
			"_" + Utilities::StringTools::real2String(m_mu) + 
			"_" + Utilities::StringTools::real2String(m_alpha) +
			".bin";

		// Fluid block
		if (fluidInfos[objIndex].type == 0)
		{
			Vector3r diag = (fluidInfos[objIndex].box.max() - fluidInfos[objIndex].box.min());
			cacheFileName = sim->getCachePath() + "/Block_" + 
				Utilities::StringTools::real2String(diag.x()) + "_" + Utilities::StringTools::real2String(diag.y()) + "_" + Utilities::StringTools::real2String(diag.z()) + 
				"_" + std::to_string(fluidInfos[objIndex].mode) +
				"_" + md5 + "_" + Utilities::StringTools::real2String(dt) +
				"_" + Utilities::StringTools::real2String(m_mu) +
				"_" + Utilities::StringTools::real2String(m_alpha) +
				".bin";
		}

		// check if cache file exists
		const bool foundCacheFile = Utilities::FileSystem::fileExists(cacheFileName);

		if (sim->getUseCache() && foundCacheFile)
		{
			// if factorization cannot be reused from another object and a cache file was found, load the cache file
			if (!foundFactorization)
			{
				LOG_INFO << "Read cached factorization: " << cacheFileName;
				BinaryFileReader binReader;
				binReader.openFile(cacheFileName);
				binReader.readSparseMatrix(obj->m_factorization->m_D);
				binReader.readSparseMatrix(obj->m_factorization->m_DT_K);
				binReader.read(obj->m_factorization->m_dt);
				binReader.read(obj->m_factorization->m_mu);
				binReader.readSparseMatrix(obj->m_factorization->m_matHTH);
#ifdef USE_AVX
				delete obj->m_factorization->m_cholesky;
				obj->m_factorization->m_cholesky = new CholeskyAVXSolver();
				obj->m_factorization->m_cholesky->load(binReader);
#else
				binReader.readSparseMatrix(obj->m_factorization->m_matL);
				binReader.readSparseMatrix(obj->m_factorization->m_matLT);
				binReader.readMatrixX(obj->m_factorization->m_permInd);
				binReader.readMatrixX(obj->m_factorization->m_permInvInd);
#endif
				binReader.closeFile();
			}

			// init vectors
			int numParticles = (int)obj->m_particleIndices.size();
#ifdef USE_AVX
			obj->m_dx_avx.resize(numParticles - obj->m_nFixed);
			obj->m_f_avx.resize(3 * numParticles);
			obj->m_xk_avx.resize(numParticles);
			obj->m_xTilde_avx.resize(numParticles);
#else
			obj->m_dx.resize(numParticles - obj->m_nFixed);
			obj->m_f.resize(3 * numParticles);
			obj->m_xk.resize(numParticles);
			obj->m_xTilde.resize(numParticles);
#endif
			int nFree = numParticles - obj->m_nFixed;
#ifdef USE_AVX
			obj->m_gradient_avx.resize(numParticles, Scalarf8(0.0f));
			// L-BFGS buffers (Scalarf8 domain)
			obj->m_lbfgs_s_avx.resize(m_lbfgsWindowSize);
			obj->m_lbfgs_y_avx.resize(m_lbfgsWindowSize);
			for (int w = 0; w < m_lbfgsWindowSize; w++)
			{
				obj->m_lbfgs_s_avx[w].resize(nFree, Scalarf8(0.0f));
				obj->m_lbfgs_y_avx[w].resize(nFree, Scalarf8(0.0f));
			}
			obj->m_lbfgs_rho.resize(m_lbfgsWindowSize, 0);
			obj->m_lbfgs_alpha.resize(m_lbfgsWindowSize, 0);
			obj->m_lbfgs_last_sol_avx.resize(nFree, Scalarf8(0.0f));
			obj->m_lbfgs_q_avx.resize(nFree, Scalarf8(0.0f));
			obj->m_lbfgs_count = 0;
			// Newton CG workspace (Scalarf8, coord-packed)
			obj->m_pcg_r_avx.resize(nFree, Scalarf8(0.0f));
			obj->m_pcg_p_avx.resize(nFree, Scalarf8(0.0f));
			obj->m_pcg_Ap_avx.resize(nFree, Scalarf8(0.0f));
			obj->m_pcg_z_avx.resize(nFree, Scalarf8(0.0f));
#else
			obj->m_dx_perm.resize(numParticles - obj->m_nFixed);
			obj->m_gradient.resize(numParticles, Vector3r::Zero());
			// L-BFGS buffers
			obj->m_lbfgs_s.resize(m_lbfgsWindowSize);
			obj->m_lbfgs_y.resize(m_lbfgsWindowSize);
			for (int w = 0; w < m_lbfgsWindowSize; w++)
			{
				obj->m_lbfgs_s[w].resize(nFree, Vector3r::Zero());
				obj->m_lbfgs_y[w].resize(nFree, Vector3r::Zero());
			}
			obj->m_lbfgs_rho.resize(m_lbfgsWindowSize, 0);
			obj->m_lbfgs_alpha.resize(m_lbfgsWindowSize, 0);
			obj->m_lbfgs_last_sol.resize(nFree, Vector3r::Zero());
			obj->m_lbfgs_q.resize(nFree, Vector3r::Zero());
			obj->m_lbfgs_count = 0;

			// Newton buffers (scalar PCG workspace)
			obj->m_pcg_r.resize(nFree, Vector3r::Zero());
			obj->m_pcg_p.resize(nFree, Vector3r::Zero());
			obj->m_pcg_Ap.resize(nFree, Vector3r::Zero());
			obj->m_pcg_z.resize(nFree, Vector3r::Zero());
#endif
			// Shared Newton buffers (both builds)
			obj->m_hessian9x9.resize(numParticles);
			obj->m_pcg_precond.resize(nFree, Matrix3r::Identity());
		}
		else    // no cache found
		{
			ElasticObject* obj = m_objects[objIndex];

			// compute new factorization if no factorization can be reused
			if (!foundFactorization)
				initFactorization(obj->m_factorization, obj->m_particleIndices, obj->m_nFixed, dt, m_mu);

			// init vectors
			int numParticles = (int)obj->m_particleIndices.size();
#ifdef USE_AVX
			obj->m_dx_avx.resize(numParticles - obj->m_nFixed);
			obj->m_f_avx.resize(3 * numParticles);
			obj->m_xk_avx.resize(numParticles);
			obj->m_xTilde_avx.resize(numParticles);
#else
			obj->m_dx.resize(numParticles - obj->m_nFixed);
			obj->m_f.resize(3 * numParticles);
			obj->m_xk.resize(numParticles);
			obj->m_xTilde.resize(numParticles);
#endif
			int nFree = numParticles - obj->m_nFixed;
#ifdef USE_AVX
			obj->m_gradient_avx.resize(numParticles, Scalarf8(0.0f));
			// L-BFGS buffers (Scalarf8 domain)
			obj->m_lbfgs_s_avx.resize(m_lbfgsWindowSize);
			obj->m_lbfgs_y_avx.resize(m_lbfgsWindowSize);
			for (int w = 0; w < m_lbfgsWindowSize; w++)
			{
				obj->m_lbfgs_s_avx[w].resize(nFree, Scalarf8(0.0f));
				obj->m_lbfgs_y_avx[w].resize(nFree, Scalarf8(0.0f));
			}
			obj->m_lbfgs_rho.resize(m_lbfgsWindowSize, 0);
			obj->m_lbfgs_alpha.resize(m_lbfgsWindowSize, 0);
			obj->m_lbfgs_last_sol_avx.resize(nFree, Scalarf8(0.0f));
			obj->m_lbfgs_q_avx.resize(nFree, Scalarf8(0.0f));
			obj->m_lbfgs_count = 0;
			// Newton CG workspace (Scalarf8, coord-packed)
			obj->m_pcg_r_avx.resize(nFree, Scalarf8(0.0f));
			obj->m_pcg_p_avx.resize(nFree, Scalarf8(0.0f));
			obj->m_pcg_Ap_avx.resize(nFree, Scalarf8(0.0f));
			obj->m_pcg_z_avx.resize(nFree, Scalarf8(0.0f));
#else
			obj->m_dx_perm.resize(numParticles - obj->m_nFixed);
			obj->m_gradient.resize(numParticles, Vector3r::Zero());
			// L-BFGS buffers
			obj->m_lbfgs_s.resize(m_lbfgsWindowSize);
			obj->m_lbfgs_y.resize(m_lbfgsWindowSize);
			for (int w = 0; w < m_lbfgsWindowSize; w++)
			{
				obj->m_lbfgs_s[w].resize(nFree, Vector3r::Zero());
				obj->m_lbfgs_y[w].resize(nFree, Vector3r::Zero());
			}
			obj->m_lbfgs_rho.resize(m_lbfgsWindowSize, 0);
			obj->m_lbfgs_alpha.resize(m_lbfgsWindowSize, 0);
			obj->m_lbfgs_last_sol.resize(nFree, Vector3r::Zero());
			obj->m_lbfgs_q.resize(nFree, Vector3r::Zero());
			obj->m_lbfgs_count = 0;

			// Newton buffers (scalar PCG workspace)
			obj->m_pcg_r.resize(nFree, Vector3r::Zero());
			obj->m_pcg_p.resize(nFree, Vector3r::Zero());
			obj->m_pcg_Ap.resize(nFree, Vector3r::Zero());
			obj->m_pcg_z.resize(nFree, Vector3r::Zero());
#endif
			// Shared Newton buffers (both builds)
			obj->m_hessian9x9.resize(numParticles);
			obj->m_pcg_precond.resize(nFree, Matrix3r::Identity());

			// write cache file
			if (sim->getUseCache() && (Utilities::FileSystem::makeDir(sim->getCachePath()) == 0))
			{
				BinaryFileWriter binWriter;
				binWriter.openFile(cacheFileName);
				binWriter.writeSparseMatrix(obj->m_factorization->m_D);
				binWriter.writeSparseMatrix(obj->m_factorization->m_DT_K);
				binWriter.write(obj->m_factorization->m_dt);
				binWriter.write(obj->m_factorization->m_mu);
				binWriter.writeSparseMatrix(obj->m_factorization->m_matHTH);
#ifdef USE_AVX
				obj->m_factorization->m_cholesky->save(binWriter);
#else 
				binWriter.writeSparseMatrix(obj->m_factorization->m_matL);
				binWriter.writeSparseMatrix(obj->m_factorization->m_matLT);
				binWriter.writeMatrixX(obj->m_factorization->m_permInd);
				binWriter.writeMatrixX(obj->m_factorization->m_permInvInd);
#endif
				binWriter.closeFile();
			}
		}
		obj->m_md5 = md5;
	}
}

/** Compute the factorization of the linear system matrix. 
* This is only done once at the beginning of the simulation. 
*/
void Elasticity_Kee2023::initFactorization(std::shared_ptr<Factorization> factorization, std::vector<unsigned int>& particleIndices, const unsigned int nFixed, const Real dt, const Real mu)
{
	Simulation* sim = Simulation::getCurrent();
	const unsigned int fluidModelIndex = m_model->getPointSetIndex();

	factorization->m_dt = dt;
	factorization->m_mu = mu;

	// init mapping to find the particle indices in the current particle group (object)
	std::vector<unsigned int> groupInv;
	groupInv.resize(m_model->numActiveParticles());

	int numParticles = (int)particleIndices.size();
	std::vector<unsigned int>& group = particleIndices;

	for (int i = 0; i < numParticles; i++)
		groupInv[group[i]] = i;

	// determine total number of neighbors
	int totalNeighbors = 0;
	for (int i = 0; i < (int)numParticles; i++)
	{
		int particleIndex = group[i];
		const int numNeighbors = (int) m_initialNeighbors[i].size();
		totalNeighbors += numNeighbors;
	}

	// init triplets for matrices D, K, M
	std::vector<Eigen::Triplet<double>> triplets_D;
	std::vector<Eigen::Triplet<double>> triplets_K;
	std::vector<Eigen::Triplet<double>> triplets_M;

	triplets_D.reserve(2 * 3 * totalNeighbors);
	triplets_K.reserve(3 * numParticles);
	triplets_M.reserve(numParticles);

	std::vector<Eigen::Triplet<double>> triplets_H;
	triplets_H.reserve(totalNeighbors + numParticles);
	std::vector<Eigen::Triplet<double>> triplets_K2;
	triplets_K2.reserve(totalNeighbors);

	const double dtd = dt;
	const double mud = mu;
	const double alphad = m_alpha;

	// init matrices D and K
	unsigned int row_index = 0;
	for (int i = 0; i < (int)numParticles; i++)
	{
		int particleIndex0 = group[i];
		const unsigned int i0 = particleIndex0;
		unsigned int particleIndex = m_initial_to_current_index[i0];

		const double restVolumes_id = m_restVolumes[i];

		// init matrix K: constant approximation (2*mu + lambda) for L-BFGS H_0
		const double lambdad = static_cast<double>(m_lambda);
		const double Kreal = (2.0 * mud + lambdad) * dtd * dtd * restVolumes_id;
		

		// init mass matrix
		// all particles have the same mass => M = mass*I
		triplets_M.push_back(Eigen::Triplet<double>(i, i, m_model->getMass(particleIndex)));
		for (int j = 0; j < 3; j++)
			triplets_K.push_back(Eigen::Triplet<double>(3 * i + j, 3 * i + j, Kreal));

		const Eigen::Vector3d xi0d = m_model->getPosition0(i0).cast<double>();
		const size_t numNeighbors = m_initialNeighbors[i0].size();

		for (unsigned int j = 0; j < numNeighbors; j++)
		{
			const unsigned int neighborIndex0 = m_initialNeighbors[i0][j];
			const unsigned int neighborIndex = m_initial_to_current_index[neighborIndex0];
			const Eigen::Vector3d xj0d = m_model->getPosition0(neighborIndex0).cast<double>();

			const Eigen::Vector3d correctedKernelGradient = m_L[particleIndex].cast<double>() * sim->gradW((xi0d - xj0d).cast<Real>()).cast<double>();
			const double restVolumes_jd = m_restVolumes[neighborIndex];
			const Eigen::Vector3d y = restVolumes_jd * correctedKernelGradient;

			// init matrix D according to Eq. 3 and 4 in the supplemental document of the paper
			triplets_D.push_back(Eigen::Triplet<double>(3 * i + 0, i, -y[0]));
			triplets_D.push_back(Eigen::Triplet<double>(3 * i + 1, i, -y[1]));
			triplets_D.push_back(Eigen::Triplet<double>(3 * i + 2, i, -y[2]));

			triplets_D.push_back(Eigen::Triplet<double>(3 * i + 0, groupInv[neighborIndex0], y[0]));
			triplets_D.push_back(Eigen::Triplet<double>(3 * i + 1, groupInv[neighborIndex0], y[1]));
			triplets_D.push_back(Eigen::Triplet<double>(3 * i + 2, groupInv[neighborIndex0], y[2]));

			double sum = 0.0;
			const Eigen::Vector3d xi_xj_0 = xi0d - xj0d;
			const double xixj0_l2 = xi_xj_0.squaredNorm();
			const double beta = 1.0; 

			// init matrix \tilde K * dt^2 according to Eq. 2 in the paper
			const double Kreal2 = alphad * dtd * dtd * mud * static_cast<double>(m_restVolumes[i]) * static_cast<double>(m_restVolumes[neighborIndex]) * (sim->W(xi_xj_0.cast<Real>()) / xixj0_l2);
			triplets_K2.push_back(Eigen::Triplet<double>(row_index, row_index, Kreal2));

			// init matrix H according to Eq. 5, 6 and 7 in the supplemental document of the paper
			for (unsigned int k = 0; k < numNeighbors; k++)
			{
				const unsigned int kIndex0 = m_initialNeighbors[i0][k];
				const unsigned int kIndex = m_initial_to_current_index[kIndex0];
				const Eigen::Vector3d xk0d = m_model->getPosition0(kIndex0).cast<double>();
				
				const Eigen::Vector3d correctedKernelGradientd = m_L[particleIndex].cast<double>() * sim->gradW((xi0d - xk0d).cast<Real>()).cast<double>();
				const double restVolumes_kd = m_restVolumes[kIndex];
				const double H3j = beta * restVolumes_kd * correctedKernelGradientd.dot(xi0d - xj0d);
				triplets_H.push_back(Eigen::Triplet<double>(row_index, groupInv[kIndex0], H3j));
				sum -= H3j;

				if (k == j)
				{
					triplets_H.push_back(Eigen::Triplet<double>(row_index, groupInv[kIndex0], beta));
				}
			}
			triplets_H.push_back(Eigen::Triplet<double>(row_index, i, sum - beta));
			row_index++;
		}
	}

	Eigen::SparseMatrix<double> D(3 * numParticles, numParticles);
	factorization->m_D.resize(3 * numParticles, numParticles);
	D.setFromTriplets(triplets_D.begin(), triplets_D.end());
	factorization->m_D = D.cast<Real>();

	//set matrices
	Eigen::SparseMatrix<double> K(3 * numParticles, 3 * numParticles); // actually 2 * dt* dt * K	
	K.setFromTriplets(triplets_K.begin(), triplets_K.end());

	Eigen::SparseMatrix<double> M(numParticles, numParticles);
	M.setFromTriplets(triplets_M.begin(), triplets_M.end());

	Eigen::SparseMatrix<double> DT_K(numParticles, 3 * numParticles);
	factorization->m_DT_K.resize(numParticles, 3 * numParticles);
	DT_K = D.transpose() * K;
	factorization->m_DT_K = DT_K.cast<Real>();

	Eigen::SparseMatrix<double> K2(totalNeighbors, totalNeighbors);
	K2.setFromTriplets(triplets_K2.begin(), triplets_K2.end());

	Eigen::SparseMatrix<double> H(totalNeighbors, numParticles);
	H.setFromTriplets(triplets_H.begin(), triplets_H.end());
	Eigen::SparseMatrix<double> HTH(numParticles, numParticles);
	factorization->m_matHTH.resize(numParticles, numParticles);
	HTH = H.transpose() * K2 * H;
	factorization->m_matHTH = HTH.cast<Real>();
	LOG_INFO << "Non zero elements (H^T * K * H): " << factorization->m_matHTH.nonZeros();

	Eigen::SparseMatrix<double> DT_K_D = DT_K * D;


	Eigen::SparseMatrix<double> M_plus_DT_K_D;
	// init linear system matrix according to Eq. 29 in the paper
	if (m_alpha != 0.0)
		M_plus_DT_K_D = (M + DT_K_D + HTH).block(0, 0, numParticles - nFixed, numParticles - nFixed);
	else      // no zero energy mode control
		M_plus_DT_K_D = (M + DT_K_D).block(0, 0, numParticles - nFixed, numParticles - nFixed);

	M_plus_DT_K_D.makeCompressed();

#ifdef USE_AVX
	// compute factorization of the matrix
	delete factorization->m_cholesky;
	factorization->m_cholesky = new CholeskyAVXSolver(M_plus_DT_K_D);
#else
	SolverLLT* solverLLT = new SolverLLT();
	solverLLT->compute(M_plus_DT_K_D);

	if (solverLLT->info() != Eigen::Success)
	{
		LOG_ERR << "Cholesky decomposition failed.";
		LOG_ERR << solverLLT->info();
		return;
	}

	factorization->m_permInd = solverLLT->permutationP().indices();
	factorization->m_permInvInd = solverLLT->permutationPinv().indices();
	factorization->m_matL = Eigen::SparseMatrix<Real, Eigen::ColMajor>(solverLLT->matrixL().cast<Real>());
	factorization->m_matLT = Eigen::SparseMatrix<Real, Eigen::ColMajor>(solverLLT->matrixU().cast<Real>());

	delete solverLLT;
#endif

	LOG_INFO << "Non zero elements (A): " << M_plus_DT_K_D.nonZeros();
}

/** Perform a step of the elasticity solver.
 */
void Elasticity_Kee2023::step()
{
	// apply accelerations
	const unsigned int numParticles = m_model->numActiveParticles();
	const Real dt = TimeManager::getCurrent()->getTimeStepSize();
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)
		for (int i = 0; i < (int)numParticles; i++)
		{
			if (m_model->getParticleState(i) == ParticleState::Active)
			{
				Vector3r& vel = m_model->getVelocity(i);
				vel += dt * m_model->getAcceleration(i);
				m_model->getAcceleration(i).setZero();
			}
		}
	}

	START_TIMING("Elasticity")
	stepElasticitySolver();					// elasticity solver using the factorization
	STOP_TIMING_AVG
}

void Elasticity_Kee2023::reset()
{
	initValues();
}

void Elasticity_Kee2023::performNeighborhoodSearchSort()
{
	const unsigned int numPart = m_model->numActiveParticles();
	if (numPart == 0)
		return;

	Simulation *sim = Simulation::getCurrent();
	auto const& d = sim->getNeighborhoodSearch()->point_set(m_model->getPointSetIndex());
	d.sort_field(&m_rotations[0]);
	d.sort_field(&m_current_to_initial_index[0]);
	d.sort_field(&m_L[0]);
	d.sort_field(&m_restVolumes[0]);

	for (unsigned int i = 0; i < numPart; i++)
		m_initial_to_current_index[m_current_to_initial_index[i]] = i;
}

/** Compute kernel gradient correction matrices (Eq. 8).
*/
void Elasticity_Kee2023::computeMatrixL()
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
				// get initial neighbor index considering the current particle order 
				const unsigned int neighborIndex0 = m_initialNeighbors[i0][j];
				const unsigned int neighborIndex = m_initial_to_current_index[neighborIndex0];

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
			L.computeInverseWithCheck(m_L[i], invertible, static_cast<Real>(1e-9));
			if (!invertible)
			{
				LOG_INFO << "Matrix not invertible.";
				MathFunctions::pseudoInverse(L, m_L[i]);
				//m_L[i] = Matrix3r::Identity();
			}
		}
	}
}

/** Precompute V_j * gradW(xi0 - xj0) for all neighbor pairs.
*   Depends only on rest positions and volumes — constant for the simulation.
*   Called once at init (and after any neighbor sort).
*   AVX version packs 8 neighbors per Vector3f8 entry (Kugelstadt pattern).
*/
void Elasticity_Kee2023::precomputeValues()
{
	Simulation* sim = Simulation::getCurrent();
	const int numParticles = (int)m_model->numActiveParticles();

	// Scalar format (always filled — used by Newton CG in both builds)
	m_precomputed_indices.clear();
	m_precomp_V_gradW.clear();
	m_precomputed_indices.resize(numParticles);

	unsigned int sumNeighbors = 0;
	for (int i = 0; i < numParticles; i++)
	{
		const unsigned int i0 = m_current_to_initial_index[i];
		const size_t numNeighbors = m_initialNeighbors[i0].size();
		m_precomputed_indices[i] = sumNeighbors;
		sumNeighbors += numNeighbors;
	}
	m_precomp_V_gradW.resize(sumNeighbors);

	#pragma omp parallel for schedule(static)
	for (int i = 0; i < numParticles; i++)
	{
		const unsigned int i0 = m_current_to_initial_index[i];
		const Vector3r xi0 = m_model->getPosition0(i0);
		const unsigned int numNeighbors = (unsigned int)m_initialNeighbors[i0].size();
		unsigned int base = m_precomputed_indices[i];

		for (unsigned int j = 0; j < numNeighbors; j++)
		{
			const unsigned int neighborIndex0 = m_initialNeighbors[i0][j];
			const unsigned int neighborCurrent = m_initial_to_current_index[neighborIndex0];
			const Vector3r xj0 = m_model->getPosition0(neighborIndex0);
			const Real V_j = m_restVolumes[neighborCurrent];
			m_precomp_V_gradW[base + j] = V_j * sim->gradW(xi0 - xj0);
		}
	}

#ifdef USE_AVX
	// AVX-packed format (additionally filled for L-BFGS force loops)
	m_precomputed_indices8.clear();
	m_precomp_V_gradW8.clear();
	m_precomputed_indices8.resize(numParticles);

	unsigned int sumBlocks = 0;
	for (int i = 0; i < numParticles; i++)
	{
		const unsigned int i0 = m_current_to_initial_index[i];
		const size_t numNeighbors = m_initialNeighbors[i0].size();
		m_precomputed_indices8[i] = sumBlocks;
		sumBlocks += (unsigned int)numNeighbors / 8u;
		if (numNeighbors % 8 != 0)
			sumBlocks++;
	}

	m_precomp_V_gradW8.resize(sumBlocks);

	#pragma omp parallel for schedule(static)
	for (int i = 0; i < numParticles; i++)
	{
		const unsigned int i0 = m_current_to_initial_index[i];
		const Vector3r xi0 = m_model->getPosition0(i0);
		const unsigned int numNeighbors = (unsigned int)m_initialNeighbors[i0].size();

		const Vector3f8 xi0_avx(xi0.cast<float>());
		unsigned int base8 = m_precomputed_indices8[i];
		unsigned int idx = 0;

		for (unsigned int j = 0; j < numNeighbors; j += 8)
		{
			const unsigned int count = std::min(numNeighbors - j, 8u);
			unsigned int nIndices[8];
			for (auto k = 0u; k < count; k++)
				nIndices[k] = m_initial_to_current_index[m_initialNeighbors[i0][j + k]];

			const Scalarf8 Vj_avx = convert_zero(nIndices, &m_restVolumes[0], count);
			const Vector3f8 xj0_avx = convertVec_zero(&m_initialNeighbors[i0][j], &m_model->getPosition0(0), count);
			const Vector3f8 gradW_avx = CubicKernel_AVX::gradW(xi0_avx - xj0_avx);
			m_precomp_V_gradW8[base8 + idx] = gradW_avx * Vj_avx;
			idx++;
		}
	}
#endif
}

void Elasticity_Kee2023::saveState(BinaryFileWriter &binWriter)
{
	binWriter.writeBuffer((char*)m_current_to_initial_index.data(), m_current_to_initial_index.size() * sizeof(unsigned int));
	binWriter.writeBuffer((char*)m_initial_to_current_index.data(), m_initial_to_current_index.size() * sizeof(unsigned int));
	binWriter.writeBuffer((char*)m_L.data(), m_L.size() * sizeof(Matrix3r));
	binWriter.writeBuffer((char*)m_rotations.data(), m_rotations.size() * sizeof(Matrix3r));
}

void Elasticity_Kee2023::loadState(BinaryFileReader &binReader)
{
	binReader.readBuffer((char*)m_current_to_initial_index.data(), m_current_to_initial_index.size() * sizeof(unsigned int));
	binReader.readBuffer((char*)m_initial_to_current_index.data(), m_initial_to_current_index.size() * sizeof(unsigned int));
	binReader.readBuffer((char*)m_L.data(), m_L.size() * sizeof(Matrix3r));
	binReader.readBuffer((char*)m_rotations.data(), m_rotations.size() * sizeof(Matrix3r));
}


/** Compute predicted position: xTilde = x + dt * v
*/
void Elasticity_Kee2023::computeXTilde(ElasticObject* obj)
{
	const std::vector<unsigned int>& group = obj->m_particleIndices;
	const int numParticles = (int)group.size();
	const Real dt = obj->m_factorization->m_dt;
#ifdef USE_AVX
	auto& xTilde = obj->m_xTilde_avx;
#else
	auto& xTilde = obj->m_xTilde;
#endif

	#pragma omp parallel for schedule(static)
	for (int i = 0; i < numParticles; i++)
	{
		const unsigned int i0 = group[i];
		const unsigned int particleIndex = m_initial_to_current_index[i0];
		const Vector3r xT = m_model->getPosition(particleIndex) + dt * m_model->getVelocity(particleIndex);
#ifdef USE_AVX
		xTilde[i] = Scalarf8((float)xT[0], (float)xT[1], (float)xT[2], 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
#else
		xTilde[i] = xT;
#endif
	}
}

/** Update velocities from solved positions: v = (xk - x) / dt
*/
void Elasticity_Kee2023::updateVelocity(ElasticObject* obj, Real fdt)
{
	const std::vector<unsigned int>& group = obj->m_particleIndices;
	const int numParticles = (int)group.size();
#ifdef USE_AVX
	const auto& xk = obj->m_xk_avx;
#else
	const auto& xk = obj->m_xk;
#endif
	const Real damping = static_cast<Real>(0.0);
	const Real invFdt = (1 - damping) / fdt;
	#pragma omp parallel for schedule(static)
	for (int i = 0; i < numParticles; i++)
	{
		const unsigned int i0 = group[i];
		const unsigned int particleIndex = m_initial_to_current_index[i0];
		if (m_model->getParticleState(particleIndex) == ParticleState::Active)
		{
			const Vector3r& x = m_model->getPosition(particleIndex);
#ifdef USE_AVX
			float v[8];
			xk[i].store(v);
			const Vector3r xk_i((Real)v[0], (Real)v[1], (Real)v[2]);
			m_model->getVelocity(particleIndex) = invFdt * (xk_i - x);
#else
			m_model->getVelocity(particleIndex) = invFdt * (xk[i] - x);
#endif
		}
	}
}

/** Evaluate the objective energy at the current xk (no gradient computation).
*   Used for backtracking line search.
*/
#ifdef USE_AVX
Real Elasticity_Kee2023::computeEnergy(ElasticObject* obj)
{
	const std::vector<unsigned int>& group = obj->m_particleIndices;
	int numParticles = (int)group.size();

	auto& D = obj->m_factorization->m_D;
	auto& HT_K_H = obj->m_factorization->m_matHTH;
	auto& f = obj->m_f_avx;          // Scalarf8 (3*numParticles entries)
	auto& xk = obj->m_xk_avx;        // Scalarf8 (numParticles entries)
	auto& xTilde = obj->m_xTilde_avx;
	const Real fdt = obj->m_factorization->m_dt;

	Real elasticEnergy = 0;
	Real massEnergy = 0;

	// F = D * xk  (coord-packed AVX sparse matvec)
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)
		for (int i = 0; i < 3 * numParticles; i++)
			f[i].setZero();

		#pragma omp for schedule(static)
		for (int k = 0; k < D.outerSize(); ++k)
		{
			for (Eigen::SparseMatrix<Real, Eigen::RowMajor>::InnerIterator it(D, k); it; ++it)
				f[it.row()] += Scalarf8((float)it.value()) * xk[it.col()];
		}

		// Elastic energy: sum V_i * Psi(F_i)
		#pragma omp for reduction(+:elasticEnergy) schedule(static)
		for (int i = 0; i < numParticles; i++)
		{
			const unsigned int i0 = group[i];
			const unsigned int particleIndex = m_initial_to_current_index[i0];

			float x0[8], x1[8], x2[8];
			f[3 * i].store(x0);
			f[3 * i + 1].store(x1);
			f[3 * i + 2].store(x2);
			Matrix3r Fi;
			Fi <<	x0[0], x0[1], x0[2],
					x1[0], x1[1], x1[2],
					x2[0], x2[1], x2[2];

			Eigen::JacobiSVD<Matrix3r> svd_i(Fi, Eigen::ComputeFullU | Eigen::ComputeFullV);
			Matrix3r U = svd_i.matrixU();
			if ((U * svd_i.matrixV().transpose()).determinant() < 0)
				U.col(2) = -U.col(2);
			Matrix3r Ri = U * svd_i.matrixV().transpose();

			elasticEnergy += m_restVolumes[particleIndex] * computePsi(Fi, Ri);
		}
	}

	// Inertial energy: (1/2) m_i ||xk_i - xTilde_i||^2
	// Vertical accumulation: each thread keeps a Scalarf8 acc, reduces once at the end.
	#pragma omp parallel reduction(+:massEnergy)
	{
		Scalarf8 acc; acc.setZero();
		#pragma omp for schedule(static) nowait
		for (int i = 0; i < numParticles; i++)
		{
			const unsigned int i0 = group[i];
			const unsigned int particleIndex = m_initial_to_current_index[i0];
			const Real mass = m_model->getMass(particleIndex);
			const Scalarf8 diff = xk[i] - xTilde[i];
			acc += Scalarf8(static_cast<float>(0.5 * mass)) * diff * diff;
		}
		massEnergy += (Real)acc.reduce();
	}

	// Zero-energy mode: (1/2) xk^T HTH xk
	Real zeroEnergy = 0;
	if (m_alpha != 0.0)
	{
		#pragma omp parallel for reduction(+:zeroEnergy) schedule(static)
		for (int k = 0; k < HT_K_H.outerSize(); ++k)
		{
			for (Eigen::SparseMatrix<Real, Eigen::ColMajor>::InnerIterator it(HT_K_H, k); it; ++it)
				zeroEnergy += it.value() * (Real)(xk[it.row()] * xk[it.col()]).reduce();
		}
		zeroEnergy *= static_cast<Real>(0.5);
	}

	return massEnergy + fdt * fdt * elasticEnergy + zeroEnergy;
}
#else
Real Elasticity_Kee2023::computeEnergy(ElasticObject* obj)
{
	const std::vector<unsigned int>& group = obj->m_particleIndices;
	int numParticles = (int)group.size();

	auto& D = obj->m_factorization->m_D;
	auto& HT_K_H = obj->m_factorization->m_matHTH;
	auto& f = obj->m_f;
	auto& xk = obj->m_xk;
	const Real fdt = obj->m_factorization->m_dt;

	Real elasticEnergy = 0;
	Real massEnergy = 0;

	// F = D * xk
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)
		for (int i = 0; i < (int)numParticles; i++)
		{
			f[3 * i].setZero();
			f[3 * i + 1].setZero();
			f[3 * i + 2].setZero();
		}

		#pragma omp for schedule(static)
		for (int k = 0; k < D.outerSize(); ++k)
		{
			for (Eigen::SparseMatrix<Real, Eigen::RowMajor>::InnerIterator it(D, k); it; ++it)
				f[it.row()] += it.value() * xk[it.col()];
		}

		// Elastic energy: sum V_i * Psi(F_i)
		#pragma omp for reduction(+:elasticEnergy) schedule(static)
		for (int i = 0; i < (int)numParticles; i++)
		{
			const unsigned int i0 = group[i];
			const unsigned int particleIndex = m_initial_to_current_index[i0];

			Matrix3r Fi;
			Fi <<	f[3 * i][0], f[3 * i][1], f[3 * i][2],
					f[3 * i + 1][0], f[3 * i + 1][1], f[3 * i + 1][2],
					f[3 * i + 2][0], f[3 * i + 2][1], f[3 * i + 2][2];

			Eigen::JacobiSVD<Matrix3r> svd_i(Fi, Eigen::ComputeFullU | Eigen::ComputeFullV);
			Matrix3r U = svd_i.matrixU();
			if ((U * svd_i.matrixV().transpose()).determinant() < 0)
				U.col(2) = -U.col(2);
			Matrix3r Ri = U * svd_i.matrixV().transpose();

			elasticEnergy += m_restVolumes[particleIndex] * computePsi(Fi, Ri);
		}
	}

	// Inertial energy: (1/2) m_i ||xk_i - xTilde_i||^2
	auto& xTilde = obj->m_xTilde;
	#pragma omp parallel for reduction(+:massEnergy) schedule(static)
	for (int i = 0; i < (int)numParticles; i++)
	{
		const unsigned int i0 = group[i];
		const unsigned int particleIndex = m_initial_to_current_index[i0];
		const Real mass = m_model->getMass(particleIndex);
		massEnergy += static_cast<Real>(0.5) * mass * (xk[i] - xTilde[i]).squaredNorm();
	}

	// Zero-energy mode: (1/2) xk^T HTH xk
	Real zeroEnergy = 0;
	if (m_alpha != 0.0)
	{
		#pragma omp parallel for reduction(+:zeroEnergy) schedule(static)
		for (int k = 0; k < HT_K_H.outerSize(); ++k)
		{
			for (Eigen::SparseMatrix<Real, Eigen::ColMajor>::InnerIterator it(HT_K_H, k); it; ++it)
				zeroEnergy += it.value() * xk[it.row()].dot(xk[it.col()]);
		}
		zeroEnergy *= static_cast<Real>(0.5);
	}

	return massEnergy + fdt * fdt * elasticEnergy + zeroEnergy;
}
#endif

/** Compute the energy and gradient of the objective function.
*
* Energy:
*   E(x) = (1/2)(x - x_tilde)^T M (x - x_tilde) + dt^2 * sum_i V_i * Psi(F_i) + (1/2) x^T H^T K2 H x
*
* Energy density Psi depends on the selected material type:
*   - Stable Neo-Hookean (Smith et al. 2018):
*       Psi = (mu/2)(I_C - 3) + (lambda/2)(J - alpha)^2 - (mu/2)log(I_C + 1)
*   - Co-rotated:
*       Psi = mu * ||F - R||^2_F + (lambda/2) * (J - 1)^2
*
* Gradient (stored in m_gradient) of the objective:
*   g_i = m_i * (sol_i - x_tilde_i)
*       - dt^2 * V_i * sum_j V_j * (P_i^T * L_i + P_j^T * L_j) * gradW(x_i0 - x_j0)
*       + (H^T K2 H xk)_i
*
* Returns the total energy E(x).
*/

#ifdef USE_AVX
Real Elasticity_Kee2023::computeEnergyAndGradient(ElasticObject* obj)
{
	Simulation *sim = Simulation::getCurrent();
	const std::vector<unsigned int>& group = obj->m_particleIndices;
	int numParticles = (int)group.size();

	auto& D = obj->m_factorization->m_D;
	auto& HT_K_H = obj->m_factorization->m_matHTH;

	auto& gradient = obj->m_gradient_avx;   // Scalarf8
	auto& f = obj->m_f_avx;                 // Scalarf8 (3*numParticles)
	auto& xk = obj->m_xk_avx;               // Scalarf8
	auto& xTilde = obj->m_xTilde_avx;       // Scalarf8

	const Real fdt = obj->m_factorization->m_dt;

	Real elasticEnergy = 0;
	Real massEnergy = 0;

	//////////////////////////////////////////////////////////////////////////
	// 1. F = D * xk  (coord-packed AVX sparse matvec, Kugelstadt Pattern A)
	//////////////////////////////////////////////////////////////////////////
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)
		for (int i = 0; i < 3 * numParticles; i++)
			f[i].setZero();

		#pragma omp for schedule(static)
		for (int k = 0; k < D.outerSize(); ++k)
		{
			for (Eigen::SparseMatrix<Real, Eigen::RowMajor>::InnerIterator it(D, k); it; ++it)
				f[it.row()] += Scalarf8((float)it.value()) * xk[it.col()];
		}

		#pragma omp for schedule(static)
		for (int i = 0; i < numParticles; i++)
		{
			const unsigned int i0 = group[i];
			const unsigned int particleIndex = m_initial_to_current_index[i0];

			float x0[8], x1[8], x2[8];
			f[3 * i].store(x0);
			f[3 * i + 1].store(x1);
			f[3 * i + 2].store(x2);

			m_F[particleIndex] <<	x0[0], x0[1], x0[2],
									x1[0], x1[1], x1[2],
									x2[0], x2[1], x2[2];
		}
	}

	//////////////////////////////////////////////////////////////////////////
	// 2. Polar decomposition via SVD (scalar — no AVX SVD available).
	//    Precompute P^T * L for each particle (stored in m_PL).
	//    Accumulate elastic energy density Psi(F_i) * V_i.
	//////////////////////////////////////////////////////////////////////////
	#pragma omp parallel default(shared)
	{
		#pragma omp for reduction(+:elasticEnergy) schedule(static)
		for (int i = 0; i < numParticles; i++)
		{
			const unsigned int i0 = group[i];
			const unsigned int particleIndex = m_initial_to_current_index[i0];

			Eigen::JacobiSVD<Matrix3r> svd_i(m_F[particleIndex], Eigen::ComputeFullU | Eigen::ComputeFullV);
			Matrix3r U = svd_i.matrixU();
			if ((U * svd_i.matrixV().transpose()).determinant() < 0)
				U.col(2) = -U.col(2);
			m_rotations[particleIndex] = U * svd_i.matrixV().transpose();

			const Matrix3r P_i = computeP(m_F[particleIndex], m_rotations[particleIndex]);
			m_PL[particleIndex] = P_i.transpose() * m_L[particleIndex];

			elasticEnergy += m_restVolumes[particleIndex] * computePsi(m_F[particleIndex], m_rotations[particleIndex]);

			gradient[i].setZero();
		}
	}

	//////////////////////////////////////////////////////////////////////////
	// 3. Elastic gradient (scalar neighbor loop; write Scalarf8 gradient).
	//    g_i = -dt^2 * V_i * sum_j V_j * (P_i^T L_i + P_j^T L_j) * gradW
	//////////////////////////////////////////////////////////////////////////
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)
		for (int i = 0; i < numParticles; i++)
		{
			const unsigned int i0 = group[i];
			const unsigned int particleIndex = m_initial_to_current_index[i0];
			const Vector3r xi0 = m_model->getPosition0(i0);
			const Real V_i = m_restVolumes[particleIndex];
			const Matrix3r& PtL_i = m_PL[particleIndex];

			const unsigned int numNeighbors = (unsigned int)m_initialNeighbors[i0].size();
			const Matrix3f8 PtLi_avx(m_PL[particleIndex].cast<float>());
			Vector3f8 force_avx;
			force_avx.setZero();

			for (unsigned int j = 0; j < numNeighbors; j += 8)
			{
				const unsigned int count = std::min(numNeighbors - j, 8u);
				unsigned int nIndices[8];
				for (auto k = 0u; k < count; k++)
					nIndices[k] = m_initial_to_current_index[m_initialNeighbors[i0][j + k]];

				const Matrix3f8 PtLj_avx = convertMat_zero(nIndices, &m_PL[0], count);
				const Vector3f8& V_gradW = m_precomp_V_gradW8[m_precomputed_indices8[particleIndex] + j / 8];
				force_avx += PtLi_avx * V_gradW + PtLj_avx * V_gradW;
			}

			Vector3r force;
			force[0] = force_avx.x().reduce();
			force[1] = force_avx.y().reduce();
			force[2] = force_avx.z().reduce();

			const Vector3r g = -(fdt * fdt * V_i * force);
			gradient[i] = Scalarf8((float)g[0], (float)g[1], (float)g[2], 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
		}
	}

	//////////////////////////////////////////////////////////////////////////
	// 4. Mass term: g_i += m_i * (xk_i - xTilde_i)
	//    Vertical accumulation for massEnergy (one .reduce() per thread).
	//////////////////////////////////////////////////////////////////////////
	#pragma omp parallel reduction(+:massEnergy)
	{
		Scalarf8 acc; acc.setZero();
		#pragma omp for schedule(static) nowait
		for (int i = 0; i < numParticles; i++)
		{
			const unsigned int i0 = group[i];
			const unsigned int particleIndex = m_initial_to_current_index[i0];
			const Real mass = m_model->getMass(particleIndex);
			const Scalarf8 diff = xk[i] - xTilde[i];
			gradient[i] += Scalarf8((float)mass) * diff;
			acc += Scalarf8(static_cast<float>(0.5 * mass)) * diff * diff;
		}
		massEnergy += (Real)acc.reduce();
	}

	//////////////////////////////////////////////////////////////////////////
	// 5. Zero-energy mode: g += H^T * K2 * H * xk  (coord-packed AVX sparse matvec)
	//////////////////////////////////////////////////////////////////////////
	Real zeroEnergy = 0;
	if (m_alpha != 0.0)
	{
		#pragma omp parallel default(shared)
		{
			#pragma omp for reduction(+:zeroEnergy) schedule(static)
			for (int k = 0; k < HT_K_H.outerSize(); ++k)
			{
				for (Eigen::SparseMatrix<Real, Eigen::ColMajor>::InnerIterator it(HT_K_H, k); it; ++it)
				{
					zeroEnergy += it.value() * (Real)(xk[it.row()] * xk[it.col()]).reduce();
					gradient[it.col()] += Scalarf8((float)it.value()) * xk[it.row()];
				}
			}
		}
		zeroEnergy *= static_cast<Real>(0.5);
	}

	return massEnergy + fdt * fdt * elasticEnergy + zeroEnergy;
}
#else
Real Elasticity_Kee2023::computeEnergyAndGradient(ElasticObject* obj)
{
	Simulation *sim = Simulation::getCurrent();
	const std::vector<unsigned int>& group = obj->m_particleIndices;
	int numParticles = (int)group.size();

	auto& D = obj->m_factorization->m_D;
	auto& HT_K_H = obj->m_factorization->m_matHTH;

	auto& gradient = obj->m_gradient;
	auto& f = obj->m_f;
	auto& xk = obj->m_xk;

	const Real fdt = obj->m_factorization->m_dt;

	Real elasticEnergy = 0;
	Real massEnergy = 0;

	//////////////////////////////////////////////////////////////////////////
	// 1. Compute deformation gradient: F = D * xk  ->  m_F
	//////////////////////////////////////////////////////////////////////////
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)
		for (int i = 0; i < numParticles; i++)
		{
			f[3 * i].setZero();
			f[3 * i + 1].setZero();
			f[3 * i + 2].setZero();
		}

		#pragma omp for schedule(static)
		for (int k = 0; k < D.outerSize(); ++k)
		{
			for (Eigen::SparseMatrix<Real, Eigen::RowMajor>::InnerIterator it(D, k); it; ++it)
				f[it.row()] += it.value() * xk[it.col()];
		}

		#pragma omp for schedule(static)
		for (int i = 0; i < numParticles; i++)
		{
			const unsigned int i0 = group[i];
			const unsigned int particleIndex = m_initial_to_current_index[i0];

			m_F[particleIndex] <<	f[3 * i][0], f[3 * i][1], f[3 * i][2],
									f[3 * i + 1][0], f[3 * i + 1][1], f[3 * i + 1][2],
									f[3 * i + 2][0], f[3 * i + 2][1], f[3 * i + 2][2];
		}
	}

	//////////////////////////////////////////////////////////////////////////
	// 2. Extract rotation R from F via polar decomposition
	//    Precompute P^T * L for each particle (stored in m_PL)
	//    Accumulate elastic energy density Psi(F_i) * V_i
	//////////////////////////////////////////////////////////////////////////
	#pragma omp parallel default(shared)
	{
		#pragma omp for reduction(+:elasticEnergy) schedule(static)
		for (int i = 0; i < (int)numParticles; i++)
		{
			const unsigned int i0 = group[i];
			const unsigned int particleIndex = m_initial_to_current_index[i0];

			Eigen::JacobiSVD<Matrix3r> svd_i(m_F[particleIndex], Eigen::ComputeFullU | Eigen::ComputeFullV);
			Matrix3r U = svd_i.matrixU();
			if ((U * svd_i.matrixV().transpose()).determinant() < 0)
				U.col(2) = -U.col(2);
			m_rotations[particleIndex] = U * svd_i.matrixV().transpose();

			const Matrix3r P_i = computeP(m_F[particleIndex], m_rotations[particleIndex]);
			m_PL[particleIndex] = P_i.transpose() * m_L[particleIndex];

			elasticEnergy += m_restVolumes[particleIndex] * computePsi(m_F[particleIndex], m_rotations[particleIndex]);

			gradient[i].setZero();
		}
	}

	//////////////////////////////////////////////////////////////////////////
	// 3. Elastic gradient: g_i = -dt^2 * V_i * sum_j V_j * (P_i^T L_i + P_j^T L_j) * gradW
	//////////////////////////////////////////////////////////////////////////
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)
		for (int i = 0; i < (int)numParticles; i++)
		{
			const unsigned int i0 = group[i];
			const unsigned int particleIndex = m_initial_to_current_index[i0];
			const Vector3r xi0 = m_model->getPosition0(i0);
			const Real V_i = m_restVolumes[particleIndex];
			const Matrix3r& PtL_i = m_PL[particleIndex];

			const size_t numNeighbors = m_initialNeighbors[i0].size();
			Vector3r force;
			force.setZero();

			for (unsigned int j = 0; j < numNeighbors; j++)
			{
				const unsigned int neighborIndex0 = m_initialNeighbors[i0][j];
				const unsigned int neighborCurrent = m_initial_to_current_index[neighborIndex0];
				const Matrix3r& PtL_j = m_PL[neighborCurrent];
				force += (PtL_i + PtL_j) * m_precomp_V_gradW[m_precomputed_indices[particleIndex] + j];
			}

			gradient[i] = -(fdt * fdt * V_i * force);
		}
	}

	//////////////////////////////////////////////////////////////////////////
	// 4. Mass term: g_i += m_i * (xk_i - xTilde_i)
	//////////////////////////////////////////////////////////////////////////
	auto& xTilde = obj->m_xTilde;
	#pragma omp parallel default(shared)
	{
		#pragma omp for reduction(+:massEnergy) schedule(static)
		for (int i = 0; i < (int)numParticles; i++)
		{
			const unsigned int i0 = group[i];
			const unsigned int particleIndex = m_initial_to_current_index[i0];
			const Real mass = m_model->getMass(particleIndex);
			gradient[i] += mass * (xk[i] - xTilde[i]);
			massEnergy += static_cast<Real>(0.5) * mass * (xk[i] - xTilde[i]).squaredNorm();
		}
	}

	//////////////////////////////////////////////////////////////////////////
	// 5. Zero-energy mode control: g += H^T * K2 * H * xk
	//////////////////////////////////////////////////////////////////////////
	Real zeroEnergy = 0;
	std::vector<Vector3r, Eigen::aligned_allocator<Vector3r>> gradient_ze(numParticles, Vector3r::Zero());
	if (m_alpha != 0.0)
	{
		#pragma omp parallel default(shared)
		{
			#pragma omp for reduction(+:zeroEnergy) schedule(static)
			for (int k = 0; k < HT_K_H.outerSize(); ++k)
			{
				for (Eigen::SparseMatrix<Real, Eigen::ColMajor>::InnerIterator it(HT_K_H, k); it; ++it)
				{
					zeroEnergy += it.value() * xk[it.row()].dot(xk[it.col()]);
					gradient_ze[it.col()] += it.value() * xk[it.row()];
				}
			}
		}
		zeroEnergy *= static_cast<Real>(0.5);

		for (int i = 0; i < (int)numParticles; i++)
			gradient[i] += gradient_ze[i];
	}

	return massEnergy + fdt * fdt * elasticEnergy + zeroEnergy;
}
#endif

/** Compute the Hessian of the elastic energy w.r.t. positions.
*
* For Newton: should be updated every iteration.
* For L-BFGS: constant proxy K = (2*mu + lambda) * V_i
* is prefactored in initFactorization() and used as H_0.
*/
void Elasticity_Kee2023::computeHessian(ElasticObject* obj)
{
	if (m_materialType == ENUM_MATERIAL_STABLE_NEOHOOKEAN)
		computeStableNeoHookeanHessian9x9(obj);
	else  // Co-rotated
		computeCorotatedHessian9x9(obj);
}

/** Compute the per-particle 9x9 Hessian d²ψ/d(vecF)² for the co-rotated model.
*
* Uses iARAP approach (Lin et al. 2022) for full decomposition:
* - Singular values extracted from quartic roots
* - R computed from Cauchy-Green invariants
* - V computed analytically from S = R^T F
*
* Co-rotated energy density: ψ = μ||F-R||² + (λ/2)(tr(R^T F) - 3)²
*
* The 9x9 Hessian decomposes into:
*   H = 2μ * H_ARAP + λ * H_vol
*
* ARAP part:
*   H_ARAP = I₉ - Σₖ λₖ (tₖ tₖᵀ)
*   λₖ = 2/(σₐ+σᵦ), clamped to 1 when σₐ+σᵦ < 2
*
* PSD projection: eigendecompose 9x9, clamp negative eigenvalues to 0.
*/
void Elasticity_Kee2023::computeCorotatedHessian9x9(ElasticObject* obj)
{
	const std::vector<unsigned int>& group = obj->m_particleIndices;
	const int numParticles = (int)group.size();
	const Real dt = obj->m_factorization->m_dt;
	const Real dt2 = dt * dt;

	const Real sqrt2inv = static_cast<Real>(1.0) / std::sqrt(static_cast<Real>(2.0));

	#pragma omp parallel for schedule(static)
	for (int i = 0; i < numParticles; i++)
	{
		const unsigned int i0 = group[i];
		const unsigned int particleIndex = m_initial_to_current_index[i0];

		const Matrix3r& Fi = m_F[particleIndex];

		// SVD: F = U * Σ * V^T
		Eigen::JacobiSVD<Matrix3r> svd_h(Fi, Eigen::ComputeFullU | Eigen::ComputeFullV);
		Matrix3r U = svd_h.matrixU();
		Matrix3r V = svd_h.matrixV();
		Vector3r sigma = svd_h.singularValues();
		// Ensure proper rotation (det = +1)
		if ((U * V.transpose()).determinant() < 0)
		{
			U.col(2) = -U.col(2);
			sigma[2] = -sigma[2];
		}

		// ======== ARAP Hessian: H = Σₖ ⟨λₖ⟩₊ qₖqₖᵀ ========
		// Eigenvalues:
		//   twist modes:   λ₀ = 1 - 2/(σ₀+σ₁), λ₁ = 1 - 2/(σ₁+σ₂), λ₂ = 1 - 2/(σ₀+σ₂)
		//   flip/scaling:  λ₃..λ₈ = 1
		Real lambda[3];
		lambda[0] = 1 - 2 / (sigma[0] + sigma[1]);
		lambda[1] = 1 - 2 / (sigma[1] + sigma[2]);
		lambda[2] = 1 - 2 / (sigma[0] + sigma[2]);

		// Clamp negative eigenvalues to 0 (PSD projection)
		for (int k = 0; k < 3; k++)
			if (lambda[k] < 0) lambda[k] = 0;

		// Eigenvectors as vec(uᵢvⱼᵀ ± uⱼvᵢᵀ) / sqrt(2)
		// Twist modes (skew-symmetric)
		Matrix3r T0 = sqrt2inv * (U.col(0) * V.col(1).transpose() - U.col(1) * V.col(0).transpose());
		Matrix3r T1 = sqrt2inv * (U.col(1) * V.col(2).transpose() - U.col(2) * V.col(1).transpose());
		Matrix3r T2 = sqrt2inv * (U.col(0) * V.col(2).transpose() - U.col(2) * V.col(0).transpose());

		// Flip modes (symmetric off-diagonal)
		Matrix3r F0 = sqrt2inv * (U.col(0) * V.col(1).transpose() + U.col(1) * V.col(0).transpose());
		Matrix3r F1 = sqrt2inv * (U.col(1) * V.col(2).transpose() + U.col(2) * V.col(1).transpose());
		Matrix3r F2 = sqrt2inv * (U.col(0) * V.col(2).transpose() + U.col(2) * V.col(0).transpose());

		// Scaling modes (diagonal)
		Matrix3r S0 = U.col(0) * V.col(0).transpose();
		Matrix3r S1 = U.col(1) * V.col(1).transpose();
		Matrix3r S2 = U.col(2) * V.col(2).transpose();

		// Vectorize (column-major)
		auto vecMap = [](const Matrix3r& M) {
			Eigen::Matrix<Real, 9, 1> v;
			v << M(0,0), M(1,0), M(2,0), M(0,1), M(1,1), M(2,1), M(0,2), M(1,2), M(2,2);
			return v;
		};

		Eigen::Matrix<Real, 9, 1> q[9];
		q[0] = vecMap(T0); q[1] = vecMap(T1); q[2] = vecMap(T2);
		q[3] = vecMap(F0); q[4] = vecMap(F1); q[5] = vecMap(F2);
		q[6] = vecMap(S0); q[7] = vecMap(S1); q[8] = vecMap(S2);

		// Reconstruct ARAP: H = 2μ Σₖ λₖ qₖqₖᵀ  (flip/scaling have λ=1)
		const Real mu2 = static_cast<Real>(2.0) * m_mu;
		Eigen::Matrix<Real, 9, 9> K = Eigen::Matrix<Real, 9, 9>::Zero();
		for (int k = 0; k < 3; k++)
			K += (mu2 * lambda[k]) * q[k] * q[k].transpose();
		for (int k = 3; k < 9; k++)
			K += mu2 * q[k] * q[k].transpose();

		// ======== Volume Hessian: λ * vec(R) * vec(R)^T ========
		// d²/dvecF² of (λ/2)(tr(R^T F) - 3)² = λ * R ⊗ R  (always PSD)
		const Matrix3r Ri = U * V.transpose();
		Eigen::Matrix<Real, 9, 1> vecR;
		vecR << Ri(0,0), Ri(1,0), Ri(2,0), Ri(0,1), Ri(1,1), Ri(2,1), Ri(0,2), Ri(1,2), Ri(2,2);
		K += m_lambda * vecR * vecR.transpose();

		// Store unscaled material Hessian
		// The dt² * Vi scaling is applied in the symmetric matvec
		obj->m_hessian9x9[i] = K;
	}
}

/** Compute the per-particle 9x9 Hessian for Stable Neo-Hookean (Smith et al. 2018).
*
* Uses iARAP decomposition (Lin et al. 2022) for U, V, sigma.
*
* Energy: ψ = μ/2 (I_C - 3) + λ/2 (J - α)² - μ/2 log(I_C + 1)
*   where α = 1 + μ/λ - μ/(4λ)
*
* Hessian decomposes into 4 terms (Eq. 22):
*   A = μ_T I + μ_F f⊗f + λ g⊗g + γ_H H_vol
*
* Direct approach: add each term separately, use volume eigendecomposition for H_vol.
*/
void Elasticity_Kee2023::computeStableNeoHookeanHessian9x9(ElasticObject* obj)
{
	const std::vector<unsigned int>& group = obj->m_particleIndices;
	const int numParticles = (int)group.size();
	const Real sqrt2inv = static_cast<Real>(1.0) / std::sqrt(static_cast<Real>(2.0));
	const Real eps_tol = static_cast<Real>(1e-12);

	#pragma omp parallel for schedule(static)
	for (int i = 0; i < numParticles; i++)
	{
		const unsigned int i0 = group[i];
		const unsigned int particleIndex = m_initial_to_current_index[i0];
		const Matrix3r& Fi = m_F[particleIndex];

		// SVD: F = U * Σ * V^T
		Eigen::JacobiSVD<Matrix3r> svd_h(Fi, Eigen::ComputeFullU | Eigen::ComputeFullV);
		Matrix3r U = svd_h.matrixU();
		Matrix3r V = svd_h.matrixV();
		Vector3r sigma = svd_h.singularValues();
		// Ensure proper rotation (det = +1)
		if ((U * V.transpose()).determinant() < 0)
		{
			U.col(2) = -U.col(2);
			sigma[2] = -sigma[2];
		}

		const Real s0 = sigma[0], s1 = sigma[1], s2 = sigma[2];
		const Real I_C = s0*s0 + s1*s1 + s2*s2;
		const Real J = s0 * s1 * s2;

		// Scalars
		const Real IC1 = I_C + 1;
		const Real mu_T = m_mu * (1 - 1 / IC1);           // Tikhonov coefficient
		const Real mu_F = 2 * m_mu / (IC1 * IC1);         // F rank-1 coefficient
		const Real gamma_H = m_lambda * (J - 1) - static_cast<Real>(0.75) * m_mu;  // volume coefficient

		// Vectorize helper
		auto vecMap = [](const Matrix3r& M) {
			Eigen::Matrix<Real, 9, 1> v;
			v << M(0,0), M(1,0), M(2,0), M(0,1), M(1,1), M(2,1), M(0,2), M(1,2), M(2,2);
			return v;
		};

		// ======== Term 1: μ_T * I₉ (always PSD) ========
		Eigen::Matrix<Real, 9, 9> K = mu_T * Eigen::Matrix<Real, 9, 9>::Identity();

		// ======== Term 2: μ_F * vec(F) ⊗ vec(F) (rank-1, always PSD) ========
		Eigen::Matrix<Real, 9, 1> vecF = vecMap(Fi);
		K += mu_F * vecF * vecF.transpose();

		// ======== Term 3: λ * vec(cof) ⊗ vec(cof) (rank-1, always PSD) ========
		Vector3r c0 = Fi.col(1).cross(Fi.col(2));
		Vector3r c1 = Fi.col(2).cross(Fi.col(0));
		Vector3r c2 = Fi.col(0).cross(Fi.col(1));
		Eigen::Matrix<Real, 9, 1> vecCof;
		vecCof << c0[0], c0[1], c0[2], c1[0], c1[1], c1[2], c2[0], c2[1], c2[2];
		K += m_lambda * vecCof * vecCof.transpose();

		// ======== Term 4: γ_H * H_vol (volume Hessian with PSD projection) ========
		// Same eigendecomposition as Co-rotated volume term

		// Group 1: 6 eigenvalues ±σᵢ
		Real volEig[6] = { s0, -s0, s1, -s1, s2, -s2 };
		Matrix3r volD[6];
		volD[0] << 0,0,0, 0,0,1, 0,-1,0;
		volD[1] << 0,0,0, 0,0,1, 0,1,0;
		volD[2] << 0,0,1, 0,0,0, -1,0,0;
		volD[3] << 0,0,1, 0,0,0, 1,0,0;
		volD[4] << 0,1,0, -1,0,0, 0,0,0;
		volD[5] << 0,1,0, 1,0,0, 0,0,0;

		for (int k = 0; k < 6; k++)
		{
			Real scaled = gamma_H * volEig[k];
			if (scaled > 0)
			{
				Matrix3r Q = sqrt2inv * U * volD[k] * V.transpose();
				Eigen::Matrix<Real, 9, 1> qv = vecMap(Q);
				K += scaled * qv * qv.transpose();
			}
		}

		// Group 2: 3 eigenvalues from depressed cubic
		if (I_C > eps_tol)
		{
			const Real sqrtIC3 = std::sqrt(I_C / static_cast<Real>(3.0));
			const Real cosArg = static_cast<Real>(3.0) * J / I_C * std::sqrt(static_cast<Real>(3.0) / I_C);
			const Real clampedCos = std::max(static_cast<Real>(-1.0), std::min(static_cast<Real>(1.0), cosArg));
			const Real theta = std::acos(clampedCos);
			const Real pi = static_cast<Real>(3.14159265358979323846);

			for (int k = 0; k < 3; k++)
			{
				Real eps_k = static_cast<Real>(2.0) * sqrtIC3 * std::cos((theta + static_cast<Real>(2.0) * pi * k) / static_cast<Real>(3.0));
				Real scaled = gamma_H * eps_k;
				if (scaled > 0)
				{
					Vector3r diag_k(s0*s2 + s1*eps_k, s1*s2 + s0*eps_k, eps_k*eps_k - s2*s2);
					Real norm_k = diag_k.norm();
					if (norm_k > eps_tol)
					{
						diag_k /= norm_k;
						Matrix3r Q = U * diag_k.asDiagonal() * V.transpose();
						Eigen::Matrix<Real, 9, 1> qv = vecMap(Q);
						K += scaled * qv * qv.transpose();
					}
				}
			}
		}

		obj->m_hessian9x9[i] = K;
	}
}

/** Assemble the 3N×3N Newton system matrix and Cholesky factorize.
*
* A[3j+b, 3k+c] = m_j δ_{jk}δ_{bc}
*                + dt² Σᵢ Vᵢ Σ_{a,a'} D[3i+a,j] H_i[3b+a,3c+a'] D[3i+a',k]
*                + HTH[j,k] δ_{bc}
*
* where H_i is the per-particle 9×9 Hessian, D is the deformation gradient operator,
* and HTH is the zero-energy mode control matrix.
*/
/** Compute block-diagonal preconditioner for Newton PCG.
*   Extracts and inverts 3×3 diagonal blocks of A = M + D^T·K·D + H^T·K_ze·H.
*   Note: K already includes dt² * Vi (pre-multiplied in computeCorotatedHessian9x9)
*/
void Elasticity_Kee2023::computeNewtonPreconditioner(ElasticObject* obj)
{
	const std::vector<unsigned int>& group = obj->m_particleIndices;
	const int numParticles = (int)group.size();
	const int nFree = numParticles - (int)obj->m_nFixed;

	auto& D = obj->m_factorization->m_D;
	auto& HTH = obj->m_factorization->m_matHTH;
	auto& precond = obj->m_pcg_precond;

	// Initialize diagonal blocks with mass
	#pragma omp parallel for schedule(static)
	for (int j = 0; j < nFree; j++)
	{
		const unsigned int j0 = group[j];
		const unsigned int pj = m_initial_to_current_index[j0];
		const Real mj = m_model->getMass(pj);
		precond[j] = Matrix3r::Identity() * mj;
	}

	// Add HTH diagonal (zero-energy mode control)
	if (m_alpha != static_cast<Real>(0.0))
	{
		for (int outer = 0; outer < HTH.outerSize(); outer++)
			for (Eigen::SparseMatrix<Real, Eigen::ColMajor>::InnerIterator it(HTH, outer); it; ++it)
			{
				const int j = (int)it.row();
				if (j >= nFree) continue;
				if (it.row() == it.col())  // diagonal only
				{
					for (int b = 0; b < 3; b++)
						precond[j](b, b) += it.value();
				}
			}
	}

	// Add elastic contribution: for each particle i, accumulate dt² * V_i * d_j^T K_i d_j to the j-th diagonal block
	// K_i is the unscaled material Hessian - dt² * V_i scaling is applied here
	const Real dt = obj->m_factorization->m_dt;
	const Real dt2 = dt * dt;

	for (int i = 0; i < numParticles; i++)
	{
		const unsigned int i0 = group[i];
		const unsigned int particleIndex = m_initial_to_current_index[i0];
		const Real V_i = m_restVolumes[particleIndex];
		const Real scale = dt2 * V_i;  // Scale factor for this particle's contribution

		const Eigen::Matrix<Real, 9, 9>& Ki = obj->m_hessian9x9[i];

		// Gather D values for particle i
		int nCols = 0;
		for (Eigen::SparseMatrix<Real, Eigen::RowMajor>::InnerIterator it(D, 3 * i); it; ++it)
			nCols++;

		std::vector<int> cols(nCols);
		std::vector<Real> dv0(nCols), dv1(nCols), dv2(nCols);

		int idx = 0;
		for (Eigen::SparseMatrix<Real, Eigen::RowMajor>::InnerIterator it(D, 3 * i); it; ++it)
		{ cols[idx] = (int)it.col(); dv0[idx] = it.value(); idx++; }
		
		idx = 0;
		for (Eigen::SparseMatrix<Real, Eigen::RowMajor>::InnerIterator it(D, 3 * i + 1); it; ++it)
		{ dv1[idx] = it.value(); idx++; }
		
		idx = 0;
		for (Eigen::SparseMatrix<Real, Eigen::RowMajor>::InnerIterator it(D, 3 * i + 2); it; ++it)
		{ dv2[idx] = it.value(); idx++; }

		// For each neighbor j, compute contribution to precond[j]
		for (int ji = 0; ji < nCols; ji++)
		{
			const int j = cols[ji];
			if (j >= nFree) continue;

			// d_j is a 9-vector: d_j[3b+a] = D[3i+b, j] * delta(a matches column index)
			// Simplified: d_j^T K_i d_j is a 3x3 block
			const Real dv[3] = { dv0[ji], dv1[ji], dv2[ji] };

			// Compute 3x3 block: dt² * V_i * Σ_{b,c} dv[b] * K_i[3b+a, 3c+a'] * dv[c] for each (a,a')
			for (int a = 0; a < 3; a++)
				for (int ap = 0; ap < 3; ap++)
				{
					Real val = 0;
					for (int b = 0; b < 3; b++)
						for (int c = 0; c < 3; c++)
							val += dv[b] * Ki(3 * b + a, 3 * c + ap) * dv[c];
					precond[j](a, ap) += scale * val;
				}
		}
	}

	// Invert each 3x3 block
	#pragma omp parallel for schedule(static)
	for (int j = 0; j < nFree; j++)
		precond[j] = precond[j].inverse();
}

#ifdef USE_AVX
/** AVX Newton matvec: AVX for D*p and HTH*p sparse matvecs only.
*   K*Dx and force loop stay scalar (CG-sensitive, proven correct scalar).
*   Pack/unpack at boundaries.
*/
void Elasticity_Kee2023::newtonMatvec(ElasticObject* obj)
{
	const std::vector<unsigned int>& group = obj->m_particleIndices;
	const int numParticles = (int)group.size();
	const int nFree = numParticles - (int)obj->m_nFixed;
	auto& p_avx = obj->m_pcg_p_avx;
	auto& Ap_avx = obj->m_pcg_Ap_avx;

	const Real dt = obj->m_factorization->m_dt;
	const Real dt2 = dt * dt;
	auto& D = obj->m_factorization->m_D;
	auto& HTH = obj->m_factorization->m_matHTH;
	auto& f_avx = obj->m_f_avx;

	// Step 1: Ap_avx = M * p_avx
	#pragma omp parallel for schedule(static)
	for (int j = 0; j < nFree; j++)
	{
		const unsigned int j0 = group[j];
		const unsigned int pj = m_initial_to_current_index[j0];
		Ap_avx[j] = Scalarf8((float)m_model->getMass(pj)) * p_avx[j];
	}

	// Step 2: Ap_avx += HTH * p_avx
	if (m_alpha != static_cast<Real>(0.0))
	{
		#pragma omp parallel default(shared)
		{
			#pragma omp for schedule(static)
			for (int k = 0; k < HTH.outerSize(); ++k)
			{
				for (Eigen::SparseMatrix<Real, Eigen::ColMajor>::InnerIterator it(HTH, k); it; ++it)
				{
					const int row = (int)it.row();
					const int col = (int)it.col();
					if (row < nFree && col < nFree)
						Ap_avx[col] = Ap_avx[col] + Scalarf8((float)it.value()) * p_avx[row];
				}
			}
		}
	}

	// Step 3: Dx = D * p  (coord-packed AVX sparse matvec)
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)
		for (int i = 0; i < 3 * numParticles; i++)
			f_avx[i].setZero();

		#pragma omp for schedule(static)
		for (int k = 0; k < D.outerSize(); ++k)
		{
			for (Eigen::SparseMatrix<Real, Eigen::RowMajor>::InnerIterator it(D, k); it; ++it)
			{
				const int col = (int)it.col();
				if (col < nFree)
					f_avx[it.row()] = f_avx[it.row()] + Scalarf8((float)it.value()) * p_avx[col];
			}
		}
	}

	// Step 4: K*Dx → dPtL (scalar — 9×9 Hessian multiply)
	const unsigned int numActiveParticles = m_model->numActiveParticles();
	std::vector<Matrix3r, Eigen::aligned_allocator<Matrix3r>> dPtL(numActiveParticles, Matrix3r::Zero());

	#pragma omp parallel for schedule(static)
	for (int i = 0; i < numParticles; i++)
	{
		const unsigned int i0 = group[i];
		const unsigned int particleIndex = m_initial_to_current_index[i0];

		float x0[8], x1[8], x2[8];
		f_avx[3 * i].store(x0);
		f_avx[3 * i + 1].store(x1);
		f_avx[3 * i + 2].store(x2);

		Eigen::Matrix<Real, 9, 1> Dxi;
		Dxi << (Real)x0[0], (Real)x1[0], (Real)x2[0],
		       (Real)x0[1], (Real)x1[1], (Real)x2[1],
		       (Real)x0[2], (Real)x1[2], (Real)x2[2];

		Eigen::Matrix<Real, 9, 1> KDx_vec = obj->m_hessian9x9[i] * Dxi;
		Matrix3r dP;
		dP.col(0) = KDx_vec.segment<3>(0);
		dP.col(1) = KDx_vec.segment<3>(3);
		dP.col(2) = KDx_vec.segment<3>(6);
		dPtL[particleIndex] = dP.transpose() * m_L[particleIndex];
	}

	// Step 5: Force loop (AVX gather with precomputed V_gradW8), accumulate into Ap_avx
	#pragma omp parallel for schedule(static)
	for (int i = 0; i < nFree; i++)
	{
		const unsigned int i0 = group[i];
		const unsigned int particleIndex = m_initial_to_current_index[i0];
		const float V_i = (float)m_restVolumes[particleIndex];

		const unsigned int numNeighbors = (unsigned int)m_initialNeighbors[i0].size();
		const Matrix3f8 dPtLi_avx(dPtL[particleIndex].cast<float>());
		Vector3f8 force_avx;
		force_avx.setZero();

		for (unsigned int j = 0; j < numNeighbors; j += 8)
		{
			const unsigned int count = std::min(numNeighbors - j, 8u);
			unsigned int nIndices[8];
			for (auto k = 0u; k < count; k++)
				nIndices[k] = m_initial_to_current_index[m_initialNeighbors[i0][j + k]];

			const Matrix3f8 dPtLj_avx = convertMat_zero(nIndices, &dPtL[0], count);
			const Vector3f8& V_gradW = m_precomp_V_gradW8[m_precomputed_indices8[particleIndex] + j / 8];
			force_avx += dPtLi_avx * V_gradW + dPtLj_avx * V_gradW;
		}

		const Scalarf8 fx(force_avx.x().reduce(), force_avx.y().reduce(), force_avx.z().reduce(), 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
		Ap_avx[i] = Ap_avx[i] - Scalarf8((float)dt2 * V_i) * fx;
	}

}
#else
/** Compute A*x = M*x + H^T*K_ze*H*x + dt²*D^T*K*D*x for Newton PCG. */
void Elasticity_Kee2023::newtonMatvec(ElasticObject* obj)
{
	const std::vector<unsigned int>& group = obj->m_particleIndices;
	const int numParticles = (int)group.size();
	const int nFree = numParticles - (int)obj->m_nFixed;
	auto& x = obj->m_pcg_p;
	auto& Ax = obj->m_pcg_Ap;

	const Real dt = obj->m_factorization->m_dt;
	const Real dt2 = dt * dt;
	auto& D = obj->m_factorization->m_D;
	auto& HTH = obj->m_factorization->m_matHTH;

	// Initialize Ax = M*x
	#pragma omp parallel for schedule(static)
	for (int j = 0; j < nFree; j++)
	{
		const unsigned int j0 = group[j];
		const unsigned int pj = m_initial_to_current_index[j0];
		const Real mj = m_model->getMass(pj);
		Ax[j] = mj * x[j];
	}

	// Add H^T*K_ze*H*x (zero-energy mode control)
	// Col-write pattern: Ax[col] += val * x[row] (HTH is symmetric, avoids data race)
	if (m_alpha != static_cast<Real>(0.0))
	{
		#pragma omp parallel default(shared)
		{
			#pragma omp for schedule(static)
			for (int k = 0; k < HTH.outerSize(); ++k)
			{
				for (Eigen::SparseMatrix<Real, Eigen::ColMajor>::InnerIterator it(HTH, k); it; ++it)
				{
					const int row = (int)it.row();
					const int col = (int)it.col();
					if (row < nFree && col < nFree)
						Ax[col] += it.value() * x[row];
				}
			}
		}
	}

	// D * x
	std::vector<Eigen::Matrix<Real, 9, 1>> Dx(numParticles);
	#pragma omp parallel for schedule(static)
	for (int i = 0; i < numParticles; i++)
	{
		Eigen::Matrix<Real, 9, 1> yi = Eigen::Matrix<Real, 9, 1>::Zero();
		for (int a = 0; a < 3; a++)
		{
			for (Eigen::SparseMatrix<Real, Eigen::RowMajor>::InnerIterator it(D, 3 * i + a); it; ++it)
			{
				const int j = (int)it.col();
				if (j < nFree)
				{
					for (int b = 0; b < 3; b++)
						yi(3 * b + a) += it.value() * x[j][b];
				}
			}
		}
		Dx[i] = yi;
	}

	// K * Dx → dPtL
	const unsigned int numActiveParticles = m_model->numActiveParticles();
	std::vector<Matrix3r, Eigen::aligned_allocator<Matrix3r>> dPtL(numActiveParticles, Matrix3r::Zero());

	#pragma omp parallel for schedule(static)
	for (int i = 0; i < numParticles; i++)
	{
		const unsigned int i0 = group[i];
		const unsigned int particleIndex = m_initial_to_current_index[i0];
		Eigen::Matrix<Real, 9, 1> KDx_vec = obj->m_hessian9x9[i] * Dx[i];
		Matrix3r dP;
		dP.col(0) = KDx_vec.segment<3>(0);
		dP.col(1) = KDx_vec.segment<3>(3);
		dP.col(2) = KDx_vec.segment<3>(6);
		dPtL[particleIndex] = dP.transpose() * m_L[particleIndex];
	}

	// Force loop
	Simulation* sim = Simulation::getCurrent();
	#pragma omp parallel for schedule(static)
	for (int i = 0; i < nFree; i++)
	{
		const unsigned int i0 = group[i];
		const unsigned int particleIndex = m_initial_to_current_index[i0];
		const Vector3r xi0 = m_model->getPosition0(i0);
		const Real V_i = m_restVolumes[particleIndex];
		const Matrix3r& dPtL_i = dPtL[particleIndex];

		const size_t numNeighbors = m_initialNeighbors[i0].size();
		Vector3r force;
		force.setZero();

		for (unsigned int jn = 0; jn < numNeighbors; jn++)
		{
			const unsigned int j0 = m_initialNeighbors[i0][jn];
			const unsigned int jCurrent = m_initial_to_current_index[j0];
			const Matrix3r& dPtL_j = dPtL[jCurrent];
			force += (dPtL_i + dPtL_j) * m_precomp_V_gradW[m_precomputed_indices[particleIndex] + jn];
		}

		Ax[i] -= dt2 * V_i * force;
	}
}
#endif

/** Solve A*dx = -gradient using Preconditioned Conjugate Gradient.
*   Uses block-diagonal preconditioner (3×3 blocks).
*   Writes result into obj->m_dx_avx / obj->m_dx.
*/
#ifdef USE_AVX
int Elasticity_Kee2023::matFreePCG(ElasticObject* obj)
{
	computeNewtonPreconditioner(obj);

	const int nFree = (int)obj->m_dx_avx.size();

	auto& gradient = obj->m_gradient_avx;
	auto& dx = obj->m_dx_avx;
	auto& r = obj->m_pcg_r_avx;
	auto& p = obj->m_pcg_p_avx;
	auto& Ap = obj->m_pcg_Ap_avx;
	auto& z = obj->m_pcg_z_avx;
	auto& precond = obj->m_pcg_precond;

	// Initialize: x_0 = 0, r_0 = b = -gradient
	#pragma omp parallel for schedule(static)
	for (int i = 0; i < nFree; i++)
	{
		dx[i].setZero();
		r[i] = Scalarf8(0.0f) - gradient[i];
	}

	// z_0 = M^{-1} r_0 (apply block-diagonal preconditioner)
	#pragma omp parallel for schedule(static)
	for (int i = 0; i < nFree; i++)
	{
		float rv[8]; r[i].store(rv);
		Vector3r ri((Real)rv[0], (Real)rv[1], (Real)rv[2]);
		Vector3r zi = precond[i] * ri;
		z[i] = Scalarf8((float)zi[0], (float)zi[1], (float)zi[2], 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
	}

	// p_0 = z_0
	#pragma omp parallel for schedule(static)
	for (int i = 0; i < nFree; i++)
		p[i] = z[i];

	// rz_old = r_0^T z_0  (dot product via Scalarf8: (r*z).reduce() sums x*x + y*y + z*z + 0+...)
	Real rz_old = 0;
	#pragma omp parallel for reduction(+:rz_old) schedule(static)
	for (int i = 0; i < nFree; i++)
		rz_old += (Real)(r[i] * z[i]).reduce();

	// Compute initial residual norm for convergence check
	Real r0_norm = 0;
	#pragma omp parallel for reduction(+:r0_norm) schedule(static)
	for (int i = 0; i < nFree; i++)
		r0_norm += (Real)(r[i] * r[i]).reduce();
	r0_norm = std::sqrt(r0_norm);

	const Real tol = m_tolCG * r0_norm;

	int iter = 0;
	for (; iter < m_maxIterCG; iter++)
	{
		// Ap = A * p
		newtonMatvec(obj);

		// alpha = rz_old / (p^T Ap)
		Real pAp = 0;
		#pragma omp parallel for reduction(+:pAp) schedule(static)
		for (int i = 0; i < nFree; i++)
			pAp += (Real)(p[i] * Ap[i]).reduce();

		if (pAp <= static_cast<Real>(0))
		{
			LOG_ERR << "PCG: pAp <= 0 at iter " << iter << ", pAp = " << pAp
			        << " (matrix not positive definite)";
			break;
		}

		const Scalarf8 alpha_avx((float)(rz_old / pAp));

		// x_{k+1} = x_k + alpha * p_k
		// r_{k+1} = r_k - alpha * Ap
		#pragma omp parallel for schedule(static)
		for (int i = 0; i < nFree; i++)
		{
			dx[i] = dx[i] + alpha_avx * p[i];
			r[i] = r[i] - alpha_avx * Ap[i];
		}

		// Check convergence
		Real r_norm = 0;
		#pragma omp parallel for reduction(+:r_norm) schedule(static)
		for (int i = 0; i < nFree; i++)
			r_norm += (Real)(r[i] * r[i]).reduce();
		r_norm = std::sqrt(r_norm);

		if (r_norm < tol)
		{
			iter++;  // count this iteration
			break;
		}

		// z_{k+1} = M^{-1} r_{k+1} (apply block-diagonal preconditioner)
		#pragma omp parallel for schedule(static)
		for (int i = 0; i < nFree; i++)
		{
			float rv[8]; r[i].store(rv);
			Vector3r ri((Real)rv[0], (Real)rv[1], (Real)rv[2]);
			Vector3r zi = precond[i] * ri;
			z[i] = Scalarf8((float)zi[0], (float)zi[1], (float)zi[2], 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
		}

		// rz_new = r_{k+1}^T z_{k+1}
		Real rz_new = 0;
		#pragma omp parallel for reduction(+:rz_new) schedule(static)
		for (int i = 0; i < nFree; i++)
			rz_new += (Real)(r[i] * z[i]).reduce();

		// beta = rz_new / rz_old
		const Scalarf8 beta_avx((float)(rz_new / rz_old));

		// p_{k+1} = z_{k+1} + beta * p_k
		#pragma omp parallel for schedule(static)
		for (int i = 0; i < nFree; i++)
			p[i] = z[i] + beta_avx * p[i];

		rz_old = rz_new;
	}

	return iter;
}
#else
int Elasticity_Kee2023::matFreePCG(ElasticObject* obj)
{
	computeNewtonPreconditioner(obj);

	const int nFree = (int)obj->m_dx.size();

	auto& gradient = obj->m_gradient;
	auto& dx = obj->m_dx;
	auto& r = obj->m_pcg_r;
	auto& p = obj->m_pcg_p;
	auto& Ap = obj->m_pcg_Ap;
	auto& z = obj->m_pcg_z;
	auto& precond = obj->m_pcg_precond;

	// Initialize: x_0 = 0, r_0 = b = -gradient
	#pragma omp parallel for schedule(static)
	for (int i = 0; i < nFree; i++)
	{
		dx[i].setZero();
		r[i] = -gradient[i];
	}

	// z_0 = M^{-1} r_0 (apply block-diagonal preconditioner)
	#pragma omp parallel for schedule(static)
	for (int i = 0; i < nFree; i++)
		z[i] = precond[i] * r[i];

	// p_0 = z_0
	#pragma omp parallel for schedule(static)
	for (int i = 0; i < nFree; i++)
		p[i] = z[i];

	// rz_old = r_0^T z_0
	Real rz_old = 0;
	#pragma omp parallel for reduction(+:rz_old) schedule(static)
	for (int i = 0; i < nFree; i++)
		rz_old += r[i].dot(z[i]);

	// Compute initial residual norm for convergence check
	Real r0_norm = 0;
	#pragma omp parallel for reduction(+:r0_norm) schedule(static)
	for (int i = 0; i < nFree; i++)
		r0_norm += r[i].squaredNorm();
	r0_norm = std::sqrt(r0_norm);

	const Real tol = m_tolCG * r0_norm;

	int iter = 0;
	for (; iter < m_maxIterCG; iter++)
	{
		// Ap = A * p
		newtonMatvec(obj);

		// alpha = rz_old / (p^T Ap)
		Real pAp = 0;
		#pragma omp parallel for reduction(+:pAp) schedule(static)
		for (int i = 0; i < nFree; i++)
			pAp += p[i].dot(Ap[i]);

		if (pAp <= static_cast<Real>(0))
		{
			LOG_ERR << "PCG: pAp <= 0 at iter " << iter << ", pAp = " << pAp
			        << " (matrix not positive definite)";
			break;
		}

		const Real alpha = rz_old / pAp;

		// x_{k+1} = x_k + alpha * p_k
		// r_{k+1} = r_k - alpha * Ap
		#pragma omp parallel for schedule(static)
		for (int i = 0; i < nFree; i++)
		{
			dx[i] += alpha * p[i];
			r[i] -= alpha * Ap[i];
		}

		// Check convergence
		Real r_norm = 0;
		#pragma omp parallel for reduction(+:r_norm) schedule(static)
		for (int i = 0; i < nFree; i++)
			r_norm += r[i].squaredNorm();
		r_norm = std::sqrt(r_norm);

		if (r_norm < tol)
		{
			iter++;  // count this iteration
			break;
		}

		// z_{k+1} = M^{-1} r_{k+1} (apply block-diagonal preconditioner)
		#pragma omp parallel for schedule(static)
		for (int i = 0; i < nFree; i++)
			z[i] = precond[i] * r[i];

		// rz_new = r_{k+1}^T z_{k+1}
		Real rz_new = 0;
		#pragma omp parallel for reduction(+:rz_new) schedule(static)
		for (int i = 0; i < nFree; i++)
			rz_new += r[i].dot(z[i]);

		// beta = rz_new / rz_old
		const Real beta = rz_new / rz_old;

		// p_{k+1} = z_{k+1} + beta * p_k
		#pragma omp parallel for schedule(static)
		for (int i = 0; i < nFree; i++)
			p[i] = z[i] + beta * p[i];

		rz_old = rz_new;
	}

	return iter;
}
#endif

/** Solve the linear system A * x = b using the prefactored Cholesky decomposition.
*
* Input: dx contains b (right-hand side)
* Output: dx contains the solution
*
* Forward substitution:  L * y  = P * b
* Backward substitution: L^T * z = y
* Inverse permutation:   dx = P^T * z
*/
void Elasticity_Kee2023::prefactorizedLLTSolve(ElasticObject* obj)
{
#ifdef USE_AVX
	auto& dx = obj->m_dx_avx;
#else
	auto& dx = obj->m_dx;
#endif
	const int n = (int)dx.size();

#ifdef USE_AVX
	// dx is std::vector<Scalarf8> with 3 active lanes per element (x, y, z, 0, 0, 0, 0, 0).
	// CholeskyAVXSolver::solve takes float* with stride, so we pack/unpack at this boundary.
	VectorXr b(3 * n), x(3 * n);

	#pragma omp parallel for schedule(static)
	for (int i = 0; i < n; i++)
	{
		float vals[8];
		dx[i].store(vals);
		b[3 * i]     = vals[0];
		b[3 * i + 1] = vals[1];
		b[3 * i + 2] = vals[2];
	}

	// Solve 3 independent n-systems (one per coordinate, stride 3)
	#pragma omp parallel for schedule(static)
	for (int c = 0; c < 3; c++)
		obj->m_factorization->m_cholesky->solve(&x[c], &b[c], 3);

	#pragma omp parallel for schedule(static)
	for (int i = 0; i < n; i++)
		dx[i] = Scalarf8(x[3 * i], x[3 * i + 1], x[3 * i + 2], 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
#else
	auto& dx_perm = obj->m_dx_perm;
	auto& permInd = obj->m_factorization->m_permInd;
	auto& permInvInd = obj->m_factorization->m_permInvInd;
	auto& matL = obj->m_factorization->m_matL;
	auto& matLT = obj->m_factorization->m_matLT;

	// Permute dx (fill-in reduction ordering)
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)
		for (int i = 0; i < n; i++)
			dx_perm[permInd[i]] = dx[i];
	}

	// Forward substitution: L * y = RHS_perm
	for (int k = 0; k < matL.outerSize(); ++k)
		for (Eigen::SparseMatrix<Real, Eigen::ColMajor>::InnerIterator it(matL, k); it; ++it)
			if (it.row() == it.col())
				dx_perm[it.row()] /= it.value();
			else
				dx_perm[it.row()] -= it.value() * dx_perm[it.col()];

	// Backward substitution: L^T * dx = y
	for (int k = (int) matLT.outerSize() - 1; k >= 0; --k)
		for (Eigen::SparseMatrix<Real, Eigen::ColMajor>::ReverseInnerIterator it(matLT, k); it; --it)
			if (it.row() == it.col())
				dx_perm[it.row()] /= it.value();
			else
				dx_perm[it.row()] -= it.value() * dx_perm[it.col()];

	// Inverse permutation
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)
		for (int i = 0; i < n; i++)
			dx[permInvInd[i]] = dx_perm[i];
	}
#endif
}

#ifdef USE_AVX
/** AVX Newton solve: energy+gradient+CG all in Scalarf8 domain.
*/
Real Elasticity_Kee2023::newtonSolve(ElasticObject* obj, int& cgIter)
{
	// Energy + gradient in AVX (fills gradient_avx, m_F, m_rotations, m_PL)
	Real energy = computeEnergyAndGradient(obj);

	// Hessian (scalar, reads m_F, writes m_hessian9x9)
	computeHessian(obj);

	// CG in AVX (reads gradient_avx, writes dx_avx)
	cgIter = matFreePCG(obj);

	return energy;
}
#else
/** Newton solve: compute energy+gradient, assemble true Hessian, solve 3N×3N system.
*   Returns dx in obj->m_dx, gradient in obj->m_gradient.
*   Returns energy at current xk.
*/
Real Elasticity_Kee2023::newtonSolve(ElasticObject* obj, int& cgIter)
{
	Real energy = computeEnergyAndGradient(obj);  // sets gradient, m_F, m_rotations
	computeHessian(obj);                          // computes 9×9 per particle, assembles 3N system
	cgIter = matFreePCG(obj);              // solves 3N scalar system, writes solution into dx

	// DEBUG: check descent direction and Hessian eigenvalues
	{
		const auto& gradient = obj->m_gradient;
		const auto& dx = obj->m_dx;  // dx
		const int n = (int)dx.size();

		// Gradient norm
		Real gradNorm = 0;
		for (int i = 0; i < n; i++)
			gradNorm += gradient[i].squaredNorm();
		gradNorm = std::sqrt(gradNorm);

		// Descent: g · dx (should be negative for descent)
		Real gDotDx = 0;
		for (int i = 0; i < n; i++)
			gDotDx += gradient[i].dot(dx[i]);

		// Sample 9x9 Hessian eigenvalues for particle 0
		if (!obj->m_hessian9x9.empty())
		{
			Eigen::SelfAdjointEigenSolver<Eigen::Matrix<Real, 9, 9>> eig(obj->m_hessian9x9[0]);
			const auto& evals = eig.eigenvalues();
			std::cout << "DEBUG Newton: E=" << energy << " |g|=" << gradNorm
			          << " g·dx=" << gDotDx << " H[0] eigs=[" << evals.transpose() << "]" << std::endl;
		}
		else
		{
			std::cout << "DEBUG Newton: E=" << energy << " |g|=" << gradNorm
			          << " g·dx=" << gDotDx << " (no Hessian)" << std::endl;
		}
	}

	return energy;
}
#endif

/** L-BFGS solve: quasi-Newton with prefactored Cholesky as H_0.
*   Uses two-loop recursion with secant pairs in a circular queue.
*   Returns dx in obj->m_dx, gradient in obj->m_gradient.
*   Returns energy at current xk.
*
*   Secant pairs are computed from actual position/gradient differences:
*     s_{k-1} = xk - xk_prev  (includes line search alpha)
*     y_{k-1} = gradient_current - gradient_prev  (via m_gradient)
*/
#ifdef USE_AVX
/** AVX L-BFGS step (full two-loop recursion).
*   Mirrors the scalar algorithm, operating on Scalarf8 vectors with 3 active
*   lanes (x, y, z, 0...0). Inner products use (a * b).reduce() which sums
*   all 8 lanes — lanes 3..7 are zero, so this gives the Vector3r dot product.
*   No line search (caller applies step size 1).
*
*   Sets dx = -H^-1 * gradient (descent direction), gradient in m_gradient_avx.
*   Returns energy at current xk.
*/
Real Elasticity_Kee2023::lbfgsSolve(ElasticObject* obj)
{
	auto& dx = obj->m_dx_avx;           // Scalarf8
	const int windowSize = m_lbfgsWindowSize;
	const int n = (int)dx.size();

	int& count = obj->m_lbfgs_count;
	auto& lbfgs_s = obj->m_lbfgs_s_avx;     // Scalarf8 secant history
	auto& lbfgs_y = obj->m_lbfgs_y_avx;     // Scalarf8 secant history
	auto& lbfgs_rho = obj->m_lbfgs_rho;         // Real
	auto& lbfgs_alpha = obj->m_lbfgs_alpha;     // Real
	auto& last_sol = obj->m_lbfgs_last_sol_avx;     // Scalarf8
	auto& gradient = obj->m_gradient_avx;           // Scalarf8
	auto& q = obj->m_lbfgs_q_avx;                   // Scalarf8

	// Save g_{k-1} into q before computeEnergyAndGradient overwrites gradient
	if (count > 0)
	{
		#pragma omp parallel for schedule(static)
		for (int i = 0; i < n; i++)
			q[i] = gradient[i];
	}

	Real energy = computeEnergyAndGradient(obj);  // gradient = g_k

	// Compute secant pairs from actual position/gradient differences
	if (count > 0 && windowSize > 0)
	{
		int slot = (count - 1) % windowSize;

		// s_{k-1} = x_k - x_{k-1}
		#pragma omp parallel for schedule(static)
		for (int i = 0; i < n; i++)
			lbfgs_s[slot][i] = obj->m_xk_avx[i] - last_sol[i];

		// y_{k-1} = g_k - g_{k-1}
		#pragma omp parallel for schedule(static)
		for (int i = 0; i < n; i++)
			lbfgs_y[slot][i] = gradient[i] - q[i];

		Real ys = 0;
		#pragma omp parallel reduction(+:ys)
		{
			Scalarf8 acc; acc.setZero();
			#pragma omp for schedule(static) nowait
			for (int i = 0; i < n; i++)
				acc += lbfgs_y[slot][i] * lbfgs_s[slot][i];
			ys += (Real)acc.reduce();
		}

		lbfgs_rho[slot] = (std::abs(ys) > static_cast<Real>(1e-10))
			? static_cast<Real>(1.0) / ys : 0;
	}

	// Save current xk for next iteration, set q = g_k for two-loop
	#pragma omp parallel for schedule(static)
	for (int i = 0; i < n; i++)
	{
		q[i] = gradient[i];
		last_sol[i] = obj->m_xk_avx[i];
	}

	int nPairs = std::min(count, windowSize);

	// Backward loop (newest -> oldest secant pair)
	for (int j = 0; j < nPairs; j++)
	{
		int idx = ((count - 1 - j) % windowSize + windowSize) % windowSize;

		Real sq = 0;
		#pragma omp parallel reduction(+:sq)
		{
			Scalarf8 acc; acc.setZero();
			#pragma omp for schedule(static) nowait
			for (int i = 0; i < n; i++)
				acc += lbfgs_s[idx][i] * q[i];
			sq += (Real)acc.reduce();
		}
		lbfgs_alpha[j] = lbfgs_rho[idx] * sq;

		const Scalarf8 alpha_avx((float)lbfgs_alpha[j]);
		#pragma omp parallel for schedule(static)
		for (int i = 0; i < n; i++)
			q[i] = q[i] - alpha_avx * lbfgs_y[idx][i];
	}

	// Apply initial inverse Hessian: dx = H_0^-1 * q
	#pragma omp parallel for schedule(static)
	for (int i = 0; i < n; i++)
		dx[i] = q[i];

	prefactorizedLLTSolve(obj);

	// Forward loop (oldest -> newest secant pair)
	for (int j = nPairs - 1; j >= 0; j--)
	{
		int idx = ((count - 1 - j) % windowSize + windowSize) % windowSize;

		Real yr = 0;
		#pragma omp parallel reduction(+:yr)
		{
			Scalarf8 acc; acc.setZero();
			#pragma omp for schedule(static) nowait
			for (int i = 0; i < n; i++)
				acc += lbfgs_y[idx][i] * dx[i];
			yr += (Real)acc.reduce();
		}
		Real beta = lbfgs_rho[idx] * yr;

		const Scalarf8 ab_avx((float)(lbfgs_alpha[j] - beta));
		#pragma omp parallel for schedule(static)
		for (int i = 0; i < n; i++)
			dx[i] = dx[i] + ab_avx * lbfgs_s[idx][i];
	}

	// Search direction: dx = -H*g
	const Scalarf8 negOne(-1.0f);
	#pragma omp parallel for schedule(static)
	for (int i = 0; i < n; i++)
		dx[i] = negOne * dx[i];

	count++;
	return energy;
}
#else
Real Elasticity_Kee2023::lbfgsSolve(ElasticObject* obj)
{
	auto& dx = obj->m_dx;
	const int windowSize = m_lbfgsWindowSize;
	const int n = (int)dx.size();

	int& count = obj->m_lbfgs_count;
	auto& lbfgs_s = obj->m_lbfgs_s;
	auto& lbfgs_y = obj->m_lbfgs_y;
	auto& lbfgs_rho = obj->m_lbfgs_rho;
	auto& lbfgs_alpha = obj->m_lbfgs_alpha;
	auto& last_sol = obj->m_lbfgs_last_sol;
	auto& gradient = obj->m_gradient;
	auto& q = obj->m_lbfgs_q;

	// Save g_{k-1} into q before computeEnergyAndGradient overwrites gradient
	if (count > 0)
	{
		#pragma omp parallel for schedule(static)
		for (int i = 0; i < n; i++)
			q[i] = gradient[i];
	}

	// Compute energy and gradient at current xk
	Real energy = computeEnergyAndGradient(obj);  // gradient = g_k

	// Compute secant pairs from actual position/gradient differences
	if (count > 0 && windowSize > 0)
	{
		int slot = (count - 1) % windowSize;

		// s_{k-1} = x_k - x_{k-1} (actual step taken, includes line search alpha)
		#pragma omp parallel for schedule(static)
		for (int i = 0; i < n; i++)
			lbfgs_s[slot][i] = obj->m_xk[i] - last_sol[i];

		// y_{k-1} = g_k - g_{k-1}  (q holds g_{k-1})
		#pragma omp parallel for schedule(static)
		for (int i = 0; i < n; i++)
			lbfgs_y[slot][i] = gradient[i] - q[i];

		Real ys = 0;
		#pragma omp parallel for reduction(+:ys) schedule(static)
		for (int i = 0; i < n; i++)
			ys += lbfgs_y[slot][i].dot(lbfgs_s[slot][i]);

		lbfgs_rho[slot] = (std::abs(ys) > static_cast<Real>(1e-10))
			? static_cast<Real>(1.0) / ys : 0;
	}

	// Save current xk for next iteration, set q = g_k for two-loop
	#pragma omp parallel for schedule(static)
	for (int i = 0; i < n; i++)
	{
		q[i] = gradient[i];
		last_sol[i] = obj->m_xk[i];
	}

	// Number of usable secant pairs
	int nPairs = std::min(count, windowSize);

	//////////////////////////////////////////////////////////////
	// L-BFGS two-loop recursion
	// q contains g_k (gradient)
	//////////////////////////////////////////////////////////////

	// Backward loop (newest to oldest secant pair)
	for (int j = 0; j < nPairs; j++)
	{
		int idx = ((count - 1 - j) % windowSize + windowSize) % windowSize;

		Real sq = 0;
		#pragma omp parallel for reduction(+:sq) schedule(static)
		for (int i = 0; i < n; i++)
			sq += lbfgs_s[idx][i].dot(q[i]);
		lbfgs_alpha[j] = lbfgs_rho[idx] * sq;

		#pragma omp parallel for schedule(static)
		for (int i = 0; i < n; i++)
			q[i] -= lbfgs_alpha[j] * lbfgs_y[idx][i];
	}

	// Apply initial inverse Hessian: r = A^{-1} * q
	#pragma omp parallel for schedule(static)
	for (int i = 0; i < n; i++)
		dx[i] = q[i];

	prefactorizedLLTSolve(obj);
	// dx now contains r = H_0^{-1} * q

	// Forward loop (oldest to newest secant pair)
	for (int j = nPairs - 1; j >= 0; j--)
	{
		int idx = ((count - 1 - j) % windowSize + windowSize) % windowSize;

		Real yr = 0;
		#pragma omp parallel for reduction(+:yr) schedule(static)
		for (int i = 0; i < n; i++)
			yr += lbfgs_y[idx][i].dot(dx[i]);
		Real beta = lbfgs_rho[idx] * yr;

		#pragma omp parallel for schedule(static)
		for (int i = 0; i < n; i++)
			dx[i] += (lbfgs_alpha[j] - beta) * lbfgs_s[idx][i];
	}

	// Search direction: dx = -H*g
	#pragma omp parallel for schedule(static)
	for (int i = 0; i < n; i++)
		dx[i] = -dx[i];

	count++;
	return energy;
}
#endif

#ifndef USE_AVX
/** Backtracking line search with Armijo condition.
*   Reads search direction dx from obj->m_dx, gradient from obj->m_gradient.
*   Returns the step size alpha (0 if line search failed).
*   Does NOT modify obj->m_xk — caller applies the step.
*/
Real Elasticity_Kee2023::lineSearch(ElasticObject* obj, Real energy, int& lsIter)
{
	lsIter = 0;
	auto& x = obj->m_xk;
	auto& dx = obj->m_dx;
	auto& gradient = obj->m_gradient;
	const int n = (int)dx.size();

	// Directional derivative: g^T * dx  (must be < 0 for descent)
	Real dirDeriv = 0;
	#pragma omp parallel for reduction(+:dirDeriv) schedule(static)
	for (int i = 0; i < n; i++)
		dirDeriv += gradient[i].dot(dx[i]);

	if (dirDeriv >= 0)
		LOG_WARN << "LS: dirDeriv=" << dirDeriv << " (not a descent direction)";

	// Backup current xk
	std::vector<Vector3r, Eigen::aligned_allocator<Vector3r>> xk_backup(n);
	#pragma omp parallel for schedule(static)
	for (int i = 0; i < n; i++)
		xk_backup[i] = x[i];

	Real alpha = static_cast<Real>(1.0);
	const Real eps = static_cast<Real>(1e-4);
	for (int ls = 0; ls < m_maxLSIter; ls++)
	{
		// Set trial point: xk = xk_backup + alpha * dx
		#pragma omp parallel for schedule(static)
		for (int i = 0; i < n; i++)
			x[i] = xk_backup[i] + alpha * dx[i];

		Real trialEnergy = computeEnergy(obj);
		lsIter++;

		// Armijo: f(x + alpha*dx) <= f(x) + c1*alpha*(g^T*dx)
		if (trialEnergy <= energy + m_lsArmijoParam * alpha * dirDeriv)
		{
			#pragma omp parallel for schedule(static)
			for (int i = 0; i < n; i++)
				x[i] = xk_backup[i];
			return alpha;
		}

		// Energy difference negligible — accept alpha
		if (std::abs(trialEnergy - energy) < eps)
		{
			#pragma omp parallel for schedule(static)
			for (int i = 0; i < n; i++)
				x[i] = xk_backup[i];
			return alpha;
		}

		alpha *= m_lsBeta;
	}

	// LS exhausted
	LOG_WARN << "LS exhausted: dirDeriv=" << dirDeriv << " E0=" << energy;

	// Restore original position
	#pragma omp parallel for schedule(static)
	for (int i = 0; i < n; i++)
		x[i] = xk_backup[i];

	return static_cast<Real>(0);
}
#endif

/** Solve the optimization problem for elastic forces.
*
* min_x  (1/2)(x - x_tilde)^T M (x - x_tilde) + dt^2 * Psi(x)
*
* Newton: A * dx = -gradient, where A = M + dt^2 * D^T * K * D (prefactored Cholesky)
* L-BFGS: quasi-Newton with prefactored A as initial inverse Hessian H_0,
*         secant pairs (s_k, y_k) in a circular queue of size m (window size)
*/
#ifdef USE_AVX
/** AVX step elasticity solver. L-BFGS only (solverType=1).
*
*   All state (xk, xTilde, dx, gradient, lbfgs workspace) lives in Scalarf8.
*   We pack once from FluidModel positions at entry, iterate preconditioned
*   gradient descent (computeEnergyAndGradient -> prefactorizedLLTSolve),
*   and unpack to velocities at the end via updateVelocity.
*   No line search in this path — step size fixed at 1.
*/
void Elasticity_Kee2023::stepElasticitySolver()
{
	START_TIMING("Elasticity_Kee2023")
	const unsigned int numActiveParticles = m_model->numActiveParticles();
	if (numActiveParticles == 0)
		return;

	precomputeValues();

	// AVX build supports both Newton (solverType=0) and L-BFGS (solverType=1).

	size_t numObjects = m_objects.size();

	for (auto objIndex = 0; objIndex < numObjects; objIndex++)
	{
		START_TIMING("objSolve")
		ElasticObject* obj = m_objects[objIndex];
		auto& xk = obj->m_xk_avx;         // Scalarf8
		auto& dx = obj->m_dx_avx;         // Scalarf8
		const std::vector<unsigned int>& group = obj->m_particleIndices;
		int numParticles = (int)group.size();
		const Real fdt = obj->m_factorization->m_dt;

		// xTilde = x + dt*v  (uses externally-set velocity, e.g. from AnimationField)
		computeXTilde(obj);

		// Fixed particles: DFSPH integrator skips them, so write xTilde back to model position
		// so they advance by the externally-prescribed velocity.
		const int nFixed = (int)obj->m_nFixed;
		const int firstFixed = numParticles - nFixed;
		for (int i = firstFixed; i < numParticles; i++)
		{
			const unsigned int i0 = group[i];
			const unsigned int particleIndex = m_initial_to_current_index[i0];
			float v[8];
			obj->m_xTilde_avx[i].store(v);
			m_model->getPosition(particleIndex) = Vector3r((Real)v[0], (Real)v[1], (Real)v[2]);
		}

		// xk = xTilde  (initial iterate)
		#pragma omp parallel for schedule(static)
		for (int i = 0; i < numParticles; i++)
			xk[i] = obj->m_xTilde_avx[i];

		const int maxIter = m_maxIter;
		const int nFree = (int)dx.size();

		obj->m_lbfgs_count = 0;

		int totalCGIter = 0;

		int iter = 0;
		for (; iter < maxIter; iter++)
		{
			// Compute search direction (stored in dx_avx)
			Real energy;
			if (m_solverType == 0)  // Newton
			{
				int cgIter;
				energy = newtonSolve(obj, cgIter);
				totalCGIter += cgIter;
			}
			else  // L-BFGS
			{
				energy = lbfgsSolve(obj);
			}

			// Convergence check on ||dx||_inf (reduce from Scalarf8: max of lanes 0..2)
			Real dxNorm = 0;
			#pragma omp parallel for schedule(static)
			for (int i = 0; i < nFree; i++)
			{
				float v[8];
				dx[i].store(v);
				const Real localMax = std::max(std::max(std::abs((Real)v[0]), std::abs((Real)v[1])), std::abs((Real)v[2]));
				#pragma omp critical
				if (localMax > dxNorm) dxNorm = localMax;
			}

			if (dxNorm < m_maxError)
				break;

			// xk += dx  (step size fixed at 1 — no line search in AVX path)
			#pragma omp parallel for schedule(static)
			for (int i = 0; i < nFree; i++)
				xk[i] += dx[i];
		}

		if (iter == maxIter)
			LOG_WARN << "Elasticity_Kee2023: not converged after " << maxIter << " iterations";

		updateVelocity(obj, fdt);

		double elapsedMs = STOP_TIMING
		LOG_INFO << "STATS frame obj=" << objIndex << " solverIter=" << iter + 1 << " lsIter=0 elapsed_ms=" << elapsedMs;
	}

	STOP_TIMING_AVG
}
#else
void Elasticity_Kee2023::stepElasticitySolver()
{
	START_TIMING("Elasticity_Kee2023")
	const unsigned int numActiveParticles = m_model->numActiveParticles();
	if (numActiveParticles == 0)
		return;

	precomputeValues();

	size_t numObjects = m_objects.size();

	for (auto objIndex = 0; objIndex < numObjects; objIndex++)
	{
		START_TIMING("objSolve")
		ElasticObject* obj = m_objects[objIndex];
		auto& xk = obj->m_xk;
		auto& dx = obj->m_dx;
		const std::vector<unsigned int>& group = obj->m_particleIndices;
		int numParticles = (int)group.size();
		const Real fdt = obj->m_factorization->m_dt;

		// Initialize: compute x_tilde = x + dt * v  (uses externally-set velocity, e.g. from AnimationField)
		computeXTilde(obj);

		// Fixed particles: DFSPH integrator skips them, so write xTilde back to model position
		// so they advance by the externally-prescribed velocity.
		const int nFixed = (int)obj->m_nFixed;
		const int firstFixed = numParticles - nFixed;
		for (int i = firstFixed; i < numParticles; i++)
		{
			const unsigned int i0 = group[i];
			const unsigned int particleIndex = m_initial_to_current_index[i0];
			m_model->getPosition(particleIndex) = obj->m_xTilde[i];
		}

		// x0 = x_tilde
		#pragma omp parallel for schedule(static)
		for (int i = 0; i < numParticles; i++)
			xk[i] = obj->m_xTilde[i];
	
		// Iterative solve
		const int maxIter = m_maxIter;
		const int nFree = (int)dx.size();

		// Renew L-BFGS buffers if window size changed at runtime
		obj->m_lbfgs_count = 0;
		if ((int)obj->m_lbfgs_s.size() != m_lbfgsWindowSize)
		{
			obj->m_lbfgs_s.resize(m_lbfgsWindowSize);
			obj->m_lbfgs_y.resize(m_lbfgsWindowSize);
			for (int w = 0; w < m_lbfgsWindowSize; w++)
			{
				obj->m_lbfgs_s[w].resize(nFree, Vector3r::Zero());
				obj->m_lbfgs_y[w].resize(nFree, Vector3r::Zero());
			}
			obj->m_lbfgs_rho.resize(m_lbfgsWindowSize, 0);
			obj->m_lbfgs_alpha.resize(m_lbfgsWindowSize, 0);
		}
		int totalLSIter = 0;
		int totalCGIter = 0;

		int iter = 0;
		for (; iter < maxIter; iter++)
		{
			// Compute search direction (stored in dx) and gradient (stored in m_gradient)
			Real energy;
			if (m_solverType == 0)  // Newton
			{
				int cgIter;
				energy = newtonSolve(obj, cgIter);
				totalCGIter += cgIter;
			}
			else  // L-BFGS
			{
				energy = lbfgsSolve(obj);
			}

			// Convergence check on ||dx||_inf before line search
			Real dxNorm = 0;
			#pragma omp parallel for schedule(static)
			for (int i = 0; i < nFree; i++)
			{
				const Real localMax = dx[i].cwiseAbs().maxCoeff();
				#pragma omp critical
				if (localMax > dxNorm) dxNorm = localMax;
			}

			if (dxNorm < m_maxError)
				break;

			// Line search
			Real alpha = static_cast<Real>(1.0);
			if (m_useLineSearch)
			{
				int lsIter;
				alpha = lineSearch(obj, energy, lsIter);
				totalLSIter += lsIter;
				if (alpha == static_cast<Real>(0))
				{
					LOG_INFO << "Elasticity_Kee2023: line search failed at iter " << iter;
					break;
				}
			}

			// xk += alpha * dx
			#pragma omp parallel for schedule(static)
			for (int i = 0; i < nFree; i++)
				xk[i] += alpha * dx[i];
		}

		if (iter == maxIter)
			LOG_WARN << "Elasticity_Kee2023: not converged after " << maxIter << " iterations";

		updateVelocity(obj, fdt);
		
		double elapsedMs = STOP_TIMING
		if (m_solverType == 0)  // Newton
			LOG_INFO << "STATS frame obj=" << objIndex << " solverIter=" << iter + 1 << " cgIter=" << totalCGIter << " avgCG=" << (iter > 0 ? totalCGIter / (iter + 1) : 0) << " lsIter=" << totalLSIter << " elapsed_ms=" << elapsedMs;
		else  // L-BFGS
			LOG_INFO << "STATS frame obj=" << objIndex << " solverIter=" << iter + 1 << " lsIter=" << totalLSIter << " elapsed_ms=" << elapsedMs;
	}

	STOP_TIMING_AVG
}
#endif


/** Compute energy density Psi based on material type.
*
*   - Stable Neo-Hookean (Smith et al. 2018 Eq. 14):
*       Psi = (μ/2)(I_C - 3) + (λ/2)(J - α)² - (μ/2)log(I_C + 1)
*   - Co-rotated:
*       Psi = μ||F - R||² + (λ/2)(tr(R^T F) - 3)²
*/
Real Elasticity_Kee2023::computePsi(const Matrix3r& F, const Matrix3r& R) const
{
	const Real J = F.determinant();

	if (m_materialType == ENUM_MATERIAL_STABLE_NEOHOOKEAN)
	{
		// Smith et al. 2018 Eq. 14
		const Real I_C = (F.transpose() * F).trace();
		const Real alpha = 1 + m_mu / m_lambda - m_mu / (4 * m_lambda);
		const Real Jma = J - alpha;
		return static_cast<Real>(0.5) * m_mu * (I_C - static_cast<Real>(3.0))
			+ static_cast<Real>(0.5) * m_lambda * Jma * Jma
			- static_cast<Real>(0.5) * m_mu * std::log(I_C + 1);
	}
	// Co-rotated: Psi = μ||F - R||² + (λ/2)(tr(R^T F) - 3)²
	const Real trRtF = (R.transpose() * F).trace();
	return m_mu * (F - R).squaredNorm() + static_cast<Real>(0.5) * m_lambda * (trRtF - static_cast<Real>(3.0)) * (trRtF - static_cast<Real>(3.0));
}

/** Compute first Piola-Kirchhoff stress P based on material type.
*
*   - Stable Neo-Hookean (Smith et al. 2018 Eq. 18):
*       P = μ(1 - 1/(I_C+1)) F + λ(J - α) cof(F)
*   - Co-rotated:
*       P = 2μ(F - R) + λ(tr(R^T F) - 3) R
*/
Matrix3r Elasticity_Kee2023::computeP(const Matrix3r& F, const Matrix3r& R) const
{
	const Real J = F.determinant();
	Matrix3r cofF;
	cofF.col(0) = F.col(1).cross(F.col(2));
	cofF.col(1) = F.col(2).cross(F.col(0));
	cofF.col(2) = F.col(0).cross(F.col(1));

	if (m_materialType == ENUM_MATERIAL_STABLE_NEOHOOKEAN)
	{
		// Smith et al. 2018 Eq. 18
		const Real I_C = (F.transpose() * F).trace();
		const Real alpha = 1 + m_mu / m_lambda - m_mu / (4 * m_lambda);
		return m_mu * (1 - 1 / (I_C + 1)) * F + m_lambda * (J - alpha) * cofF;
	}
	// Co-rotated: P = 2μ(F - R) + λ(tr(R^T F) - 3) R
	const Real trRtF = (R.transpose() * F).trace();
	return static_cast<Real>(2.0) * m_mu * (F - R) + m_lambda * (trRtF - static_cast<Real>(3.0)) * R;
}

