#include "Elasticity_Kugelstadt2021.h"
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

using namespace SPH;
using namespace GenParam;


int Elasticity_Kugelstadt2021::ITERATIONS_V = -1;
int Elasticity_Kugelstadt2021::MAX_ITERATIONS_V = -1;
int Elasticity_Kugelstadt2021::MAX_ERROR_V = -1;
int Elasticity_Kugelstadt2021::ALPHA = -1;
int Elasticity_Kugelstadt2021::MASS_DAMPING_COEFF = -1;
int Elasticity_Kugelstadt2021::STIFFNESS_DAMPING_COEFF = -1;
int Elasticity_Kugelstadt2021::MAX_NEIGHBORS = -1;


Elasticity_Kugelstadt2021::Elasticity_Kugelstadt2021(FluidModel *model) :
	ElasticityBase(model)
{
	const unsigned int numParticles = model->numActiveParticles();
	m_restVolumes.resize(numParticles);
	m_current_to_initial_index.resize(numParticles);
	m_initial_to_current_index.resize(numParticles);
	m_initialNeighbors.resize(numParticles);
	m_rotations.resize(numParticles, Matrix3r::Identity());
	m_stress.resize(numParticles);
	m_L.resize(numParticles);								// kernel gradient correction matrix L
	m_F.resize(numParticles);								// deformation gradient
	m_RL.resize(numParticles);								// stores the rotation matrix times the matrix L
	m_vDiff.resize(numParticles, Vector3r::Zero());			// velocity difference used for the warm start of the volume cg solver

	m_iterationsV = 0;
	m_maxIterV = 100;
	m_maxErrorV = static_cast<Real>(1.0e-4);
	m_alpha = 0.0;
	m_massDampingCoeff = 0.0;
	m_stiffnessDampingCoeff = 0.0;
	m_youngsModulus = 5000000;
	m_maxNeighbors = -1;

	model->addField({ "rest volume", FieldType::Scalar, [&](const unsigned int i) -> Real* { return &m_restVolumes[i]; }, true });
	model->addField({ "rotation", FieldType::Matrix3, [&](const unsigned int i) -> Real* { return &m_rotations[i](0,0); } });
	model->addField({ "stress", FieldType::Scalar, [&](const unsigned int i) -> Real* { return &m_stress[i]; } });
	model->addField({ "deformation gradient", FieldType::Matrix3, [&](const unsigned int i) -> Real* { return &m_F[i](0,0); } });
	model->addField({ "correction matrix", FieldType::Matrix3, [&](const unsigned int i) -> Real* { return &m_L[i](0,0); } });
}

Elasticity_Kugelstadt2021::~Elasticity_Kugelstadt2021(void)
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

void Elasticity_Kugelstadt2021::deferredInit()
{
	initValues();
}


void Elasticity_Kugelstadt2021::initParameters()
{
	ParameterBase::GetFunc<Real> getFct = [&]()-> Real { return m_youngsModulus; };
	ParameterBase::SetFunc<Real> setFct = [&](Real val) 
	{ 
		m_youngsModulus = val; 
		m_mu = m_youngsModulus / (static_cast<Real>(2.0) * (static_cast<Real>(1.0) + m_poissonRatio));
		m_lambda = m_youngsModulus * m_poissonRatio / ((static_cast<Real>(1.0) + m_poissonRatio) * (static_cast<Real>(1.0) - static_cast<Real>(2.0) * m_poissonRatio));

		if (Simulation::getCurrent()->isSimulationInitialized())		// if Young's modulus has changed, recompute the factorization
			Simulation::getCurrent()->reset();
	};
	YOUNGS_MODULUS = createNumericParameter("youngsModulus", "Young`s modulus", getFct, setFct);
	setGroup(YOUNGS_MODULUS, "Elasticity");
	setDescription(YOUNGS_MODULUS, "Stiffness of the elastic material");
	RealParameter* rparam = static_cast<RealParameter*>(getParameter(YOUNGS_MODULUS));
	rparam->setMinValue(0.0);

	ParameterBase::GetFunc<Real> getFct2 = [&]()-> Real { return m_poissonRatio; };
	ParameterBase::SetFunc<Real> setFct2 = [&](Real val)
	{
		m_poissonRatio = val;
		m_mu = m_youngsModulus / (static_cast<Real>(2.0) * (static_cast<Real>(1.0) + m_poissonRatio));
		m_lambda = m_youngsModulus * m_poissonRatio / ((static_cast<Real>(1.0) + m_poissonRatio) * (static_cast<Real>(1.0) - static_cast<Real>(2.0) * m_poissonRatio));

		if (Simulation::getCurrent()->isSimulationInitialized())		// if Poisson ration has changed, recompute the factorization
			Simulation::getCurrent()->reset();
	};
	POISSON_RATIO = createNumericParameter("poissonsRatio", "Poisson`s ratio", getFct2, setFct2);
	setGroup(POISSON_RATIO, "Elasticity");
	setDescription(POISSON_RATIO, "Ratio of transversal expansion and axial compression");
	rparam = static_cast<RealParameter*>(getParameter(POISSON_RATIO));
	rparam->setMinValue(static_cast<Real>(-1.0 + 1e-4));
	rparam->setMaxValue(static_cast<Real>(0.5 - 1e-4));

	ITERATIONS_V = createNumericParameter("volumeIterations", "Iterations", &m_iterationsV);
	setGroup(ITERATIONS_V, "Elasticity");
	setDescription(ITERATIONS_V, "Iterations required by the volume solver.");
	getParameter(ITERATIONS_V)->setReadOnly(true);

	MAX_ITERATIONS_V = createNumericParameter("volumeMaxIter", "Max. iterations (volume solver)", &m_maxIterV);
	setGroup(MAX_ITERATIONS_V, "Elasticity");
	setDescription(MAX_ITERATIONS_V, "Max. iterations of the volume solver");
	static_cast<NumericParameter<unsigned int>*>(getParameter(MAX_ITERATIONS_V))->setMinValue(0);

	MAX_ERROR_V = createNumericParameter("volumeMaxError", "Max. volume error", &m_maxErrorV);
	setGroup(MAX_ERROR_V, "Elasticity");
	setDescription(MAX_ERROR_V, "Max. error of the volume solver");
	rparam = static_cast<RealParameter*>(getParameter(MAX_ERROR_V));
	rparam->setMinValue(1e-9);

	ALPHA = createNumericParameter("alpha", "Zero-energy modes suppression", &m_alpha);
	setGroup(ALPHA, "Elasticity");
	setDescription(ALPHA, "Coefficent for zero-energy modes suppression method");
	rparam = static_cast<RealParameter*>(getParameter(ALPHA));
	rparam->setMinValue(0.0);

	ParameterBase::GetFunc<Real> getFct3 = [&]()-> Real { return m_massDampingCoeff; };
	ParameterBase::SetFunc<Real> setFct3 = [&](Real val)
	{
		m_massDampingCoeff = val;
		if (Simulation::getCurrent()->isSimulationInitialized())		// if value has changed, recompute the factorization
			Simulation::getCurrent()->reset();
	};
	MASS_DAMPING_COEFF = createNumericParameter("massDamping", "Mass damping coeff.", getFct3, setFct3);
	setGroup(MASS_DAMPING_COEFF, "Elasticity");
	setDescription(MASS_DAMPING_COEFF, "Mass damping coefficient (Rayleigh damping)");
	rparam = static_cast<RealParameter*>(getParameter(MASS_DAMPING_COEFF));
	rparam->setMinValue(0.0);

	ParameterBase::GetFunc<Real> getFct4 = [&]()-> Real { return m_stiffnessDampingCoeff; };
	ParameterBase::SetFunc<Real> setFct4 = [&](Real val)
	{
		m_stiffnessDampingCoeff = val;
		if (Simulation::getCurrent()->isSimulationInitialized())		// if value has changed, recompute the factorization
			Simulation::getCurrent()->reset();
	};
	STIFFNESS_DAMPING_COEFF = createNumericParameter("stiffnessDamping", "Stiffness damping coeff.", getFct4, setFct4);
	setGroup(STIFFNESS_DAMPING_COEFF, "Elasticity");
	setDescription(STIFFNESS_DAMPING_COEFF, "Stiffness damping coefficient (Rayleigh damping)");
	rparam = static_cast<RealParameter*>(getParameter(STIFFNESS_DAMPING_COEFF));
	rparam->setMinValue(0.0);

	ParameterBase::GetFunc<int> getFct5 = [&]()-> int { return m_maxNeighbors; };
	ParameterBase::SetFunc<int> setFct5 = [&](int val)
	{
		m_maxNeighbors = val;
		if (Simulation::getCurrent()->isSimulationInitialized())		// if value has changed, recompute the factorization
			Simulation::getCurrent()->reset();
	};
	MAX_NEIGHBORS = createNumericParameter("maxNeighbors", "Max. neighbors", getFct5, setFct5);
	setGroup(MAX_NEIGHBORS, "Elasticity");
	setDescription(MAX_NEIGHBORS, "Maximum number of neighbors that are considered.");

	ParameterBase::GetVecFunc<Real> getFct6 = [&]()-> Real* { return m_fixedBoxMin.data(); };
	ParameterBase::SetVecFunc<Real> setFct6 = [&](Real* val)
	{
		m_fixedBoxMin = Vector3r(val[0], val[1], val[2]);
		determineFixedParticles();
	};
	FIXED_BOX_MIN = createVectorParameter("fixedBoxMin", "Fixed box min", 3u, getFct6, setFct6);
	setGroup(FIXED_BOX_MIN, "Elasticity");
	setDescription(FIXED_BOX_MIN, "Minimum point of box of which contains the fixed particles.");
	getParameter(FIXED_BOX_MIN)->setReadOnly(true);


	ParameterBase::GetVecFunc<Real> getFct7 = [&]()-> Real* { return m_fixedBoxMax.data(); };
	ParameterBase::SetVecFunc<Real> setFct7 = [&](Real* val)
	{
		m_fixedBoxMax = Vector3r(val[0], val[1], val[2]);
		determineFixedParticles();
	};
	FIXED_BOX_MAX = createVectorParameter("fixedBoxMax", "Fixed box max", 3u, getFct7, setFct7);
	setGroup(FIXED_BOX_MAX, "Elasticity");
	setDescription(FIXED_BOX_MAX, "Maximum point of box of which contains the fixed particles.");
	getParameter(FIXED_BOX_MAX)->setReadOnly(true);
}

/** Compute an MD4 check sum using the neighborhood structure in order to 
* recognize known particle models (cache).
*/
std::string Elasticity_Kugelstadt2021::computeMD5(const unsigned int objIndex)
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
void Elasticity_Kugelstadt2021::initValues()
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

			// compute rest volume
			Real density = model->getMass(i) * sim->W_zero();
			const Vector3r &xi0 = model->getPosition0(i);
			for (size_t j = 0; j < m_initialNeighbors[i].size(); j++)
			{
				const unsigned int neighborIndex0 = m_initialNeighbors[i][j];
				const Vector3r& xj0 = model->getPosition0(neighborIndex0);
				density += model->getMass(neighborIndex0) * sim->W(xi0 - xj0);
			}
			m_restVolumes[i] = model->getMass(i) / density;
			m_rotations[i].setIdentity();
			m_stress[i] = 0.0;
			m_F[i].setIdentity();
			m_vDiff[i].setZero();
			m_RL[i].setIdentity();
		}
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
void Elasticity_Kugelstadt2021::findObjects()
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
			Comparator(Elasticity_Kugelstadt2021* _this) : m_this(_this) {};
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

			Elasticity_Kugelstadt2021* m_this;
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
void Elasticity_Kugelstadt2021::initSystem()
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
			"_" + Utilities::StringTools::real2String(m_massDampingCoeff) +
			"_" + Utilities::StringTools::real2String(m_stiffnessDampingCoeff) +
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
				"_" + Utilities::StringTools::real2String(m_massDampingCoeff) +
				"_" + Utilities::StringTools::real2String(m_stiffnessDampingCoeff) +
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
				binReader.readSparseMatrix(obj->m_factorization->m_dampingMatrix);
				binReader.read(obj->m_factorization->m_dt);
				binReader.read(obj->m_factorization->m_mu);
				binReader.readSparseMatrix(obj->m_factorization->m_matHTH);
				delete obj->m_factorization->m_cholesky;
				obj->m_factorization->m_cholesky = new CholeskyAVXSolver();
				obj->m_factorization->m_cholesky->load(binReader);
				binReader.closeFile();
			}

			// init vectors
			int numParticles = (int)obj->m_particleIndices.size();
			obj->m_f_avx.resize(3 * numParticles);
			obj->m_sol_avx.resize(numParticles);
			obj->m_v_avx.resize(numParticles);
			obj->m_RHS.resize(numParticles - obj->m_nFixed);
			obj->m_rhs.resize(3 * numParticles - obj->m_nFixed);
			obj->m_sol.resize(3 * numParticles - obj->m_nFixed);

			int vecSize;
			if (numParticles % 8 == 0) vecSize = numParticles / 8;
			else vecSize = numParticles / 8 + 1;
			obj->m_quats.resize(vecSize);
			for (int i = 0; i < vecSize; i++)
				obj->m_quats[i] = Quaternion8f();
		}
		else    // no cache found
		{
			ElasticObject* obj = m_objects[objIndex];

			// compute new factorization if no factorization can be reused
			if (!foundFactorization)
				initFactorization(obj->m_factorization, obj->m_particleIndices, obj->m_nFixed, dt, m_mu);

			// init vectors
			int numParticles = (int)obj->m_particleIndices.size();
			obj->m_f_avx.resize(3 * numParticles);
			obj->m_sol_avx.resize(numParticles);
			obj->m_v_avx.resize(numParticles);
			obj->m_RHS.resize(numParticles - obj->m_nFixed);
			obj->m_rhs.resize(3 * numParticles - obj->m_nFixed);
			obj->m_sol.resize(3 * numParticles - obj->m_nFixed);

			int vecSize;
			if (numParticles % 8 == 0) vecSize = numParticles / 8;
			else vecSize = numParticles / 8 + 1;
			obj->m_quats.resize(vecSize);
			for (int i = 0; i < vecSize; i++)
				obj->m_quats[i] = Quaternion8f();

			// write cache file
			if (sim->getUseCache() && (Utilities::FileSystem::makeDir(sim->getCachePath()) == 0))
			{
				BinaryFileWriter binWriter;
				binWriter.openFile(cacheFileName);
				binWriter.writeSparseMatrix(obj->m_factorization->m_D);
				binWriter.writeSparseMatrix(obj->m_factorization->m_DT_K);
				binWriter.writeSparseMatrix(obj->m_factorization->m_dampingMatrix);
				binWriter.write(obj->m_factorization->m_dt);
				binWriter.write(obj->m_factorization->m_mu);
				binWriter.writeSparseMatrix(obj->m_factorization->m_matHTH);
				obj->m_factorization->m_cholesky->save(binWriter);
				binWriter.closeFile();
			}
		}
		obj->m_md5 = md5;
	}
}

/** Compute the factorization of the linear system matrix. 
* This is only done once at the beginning of the simulation. 
*/
void Elasticity_Kugelstadt2021::initFactorization(std::shared_ptr<Factorization> factorization, std::vector<unsigned int>& particleIndices, const unsigned int nFixed, const Real dt, const Real mu)
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
		const unsigned int numNeighbors = m_initialNeighbors[i].size();
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
		int particleIndex = group[i];
		const unsigned int i0 = m_current_to_initial_index[particleIndex];

		const double restVolumes_id = m_restVolumes[i];

		// init matrix K according to Eq. 13 in the paper
		// Set triplets for the matrix K. Directly multiply the factor 2*dt*dt into K
		const double Kreal = 2.0 * dtd * dtd * mud * restVolumes_id;
		

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
			const unsigned int neighborIndex = m_initial_to_current_index[m_initialNeighbors[i0][j]];
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
				const unsigned int kIndex = m_initial_to_current_index[m_initialNeighbors[i0][k]];
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
	if ((m_massDampingCoeff != 0.0) || (m_stiffnessDampingCoeff != 0.0))
	{
		Eigen::SparseMatrix<double> dampingMatrix(numParticles, numParticles);
		factorization->m_dampingMatrix.resize(numParticles, numParticles);
		dampingMatrix = dtd * (static_cast<double>(m_massDampingCoeff) * M) + static_cast<double>(m_stiffnessDampingCoeff) / dtd * (DT_K_D + HTH);
		factorization->m_dampingMatrix = dampingMatrix.cast<Real>();

		// init linear system matrix according to Eq. 29 in the paper (+ the Rayleigh damping matrix)
		if (m_alpha != 0.0)
			M_plus_DT_K_D = (M + DT_K_D + HTH + dampingMatrix).block(0, 0, numParticles - nFixed, numParticles - nFixed);
		else      // no zero energy mode control
			M_plus_DT_K_D = (M + DT_K_D + dampingMatrix).block(0, 0, numParticles - nFixed, numParticles - nFixed);
	}
	else
	{
		// init linear system matrix according to Eq. 29 in the paper
		if (m_alpha != 0.0)
			M_plus_DT_K_D = (M + DT_K_D + HTH).block(0, 0, numParticles - nFixed, numParticles - nFixed);
		else      // no zero energy mode control
			M_plus_DT_K_D = (M + DT_K_D).block(0, 0, numParticles - nFixed, numParticles - nFixed);
	}

	M_plus_DT_K_D.makeCompressed();

	// compute factorization of the matrix
	delete factorization->m_cholesky;
	factorization->m_cholesky = new CholeskyAVXSolver(M_plus_DT_K_D);

	LOG_INFO << "Non zero elements (A): " << M_plus_DT_K_D.nonZeros();
}

/** Perform a step of the elasticity solver.
 */
void Elasticity_Kugelstadt2021::step()
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
	
	stepVolumeSolver();						// volume solver using a cg method
	STOP_TIMING_AVG
}

/** Solve the linear system for the stretching forces including zero energy mode control
* using the precomputed matrix factorization.
*/
void Elasticity_Kugelstadt2021::stepElasticitySolver()
{
	START_TIMING("Elasticity_Kugelstadt2021")
	const unsigned int numActiveParticles = m_model->numActiveParticles();
	if (numActiveParticles == 0)
		return;
	const Real dt = TimeManager::getCurrent()->getTimeStepSize();
	Simulation* sim = Simulation::getCurrent();

	size_t numObjects = m_objects.size();

	// solve the systems for each object separately
	for (auto objIndex = 0; objIndex < numObjects; objIndex++)
	{
		ElasticObject* obj = m_objects[objIndex];
		const std::vector<unsigned int>& group = obj->m_particleIndices;
		int numParticles = (int)group.size();

		auto& D = obj->m_factorization->m_D;
		auto& DT_K = obj->m_factorization->m_DT_K;
		auto& HT_K_H = obj->m_factorization->m_matHTH;

		auto& dampingMatrix = obj->m_factorization->m_dampingMatrix;
		auto& RHS = obj->m_RHS;
		auto& f_avx = obj->m_f_avx;
		auto& sol_avx = obj->m_sol_avx;
		auto& v_avx = obj->m_v_avx;
		auto& quats = obj->m_quats;

		START_TIMING("advect x & Dx")
		#pragma omp parallel default(shared)
		{
			//////////////////////////////////////////////////////////////////////////
			// advect particles to get \tilde x in Eq. 29:
			// store the 3 components of the advected positions in f_avx
			// store the 3 components of dt*velocity in v_avx
			//////////////////////////////////////////////////////////////////////////
			#pragma omp for schedule(static)  
			for (int i = 0; i < (int)numParticles; i++)
			{
				const unsigned int i0 = group[i];
				const unsigned int particleIndex = m_initial_to_current_index[i0];
				const Vector3r& xi0 = m_model->getPosition0(i0);
				const size_t numNeighbors = m_initialNeighbors[i0].size();

				const Real fdt = obj->m_factorization->m_dt;
				const Vector3r x = m_model->getPosition(particleIndex);
				const Vector3r dv = fdt * m_model->getVelocity(particleIndex);
				const Vector3r xNew = x + dv;
				// copy the 3 coordinates to sol_avx
				v_avx[i] = Scalarf8(dv[0], dv[1], dv[2], 0, 0, 0, 0, 0);
				sol_avx[i] = Scalarf8(xNew[0], xNew[1], xNew[2], 0, 0, 0, 0, 0);
				f_avx[3 * i].setZero();
				f_avx[3 * i + 1].setZero();
				f_avx[3 * i + 2].setZero();
			}

			//////////////////////////////////////////////////////////////////////////
			// compute deformation gradient F by sparse matrix-vector product in parallel:
			//
			// f_avx = D * x_advected		(Eq. 12 and Eq. 29)
			//////////////////////////////////////////////////////////////////////////
			#pragma omp for schedule(static)  
			for (int k = 0; k < D.outerSize(); ++k)
			{
				for (Eigen::SparseMatrix<float, Eigen::RowMajor>::InnerIterator it(D, k); it; ++it)
				{
					f_avx[it.row()] += Scalarf8(it.value()) * sol_avx[it.col()];
				}
			}
		}
		STOP_TIMING_AVG;

		int vecSize;
		if (numParticles % 8 == 0) 
			vecSize = numParticles / 8;
		else 
			vecSize = numParticles / 8 + 1;

		START_TIMING("extract rot");
		#pragma omp parallel default(shared)
		{
			//////////////////////////////////////////////////////////////////////////
			// extract deformation gradient from avx values f_avx
			//////////////////////////////////////////////////////////////////////////
			#pragma omp for schedule(static)  
			for (int i = 0; i < (int)numParticles; i++)
			{
				const unsigned int i0 = group[i];
				const unsigned int particleIndex = m_initial_to_current_index[i0];

				// copy data from f_avx to m_F
				float x0[8], x1[8], x2[8];
				f_avx[3*i].store(x0);
				f_avx[3*i+1].store(x1);
				f_avx[3*i+2].store(x2);

				m_F[particleIndex] <<	x0[0], x0[1], x0[2],
										x1[0], x1[1], x1[2],
										x2[0], x2[1], x2[2];
			}

			//////////////////////////////////////////////////////////////////////////
			// extract rotation from F by analytic polar decomposition (Kugelstadt et al. 2018)
			//////////////////////////////////////////////////////////////////////////
			#pragma omp for schedule(static)  
			for (int i = 0; i < vecSize; i++)
			{
				const int count = std::min(numParticles - i*8, 8);

				// store the deformation gradient of 8 particles in avx vectors
				int idx[8];
				for (int j=0; j < count; j++)
					idx[j] = m_initial_to_current_index[group[8*i + j]];
				for (int j = count; j < 8; j++)
					idx[j] = 0;

				Vector3f8 F1, F2, F3;	//columns of the deformation gradient
				F1 = Vector3f8(m_F[idx[0]].col(0), m_F[idx[1]].col(0), m_F[idx[2]].col(0), m_F[idx[3]].col(0), m_F[idx[4]].col(0), m_F[idx[5]].col(0), m_F[idx[6]].col(0), m_F[idx[7]].col(0));
				F2 = Vector3f8(m_F[idx[0]].col(1), m_F[idx[1]].col(1), m_F[idx[2]].col(1), m_F[idx[3]].col(1), m_F[idx[4]].col(1), m_F[idx[5]].col(1), m_F[idx[6]].col(1), m_F[idx[7]].col(1));
				F3 = Vector3f8(m_F[idx[0]].col(2), m_F[idx[1]].col(2), m_F[idx[2]].col(2), m_F[idx[3]].col(2), m_F[idx[4]].col(2), m_F[idx[5]].col(2), m_F[idx[6]].col(2), m_F[idx[7]].col(2));
				
				// perform polar decomposition
				Quaternion8f& q = quats[i];
				APD_Newton_AVX(F1, F2, F3, q);

				//transform quaternion to rotation matrix
				Vector3f8 R1, R2, R3;	//columns of the rotation matrix
				quats[i].toRotationMatrix(R1, R2, R3);

				//////////////////////////////////////////////////////////////////////////
				// R := R-F
				//////////////////////////////////////////////////////////////////////////
				R1 -= F1;
				R2 -= F2;
				R3 -= F3;

				//////////////////////////////////////////////////////////////////////////
				// store result in f_avx
				// f_avx has size 3*n and 3 rows contain F-R
				//////////////////////////////////////////////////////////////////////////
				std::array<Vector3r, 8> v0, v1, v2;
				R1.store(v0.data());
				R2.store(v1.data());
				R3.store(v2.data());
				for (auto j = 0; j < count; j++)
				{
					f_avx[24*i + 3*j] =     Scalarf8(v0[j][0], v1[j][0], v2[j][0], 0, 0, 0, 0, 0);
					f_avx[24*i + 3*j + 1] = Scalarf8(v0[j][1], v1[j][1], v2[j][1], 0, 0, 0, 0, 0);
					f_avx[24*i + 3*j + 2] = Scalarf8(v0[j][2], v1[j][2], v2[j][2], 0, 0, 0, 0, 0);
					if (8*i+j < RHS.size())
						RHS[8*i+j] = Scalarf8(0.0f);
				}
			}
		}
		STOP_TIMING_AVG

		//////////////////////////////////////////////////////////////////////////
		// Compute right hand side
		//////////////////////////////////////////////////////////////////////////		
	
		START_TIMING("rhs")

		#pragma omp parallel default(shared)
		{
			//////////////////////////////////////////////////////////////////////////
			// Compute D^T K * (R-F)		(Eq. 29)
			// Note: K already contains the factor: 2 dt^2
			//////////////////////////////////////////////////////////////////////////		
			#pragma omp for schedule(static)  
			for (int k = 0; k < DT_K.outerSize(); ++k)
			{
				for (Eigen::SparseMatrix<float, Eigen::RowMajor>::InnerIterator it(DT_K, k); it; ++it)
				{
					if (it.row() < (int) RHS.size())
						RHS[it.row()] += Scalarf8(it.value()) * f_avx[it.col()];
				}
			}
		}

		//////////////////////////////////////////////////////////////////////////
		// If zero energy model control is turned on, 
		// add the following term to the right hand side (Eq. 29): 
		// H^T * K2 * H * x_advected
		// Note: K2 already contains the factor: dt^2
		//////////////////////////////////////////////////////////////////////////		
		if (m_alpha != 0.0)
		{
			#pragma omp parallel default(shared)
			{
				#pragma omp for schedule(static)  
				for (int k = 0; k < HT_K_H.outerSize(); ++k)
				{
					for (Eigen::SparseMatrix<Real, Eigen::ColMajor>::InnerIterator it(HT_K_H, k); it; ++it)
					{
						if (it.col() < (int)RHS.size())
							RHS[it.col()] -= Scalarf8(static_cast<float>(it.value())) * sol_avx[it.row()];
					}
				}

			}
		}

		//////////////////////////////////////////////////////////////////////////
		// If damping is turned on, 
		// add the following term to the right hand side: 
		// DampingMatrix * x_advected
		// Note: DampingMatrix already contains the factor: dt
		//////////////////////////////////////////////////////////////////////////		
		if ((m_massDampingCoeff != 0.0) || (m_stiffnessDampingCoeff != 0.0))
		{
			#pragma omp parallel default(shared)
			{
				#pragma omp for schedule(static)  
				for (int k = 0; k < dampingMatrix.outerSize(); ++k)
				{
					for (Eigen::SparseMatrix<float, Eigen::ColMajor>::InnerIterator it(dampingMatrix, k); it; ++it)
					{
						if (it.col() < (int)RHS.size())
							RHS[it.col()] -= Scalarf8(it.value()) * v_avx[it.row()];
					}
				}
			}
		}

		//////////////////////////////////////////////////////////////////////////
		// Copy the right hand side in a vector for the cholesky solver
		//////////////////////////////////////////////////////////////////////////
		#pragma omp parallel default(shared)
		{
			#pragma omp for schedule(static)  
			for (int i = 0; i < (int)RHS.size(); i++)
			{
				float x[8];
				RHS[i].store(x);
				obj->m_rhs[3 * i] = x[0];
				obj->m_rhs[3 * i + 1] = x[1];
				obj->m_rhs[3 * i + 2] = x[2];
			}
		}

		STOP_TIMING_AVG;
	}


	//////////////////////////////////////////////////////////////////////////
	// Solve linear system 
	//////////////////////////////////////////////////////////////////////////		
	START_TIMING("solve SLE")

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < 3*numObjects; i++)
		{
			int objIndex = i / 3;
			int index = i % 3;
			ElasticObject* obj = m_objects[objIndex];
			obj->m_factorization->m_cholesky->solve(&obj->m_sol[index], &obj->m_rhs[index], /* stride = */ 3);
		}
	}
	for (auto objIndex = 0; objIndex < numObjects; objIndex++)
	{
		ElasticObject* obj = m_objects[objIndex];
		const std::vector<unsigned int>& group = obj->m_particleIndices;
		#pragma omp parallel default(shared)
		{
			#pragma omp for schedule(static)  
			for (int i = 0; i < (int) obj->m_sol.size()/3; i++)
			{
				const unsigned int i0 = group[i];
				const unsigned int particleIndex = m_initial_to_current_index[i0];
				if (m_model->getParticleState(particleIndex) == ParticleState::Active)
				{
					Vector3r& vi = m_model->getVelocity(particleIndex);
					const Vector3r &dx = obj->m_sol.segment<3>(3*i);
					const Real fdt = obj->m_factorization->m_dt;
					vi += (1.0 / fdt) * dx;
				}
			}
		}		
	}
	STOP_TIMING_AVG
	
	STOP_TIMING_AVG
}

void Elasticity_Kugelstadt2021::reset()
{
	initValues();
}

void Elasticity_Kugelstadt2021::performNeighborhoodSearchSort()
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

	// update quaterions which are needed for warmstart
	rotationMatricesToAVXQuaternions();
}

/** Extract rotation matrices from deformation gradients.
*/
void Elasticity_Kugelstadt2021::computeRotations()
{
	const Real dt = TimeManager::getCurrent()->getTimeStepSize();

	size_t numObjects = m_objects.size();
	for (auto objIndex = 0; objIndex < numObjects; objIndex++)
	{
		ElasticObject* obj = m_objects[objIndex];
		const std::vector<unsigned int>& group = obj->m_particleIndices;
		int numParticles = (int)group.size();

		int vecSize;
		if (numParticles % 8 == 0)
			vecSize = numParticles / 8;
		else
			vecSize = numParticles / 8 + 1;

		auto& D = obj->m_factorization->m_D;
		auto& f_avx = obj->m_f_avx;
		auto& sol_avx = obj->m_sol_avx;
		auto& quats = obj->m_quats;

		//////////////////////////////////////////////////////////////////////////
		// advect particles
		//////////////////////////////////////////////////////////////////////////
		#pragma omp parallel default(shared)
		{
			#pragma omp for schedule(static)  
			for (int i = 0; i < (int)numParticles; i++)
			{
				const unsigned int i0 = group[i];
				const unsigned int particleIndex = m_initial_to_current_index[i0];

				const Vector3r& xi0 = m_model->getPosition0(i0);
				const size_t numNeighbors = m_initialNeighbors[i0].size();

				Vector3r xNew = m_model->getPosition(particleIndex);
				xNew += dt * m_model->getVelocity(particleIndex);

				// copy the 3 coordinates to sol_avx
				sol_avx[i] = Scalarf8(xNew[0], xNew[1], xNew[2], 0, 0, 0, 0, 0);
				f_avx[3 * i].setZero();
				f_avx[3 * i + 1].setZero();
				f_avx[3 * i + 2].setZero();
			}

			// compute sparse matrix-vector product in parallel
			#pragma omp for schedule(static)  
			for (int k = 0; k < D.outerSize(); ++k)
			{
				for (Eigen::SparseMatrix<float, Eigen::RowMajor>::InnerIterator it(D, k); it; ++it)
				{
					f_avx[it.row()] += Scalarf8(it.value()) * sol_avx[it.col()];
				}
			}
		
			#pragma omp for schedule(static)  
			for (int i = 0; i < (int)numParticles; i++)
			{
				const unsigned int i0 = group[i];
				const unsigned int particleIndex = m_initial_to_current_index[i0];

				// copy data from f_avx to m_F
				float x0[8], x1[8], x2[8];
				f_avx[3 * i].store(x0);
				f_avx[3 * i + 1].store(x1);
				f_avx[3 * i + 2].store(x2);

				m_F[particleIndex] << x0[0], x0[1], x0[2],
										x1[0], x1[1], x1[2],
										x2[0], x2[1], x2[2];
			}

			#pragma omp for schedule(static)  
			for (int i = 0; i < vecSize; i++)
			{
				const int count = std::min(numParticles - i * 8, 8);

				// store the deformation gradient of 8 particles in avx vectors
				int idx[8];
				for (int j = 0; j < count; j++)
					idx[j] = m_initial_to_current_index[group[8 * i + j]];
				for (int j = count; j < 8; j++)
					idx[j] = 0;

				Vector3f8 F1, F2, F3;	//columns of the deformation gradient
				F1 = Vector3f8(m_F[idx[0]].col(0), m_F[idx[1]].col(0), m_F[idx[2]].col(0), m_F[idx[3]].col(0), m_F[idx[4]].col(0), m_F[idx[5]].col(0), m_F[idx[6]].col(0), m_F[idx[7]].col(0));
				F2 = Vector3f8(m_F[idx[0]].col(1), m_F[idx[1]].col(1), m_F[idx[2]].col(1), m_F[idx[3]].col(1), m_F[idx[4]].col(1), m_F[idx[5]].col(1), m_F[idx[6]].col(1), m_F[idx[7]].col(1));
				F3 = Vector3f8(m_F[idx[0]].col(2), m_F[idx[1]].col(2), m_F[idx[2]].col(2), m_F[idx[3]].col(2), m_F[idx[4]].col(2), m_F[idx[5]].col(2), m_F[idx[6]].col(2), m_F[idx[7]].col(2));


				// perform analytic polar decomposition
				Quaternion8f& q = quats[i];
				APD_Newton_AVX(F1, F2, F3, q);
				//transform quaternion to rotation matrix
				Vector3f8 R1, R2, R3;	//columns of the rotation matrix
				quats[i].toRotationMatrix(R1, R2, R3);

				std::array<Vector3r, 8> r0, r1, r2;
				R1.store(r0.data());
				R2.store(r1.data());
				R3.store(r2.data());
				for (auto j = 0; j < count; j++)
				{
					m_rotations[idx[j]].row(0) = r0[j];
					m_rotations[idx[j]].row(1) = r1[j];
					m_rotations[idx[j]].row(2) = r2[j];
				}
			}
		}
	}
}

/** Convert all rotation matrices to AVX quaternions.
*/
void Elasticity_Kugelstadt2021::rotationMatricesToAVXQuaternions()
{
	size_t numObjects = m_objects.size();
	for (auto objIndex = 0; objIndex < numObjects; objIndex++)
	{
		ElasticObject* obj = m_objects[objIndex];
		const std::vector<unsigned int>& group = obj->m_particleIndices;
		int numParticles = (int)group.size();

		int vecSize;
		if (numParticles % 8 == 0)
			vecSize = numParticles / 8;
		else
			vecSize = numParticles / 8 + 1;

		auto& quats = obj->m_quats;

		#pragma omp parallel default(shared)
		{
			#pragma omp for schedule(static)  
			for (int i = 0; i < vecSize; i++)
			{
				const int count = std::min(numParticles - i * 8, 8);

				// store the deformation gradient of 8 particles in avx vectors
				int idx[8];
				for (int j = 0; j < count; j++)
					idx[j] = m_initial_to_current_index[group[8 * i + j]];
				for (int j = count; j < 8; j++)
					idx[j] = 0;

				Quaternionr q[8];
				for (auto j = 0; j < count; j++)
					q[j] = Quaternionr(m_rotations[idx[j]].transpose());
				for (auto j = count; j < 8; j++)
					q[j] = Quaternionr();

				Quaternion8f& q_avx = quats[i];
				q_avx.set(q);
			}
		}
	}
}

/** Compute kernel gradient correction matrices (Eq. 8).
*/
void Elasticity_Kugelstadt2021::computeMatrixL()
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

// ----------------------------------------------------------------------------------------------
//computes the APD of 8 deformation gradients. (Alg. 3 from the paper: Kugelstadt et al. "Fast Corotated FEM using Operator Splitting", CGF 2018)
inline void Elasticity_Kugelstadt2021::APD_Newton_AVX(const Vector3f8& F1, const Vector3f8& F2, const Vector3f8& F3, Quaternion8f& q)
{
	//one iteration is sufficient for plausible results
	for (int it = 0; it < 1; it++)
	{
		//transform quaternion to rotation matrix
		Matrix3f8 R;
		q.toRotationMatrix(R);

		//columns of B = RT * F
		Vector3f8 B0 = R.transpose() * F1;
		Vector3f8 B1 = R.transpose() * F2;
		Vector3f8 B2 = R.transpose() * F3;

		Vector3f8 gradient(B2[1] - B1[2], B0[2] - B2[0], B1[0] - B0[1]);

		//compute Hessian, use the fact that it is symmetric
		Scalarf8 h00 = B1[1] + B2[2];
		Scalarf8 h11 = B0[0] + B2[2];
		Scalarf8 h22 = B0[0] + B1[1];
		Scalarf8 h01 = Scalarf8(-0.5) * (B1[0] + B0[1]);
		Scalarf8 h02 = Scalarf8(-0.5) * (B2[0] + B0[2]);
		Scalarf8 h12 = Scalarf8(-0.5) * (B2[1] + B1[2]);

		Scalarf8 detH = Scalarf8(-1.0) * h02 * h02 * h11 + Scalarf8(2.0) * h01 * h02 * h12 - h00 * h12 * h12 - h01 * h01 * h22 + h00 * h11 * h22;

		Vector3f8 omega;
		//compute symmetric inverse
		const Scalarf8 factor = Scalarf8(-0.25) / detH;
		omega[0] = (h11 * h22 - h12 * h12) * gradient[0]
			+ (h02 * h12 - h01 * h22) * gradient[1]
			+ (h01 * h12 - h02 * h11) * gradient[2];
		omega[0] *= factor;

		omega[1] = (h02 * h12 - h01 * h22) * gradient[0]
			+ (h00 * h22 - h02 * h02) * gradient[1]
			+ (h01 * h02 - h00 * h12) * gradient[2];
		omega[1] *= factor;

		omega[2] = (h01 * h12 - h02 * h11) * gradient[0]
			+ (h01 * h02 - h00 * h12) * gradient[1]
			+ (h00 * h11 - h01 * h01) * gradient[2];
		omega[2] *= factor;

		omega = Vector3f8::blend(abs(detH) < 1.0e-9f, gradient * Scalarf8(-1.0), omega);	//if det(H) = 0 use gradient descent, never happened in our tests, could also be removed 

		//instead of clamping just use gradient descent. also works fine and does not require the norm
		Scalarf8 useGD = blend(omega * gradient > Scalarf8(0.0), Scalarf8(1.0), Scalarf8(-1.0));
		omega = Vector3f8::blend(useGD > Scalarf8(0.0), gradient * Scalarf8(-0.125), omega);

		Scalarf8 l_omega2 = omega.squaredNorm();
		const Scalarf8 w = (1.0 - l_omega2) / (1.0 + l_omega2);
		const Vector3f8 vec = omega * (2.0 / (1.0 + l_omega2));
		q = q * Quaternion8f(vec.x(), vec.y(), vec.z(), w);		//no normalization needed because the Cayley map returs a unit quaternion
	}
}


void Elasticity_Kugelstadt2021::saveState(BinaryFileWriter &binWriter)
{
	binWriter.writeBuffer((char*)m_current_to_initial_index.data(), m_current_to_initial_index.size() * sizeof(unsigned int));
	binWriter.writeBuffer((char*)m_initial_to_current_index.data(), m_initial_to_current_index.size() * sizeof(unsigned int));
	binWriter.writeBuffer((char*)m_L.data(), m_L.size() * sizeof(Matrix3r));
	binWriter.writeBuffer((char*)m_rotations.data(), m_rotations.size() * sizeof(Matrix3r));
	binWriter.writeBuffer((char*)m_vDiff.data(), m_vDiff.size() * sizeof(Vector3r));
}

void Elasticity_Kugelstadt2021::loadState(BinaryFileReader &binReader)
{
	binReader.readBuffer((char*)m_current_to_initial_index.data(), m_current_to_initial_index.size() * sizeof(unsigned int));
	binReader.readBuffer((char*)m_initial_to_current_index.data(), m_initial_to_current_index.size() * sizeof(unsigned int));
	binReader.readBuffer((char*)m_L.data(), m_L.size() * sizeof(Matrix3r));
	binReader.readBuffer((char*)m_rotations.data(), m_rotations.size() * sizeof(Matrix3r));
	binReader.readBuffer((char*)m_vDiff.data(), m_vDiff.size() * sizeof(Vector3r));

	// update quaternions which are needed for warmstart
	rotationMatricesToAVXQuaternions();
}

/** Solver for the volume conservation forces (Eq. 30).
 */
void Elasticity_Kugelstadt2021::stepVolumeSolver()
{
	const unsigned int numParticles = m_model->numActiveParticles();
	if (numParticles == 0)
		return;
	const Real dt = TimeManager::getCurrent()->getTimeStepSize();

	// extract current rotations
	START_TIMING("computeRotations");
	computeRotations();
	STOP_TIMING_AVG;

	// precompute some values to improve performance
	START_TIMING("precomputeValues");
	precomputeValues();
	STOP_TIMING_AVG;
	

	START_TIMING("Volume solver")

	//////////////////////////////////////////////////////////////////////////
	// Init linear system solver 
	//////////////////////////////////////////////////////////////////////////
	MatrixReplacement A(3 * m_model->numActiveParticles(), matrixVecProd, (void*)this);

	m_solver.setTolerance(m_maxErrorV);
	m_solver.setMaxIterations(m_maxIterV);
	m_solver.compute(A);

	VectorXr b(3 * numParticles);
	VectorXr x(3 * numParticles);
	VectorXr g(3 * numParticles);

	// Compute right hand side of the system in Eq. 30
	computeRHS(b);

	// warmstart
	#pragma omp parallel for schedule(static) 
	for (int i = 0; i < (int)numParticles; i++)
	{		
		if (m_model->getParticleState(i) == ParticleState::Active)
			g.segment<3>(3 * i) = m_model->getVelocity(i) + m_vDiff[i];
		else
		{
			b.segment<3>(3 * i) = m_model->getVelocity(i);
			g.segment<3>(3 * i) = m_model->getVelocity(i);
		}
	}

	//////////////////////////////////////////////////////////////////////////
	// Solve linear system 
	//////////////////////////////////////////////////////////////////////////
	START_TIMING("Volume - CG solve");
	x = m_solver.solveWithGuess(b, g);
	m_iterationsV = (int)m_solver.iterations();
	STOP_TIMING_AVG;
	INCREASE_COUNTER("Elasticity - CG iterations", static_cast<Real>(m_iterationsV));

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			if (m_model->getParticleState(i) == ParticleState::Active)
			{
				Vector3r& vi = m_model->getVelocity(i);
				const Vector3r new_vi = x.segment<3>(3 * i);
				m_vDiff[i] = new_vi - vi;
				vi = new_vi;
			}
		}
	}
	STOP_TIMING_AVG;
}


#ifdef USE_AVX

/** Matrix vector product used by the matrix-free conjugate gradient 
* solver to solve the system in Eq. 30.
*/
void Elasticity_Kugelstadt2021::matrixVecProd(const Real* vec, Real *result, void *userData)
{
	Simulation *sim = Simulation::getCurrent();
	Elasticity_Kugelstadt2021* elasticity = static_cast<Elasticity_Kugelstadt2021*>(userData);
	FluidModel *model = elasticity->getModel();
	const unsigned int numParticles = model->numActiveParticles();
	const Real dt = TimeManager::getCurrent()->getTimeStepSize();

	const auto &current_to_initial_index = elasticity->m_current_to_initial_index;
	const auto &initial_to_current_index = elasticity->m_initial_to_current_index;
	const auto &initialNeighbors = elasticity->m_initialNeighbors;
	const auto &restVolumes = elasticity->m_restVolumes;
	auto &stress = elasticity->m_stress;
	const auto& precomp_RL_gradW8 = elasticity->m_precomp_RL_gradW8;
	const auto& precomp_RLj_gradW8 = elasticity->m_precomp_RLj_gradW8;
	const auto& precomputed_indices8 = elasticity->m_precomputed_indices8;

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			if (model->getParticleState(i) == ParticleState::Active)
			{
				const unsigned int i0 = current_to_initial_index[i];
				const Vector3r& pi = Eigen::Map<const Vector3r>(&vec[3 * i], 3);
				const unsigned int numNeighbors = (unsigned int)initialNeighbors[i0].size();

				const Vector3f8 pi_avx(pi);

				//////////////////////////////////////////////////////////////////////////
				// compute corotated deformation gradient 
				//////////////////////////////////////////////////////////////////////////
				Scalarf8 trace_avx;
				trace_avx.setZero();
				for (unsigned int j = 0; j < numNeighbors; j += 8)
				{
					const unsigned int count = std::min(numNeighbors - j, 8u);
					unsigned int nIndices[8];
					for (auto k = 0u; k < count; k++)
						nIndices[k] = initial_to_current_index[initialNeighbors[i0][j + k]];
					const Vector3f8 pj_avx = convertVec_zero(nIndices, &vec[0], count);
					const Vector3f8 pj_pi = pj_avx - pi_avx;

					const Vector3f8& correctedRotatedKernel = precomp_RL_gradW8[precomputed_indices8[i] + j / 8];

					// We need the trace of the strain tensor. Therefore, we are only interested in the
					// diagonal elements which are F_ii-1. However, instead of computing the trace
					// of a dyadic product, we can determine the dot product of the vectors to get the trace
					// of F. Then the trace of the strain is the result minus 3. 
					trace_avx += pj_pi.dot(correctedRotatedKernel);
				}
				Real trace = trace_avx.reduce();
				trace *= dt;

				//////////////////////////////////////////////////////////////////////////
				// First Piola Kirchhoff stress = lambda trace(epsilon) I
				//////////////////////////////////////////////////////////////////////////
				stress[i] = elasticity->m_lambda * trace;
			}
			else
				stress[i] = 0.0;
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
				const unsigned int numNeighbors = (unsigned int) initialNeighbors[i0].size();
				const Scalarf8 stress_i_avx(stress[i]);
				const Scalarf8 restVolume_i_avx(restVolumes[i]);

				//////////////////////////////////////////////////////////////////////////
				// Compute elastic force f^l_i in Eq. 30
				//////////////////////////////////////////////////////////////////////////
				Vector3f8 force_avx;
				force_avx.setZero();
				for (unsigned int j = 0; j < numNeighbors; j+=8)
				{
					const unsigned int count = std::min(numNeighbors - j, 8u);

					unsigned int nIndices[8];
					for (auto k = 0u; k < count; k++)
						nIndices[k] = initial_to_current_index[initialNeighbors[i0][j + k]];
	
					const Scalarf8 stress_j_avx = convert_zero(nIndices, &stress[0], count);
					const Scalarf8 restVolume_j_avx = convert_zero(nIndices, &restVolumes[0], count);

					const Vector3f8& V_RL_gradWi = precomp_RL_gradW8[precomputed_indices8[i] + j / 8];
					const Vector3f8& V_RL_gradWj = precomp_RLj_gradW8[precomputed_indices8[i] + j / 8];
					force_avx += (V_RL_gradWi * (restVolume_i_avx * stress_i_avx) + V_RL_gradWj * (restVolume_j_avx * stress_j_avx));
				}

				const Real factor = dt / model->getMass(i);
				result[3 * i] = vec[3 * i] - factor * force_avx.x().reduce();
				result[3 * i + 1] = vec[3 * i + 1] - factor * force_avx.y().reduce();
				result[3 * i + 2] = vec[3 * i + 2] - factor * force_avx.z().reduce();
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

/** Compute the deformation gradient for a particle i.
 */
void Elasticity_Kugelstadt2021::computeF(const unsigned int i, const Real* x, Elasticity_Kugelstadt2021 *e)
{
	const unsigned int i0 = e->m_current_to_initial_index[i];
	const Vector3r& vi = Eigen::Map<const Vector3r>(&x[3 * i], 3);
	const unsigned int numNeighbors = (unsigned int)e->m_initialNeighbors[i0].size();

	const Vector3f8 vi_avx(vi);

	//////////////////////////////////////////////////////////////////////////
	// compute corotated deformation gradient (Eq. 12)
	//////////////////////////////////////////////////////////////////////////
	Matrix3f8 F_avx;
	F_avx.setZero();

	//////////////////////////////////////////////////////////////////////////
	// Fluid
	//////////////////////////////////////////////////////////////////////////
	for (unsigned int j = 0; j < numNeighbors; j += 8)
	{
		const unsigned int count = std::min(numNeighbors - j, 8u);
		unsigned int nIndices[8];
		for (auto k = 0u; k < count; k++)
			nIndices[k] = e->m_initial_to_current_index[e->m_initialNeighbors[i0][j + k]];
		const Vector3f8 vj_avx = convertVec_zero(nIndices, &x[0], count);

		const Vector3f8 vj_vi = vj_avx - vi_avx;
		const Vector3f8& correctedKernel = e->m_precomp_L_gradW8[e->m_precomputed_indices8[i] + j / 8];
		Matrix3f8 dyad;
		dyadicProduct(vj_vi, correctedKernel, dyad);
		F_avx += dyad;
	}
	e->m_F[i] = F_avx.reduce();

	if (Simulation::getCurrent()->is2DSimulation())
		e->m_F[i](2, 2) = 1.0;
}

/** Compute right hand side of the linear system of the volume solver (Eq. 30).
*/
void Elasticity_Kugelstadt2021::computeRHS(VectorXr & rhs)
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int numParticles = m_model->numActiveParticles();
	const unsigned int fluidModelIndex = m_model->getPointSetIndex();
	FluidModel *model = m_model;
	const Real dt = TimeManager::getCurrent()->getTimeStepSize();

	#pragma omp parallel default(shared)
	{
		// update the deformation gradient
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			if (model->getParticleState(i) == ParticleState::Active)
				computeF(i, &m_model->getPosition(0)[0], this);
		}
	
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			if (m_model->getParticleState(i) == ParticleState::Active)
			{
				// trace of R^T F - I
				//const Matrix3r RTF = (m_rotations[i].transpose() * m_F[i]);
				//const Real trace = RTF(0, 0) + RTF(1, 1) + RTF(2, 2) - 3.0;
				
				// short form of: trace(R^T F - I) 
				const Real trace =	m_rotations[i].col(0).dot(m_F[i].col(0)) + 
									m_rotations[i].col(1).dot(m_F[i].col(1)) + 
									m_rotations[i].col(2).dot(m_F[i].col(2)) - static_cast<Real>(3.0);
				
				m_stress[i] = m_lambda * trace;
			}
			else
				m_stress[i] = 0.0;
		}
	}

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			if (model->getParticleState(i) == ParticleState::Active)
			{
				const unsigned int i0 = m_current_to_initial_index[i];
				const unsigned int numNeighbors = (unsigned int)m_initialNeighbors[i0].size();
				const Scalarf8 stress_i_avx(m_stress[i]);
				const Scalarf8 restVolume_i_avx(m_restVolumes[i]);

				//////////////////////////////////////////////////////////////////////////
				// Compute elastic force (Eq. 30)
				//////////////////////////////////////////////////////////////////////////
				Vector3f8 force_avx;
				force_avx.setZero();
				for (unsigned int j = 0; j < numNeighbors; j += 8)
				{
					const unsigned int count = std::min(numNeighbors - j, 8u);
					unsigned int nIndices[8];
					for (auto k = 0u; k < count; k++)
						nIndices[k] = m_initial_to_current_index[m_initialNeighbors[i0][j + k]];

					const Scalarf8 stress_j_avx = convert_zero(nIndices, &m_stress[0], count);
					const Scalarf8 restVolume_j_avx = convert_zero(nIndices, &m_restVolumes[0], count);

					const Vector3f8& correctedRotatedKernel_i = m_precomp_RL_gradW8[m_precomputed_indices8[i] + j / 8];
					const Vector3f8& correctedRotatedKernel_j = m_precomp_RLj_gradW8[m_precomputed_indices8[i] + j / 8];
					const Vector3f8 PWi = correctedRotatedKernel_i * stress_i_avx;
					const Vector3f8 PWj = correctedRotatedKernel_j * stress_j_avx;
					force_avx += (PWi * restVolume_i_avx + PWj * restVolume_j_avx);
				}
				Vector3r force;
				force[0] = force_avx.x().reduce();
				force[1] = force_avx.y().reduce();
				force[2] = force_avx.z().reduce();

				rhs.segment<3>(3 * i) = model->getVelocity(i) + dt * (1.0 / model->getMass(i) * force);
			}
			else
			{
				rhs.segment<3>(3 * i) = model->getVelocity(i);
			}
		}
	}
}

#else

/** Matrix vector product used by the matrix-free conjugate gradient
* solver to solve the system in Eq. 30.
*/
void Elasticity_Kugelstadt2021::matrixVecProd(const Real* vec, Real *result, void *userData)
{
	Simulation *sim = Simulation::getCurrent();
	Elasticity_Kugelstadt2021* elasticity = static_cast<Elasticity_Kugelstadt2021*>(userData);
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
	const auto &precomp_RL_gradW = elasticity->m_precomp_RL_gradW;
	const auto& precomp_RLj_gradW = elasticity->m_precomp_RLj_gradW;
	const auto &precomputed_indices = elasticity->m_precomputed_indices;

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			const unsigned int i0 = current_to_initial_index[i];
			const Vector3r &pi = Eigen::Map<const Vector3r>(&vec[3 * i], 3);
			const Vector3r &xi0 = model->getPosition0(i0);
			const unsigned int numNeighbors = (unsigned int) initialNeighbors[i0].size();

			const Vector3f8 pi_avx(pi);
			const Vector3f8 xi0_avx(xi0);

 			//////////////////////////////////////////////////////////////////////////
 			// compute corotated deformation gradient 
 			//////////////////////////////////////////////////////////////////////////
 
  			//////////////////////////////////////////////////////////////////////////
 			// Fluid
 			//////////////////////////////////////////////////////////////////////////
			Real trace = 0.0;
 			for (unsigned int j = 0; j < numNeighbors; j++)
 			{
 				const unsigned int neighborIndex = initial_to_current_index[initialNeighbors[i0][j]];
 				// get initial neighbor index considering the current particle order 
 				const unsigned int neighborIndex0 = initialNeighbors[i0][j];
 
 				const Vector3r &pj = Eigen::Map<const Vector3r>(&vec[3 * neighborIndex], 3);
 				const Vector3r &xj0 = model->getPosition0(neighborIndex0);
 				const Vector3r pj_pi = pj - pi;
 				const Vector3r xi_xj_0 = xi0 - xj0;
 				//const Vector3r correctedRotatedKernel = RL[i] * sim->gradW(xi_xj_0);
				const Vector3r correctedRotatedKernel = precomp_RL_gradW[precomputed_indices[i] + j];

				// We need the trace of the strain tensor. Therefore, we are only interested in the
				// diagonal elements which are F_ii-1. However, instead of computing the trace
				// of a dyadic product, we can determine the dot product of the vectors to get the trace
				// of F. Then the trace of the strain is the result minus 3. 
				trace += pj_pi.dot(correctedRotatedKernel);
				//nablaU += restVolumes[neighborIndex] * pj_pi * correctedRotatedKernel.transpose();
 			}
			trace *= dt;

			//////////////////////////////////////////////////////////////////////////
			// First Piola Kirchhoff stress = 2 mu epsilon + lambda trace(epsilon) I
			//////////////////////////////////////////////////////////////////////////
			//const Real trace = strain[0] + strain[1] + strain[2];
			stress[i] = m_lambda*trace;
		}
	}

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			if (model->getParticleState(i) != ParticleState::Fixed)
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
					const Vector3r correctedRotatedKernel_i = precomp_RL_gradW[precomputed_indices[i] + j];
					const Vector3r correctedRotatedKernel_j = -precomp_RLj_gradW[precomputed_indices[i] + j];
					const Vector3r PWi = stress[i] * correctedRotatedKernel_i;
					const Vector3r PWj = stress[neighborIndex] * correctedRotatedKernel_j;
					force += (restVolumes[i] * PWi - restVolumes[neighborIndex] * PWj);
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

/** Compute right hand side of the linear system of the volume solver (Eq. 30).
*/
void Elasticity_Kugelstadt2021::computeRHS(VectorXr & rhs)
{
	Simulation *sim = Simulation::getCurrent();
	const unsigned int numParticles = m_model->numActiveParticles();
	const unsigned int fluidModelIndex = m_model->getPointSetIndex();
	FluidModel *model = m_model;
	const Real dt = TimeManager::getCurrent()->getTimeStepSize();

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
 			// compute corotated deformation gradient (Eq. 12)
 			//////////////////////////////////////////////////////////////////////////
			//m_F[i].setZero();
 
  			//////////////////////////////////////////////////////////////////////////
 			// Fluid
 			//////////////////////////////////////////////////////////////////////////
			Real trace = 0.0;
 			for (unsigned int j = 0; j < numNeighbors; j++)
 			{
 				const unsigned int neighborIndex = m_initial_to_current_index[m_initialNeighbors[i0][j]];
 				// get initial neighbor index considering the current particle order 
 				const unsigned int neighborIndex0 = m_initialNeighbors[i0][j];
 
 				const Vector3r &xj = model->getPosition(neighborIndex);
 				//const Vector3r &xj0 = m_model->getPosition0(neighborIndex0);
 				const Vector3r xj_xi = xj - xi;
 				//const Vector3r xi_xj_0 = xi0 - xj0;
 				//const Vector3r correctedRotatedKernel = m_RL[i] * sim->gradW(xi_xj_0);
				const Vector3r correctedRotatedKernel = m_precomp_RL_gradW[m_precomputed_indices[i] + j];

 				//m_F[i] += m_restVolumes[neighborIndex] * xj_xi * correctedRotatedKernel.transpose();

				// We need the trace of the strain tensor. Therefore, we are only interested in the
				// diagonal elements which are F_ii-1. However, instead of computing the trace
				// of a dyadic product, we can determine the dot product of the vectors to get the trace
				// of F. Then the trace of the strain is the result minus 3. 
				trace += m_restVolumes[neighborIndex] * xj_xi.dot(correctedRotatedKernel);
 			}

			/*if (sim->is2DSimulation())
				m_F[i](2, 2) = 1.0;*/

 			//////////////////////////////////////////////////////////////////////////
 			// compute Cauchy strain: epsilon = 0.5 (F + F^T) - I
 			//////////////////////////////////////////////////////////////////////////
 			//Vector3r strain;
 			//strain[0] = m_F[i](0, 0) - static_cast<Real>(1.0);						// \epsilon_{00}
 			//strain[1] = m_F[i](1, 1) - static_cast<Real>(1.0);						// \epsilon_{11}
 			//strain[2] = m_F[i](2, 2) - static_cast<Real>(1.0);						// \epsilon_{22}

																					//////////////////////////////////////////////////////////////////////////
			// First Piola Kirchhoff stress = 2 mu epsilon + lambda trace(epsilon) I
			//////////////////////////////////////////////////////////////////////////
			//const Real trace = strain[0] + strain[1] + strain[2];
			trace -= 3.0;
			m_stress[i] = m_lambda*trace;
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
				const Vector3r gradW = sim->gradW(xi0 - xj0);
				//const Vector3r correctedRotatedKernel_i = m_RL[i] * gradW;
				const Vector3r correctedRotatedKernel_i = m_precomp_RL_gradW[m_precomputed_indices[i] + j];
				const Vector3r correctedRotatedKernel_j = -m_precomp_RLj_gradW[m_precomputed_indices[i] + j];
				const Vector3r PWi = m_stress[i] * correctedRotatedKernel_i;
				const Vector3r PWj = m_stress[neighborIndex] * correctedRotatedKernel_j;
				force += m_restVolumes[i] * m_restVolumes[neighborIndex] * (PWi - PWj);
			}

			rhs.segment<3>(3 * i) = model->getVelocity(i) + dt * (model->getAcceleration(i) + 1.0 / model->getMass(i) * force);
		}
	}
}

#endif

/** Precompute some values and products to improve the performance of the solvers.
*/
void Elasticity_Kugelstadt2021::precomputeValues()
{
	Simulation* sim = Simulation::getCurrent();
	const unsigned int fluidModelIndex = m_model->getPointSetIndex();

	const int numParticles = (int)m_model->numActiveParticles();

#ifdef USE_AVX
	m_precomputed_indices8.clear();
	m_precomp_L_gradW8.clear();
	m_precomp_RL_gradW8.clear();
	m_precomp_RLj_gradW8.clear();
	m_precomputed_indices8.resize(numParticles);
#else
	m_precomputed_indices.clear();
	m_precomp_RL_gradW.clear();
	m_precomp_RLj_gradW.clear();
	m_precomputed_indices.resize(numParticles);
#endif


	unsigned int sumNeighborParticles = 0;
	unsigned int sumNeighborParticles8 = 0;
	for (int i = 0; i < numParticles; i++)
	{
		const unsigned int i0 = m_current_to_initial_index[i];
		const size_t numNeighbors = m_initialNeighbors[i0].size();

#ifdef USE_AVX
		m_precomputed_indices8[i] = sumNeighborParticles8;

		// steps of 8 values due to avx
		sumNeighborParticles8 += (unsigned int) numNeighbors / 8u;
		if (numNeighbors % 8 != 0)
			sumNeighborParticles8++;
#else
		m_precomputed_indices[i] = sumNeighborParticles;
		sumNeighborParticles += numNeighbors;
#endif
	}

#ifdef USE_AVX
	m_precomp_L_gradW8.resize(sumNeighborParticles8);
	m_precomp_RL_gradW8.resize(sumNeighborParticles8);
	m_precomp_RLj_gradW8.resize(sumNeighborParticles8);
#else
	m_precomp_RL_gradW.resize(sumNeighborParticles);
	m_precomp_RLj_gradW.resize(sumNeighborParticles);
#endif

	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static) 
		for (int i = 0; i < (int)numParticles; i++)
		{
			m_RL[i] = m_rotations[i] * m_L[i];
		}
		#pragma omp for schedule(static) 
		for (int i = 0; i < (int)numParticles; i++)
		{
			const unsigned int i0 = m_current_to_initial_index[i];
			const Vector3r& xi0 = m_model->getPosition0(i0);
			const Vector3f8 xi0_avx(xi0);
			const unsigned int numNeighbors = (unsigned int)m_initialNeighbors[i0].size();
			Scalarf8 restVolume_i_avx(m_restVolumes[i]);

#ifdef USE_AVX
			unsigned int base8 = m_precomputed_indices8[i];
			unsigned int idx = 0;
			Matrix3f8 L_avx(m_L[i]);
			Matrix3f8 RL_avx(m_RL[i]);
			for (unsigned int j = 0; j < numNeighbors; j += 8)
			{
				const unsigned int count = std::min(numNeighbors - j, 8u);

				unsigned int nIndices[8];
				for (auto k = 0u; k < count; k++)
					nIndices[k] = m_initial_to_current_index[m_initialNeighbors[i0][j + k]];

				const Scalarf8 restVolume_j_avx = convert_zero(nIndices, &m_restVolumes[0], count);
				const Vector3f8 xj0_avx = convertVec_zero(&m_initialNeighbors[i0][j], &m_model->getPosition0(0), count);
				const Matrix3f8 RLj_avx = convertMat_zero(nIndices, &m_RL[0], count);

				const Vector3f8 gradW = CubicKernel_AVX::gradW(xi0_avx - xj0_avx);
				m_precomp_L_gradW8[base8 + idx] = L_avx * gradW * restVolume_j_avx;
				m_precomp_RL_gradW8[base8 + idx] = RL_avx * gradW * restVolume_j_avx;
				m_precomp_RLj_gradW8[base8 + idx] = RLj_avx * gradW * restVolume_i_avx;
				idx++;
			}
#else
			unsigned int base = m_precomputed_indices[i];

			for (unsigned int j = 0; j < numNeighbors; j++)
			{
				// get initial neighbor index considering the current particle order 
				const unsigned int neighborIndex0 = m_initialNeighbors[i0][j];
				const unsigned int neighborIndex = m_initial_to_current_index[m_initialNeighbors[i0][j]];
				const Vector3r& xj0 = m_model->getPosition0(neighborIndex0);
				Vector3r gradW = sim->gradW(xi0 - xj0);
				m_precomp_RL_gradW[base + j] = m_restVolumes[neighborIndex] * m_RL[i] * gradW;
				m_precomp_RLj_gradW[base + j] = -m_restVolumes[i] * m_RL[neighborIndex] * gradW;
			}
#endif
		}
	}
}



