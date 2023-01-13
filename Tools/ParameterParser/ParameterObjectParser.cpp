#include "ParameterObjectParser.h"
#include "Simulator/SimulatorBase.h"
#include "SPlisHSPlasH/Simulation.h"
#include "SPlisHSPlasH/TimeManager.h"
#include "SPlisHSPlasH/Viscosity/ViscosityBase.h"
#include "SPlisHSPlasH/Drag/DragBase.h"
#include "SPlisHSPlasH/SurfaceTension/SurfaceTensionBase.h"
#include "SPlisHSPlasH/Vorticity/VorticityBase.h"
#include "SPlisHSPlasH/Elasticity/ElasticityBase.h"
#include "SPlisHSPlasH/XSPH.h"
#include <vector>

using namespace GenParam;
using namespace SPH;
using namespace Eigen;
using namespace std;
using namespace Utilities;


void ParameterObjectParser::parseParameters()
{
	// generate SimulatorBase object
	SimulatorBase* base = new SimulatorBase();
	Utilities::logger.removeSink(Utilities::logger.getSinks()[0]);		// avoid output
	base->initParameters();

	// generate Simulation object
	Simulation* sim = Simulation::getCurrent();
	sim->init(0.025, false);

	// generate TimeManager object
	TimeManager* tm = TimeManager::getCurrent();
	tm->initParameters();

	// add dummy fluid model with standard ID "Fluid"
	std::vector<Vector3r> x(1, Vector3r::Zero());
	std::vector<Vector3r> v(1, Vector3r::Zero());
	std::vector<unsigned int> id(1, 0);
	sim->addFluidModel("Fluid", (unsigned int) x.size(), x.data(), v.data(), id.data(), 100u);
	FluidModel* fluidModel = sim->getFluidModel(0);
	sim->setSimulationMethod(4);

	// parse SimulatorBase object
	parseParameterObject(base);
	// parse TimeManager object
	parseParameterObject(tm);
	// parse Simulation object
	parseParameterObject(sim);
	// parse all pressure solvers
	for (unsigned int i = 0; i < static_cast<int>(SimulationMethods::NumSimulationMethods); i++)
	{
		sim->setSimulationMethod(i);
		parseParameterObject(sim->getTimeStep());
	}

	// parse FluidModel
	parseParameterObject(fluidModel);
	// parse XSPH method
	parseParameterObject(fluidModel->getXSPH());

	// Drag
	auto& dragMethods = sim->getDragMethods();
	for (unsigned int i = 0; i < dragMethods.size(); i++)
	{
		fluidModel->setDragMethod(i);
		parseParameterObject(fluidModel->getDragBase());
	}

	// Surface tension
	auto& stMethods = sim->getSurfaceTensionMethods();
	for (unsigned int i = 0; i < stMethods.size(); i++)
	{
		fluidModel->setSurfaceTensionMethod(i);
		parseParameterObject(fluidModel->getSurfaceTensionBase());
	}

	// Visco 
	auto& viscoMethods = sim->getViscosityMethods();
	for (unsigned int i = 0; i < viscoMethods.size(); i++)
	{
		fluidModel->setViscosityMethod(i);
		parseParameterObject(fluidModel->getViscosityBase());
	}

	// Vorticity 
	auto& vorticiyMethods = sim->getVorticityMethods();
	for (unsigned int i = 0; i < vorticiyMethods.size(); i++)
	{
		fluidModel->setVorticityMethod(i);
		parseParameterObject(fluidModel->getVorticityBase());
	}

	// Elasticity 
	auto& elasticityMethods = sim->getElasticityMethods();
	for (unsigned int i = 0; i < elasticityMethods.size(); i++)
	{
		fluidModel->setElasticityMethod(i);
		parseParameterObject(fluidModel->getElasticityBase());
	}

	// cleanup
	delete Simulation::getCurrent();
	delete base;

	// parse example fluid block
	FluidBlockParameterObject block;
	block.initParameters();
	parseParameterObject(&block);

	// write example fluid model
	FluidModelParameterObject fm;
	fm.initParameters();
	parseParameterObject(&fm);

	// write example emitter
	EmitterParameterObject em;
	em.initParameters();
	parseParameterObject(&em);

	// write example emitter
	AnimationFieldParameterObject af;
	af.initParameters();
	parseParameterObject(&af);

	// write example material
	MaterialParameterObject mat;
	mat.initParameters();
	parseParameterObject(&mat);

	// write example rigid body
	BoundaryParameterObject bpo;
	bpo.initParameters();
	parseParameterObject(&bpo);
}


void ParameterObjectParser::parseParameterObject(ParameterObject* paramObj)
{
	if (paramObj == nullptr)
		return;

	for (auto cb : m_paramObjCB)
		if (cb)
			cb(paramObj);

	const unsigned int numParams = paramObj->numParameters();
	for (unsigned int i = 0; i < numParams; i++)
	{
		ParameterBase* paramBase = paramObj->getParameter(i);
		if (paramBase->getType() == RealParameterType)
		{
			for (auto cb : m_realParamCB)
				if (cb)
					cb(static_cast<NumericParameter<Real>*>(paramBase));
		}
		else if (paramBase->getType() == ParameterBase::UINT32)
		{
			for (auto cb : m_uint32ParamCB)
				if (cb)
					cb(static_cast<NumericParameter<unsigned int>*>(paramBase));
		}
		else if (paramBase->getType() == ParameterBase::UINT16)
		{
			for (auto cb : m_uint16ParamCB)
				if (cb)
					cb(static_cast<NumericParameter<unsigned short>*>(paramBase));
		}
		else if (paramBase->getType() == ParameterBase::UINT8)
		{
			for (auto cb : m_uint8ParamCB)
				if (cb)
					cb(static_cast<NumericParameter<unsigned char>*>(paramBase));
		}
		else if (paramBase->getType() == ParameterBase::INT32)
		{
			for (auto cb : m_int32ParamCB)
				if (cb)
					cb(static_cast<NumericParameter<int>*>(paramBase));
		}
		else if (paramBase->getType() == ParameterBase::INT16)
		{
			for (auto cb : m_int16ParamCB)
				if (cb)
					cb(static_cast<NumericParameter<short>*>(paramBase));
		}
		else if (paramBase->getType() == ParameterBase::INT8)
		{
			for (auto cb : m_int8ParamCB)
				if (cb)
					cb(static_cast<NumericParameter<char>*>(paramBase));
		}
		else if (paramBase->getType() == ParameterBase::ENUM)
		{
			for (auto cb : m_enumParamCB)
				if (cb)
					cb(static_cast<EnumParameter*>(paramBase));
		}
		else if (paramBase->getType() == ParameterBase::BOOL)
		{
			for (auto cb : m_boolParamCB)
				if (cb)
					cb(static_cast<BoolParameter*>(paramBase));
		}
		else if (paramBase->getType() == RealVectorParameterType)
		{
			for (auto cb : m_vecRealParamCB)
				if (cb)
					cb(static_cast<RealVectorParameter*>(paramBase));
		}
		else if (paramBase->getType() == ParameterBase::STRING)
		{
			for (auto cb : m_stringParamCB)
				if (cb)
					cb(static_cast<StringParameter*>(paramBase));
		}
		else if (paramBase->getType() == ParameterBase::VEC_UINT32)
		{
			for (auto cb : m_vecUintParamCB)
				if (cb)
					cb(static_cast<VectorParameter<unsigned int>*>(paramBase));
		}
	}
	LOG_INFO << "-----------------------------------------------------";
}