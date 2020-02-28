#ifndef __PBDWrapper_h__
#define __PBDWrapper_h__

#include "SPlisHSPlasH/Common.h"
#include <string>
#include "Simulation/SimulationModel.h"
#include "Simulation/CubicSDFCollisionDetection.h"

namespace PBD
{
	class TimeStepController;
}


class PBDWrapper
{
protected:
	PBD::SimulationModel m_model;
	PBD::CubicSDFCollisionDetection m_cd;
	PBD::TimeStepController *m_timeStep;

	short m_clothSimulationMethod = 2;
	short m_solidSimulationMethod = 2;
	short m_bendingMethod = 2;
	std::string m_sceneName;
	std::string m_sceneFileName;
	bool m_enableMayaExport = false;
	Real m_dampingCoeff = 0.0;

public:
	struct RBData 
	{
		Vector3r x;
		Matrix3r R;
		Vector3r scale;
		std::string objFile;
		int collisionType;
		Real restitution;
		Real friction;
	};

	PBDWrapper();
	~PBDWrapper();

	void reset();
	void initModel(const Real timeStepSize);

	/** Read rigid body scene and create the rigid body model
	*/
	void readScene(const std::string &sceneFileName, const std::vector< RBData> &additionalRigidBodies);
	void initTriangleModelConstraints();
	void initTetModelConstraints();

	void timeStep();
	void updateVisModels();

	void loadObj(const std::string &filename, PBD::VertexData &vd, Utilities::IndexedFaceMesh &mesh, const Vector3r &scale);

	PBD::SimulationModel &getSimulationModel() { return m_model; }
	PBD::DistanceFieldCollisionDetection &getCollisionDetection() { return m_cd; }
	PBD::TimeStepController &getTimeStepController();
	
	Real getDampingCoeff() const { return m_dampingCoeff; }
	void setDampingCoeff(Real val) { m_dampingCoeff = val; }
	short getClothSimulationMethod() const { return m_clothSimulationMethod; }
	void setClothSimulationMethod(short val) { m_clothSimulationMethod = val; }
	short getSolidSimulationMethod() const { return m_solidSimulationMethod; }
	void setSolidSimulationMethod(short val) { m_solidSimulationMethod = val; }
	short getBendingMethod() const { return m_bendingMethod; }
	void setBendingMethod(short val) { m_bendingMethod = val; }
};



#endif
