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

	short m_clothSimulationMethod = 2;
	short m_solidSimulationMethod = 2;
	short m_bendingMethod = 2;
	Real m_distanceStiffness;
	Real m_xxStiffness;
	Real m_yyStiffness;
	Real m_xyStiffness;
	Real m_xyPoissonRatio;
	Real m_yxPoissonRatio;
	bool m_normalizeStretch;
	bool m_normalizeShear;
	Real m_bendingStiffness;
	Real m_solidStiffness;
	Real m_solidPoissonRatio;
	bool m_solidNormalizeStretch;
	bool m_solidNormalizeShear;
	Real m_volumeStiffness;
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
	Real getBendingStiffness() const { return m_bendingStiffness; }
	void setBendingStiffness(Real val) { m_bendingStiffness = val; }
	Real getDistanceStiffness() const { return m_distanceStiffness; }
	void setDistanceStiffness(Real val) { m_distanceStiffness = val; }
	Real getXXStiffness() const { return m_xxStiffness; }
	void setXXStiffness(Real val) { m_xxStiffness = val; }
	Real getYYStiffness() const { return m_yyStiffness; }
	void setYYStiffness(Real val) { m_yyStiffness = val; }
	Real getXYStiffness() const { return m_xyStiffness; }
	void setXYStiffness(Real val) { m_xyStiffness = val; }
	Real getXYPoissonRatio() const { return m_xyPoissonRatio; }
	void setXYPoissonRatio(Real val) { m_xyPoissonRatio = val; }
	Real getYXPoissonRatio() const { return m_yxPoissonRatio; }
	void setYXPoissonRatio(Real val) { m_yxPoissonRatio = val; }
	bool getNormalizeStretch() const { return m_normalizeStretch; }
	void setNormalizeStretch(bool val) { m_normalizeStretch = val; }
	bool getNormalizeShear() const { return m_normalizeShear; }
	void setNormalizeShear(bool val) { m_normalizeShear = val; }
	Real getSolidStiffness() const { return m_solidStiffness; }
	void setSolidStiffness(Real val) { m_solidStiffness = val; }
	Real getSolidPoissonRatio() const { return m_solidPoissonRatio; }
	void setSolidPoissonRatio(Real val) { m_solidPoissonRatio = val; }
	bool getSolidNormalizeStretch() const { return m_solidNormalizeStretch; }
	void setSolidNormalizeStretch(bool val) { m_solidNormalizeStretch = val; }
	bool getSolidNormalizeShear() const { return m_solidNormalizeShear; }
	void setSolidNormalizeShear(bool val) { m_solidNormalizeShear = val; }
	Real getVolumeStiffness() const { return m_volumeStiffness; }
	void setVolumeStiffness(Real val) { m_volumeStiffness = val; }
};



#endif
