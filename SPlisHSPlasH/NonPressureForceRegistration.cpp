#include "Simulation.h"

#include "SPlisHSPlasH/Drag/DragForce_Macklin2014.h"
#include "SPlisHSPlasH/Drag/DragForce_Gissler2017.h"

#include "SPlisHSPlasH/Viscosity/Viscosity_Standard.h"
#include "SPlisHSPlasH/Viscosity/Viscosity_Bender2017.h"
#include "SPlisHSPlasH/Viscosity/Viscosity_Peer2015.h"
#include "SPlisHSPlasH/Viscosity/Viscosity_Peer2016.h"
#include "SPlisHSPlasH/Viscosity/Viscosity_Takahashi2015.h"
#include "SPlisHSPlasH/Viscosity/Viscosity_Weiler2018.h"

#include "SPlisHSPlasH/Vorticity/VorticityConfinement.h"
#include "SPlisHSPlasH/Vorticity/MicropolarModel_Bender2017.h"
#include "SPlisHSPlasH/Vorticity/VorticityRefinement_Liu2021.h"

#include "Elasticity/Elasticity_Becker2009.h"
#include "Elasticity/Elasticity_Peer2018.h"
#include "Elasticity/Elasticity_Kugelstadt2021.h"

#include "SurfaceTension/SurfaceTension_Becker2007.h"
#include "SurfaceTension/SurfaceTension_Akinci2013.h"
#include "SurfaceTension/SurfaceTension_He2014.h"
#include "SurfaceTension/SurfaceTension_Jeske2023.h"
#ifdef USE_THIRD_PARTY_METHODS
#include "SurfaceTension/SurfaceTension_ZorillaRitter2020.h"
#endif


using namespace SPH;

void Simulation::registerNonpressureForces()
{
	addDragMethod("None", [](FluidModel*) -> NonPressureForceBase* { return nullptr; });
	addDragMethod(DragForce_Macklin2014::METHOD_NAME, DragForce_Macklin2014::creator);
	addDragMethod(DragForce_Gissler2017::METHOD_NAME, DragForce_Gissler2017::creator);

	addElasticityMethod("None", [](FluidModel*) -> NonPressureForceBase* { return nullptr; });
	addElasticityMethod(Elasticity_Becker2009::METHOD_NAME, Elasticity_Becker2009::creator);
	addElasticityMethod(Elasticity_Peer2018::METHOD_NAME, Elasticity_Peer2018::creator);
	addElasticityMethod(Elasticity_Kugelstadt2021::METHOD_NAME, Elasticity_Kugelstadt2021::creator);

	addSurfaceTensionMethod("None", [](FluidModel*) -> NonPressureForceBase* { return nullptr; });
	addSurfaceTensionMethod(SurfaceTension_Becker2007::METHOD_NAME, SurfaceTension_Becker2007::creator);
	addSurfaceTensionMethod(SurfaceTension_Akinci2013::METHOD_NAME, SurfaceTension_Akinci2013::creator);
	addSurfaceTensionMethod(SurfaceTension_He2014::METHOD_NAME, SurfaceTension_He2014::creator);
	addSurfaceTensionMethod(SurfaceTension_Jeske2023::METHOD_NAME, SurfaceTension_Jeske2023::creator);
#ifdef USE_THIRD_PARTY_METHODS
	addSurfaceTensionMethod(SurfaceTension_ZorillaRitter2020::METHOD_NAME, SurfaceTension_ZorillaRitter2020::creator);
#endif

	addViscosityMethod("None", [](FluidModel*) -> NonPressureForceBase* { return nullptr; });
	addViscosityMethod(Viscosity_Standard::METHOD_NAME, Viscosity_Standard::creator);
	addViscosityMethod(Viscosity_Bender2017::METHOD_NAME, Viscosity_Bender2017::creator);
	addViscosityMethod(Viscosity_Peer2015::METHOD_NAME, Viscosity_Peer2015::creator);
	addViscosityMethod(Viscosity_Peer2016::METHOD_NAME, Viscosity_Peer2016::creator);
	addViscosityMethod(Viscosity_Takahashi2015::METHOD_NAME, Viscosity_Takahashi2015::creator);
	addViscosityMethod(Viscosity_Weiler2018::METHOD_NAME, Viscosity_Weiler2018::creator);

	addVorticityMethod("None", [](FluidModel*) -> NonPressureForceBase* { return nullptr; });
	addVorticityMethod(MicropolarModel_Bender2017::METHOD_NAME, MicropolarModel_Bender2017::creator);
	addVorticityMethod(VorticityConfinement::METHOD_NAME, VorticityConfinement::creator);
	addVorticityMethod(VorticityRefinement_Liu2021::METHOD_NAME, VorticityRefinement_Liu2021::creator);
}
