# Implementing a new particle/rigid body data exporter

All exporters are implemented in the same way in SPlisHSPlasH. In the following we explain the implementation of such an exporter method using as example a new rigid body exporter. 

SPlisHSPlasH organizes the exporters in `/Simulator/Exporter/` and thus any changes or additions are intended to take place in this directory. The user can add new data exporters by creating new or copying and modifying existing exporter class files and registering these inside the build system and the source code.

## Creating a new class

If you want to create a new exporter class from scratch, you should take a look at existing exporters in SPlisHSPlasH. In short, every exporter inherits from the base class `ExporterBase`. A minimal working derived class would look like this:

**RigidBodyExporter_MyFormat.h**

```cpp
#ifndef __RigidBodyExporter_MyFormat_h__
#define __RigidBodyExporter_MyFormat_h__

#include "ExporterBase.h"

namespace SPH
{
	/** \brief Rigid body exporter for the OBJ format.
	*/
	class RigidBodyExporter_MyFormat : public ExporterBase
	{
	protected: 
		bool m_isFirstFrame;
		std::string m_exportPath;

	public:
		RigidBodyExporter_MyFormat(SimulatorBase* base);
		RigidBodyExporter_MyFormat(const RigidBodyExporter_MyFormat&) = delete;
         RigidBodyExporter_MyFormat& operator=(const RigidBodyExporter_MyFormat&) = delete;
		virtual ~RigidBodyExporter_MyFormat(void);

		virtual void init(const std::string& outputPath);
		virtual void step(const unsigned int frame);
		virtual void reset();
		virtual void setActive(const bool active); 
	};
}
#endif
```

**RigidBodyExporter_MyFormat.cpp**

```cpp
#include "RigidBodyExporter_MyFormat.h"
#include <Utilities/Logger.h>
#include <Utilities/FileSystem.h>
#include "SPlisHSPlasH/Simulation.h"

using namespace SPH;
using namespace Utilities;

RigidBodyExporter_MyFormat::RigidBodyExporter_MyFormat(SimulatorBase* base) :
	ExporterBase(base)
{
	m_isFirstFrame = true;
}

RigidBodyExporter_MyFormat::~RigidBodyExporter_MyFormat(void)
{
}

void RigidBodyExporter_MyFormat::init(const std::string& outputPath)
{
    // define output path for the data
	m_exportPath = FileSystem::normalizePath(outputPath + "/my_format");
}

void RigidBodyExporter_MyFormat::step(const unsigned int frame)
{
    // check if the exporter is active
	if (!m_active)
		return;

	// check if we have a static model
	bool isStatic = true;
	for (unsigned int i = 0; i < sim->numberOfBoundaryModels(); i++)
	{
		BoundaryModel* bm = sim->getBoundaryModel(i);
		if (bm->getRigidBodyObject()->isDynamic())
		{
			isStatic = false;
			break;
		}
	}

    // If we have a static model, write the data only for the first frame.
    // Otherwise each frame is exported.
	if (m_isFirstFrame || !isStatic)
	{
        [...]
    }
   	m_isFirstFrame = false;
}

void RigidBodyExporter_MyFormat::reset()
{
	m_isFirstFrame = true;
}

void RigidBodyExporter_MyFormat::setActive(const bool active)
{
	ExporterBase::setActive(active);
    // create output folder
	if (m_active)
		FileSystem::makeDirs(m_exportPath);
}
```

including the following:

* a constructor with `SimulatorBase*` as the sole parameter `RigidBodyExporter_MyFormat(SimulatorBase* base)`
* a `init(const std::string& outputPath)` method which should define the export path, `outputPath` contains the path of the current output directory of SPlisHSPlasH
* a step function `void step()` called for each frame that should be exported
* a reset function `void reset()` called on every reset of the simulation
* a function `void setActive(const bool active)` which is called when the user activates the exporter

In our example the exporter path is defined int the function `init`. When the user activates the exporter, e.g. in the GUI, the corresponding directory is created. In the function `step` the rigid body data can be written. Note that we added some code so that static rigid bodies are only exported once and not for each frame since they never change. 

## Registering the exporter

To add our new exporter, we have to integrate it into the build process and the source code. 

### Adding to the build process

Simply add the class files `RigidBodyExporter_MyFormat.h` and `RigidBodyExporter_MyFormat.cpp` to the `CMakeLists.txt` in the `/Simulator/` directory. This can be done by adding the relative file paths to the respective variables `EXPORTER_HEADRE_FILES` and `EXPORTER_SOURCE_FILES`:

```cmake
set(EXPORTER_HEADER_FILES
	[...]
	Exporter/RigidBodyExporter_MyFormat.h
)

set(EXPORTER_SOURCE_FILES
	[...]
	Exporter/RigidBodyExporter_MyFormat.cpp
)
```

### Integration in the source code

Any exporter is registered in the file `ExporterRegistration.cpp`, which can be found in the `/Simulator/` directory. Adding our new exporter is done by adding the following line to the function `void SimulatorBase::createExporters()`:

```cpp 
addRigidBodyExporter("enableRigidBodyMyFormatExport", "Rigid Body MyFormat Exporter", "Enable/disable rigid body My Format export.", new RigidBodyExporter_MyFormat(this));
```

and including `Exporter/RigidBodyExporter_MyFormat.h`. The first string defines a key which can be used in the json scene files to activate your exporter. The second string defines the name of your exporter which will appear in the GUI. This name can also be used to activate your exporter in C++ or Python. The last string contains a description of the exporter which is used as tool tip in the GUI. 

After these additions and building SPlisHSPlasH, our new exporter is available inside the simulation.

## Implementing a new exporter in Python

You can also implement a new exporter using our Python interface. You can find an example here: `pySPlisHSPlasH\examples\custom_exporter.py`.