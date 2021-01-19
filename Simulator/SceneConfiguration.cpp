#include "SceneConfiguration.h"

using namespace SPH;
using namespace std;

SceneConfiguration* SceneConfiguration::m_current = nullptr;

SceneConfiguration::SceneConfiguration () 
{
	m_sceneFile = "";
}

SceneConfiguration::~SceneConfiguration () 
{
	for (unsigned int i = 0; i < m_scene.boundaryModels.size(); i++)
		delete m_scene.boundaryModels[i];
	m_scene.boundaryModels.clear();

	for (unsigned int i = 0; i < m_scene.fluidModels.size(); i++)
		delete m_scene.fluidModels[i];
	m_scene.fluidModels.clear();

	for (unsigned int i = 0; i < m_scene.fluidBlocks.size(); i++)
		delete m_scene.fluidBlocks[i];
	m_scene.fluidBlocks.clear();

	for (unsigned int i = 0; i < m_scene.materials.size(); i++)
		delete m_scene.materials[i];
	m_scene.materials.clear();

	for (unsigned int i = 0; i < m_scene.emitters.size(); i++)
		delete m_scene.emitters[i];
	m_scene.emitters.clear();

	for (unsigned int i = 0; i < m_scene.animatedFields.size(); i++)
		delete m_scene.animatedFields[i];
	m_scene.animatedFields.clear();

	m_current = nullptr;
}

SceneConfiguration* SceneConfiguration::getCurrent ()
{
	if (m_current == nullptr)
	{
		m_current = new SceneConfiguration ();
	}
	return m_current;
}

void SceneConfiguration::setCurrent (SceneConfiguration* sc)
{
	m_current = sc;
}

bool SceneConfiguration::hasCurrent()
{
	return (m_current != nullptr);
}

