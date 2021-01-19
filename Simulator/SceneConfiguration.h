#ifndef __SceneConfiguration_h__
#define __SceneConfiguration_h__

#include "SPlisHSPlasH/Utilities/SceneLoader.h"

namespace SPH
{
	/** \brief Class to store the scene configuration that is imported from the scene file.
	*/
	class SceneConfiguration
	{
	private:
		static SceneConfiguration *m_current;

	protected:
		Utilities::SceneLoader::Scene m_scene;
		std::string m_sceneFile;

	public:
		SceneConfiguration();
		SceneConfiguration(const SceneConfiguration&) = delete;
		SceneConfiguration& operator=(const SceneConfiguration&) = delete;
		~SceneConfiguration();

		// Singleton
		static SceneConfiguration* getCurrent ();
		static void setCurrent (SceneConfiguration* sc);
		static bool hasCurrent();

		void setSceneFile(const std::string& file) { m_sceneFile = file; }
		const std::string& getSceneFile() const { return m_sceneFile; }

		Utilities::SceneLoader::Scene& getScene() { return m_scene; }

	};
}

#endif
