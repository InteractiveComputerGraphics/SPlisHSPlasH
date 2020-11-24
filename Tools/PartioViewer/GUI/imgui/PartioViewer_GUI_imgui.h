#ifndef __PartioViewer_GUI_imgui_h__
#define __PartioViewer_GUI_imgui_h__

#include "SPlisHSPlasH/Common.h"
#include "../PartioViewer_GUI_Base.h"
#include <vector>

namespace SPH 
{	
	class PartioViewer;

	class PartioViewer_GUI_imgui : public PartioViewer_GUI_Base
	{
		public:
			PartioViewer_GUI_imgui(PartioViewer *viewer);
			virtual ~PartioViewer_GUI_imgui();

		protected:
			PartioViewer *m_viewer;
			unsigned int m_currentFluidModel;			
			bool m_renderWalls;
			bool m_showBBox;
			Vector3r m_camPos;
			Vector3r m_camLookat;
			std::map<int, int> m_mapColorField2Attr;

			void initImgui();
			void initImguiParameters();
			void renderScene();
			void determineMinMaxValues();

			static void selection(const Vector2i &start, const Vector2i &end, void *clientData);				   

			void destroy();
			void createSimulationParameterGUI();

		public:
			virtual void init();
			virtual void initParameterGUI();
			virtual void addOptions(cxxopts::Options &options);
			virtual void parseOptions(cxxopts::ParseResult &result);
			virtual void render();
			virtual void update();
			virtual void cleanup();
			virtual void run();
			virtual void stop();
			virtual unsigned int getWidth() const;
			virtual unsigned int getHeight() const;
	};
}

#endif