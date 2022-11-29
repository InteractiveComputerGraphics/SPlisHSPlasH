#ifndef __Simulator_GUI_imgui_h__
#define __Simulator_GUI_imgui_h__

#include "SPlisHSPlasH/Common.h"
#include "../Simulator_GUI_Base.h"
#include <vector>

struct ImFont;
struct ImGuiContext;
struct ImGuiSettingsHandler;
struct ImGuiTextBuffer;

namespace SPH 
{	
	struct UserSettings
	{
		int scaleIndex;
		int win_x, win_y;
		int win_width, win_height;

		UserSettings() { win_x = 0; win_y = 0; win_width = 1280; win_height = 960; scaleIndex = 0; }
	};

	class Simulator_GUI_imgui : public Simulator_GUI_Base
	{
		public:
			Simulator_GUI_imgui(SimulatorBase *base);
			virtual ~Simulator_GUI_imgui();

		protected:
			unsigned int m_currentFluidModel;
			Vector3r m_oldMousePos;
			std::vector<std::string> m_colorFieldNames;
			std::vector<std::vector<unsigned int>> m_selectedParticles;			
			std::vector<ImFont*> m_fonts;
			std::vector<float> m_scales;
			unsigned int m_currentScaleIndex;
			ImGuiContext* m_context;
			UserSettings m_userSettings;

			std::vector<std::vector<unsigned int>>& getSelectedParticles() { return m_selectedParticles; }
			void initImgui();
			void initStyle();
			void initImguiParameters();
			void renderBoundary();
			void particleInfo();

			static void selection(const Vector2i &start, const Vector2i &end, void *clientData);
			static void mouseMove(int x, int y, void *clientData);

			void switchPause();
			static void switchDrawMode();
				
			void destroy();

			static void writeIni(ImGuiContext* ctx, ImGuiSettingsHandler* handler, ImGuiTextBuffer* out_buf);
			static void readIni(ImGuiContext* ctx, ImGuiSettingsHandler* handler, void* entry, const char* line);
			static void* readOpenIni(ImGuiContext* ctx, ImGuiSettingsHandler* handler, const char* name);
			static void applySettings(ImGuiContext* ctx, ImGuiSettingsHandler* handler);

		public:
			virtual void init(int argc, char **argv, const char *name);
			virtual void initParameterGUI();
			virtual void initSimulationParameterGUI();
			virtual void render();
			virtual void reset();
			virtual void update();
			virtual void cleanup();
			virtual void run();
			virtual void stop();
			virtual void addKeyFunc(char k, std::function<void()> const& func);

			void createSimulationParameterGUI();
	};
}

#endif