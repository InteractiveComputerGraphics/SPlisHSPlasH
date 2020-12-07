/*
MIT License

Copyright (c) 2020 Fernando Zorilla and Marcel Ritter

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

/*
Find supplementary material here: 
https://github.com/gileoo/MCSurfaceTension2020-Supplemental
Please get in contact for feedback/support.
*/


#ifndef __SURFACETENSION_ZORILLA2020_H__
#define __SURFACETENSION_ZORILLA2020_H__

#include "SPlisHSPlasH/Common.h"
#include "SPlisHSPlasH/FluidModel.h"
#include "Utilities/StringTools.h"
#include "SurfaceTensionBase.h"
#include <chrono>

namespace SPH
{

	/** This class implements the surface tension method introduced
	*   by Zorilla, Ritter, Sappl, Rauch, Harders:
	*   extended version 2020: https://doi.org/10.3390/computers9020023
	*   and
	*   original version 2019: https://diglib.eg.org/handle/10.2312/cgvc20191260
	*/
	class SurfaceTension_ZorillaRitter2020 : public SurfaceTensionBase
	{
	public:

		SurfaceTension_ZorillaRitter2020(FluidModel* model);
		virtual ~SurfaceTension_ZorillaRitter2020(void);

		/** Linear classifier. Divides into surface or non-surface particle. The function is equivalent
			to the network classifier. Also, inspect lines 344 to 348 in the cpp file
			for how to compute the required input.
		 * \param com       normalized center of mass / number of neighbors
		 * \param non       number of neighbors
		 * \param d_offset  constant parameter d
		 * \return          true if surface, false otherwise
		 **/
		static bool classifySurfaceParticle(double com, int non, double d_offset = 0.0);

		 void performNeighborhoodSearchSort() override;

	private:

		// -- GUI handle variables
		static int SURFACE_TENSION_MEDIAN;

		static int VERSION;
		static int ENUM_VERSION_V2019;
		static int ENUM_VERSION_V2020;

		// V2019 only
		static int CSD;
		static int R2MULT;
		static int TAU;
		static int CLASS_D;
		static int CLASS_D_OFF;
		

		// V2020 only
		static int PCA_NRM_MODE;
		static int PCA_NRM_MIX;
		static int PCA_CUR_MIX;
		static int FIX_SAMPLES;
		static int NEIGH_LIMIT;

		static int SAMPLING;
		static int SAMPLING_HALTON;
		static int SAMPLING_RND;

		static int NORMAL_MODE;
		static int NORMAL_PCA;
		static int NORMAL_MC;
		static int NORMAL_MIX;

		static int SMOOTH_PASSES;
		static int TEMPORAL_SMOOTH;


		// -- Control Parameters

		// switch step update function by version
		enum class StepVersion { V2019, V2020 };
		StepVersion m_step_version;

		// -- Both Versions:
		int  m_Csd;       // number of samples per particle per second
		Real m_tau;       // smoothing factor, default 0.5
		Real m_r2mult;    // r1 to R2 factor, default 0.8
		Real m_r1;        // radius of current particle
		Real m_r2;        // radius of neighbor particles
		Real m_class_k;   // slope of the linear classifier
		Real m_class_d;   // constant of the linear classifier
		bool m_temporal_smoothing;

		// -- containers per particle (V2019)
		//    there are more than minimally required 
		//    for analysis, visualization and debugging
		std::vector<Vector3r> m_mc_normals;          // Monte Carlo surface normals
		std::vector<Vector3r> m_mc_normals_smooth;   // smoothed normals
		std::vector<Real>     m_mc_curv;             // Monte Carlo surface curvature
		std::vector<Real>     m_mc_curv_smooth;      // smoothed curvature
		std::vector<Real>     m_classifier_output;   // outut of the surface classifier

		// -- V2020
		int    m_CsdFix;            // number of samples per computational step
		Real   m_class_d_off;       // offset of classifier d used for PCA neighbors
		Real   m_pca_N_mix;         // mixing factor of PCA normal and MC normal
		Real   m_pca_C_mix;         // mixing factor of PCA curvature and MC curvature
		int    m_neighs_limit;      // maximum nr of neighbors used in PCA computation
		int    m_CS_smooth_passes;  // nr of smoohting passes

		//     Switch between precomputed halton and random sampling mode
		enum class RandomMethod { HALTON, RND, SIZE };  // Halton23, Random 
		RandomMethod m_halton_sampling;
		
		//     Switch between 3 different normal vector variants
		enum class NormalMethod { PCA, MC, MIX, SIZE }; // PCA, Monte Carlo, Mixed
		NormalMethod m_normal_mode;

		// -- containers per particle (V2020)
		//    there are more than minimally required 
		//    for analysis, visualization and debugging
		std::vector<Vector3r> m_pca_normals;       // surface normal by PCA
		std::vector<Real>     m_pca_curv;          // curvature estimate by spherity
		std::vector<Real>     m_pca_curv_smooth;   // smoothed curvature
		std::vector<Real>     m_final_curvatures;
		std::vector<Real>     m_final_curvatures_old;
		std::vector<Real>     m_classifier_input;

		//-- simluation step function
		virtual void step();
		void stepZorilla();  // Version: https://diglib.eg.org/handle/10.2312/cgvc20191260
		void stepRitter();   // Version: https://doi.org/10.3390/computers9020023

		virtual void initParameters();
		virtual void reset();

		/** More configurable linear classifier (non fix k and d) Equivalent to 
		  * classifySurfaceParticle( ... ) above.
		 **/
		bool classifyParticleConfigurable( double com, int non, double d_offset = 0.0 );

		/** Neural network classifier. Divides into surface or non-surface particle.
			Weights are pre-trained.
		 * \param com       normalized center of mass / number of neighbors
		 * \param non       number of neighbors
		 * \param d_offset  constant parameter d
		 * \return          true if surface, false otherwise
		 **/
		bool evaluateNetwork(double com, int non);

		/** Get N random samples on a sphere with radius r (at the origin). */
		std::vector<Vector3r> getSphereSamplesRnd( int N, Real supportRadius );

		/** Get N random samples on a sphere with radius r (at the origin) by using a 
		    precomputed look-up table of vectors (vec3 serialized). */
		std::vector<Vector3r> getSphereSamplesLookUp(
			int N, Real supportRadius, int start, const std::vector<double>& vec3, int mod );

		/** Get N random samples on a sphere with radius r (at the origin) by using a
			precomputed look-up table of vectors (vec3 serialized). */
		std::vector<Vector3r> getSphereSamplesLookUp(
			int N, Real supportRadius, int start, const std::vector<float>& vec3, int mod );



		void resizeV20States( size_t N );
		void resizeV19States( size_t N );
	
		//-- Fuction to setup GUI elements and callback
		template<class T>
		void setupGUIParam(int& PARAMID, std::string name, std::string group, std::string description, T* val)
		{
			std::vector<std::string> tmp;
			Utilities::StringTools::tokenize(name, tmp, " ");

			PARAMID = createNumericParameter<T>( "surfTZR" + tmp[0], name, val);
			setGroup(PARAMID, group);
			setDescription(PARAMID, description);
			GenParam::NumericParameter<T>* rparam = static_cast<GenParam::NumericParameter<T>*>(getParameter(PARAMID));
		}

		void setupGUIEnum(int& PARAMID, std::string name, std::string group, std::string description,
			std::vector<std::pair<std::string, int*>> enum_ids, 
			const GenParam::ParameterBase::GetFunc<int>& getter,
			const GenParam::ParameterBase::SetFunc<int>& setter);


		int getVersionMethod() const
		{
			return static_cast<int>( m_step_version );
		}

		void setVersionMethod( const int a )
		{
			m_step_version = static_cast<StepVersion>(a % 2);
			
			// -- set default values, on method switch	
			if (m_step_version == StepVersion::V2020 &&
				m_normal_mode == NormalMethod::MIX )
				m_class_d = 13;
			else
				m_class_d = 28;			
		}

		int getSamplingMethod() const
		{
			return static_cast<int>(m_halton_sampling);
		}

		void setSamplingMethod(const int a)
		{
			m_halton_sampling = static_cast<RandomMethod>(a % int(RandomMethod::SIZE) );
		}


		int getNormalMethod() const
		{
			return static_cast<int>(m_normal_mode);
		}

		void setNormalMethod(const int a)
		{
			m_normal_mode = static_cast<NormalMethod>(a % int(NormalMethod::SIZE));

			// -- set default values, on method switch	
			if (m_step_version == StepVersion::V2020 &&
				m_normal_mode == NormalMethod::MIX)
				m_class_d = 13;
			else
				m_class_d = 28;
		}

		//-- Accessor functions for field hanlding
		FORCE_INLINE Real& getClassifierOutput( const unsigned int i)
		{
			return m_classifier_output[i];
		}

		FORCE_INLINE Real& getFinalCurvature( const unsigned int i)
		{
			return m_final_curvatures[i];
		}


		// -- helper functions:
		static Vector3r anglesToVec(Real theta, Real phi, Real radius);
	};
}

#endif 