#pragma once

#include <string>
#include "SPlisHSPlasH/TriangleMesh.h"
#include "Discregrid/All"
#include <memory>


namespace SPH
{
	class SamplingBase
	{
	public: 
		struct Region
		{
		public:
			typedef Real value_type;

			Region() {}
			Region(Real minx, Real miny, Real minz,
				Real maxx, Real maxy, Real maxz)
			{
				m_min = Vector3r(minx, miny, minz);
				m_max = Vector3r(maxx, maxy, maxz);
			}

			Vector3r m_min;
			Vector3r m_max;
		};

	protected:
		Real m_radius;
		Vector3r m_scale;
		Real m_diameter;
		bool m_invert;
		Eigen::Matrix<unsigned int, 3, 1> m_resolutionSDF;
		Region m_region;
		bool m_useRegion;
		std::vector<Vector3r> m_x;
		Vector3r m_bbmin, m_bbmax;
		int m_output_format; // 0: partio, 1: vtk
		std::shared_ptr<Discregrid::CubicLagrangeDiscreteGrid> m_distanceField;
		const Real m_eps = static_cast<Real>(1.0e-5);
		TriangleMesh m_mesh;
		std::string m_outputPath;
		bool m_useCache;

		void generateSDF(SPH::TriangleMesh& mesh);
		void computeBoundingBox(TriangleMesh& mesh);		
		double distance(const Vector3r& x, const Real tolerance);
		double distance(const Vector3r& x, const Real tolerance, Vector3r& normal, Vector3r& nextSurfacePoint);
		void sampleObject(const int mode);

		virtual void generateSamples() = 0;

		// VTK expects big endian
		template<typename T>
		inline static void swapByteOrder(T* v)
		{
			constexpr size_t n = sizeof(T);
			uint8_t* bytes = reinterpret_cast<uint8_t*>(v);
			for (unsigned int c = 0u; c < n / 2; c++)
				std::swap(bytes[c], bytes[n - c - 1]);
		}

	public:
		SamplingBase();
		~SamplingBase();

		virtual void generateSampling(const std::string& inputFile, const std::string& outputFile);

		static void writeParticlesVTK(const std::string& fileName, std::vector<Vector3r>& x);

		Real getRadius() const { return m_radius; }
		void setRadius(const Real radius) { this->m_radius = radius; m_diameter = static_cast<Real>(2.0)*m_radius; }
		const Vector3r& getScale() const {	return m_scale;	}
		void setScale(const Vector3r& scale) { m_scale = scale;	}
		const Eigen::Matrix<unsigned, 3, 1>& getResolutionSdf() const {	return m_resolutionSDF; }
		void setResolutionSdf(const Eigen::Matrix<unsigned, 3, 1>& resolution_sdf) { m_resolutionSDF = resolution_sdf; }
		const Region& get_m_region() const { return m_region; }
		void setRegion(const Region& region) { m_region = region; }
		bool useRegion() const { return m_useRegion; }
		void setUseRegion(const bool use_region) { m_useRegion = use_region; }
		int getOutputFormat() const { return m_output_format; }
		bool getInvert() const { return m_invert; }
		void setInvert(const bool invert) { m_invert = invert; }
		bool getUseCache() const { return m_useCache; }
		void setUseCache(const bool b) { m_useCache = b; }
		std::string& getOutputPath() { return m_outputPath;  }
		void setOutputPath(const std::string& path) { m_outputPath = path; }
	};
}
