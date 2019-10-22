#ifndef __GaussQuadrature_h__
#define __GaussQuadrature_h__

#include <Eigen/Dense>

namespace SPH
{
	class GaussQuadrature
	{
	public:

		using Integrand = std::function<double(Eigen::Vector3d const&)>;
		using Domain = Eigen::AlignedBox3d;

		static double integrate(Integrand integrand, Domain const& domain, unsigned int p);
		static void exportSamples(unsigned int p);
	};
}

#endif

