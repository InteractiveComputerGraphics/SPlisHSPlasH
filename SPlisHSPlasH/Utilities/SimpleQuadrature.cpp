#include "SimpleQuadrature.h"

#include <algorithm>
#include <numeric>
#include <iostream>
#include "Utilities/Logger.h"

using namespace SPH;

std::vector<Eigen::Vector3d> SimpleQuadrature::m_samplePoints;
double SimpleQuadrature::m_volume = 0.0;


double SimpleQuadrature::integrate(Integrand integrand)
{
	double res = 0.0;
 	for (unsigned int i = 0; i < m_samplePoints.size(); i++)
	{
		res += m_volume * integrand(m_samplePoints[i].cast<double>());
	}
    return res;
}

void SPH::SimpleQuadrature::determineSamplePointsInSphere(const double radius, unsigned int p)
{
	if (p < 1)
		p = 1;

	m_samplePoints.clear();
	m_samplePoints.reserve(p*p*p);
	const double radius2 = radius*radius;
	const double stepSize = 2.0 * radius / (double)p;
	const double start = -radius + 0.5*stepSize;
	m_volume = stepSize*stepSize*stepSize;

	Eigen::Vector3d pos;
	pos[0] = start;
	for (unsigned int i = 0; i < p; i++)
	{
		pos[1] = start;
		for (unsigned int j = 0; j < p; j++)
		{
			pos[2] = start;
			for (unsigned int k = 0; k < p; k++)
			{
				// test if sample point is in support radius and if it is not the origin
				const double pn = pos.squaredNorm();
				if (pn < radius2)
				{
					m_samplePoints.push_back(pos);
				}
				pos[2] += stepSize;
			}
			pos[1] += stepSize;
		}
		pos[0] += stepSize;
	}
}

void SPH::SimpleQuadrature::determineSamplePointsInCircle(const double radius, unsigned int p)
{
	if (p < 1)
		p = 1;

	m_samplePoints.clear();
	m_samplePoints.reserve(p*p);
	const double radius2 = radius*radius;
	const double stepSize = 2.0 * radius / (double)p;
	const double start = -radius + 0.5*stepSize;
	m_volume = stepSize*stepSize;

	Eigen::Vector3d pos;
	pos[0] = start;
	for (unsigned int i = 0; i < p; i++)
	{
		pos[1] = start;
		for (unsigned int j = 0; j < p; j++)
		{
			pos[2] = 0.0;

			// test if sample point is in support radius and if it is not the origin
			const double pn = pos.squaredNorm();
			if (pn < radius2)
			{
				m_samplePoints.push_back(pos);
			}
			pos[1] += stepSize;
		}
		pos[0] += stepSize;
	}
	LOG_INFO << "Number of sampling points: " << m_samplePoints.size() << "\n";
}
