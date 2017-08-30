#include "DragForce_Gissler2017.h"
#include "SPlisHSPlasH/TimeManager.h"

#define _USE_MATH_DEFINES
#include <math.h>

using namespace SPH;

// C.Gissler, S.Band, A.Peer, M.Ihmsen, M.Teschner,
// "Approximate Air-Fluid Interactions for SPH,"
// VRIPHYS 2017

DragForce_Gissler2017::DragForce_Gissler2017(FluidModel *model) :
	DragBase(model)
{
}

DragForce_Gissler2017::~DragForce_Gissler2017(void)
{
}

void DragForce_Gissler2017::step()
{
	const Real supportRadius = m_model->getSupportRadius();
	const Real radius = m_model->getParticleRadius();
	const Real diam = 2.0*radius;
	static const Real pi = static_cast<Real>(M_PI);
	const Real rho_l = m_model->getDensity0();

	// Air velocity.
	const Vector3r va(0, 0, 0);

	const Real L = cbrt(0.75 / pi) * diam;

	const Real inv_td = 0.5*C_d * mu_l / (rho_l * L*L);
	const Real td = 1.0 / inv_td;
	Real omegaSquare = C_k * sigma / (rho_l * L*L*L) - inv_td*inv_td;
	std::max(omegaSquare, 0.0);
	const Real omega = sqrt(omegaSquare);

	// Equation (6)
	Real val = td*td*omegaSquare;
	val = sqrt(val+1.0) + td*omega;
	val = std::max(val, -0.5 * pi);
	val = std::min(val, 0.5 * pi);
	const Real t_max = -2.0 * (atan(val) - pi) / omega;

	// Equation (7)
	const Real c_def = 1.0 - exp(-t_max / td) * (cos(omega * t_max) + 1.0/(omega*td) * sin(omega * t_max));

	// Weber number without velocity
	const Real We_i_wo_v = rho_a * L / sigma;

	// Equation (8)
	const Real y_coeff = (C_F * We_i_wo_v * c_def) / (C_k * C_b);

	const Real n_full = 38;
	const Real n_full_23 = n_full * 2.0/3.0;
	
	const unsigned int numParticles = m_model->numActiveParticles();
	#pragma omp parallel default(shared)
	{
		#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			const Vector3r &vi = m_model->getVelocity(0, i);
			Vector3r v_i_rel = va - vi;
			const Real vi_rel_square = v_i_rel.squaredNorm();
			const Real vi_rel_norm = sqrt(vi_rel_square);
			const Real We_i = We_i_wo_v * vi_rel_square;

			Vector3r v_i_rel_n = v_i_rel;
			if (vi_rel_norm <= 1.0e-6)
				continue;
			// Else.
			v_i_rel_n = v_i_rel_n * (1.0 / vi_rel_norm);
 
			// Equation (8)
			const Real y_i_max = std::min(vi_rel_square * y_coeff, 1.0);

			const Real Re_i = 2.0*std::max((rho_a * vi_rel_norm * L) / mu_a, 0.1);

			Real C_Di_sphere;
			if (Re_i <= 1000.0)
				C_Di_sphere = 24.0 / Re_i * (1.0 + 1.0/6.0 * pow(Re_i, 2.0/3.0));
			else
				C_Di_sphere = 0.424;

			// Equation (9)
			const Real C_Di_Liu = C_Di_sphere * (1.0 + 2.632 * y_i_max);

			unsigned int numNeighbors = m_model->numberOfNeighbors(0, i);
			for (unsigned int pid = 1; pid < m_model->numberOfPointSets(); pid++)
				numNeighbors += m_model->numberOfNeighbors(pid, i);

			// Equation (10)
			Real C_Di;
			const Real factor_n = std::min(n_full_23, (Real) numNeighbors) / n_full_23;
			if (numNeighbors == 0)
				C_Di = C_Di_Liu;
			else
				C_Di = (1.0 - factor_n) * C_Di_Liu + factor_n;

			// Equation (12)
			const Real h1 = (L + C_b*L*y_i_max);
			const Real A_i_droplet = pi * h1*h1;

			// Equation (13)
			const Real A_i_unoccluded = (1.0 - factor_n) * A_i_droplet + factor_n * diam*diam;

			//////////////////////////////////////////////////////////////////////////
			// Fluid
			//////////////////////////////////////////////////////////////////////////
			Real max_v_x = 0.0;
			const Vector3r &xi = m_model->getPosition(0, i);
			for (unsigned int j = 0; j < m_model->numberOfNeighbors(0, i); j++)
			{
				const unsigned int neighborIndex = m_model->getNeighbor(0, i, j);
				const Vector3r &xj = m_model->getPosition(0, neighborIndex);
				Vector3r xixj = xi - xj;
				xixj.normalize();
				const Real x_v = v_i_rel_n.dot(xixj);
				max_v_x = std::max(max_v_x, x_v);
			}
			// Equation (15)
			const Real w_i = std::max(0.0, std::min(1.0, 1.0 - max_v_x));

			// Equation (14)
			const Real A_i = w_i * A_i_unoccluded;

			// Drag force. Additionally dividing by mass to get acceleration.
			Vector3r &ai = m_model->getAcceleration(i);
			ai += m_dragCoefficient * static_cast<Real>(0.5) / m_model->getMass(i) * rho_a * (v_i_rel * vi_rel_norm) * C_Di * A_i;
		}
	}
}


void DragForce_Gissler2017::reset()
{
}

