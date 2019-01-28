#include "SPHKernels.h"

using namespace SPH;

Real CubicKernel::m_radius;
Real CubicKernel::m_k;
Real CubicKernel::m_l;
Real CubicKernel::m_W_zero;

Real Poly6Kernel::m_radius;
Real Poly6Kernel::m_k;
Real Poly6Kernel::m_l;
Real Poly6Kernel::m_m;
Real Poly6Kernel::m_W_zero;

Real SpikyKernel::m_radius;
Real SpikyKernel::m_k;
Real SpikyKernel::m_l;
Real SpikyKernel::m_W_zero;

Real WendlandQuinticC2Kernel::m_radius;
Real WendlandQuinticC2Kernel::m_k;
Real WendlandQuinticC2Kernel::m_l;
Real WendlandQuinticC2Kernel::m_W_zero;

Real CohesionKernel::m_radius;
Real CohesionKernel::m_k;
Real CohesionKernel::m_c;
Real CohesionKernel::m_W_zero;

Real AdhesionKernel::m_radius;
Real AdhesionKernel::m_k;
Real AdhesionKernel::m_W_zero;

Real CubicKernel2D::m_radius;
Real CubicKernel2D::m_k;
Real CubicKernel2D::m_l;
Real CubicKernel2D::m_W_zero;

Real WendlandQuinticC2Kernel2D::m_radius;
Real WendlandQuinticC2Kernel2D::m_k;
Real WendlandQuinticC2Kernel2D::m_l;
Real WendlandQuinticC2Kernel2D::m_W_zero;

