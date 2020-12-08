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

#include <iostream>
#include <fstream>
#include <map>

#include "SPlisHSPlasH/Simulation.h"
#include "SPlisHSPlasH/TimeManager.h"

#include "SPlisHSPlasH/Utilities/MathFunctions.h"

#include "Utilities/StringTools.h"

#include "SurfaceTension_ZorillaRitter2020_haltonVec323.h"
#include "SurfaceTension_ZorillaRitter2020.h"


using namespace SPH;
using namespace std;
using namespace chrono;


// -- GUI variable handles
int SurfaceTension_ZorillaRitter2020::SURFACE_TENSION_MEDIAN = -1;

int SurfaceTension_ZorillaRitter2020::VERSION = -1;
int SurfaceTension_ZorillaRitter2020::ENUM_VERSION_V2019 = -1;
int SurfaceTension_ZorillaRitter2020:: ENUM_VERSION_V2020 = -1;

int SurfaceTension_ZorillaRitter2020::CSD    = -1;
int SurfaceTension_ZorillaRitter2020::R2MULT = -1;
int SurfaceTension_ZorillaRitter2020::TAU = -1;
int SurfaceTension_ZorillaRitter2020::PCA_NRM_MODE = -1;
int SurfaceTension_ZorillaRitter2020::CLASS_D = -1;
int SurfaceTension_ZorillaRitter2020::CLASS_D_OFF = -1;
int SurfaceTension_ZorillaRitter2020::PCA_NRM_MIX = -1;
int SurfaceTension_ZorillaRitter2020::PCA_CUR_MIX = -1;
int SurfaceTension_ZorillaRitter2020::FIX_SAMPLES = -1; 
int SurfaceTension_ZorillaRitter2020::NEIGH_LIMIT = -1;

int SurfaceTension_ZorillaRitter2020::SAMPLING        = -1;
int SurfaceTension_ZorillaRitter2020::SAMPLING_HALTON = -1;
int SurfaceTension_ZorillaRitter2020::SAMPLING_RND    = -1;

int SurfaceTension_ZorillaRitter2020::NORMAL_MODE = -1;
int SurfaceTension_ZorillaRitter2020::NORMAL_PCA  = -1;
int SurfaceTension_ZorillaRitter2020::NORMAL_MC   = -1;
int SurfaceTension_ZorillaRitter2020::NORMAL_MIX  = -1;

int SurfaceTension_ZorillaRitter2020::SMOOTH_PASSES = -1;
int SurfaceTension_ZorillaRitter2020::TEMPORAL_SMOOTH = -1;


void SPH::SurfaceTension_ZorillaRitter2020::setupGUIEnum(int& PARAMID, string name, string group, string description,
	std::vector<pair<string, int*>> enum_ids,
	const GenParam::ParameterBase::GetFunc<int>& getter,
	const GenParam::ParameterBase::SetFunc<int>& setter)
{
	vector<string> tmp;
	Utilities::StringTools::tokenize(name, tmp, " ");

	PARAMID = createEnumParameter("surfTZR" + tmp[0], name, getter, setter);
	setGroup(PARAMID, group);
	setDescription(PARAMID, description);
	GenParam::EnumParameter* enumParam = static_cast<GenParam::EnumParameter*>(getParameter(PARAMID));
	for (auto& s : enum_ids)
		enumParam->addEnumValue(s.first, *s.second);
}

void SurfaceTension_ZorillaRitter2020::resizeV19States( size_t N )
{
	m_mc_normals.resize( N, Vector3r::Zero() );
	m_final_curvatures.resize( N, 0.0 );
}

void SurfaceTension_ZorillaRitter2020::resizeV20States( size_t N )
{
	m_pca_curv.resize( N, 0.0 );
	m_pca_curv_smooth.resize( N, 0.0 );
	m_mc_curv.resize(N, 0.0);
	m_mc_curv_smooth.resize(  N, 0.0 );

	m_mc_normals_smooth.resize( N, Vector3r::Zero() );
	m_pca_normals.resize( N, Vector3r::Zero() );

	m_final_curvatures_old.resize( N, 0.0 );

	m_classifier_input.resize( N, 0.0 );

	m_classifier_output.resize(N, 0.0);
}

SurfaceTension_ZorillaRitter2020::SurfaceTension_ZorillaRitter2020( FluidModel* model ) 
	: SurfaceTensionBase( model )
	, m_step_version(StepVersion::V2020)
	, m_Csd(10000) // 10000 // 36000 // 48000 // 60000
	, m_tau(0.5)
	, m_r2mult(0.8)
	, m_r1(Simulation::getCurrent()->getSupportRadius())
	, m_r2(m_r2mult* m_r1)
	, m_class_k( 74.688796680497925 )
	, m_class_d( 12 )
	, m_temporal_smoothing( false )
	, m_CsdFix( -1 )
	, m_class_d_off( 2 )
	, m_pca_N_mix( 0.75 )
	, m_pca_C_mix( 0.5 )
	, m_neighs_limit( 16 )	
	, m_CS_smooth_passes( 1 )
	, m_halton_sampling(RandomMethod::HALTON)
	, m_normal_mode( NormalMethod::MC )
{
	m_model->addField( { "Curv.Final (SurfZR20)", FieldType::Scalar, [this](const unsigned int i ) -> Real* { return &this->getFinalCurvature( i ); } } );
	m_model->addField( { "Surf.Class (SurfZR20)", FieldType::Scalar, [this](const unsigned int i) -> Real* { return &this->getClassifierOutput( i ); } });

	// -- Both Versions
	resizeV19States( model->numParticles() );

	// -- V2020 only
	resizeV20States( model->numParticles() );

//	initParameters();
}



// -- setup GUI control parameters
void SurfaceTension_ZorillaRitter2020::initParameters()
{
	using enumGet = GenParam::ParameterBase::GetFunc<int>;
	using enumSet = GenParam::ParameterBase::SetFunc<int>;

	SurfaceTensionBase::initParameters();

	const string grp19 = "Surface tension"; // "Surface tension ZR19";
	const string grp20 = "Surface tension"; // "Surface tension ZR20";

	enumGet getVersionFct = std::bind( &SurfaceTension_ZorillaRitter2020::getVersionMethod, this );
	enumSet setVersionFct = std::bind( &SurfaceTension_ZorillaRitter2020::setVersionMethod, this, std::placeholders::_1 );

	setupGUIEnum(VERSION, "version", "Surface tension", "Method or extended method.",
		{ { "V2019", &ENUM_VERSION_V2019 }, { "V2020", &ENUM_VERSION_V2020 } },
		getVersionFct, setVersionFct );

	// 2019 & 2020
	setupGUIParam<int> (CSD,     "Csd (19)",       grp19, "Samples per Second.",      &m_Csd);
	setupGUIParam<Real>(R2MULT,  "r-ratio (19)", grp19, "R1 to R2 ratio.",        &m_r2mult);
	setupGUIParam<Real>(TAU,     "tau (19)",     grp19, "Smoothing factor tau.",  &m_tau);
	setupGUIParam<Real>(CLASS_D, "d (19)",      grp19, "Classifier constant d.", &m_class_d);

	// 2020
	enumGet getSamplingFct = std::bind(&SurfaceTension_ZorillaRitter2020::getSamplingMethod, this);
	enumSet setSamplingFct = std::bind(&SurfaceTension_ZorillaRitter2020::setSamplingMethod, this, std::placeholders::_1);

	setupGUIEnum(SAMPLING, "sampling (20)", grp20, "Monte Carlo samping method.",
		{ { "Halton", &SAMPLING_HALTON }, { "Rnd", &SAMPLING_RND } },
		getSamplingFct, setSamplingFct);

	enumGet getNormalFct = std::bind(&SurfaceTension_ZorillaRitter2020::getNormalMethod, this);
	enumSet setNormalFct = std::bind(&SurfaceTension_ZorillaRitter2020::setNormalMethod, this, std::placeholders::_1);

	setupGUIEnum(NORMAL_MODE, "normal-mode (20)", grp20, "Normal estimation method.",
		{ { "PCA", &NORMAL_PCA}, { "Monte Carlo", &NORMAL_MC}, { "Mix", &NORMAL_MIX} },
		getNormalFct, setNormalFct);

	setupGUIParam<int> (FIX_SAMPLES, "MCSamples (20)", grp20, "Fixed nr of MC samples per step.", &m_CsdFix);
	setupGUIParam<Real>(PCA_CUR_MIX, "PcaMixCur (20)", grp20, "Factor to mix-in pca curvature.", &m_pca_C_mix);
	setupGUIParam<Real>(PCA_NRM_MIX, "PcaMixNrm (20)", grp20, "Factor to mix-in pca normal.", &m_pca_N_mix);

	TEMPORAL_SMOOTH = createBoolParameter("surfTZRtemporalSmooth" , "temporalSmoothing (20)", &m_temporal_smoothing);
	setGroup(TEMPORAL_SMOOTH, grp20);
	setDescription(TEMPORAL_SMOOTH, "Enable temporal smoothing");
}



SurfaceTension_ZorillaRitter2020::~SurfaceTension_ZorillaRitter2020( void )
{
	Simulation* sim = Simulation::getCurrent();

	m_model->removeFieldByName( "Curv.Final (SurfZR20)" );
	m_model->removeFieldByName( "Surf.Class (SurfZR20)" );

}



bool SurfaceTension_ZorillaRitter2020::evaluateNetwork( double com, int non )
{
	double W111 = 0.94356692; double W112 = -38.35266876;
	double W121 = -0.95363748; double W122 = 0.50992483;

	double B11 = 0; double B12 = -4.76320028;

	double W211 = -1.10649812; double W212 = -0.43858105;
	double W221 = 0.92051727; double W222 = -0.99366635;

	double B21 = -8.12260151; double B22 = 0.95201129;


	double L11in = com * W111 + non * W121 + B11;
	double L12in = com * W112 + non * W122 + B12;

	double L11out = std::max( 0.0, L11in );
	double L12out = std::max( 0.0, L12in );

	double L21in = L11out * W211 + L12out * W221 + B21;
	double L22in = L11out * W212 + L12out * W222 + B22;

	double L21out = 1 / (1 + exp( -L21in ));
	double L22out = 1 / (1 + exp( -L22in ));

	if (L21out < L22out) { return true; }
	else { return false; }

}


bool SurfaceTension_ZorillaRitter2020::classifySurfaceParticle( double com, int non, double d_offset)
{
	double neighborsOnTheLine = 74.688796680497925 * com + 28.0; // pre-multiplied

	if (non <= neighborsOnTheLine)
		return true;
	else
		return false;
}


bool SurfaceTension_ZorillaRitter2020::classifyParticleConfigurable( double com, int non, double d_offset )
{
	double neighborsOnTheLine = m_class_k * com + m_class_d + d_offset; // pre-multiplied

	if (non <= neighborsOnTheLine) 
		return true;
	else 
		return false;
}

void SurfaceTension_ZorillaRitter2020::performNeighborhoodSearchSort()
{
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();

	const unsigned int fluidModelIndex = m_model->getPointSetIndex();

	auto const& d = sim->getNeighborhoodSearch()->point_set(fluidModelIndex);

	// Both version state fields
	d.sort_field(&m_mc_normals[0]);
	d.sort_field(&m_final_curvatures[0]);

	// Both version 2020 state fields
	if (m_step_version == StepVersion::V2020)
	{
		d.sort_field(&m_pca_curv[0]);
		d.sort_field(&m_pca_curv_smooth[0]);
		d.sort_field(&m_mc_curv[0]);
		d.sort_field(&m_mc_curv_smooth[0]);

		d.sort_field(&m_mc_normals_smooth[0]);
		d.sort_field(&m_pca_normals[0]);

		d.sort_field(&m_final_curvatures_old[0]);

		d.sort_field(&m_classifier_input[0]);

		d.sort_field(&m_classifier_output[0]);
	}
	
}

void SurfaceTension_ZorillaRitter2020::step()
{
	// -- switch version
	if( m_step_version == StepVersion::V2019 )
		stepZorilla();
	else
		stepRitter();
}



Vector3r SurfaceTension_ZorillaRitter2020::anglesToVec( Real theta, Real phi, Real radius )
{
	const Real r_sin_phi = radius * sin( phi );

	return Vector3r(
		r_sin_phi * cos( theta ),
		r_sin_phi * sin( theta ),
		radius * cos( phi ) );
}

std::vector<Vector3r> SurfaceTension_ZorillaRitter2020::getSphereSamplesRnd( int N, Real supportRadius )
{
	std::vector<Vector3r> points(N);

	for (int i = 0; i < N; i++)
	{
		Real theta = static_cast<Real>(2.0 * M_PI) * static_cast<Real>(rand()) / static_cast<Real>( RAND_MAX );
		Real phi = acos(static_cast<Real>(1.0) - static_cast<Real>(2.0) * static_cast<Real>(rand()) / static_cast<Real>( RAND_MAX ) );

		points[i] = anglesToVec( theta, phi, supportRadius );
	}

	return points;
}

/** Step function used in the first version of the paper 
 *	https://diglib.eg.org/handle/10.2312/cgvc20191260
*/
void SurfaceTension_ZorillaRitter2020::stepZorilla()
{
	Simulation* sim = Simulation::getCurrent();

	const Real supportRadius = sim->getSupportRadius();
	const unsigned int numParticles = m_model->numActiveParticles();
	const Real k = m_surfaceTension;

	Real timeStep = SPH::TimeManager::getCurrent()->getTimeStepSize();
	Real r1 = supportRadius;
	Real r2 = m_r2mult * r1;        // m_r2mult=0.8 best results empirically
	const unsigned int NumberOfPoints = int( m_Csd * timeStep );

	const unsigned int fluidModelIndex = m_model->getPointSetIndex();

	FluidModel* model = m_model;

#pragma omp parallel default(shared)
	{
#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			m_mc_normals[i] = Vector3r::Zero();
			m_classifier_output[i] = 0;

			Vector3r centerofMasses = Vector3r::Zero();
			int numberOfNeighbours = sim->numberOfNeighbors(fluidModelIndex, fluidModelIndex, i );


			if (sim->numberOfNeighbors(fluidModelIndex, fluidModelIndex, i) == 0)
				continue;

			const Vector3r& xi = m_model->getPosition( i );

			forall_fluid_neighbors_in_same_phase(
				Vector3r xjxi = (xj - xi);
				centerofMasses += xjxi / supportRadius; )

			if (classifyParticleConfigurable( centerofMasses.norm() / double( numberOfNeighbours ), numberOfNeighbours )) //EvaluateNetwork also possible
			{
				std::vector<Vector3r> points = getSphereSamplesRnd( NumberOfPoints, supportRadius );

				forall_fluid_neighbors_in_same_phase(
					Vector3r xjxi = (xj - xi);
					for (int p = static_cast<int>(points.size()) - 1; p >= 0; --p)
					{
						Vector3r vec = (points[p] - xjxi);
						Real dist = vec.squaredNorm();
						if (dist <= pow((r2 / r1), 2) * supportRadius * supportRadius)
						{
							points.erase(points.begin() + p);
						}
					} )

				for (int p = static_cast<int>(points.size()) - 1; p >= 0; --p)
				{
					m_mc_normals[i] += points[p];
				}

				if (points.size() > 0)
				{
					m_mc_normals[i].normalize();

					m_mc_curv[i] = (static_cast<Real>(1.0) / supportRadius) * static_cast<Real>(-2.0) * pow( (static_cast<Real>(1.0) - (r2 * r2 / (r1 * r1))), static_cast<Real>(-0.5) ) *
						cos( acos(static_cast<Real>(1.0 - 2.0 * (points.size() / static_cast<Real>( NumberOfPoints ))) ) + asin( r2 / r1 ) );

					m_classifier_output[i] = 1;
				}
				else
				{
					m_mc_curv[i] = 0.0;
					m_mc_normals[i] = Vector3r::Zero();
					m_classifier_output[i] = 0.5;
				}
			}
		}
	}


#pragma omp parallel default(shared)
	{
#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			if (m_mc_normals[i] != Vector3r::Zero())
			{
				const Vector3r& xi = m_model->getPosition(i);
				Vector3r normalCorrection = Vector3r::Zero();
				Vector3r& ai = m_model->getAcceleration(i);
				Real curvatureCorrection = 0;
				Real correctionFactor = 0.0;

				forall_fluid_neighbors_in_same_phase(
					if (m_mc_normals[neighborIndex] != Vector3r::Zero())
					{
						Vector3r xjxi = (xj - xi);
						Real distanceji = xjxi.norm();

						normalCorrection += m_mc_normals[neighborIndex] * (1 - distanceji / supportRadius);
						curvatureCorrection += m_mc_curv[neighborIndex] * (1 - distanceji / supportRadius);
						correctionFactor += (1 - distanceji / supportRadius);
					})

					normalCorrection.normalize();
					normalCorrection = (1 - m_tau) * m_mc_normals[i] + m_tau * normalCorrection;
					normalCorrection.normalize();


					m_final_curvatures[i] =
						((static_cast<Real>(1.0) - m_tau) * m_mc_curv[i] + m_tau * curvatureCorrection)
						/ ((static_cast<Real>(1.0) - m_tau) + m_tau * correctionFactor);

					ai -= (static_cast<Real>(1.0) / m_model->getMass(i)) * normalCorrection * k * m_final_curvatures[i];

			}
			else
				m_final_curvatures[i] = 0.0;
		}

	}

}


// -- Helper functions for extended step function
vector<Vector3r> SurfaceTension_ZorillaRitter2020::getSphereSamplesLookUp(
	int N, Real supportRadius, int start, 
	const vector<float>& vec3, int mod )
{
	vector<Vector3r> points( N );
	int s = (start / 3) * 3; // ensure to be dividable by 3
	for (int i = 0; i < N; i++)
	{
		int i3 = s + 3 * i;
		points[i] = supportRadius * Vector3r( vec3[i3 % mod], vec3[(i3 + 1) % mod], vec3[(i3 + 2) % mod] );
	}
	return points;
}


/** Step update function used in the extended version of the paper
 *	https://doi.org/10.3390/computers9020023
 */
void SurfaceTension_ZorillaRitter2020::stepRitter()
{
	auto& tm = *SPH::TimeManager::getCurrent();
	Real timeStep = tm.getTimeStepSize();
	size_t step = static_cast<size_t>(floor( tm.getTime() / timeStep ));

	m_r2 = m_r1 * m_r2mult;

	Simulation* sim = Simulation::getCurrent();

	const Real supportRadius = sim->getSupportRadius();
	const unsigned int numParticles = m_model->numActiveParticles();
	const Real k = m_surfaceTension;

	unsigned int NrOfSamples;

	const unsigned int fluidModelIndex = m_model->getPointSetIndex();
	FluidModel* model = m_model;

	if (m_CsdFix > 0)
		NrOfSamples = m_CsdFix;
	else
		NrOfSamples = int( m_Csd * timeStep );

	// ################################################################################################
	// ## first pass, compute classification and first estimation for normal and curvature (Montecarlo)
	// ################################################################################################

#pragma omp parallel default(shared)
	{
#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			// init or reset arrays
			m_mc_normals[i] = Vector3r::Zero();

			m_pca_normals[i] = Vector3r::Zero();
			m_mc_normals_smooth[i] = Vector3r::Zero();

			m_mc_curv[i] = 0.0;
			m_mc_curv_smooth[i] = 0.0;
			m_pca_curv[i] = 0.0;
			m_pca_curv_smooth[i] = 0.0;
			m_final_curvatures[i] = 0.0;


			// -- compute center of mass of current particle

			Vector3r centerofMasses = Vector3r::Zero();
			int numberOfNeighbours = sim->numberOfNeighbors(fluidModelIndex, fluidModelIndex, i );

			if (numberOfNeighbours == 0)
			{
				m_mc_curv[i] = static_cast<Real>(1.0) / supportRadius;
				continue;
			}

			const Vector3r& xi = m_model->getPosition( i );

			forall_fluid_neighbors_in_same_phase(
				Vector3r xjxi = (xj - xi);
				centerofMasses += xjxi; )
			
			centerofMasses /= supportRadius;

			// cache classifier input, could also be recomputed later to avoid caching
			m_classifier_input[i] = centerofMasses.norm() / static_cast<Real>( numberOfNeighbours );
			

			// -- if it is a surface classified particle
			if (classifyParticleConfigurable( m_classifier_input[i], numberOfNeighbours )) //EvaluateNetwork also possible
			{

				// -- create monte carlo samples on particle
				std::vector<Vector3r> points;

				if( m_halton_sampling == RandomMethod::HALTON )
					points = getSphereSamplesLookUp(
						NrOfSamples, supportRadius, i * NrOfSamples, haltonVec323, static_cast<int>(haltonVec323.size()) ); // 8.5 // 15.0(double) // 9.0(float)
				else
					points = getSphereSamplesRnd( NrOfSamples, supportRadius );


				//  -- remove samples covered by neighbor spheres
				forall_fluid_neighbors_in_same_phase(
					Vector3r xjxi = (xj - xi);
					for (int p = static_cast<int>(points.size()) - 1; p >= 0; --p)
					{
						Vector3r vec = (points[p] - xjxi);
						Real dist = vec.squaredNorm();

						if (dist <= pow( (m_r2 / m_r1), 2 ) * supportRadius * supportRadius)
							points.erase( points.begin() + p );
					} )

				// -- estimate normal by left over sample directions
				for (int p = static_cast<int>(points.size()) - 1; p >= 0; --p)
					m_mc_normals[i] += points[p];


				// -- if surface classified and non-overlapping neighborhood spheres
				if (points.size() > 0)
				{
					m_mc_normals[i].normalize();

					// -- estimate curvature by sample ratio and particle radii
					m_mc_curv[i] = (static_cast<Real>(1.0) / supportRadius) * static_cast<Real>(-2.0) * pow( (static_cast<Real>(1.0) - (m_r2 * m_r2 / (m_r1 * m_r1))), static_cast<Real>(-0.5) ) *
						cos( acos(static_cast<Real>(1.0) - static_cast<Real>(2.0) * (static_cast<Real>( points.size() ) / static_cast<Real>( NrOfSamples )) ) + asin( m_r2 / m_r1 ) );

					m_classifier_output[i] = 1.0; // -- used to visualize surface points (blue in the paper)
				}
				else
				{
					// -- correct false positives to inner points
					m_mc_normals[i] = Vector3r::Zero();
					m_mc_curv[i] = 0.0;
					m_classifier_output[i] = 0.5; // -- used for visualize post-correction points (white in the paper)
				}
			}
			else
			{
				// -- used to visualize inner points (green in the paper)
				m_classifier_output[i] = 0.0;
			}

		}
	}

	// ################################################################################################
	// ## second pass, compute normals and curvature and compute PCA normal 
	// ################################################################################################

#pragma omp parallel default(shared)
	{
#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			if (m_mc_normals[i] != Vector3r::Zero())
			{
				const Vector3r& xi = m_model->getPosition( i );
				Vector3r normalCorrection = Vector3r::Zero();
				Vector3r& ai = m_model->getAcceleration( i );


				Real correctionForCurvature = 0;
				Real correctionFactor = 0.0;

				Vector3r centroid = xi;
				Vector3r surfCentDir = Vector3r::Zero();

				// collect neighbors
				multimap<Real, size_t> neighs;

				Eigen::Matrix<Real, 3, 3> t = Eigen::Matrix<Real, 3, 3>::Zero();
				int t_count = 0;
				Vector3r neighCent = Vector3r::Zero();

				int nrNeighhbors = sim->numberOfNeighbors(fluidModelIndex, fluidModelIndex, i );

				forall_fluid_neighbors_in_same_phase(
					if (m_mc_normals[neighborIndex] != Vector3r::Zero())
					{
						Vector3r& xj = m_model->getPosition( neighborIndex );
						Vector3r xjxi = (xj - xi);

						surfCentDir += xjxi;
						centroid += xj;
						t_count++;

						Real distanceji = xjxi.norm();

						// collect neighbors sorted by distance relative to xi as estimate
						if ( m_normal_mode != NormalMethod::MC )
						{
							Real distanceji = xjxi.norm();
							neighs.insert( { distanceji, neighborIndex } );
						}

						normalCorrection += m_mc_normals[neighborIndex] * (1 - distanceji / supportRadius);
						correctionForCurvature += m_mc_curv[neighborIndex] * (1 - distanceji / supportRadius);
						correctionFactor += (1 - distanceji / supportRadius);
					}

					else if( m_normal_mode != NormalMethod::MC
						&& classifyParticleConfigurable( m_classifier_input[neighborIndex], nrNeighhbors, m_class_d_off ))
					{
						Vector3r& xj = m_model->getPosition( neighborIndex );
						Vector3r xjxi = (xj - xi);
						surfCentDir += xjxi;
						centroid += xj;
						t_count++;
						Real distanceji = xjxi.norm();
						neighs.insert( { distanceji, neighborIndex } );
					} )

				bool do_pca = false;

				if (m_normal_mode == NormalMethod::PCA)
				{
					do_pca = true;
				}
				else if (m_normal_mode == NormalMethod::MIX)
				{
					if (m_pca_N_mix >= 0.0 && m_pca_C_mix >= 0.0)
						do_pca = true;
				}

				if (do_pca)
				{
					centroid /= static_cast<Real>(t_count + 1);
					surfCentDir.normalize();

					const size_t neigh_limit = m_neighs_limit;
					size_t neigh_count = 0;

					for (pair<Real, size_t> nn : neighs)
					{
						size_t j = nn.second;
						Vector3r& xj = m_model->getPosition( j );
						Vector3r xjxi = (xj - centroid);

						Real distanceji = xjxi.norm();

						Real w = CubicKernel::W( distanceji ) / CubicKernel::W_zero();
						neighCent += w * xj;

						Eigen::Matrix<Real, 3, 3> dy = xjxi * xjxi.transpose();

						t = t + w * dy;
						t_count++;

						if (neigh_count >= neigh_limit - 1)
							break;

						neigh_count++;
					}

					if (neighs.size() >= 4)
					{
						Vector3r pdt_ev;
						Matrix3r pdt_ec;
						MathFunctions::eigenDecomposition(t, pdt_ec, pdt_ev);

						// sort values smallest to greatest
						if (pdt_ev(0) > pdt_ev(1))
						{
							swap(pdt_ev(0), pdt_ev(1));
							swap(pdt_ec(0, 0), pdt_ec(0, 1));
							swap(pdt_ec(1, 0), pdt_ec(1, 1));
							swap(pdt_ec(2, 0), pdt_ec(2, 1));

						}
						if (pdt_ev(0) > pdt_ev(2))
						{
							swap(pdt_ev(0), pdt_ev(2));
							swap(pdt_ec(0, 0), pdt_ec(0, 2));
							swap(pdt_ec(1, 0), pdt_ec(1, 2));
							swap(pdt_ec(2, 0), pdt_ec(2, 2));
						}
						if (pdt_ev(1) > pdt_ev(2))
						{
							swap(pdt_ev(1), pdt_ev(2));
							swap(pdt_ec(0, 1), pdt_ec(0, 2));
							swap(pdt_ec(1, 1), pdt_ec(1, 2));
							swap(pdt_ec(2, 1), pdt_ec(2, 2));
						}

						Vector3r minor = Vector3r(pdt_ec(0, 0), pdt_ec(1, 0), pdt_ec(2, 0));
						
						if (minor.dot(m_mc_normals[i]) < 0.0)
							m_pca_normals[i] = -1.0 * minor;
						else
							m_pca_normals[i] = minor;

						m_pca_curv[i] = static_cast<Real>(3.0) * pdt_ev(0) / (pdt_ev(0) + pdt_ev(1) + pdt_ev(2) );

						if (surfCentDir.dot(m_pca_normals[i]) > 0.0)
							m_pca_curv[i] *= -1;

						m_pca_normals[i].normalize();
					}
					else
					{
						m_pca_normals[i] = Vector3r::Zero();
						m_pca_curv[i] = 0.0;
					}
				}

				normalCorrection.normalize();
				m_mc_normals_smooth[i] = (1 - m_tau) * m_mc_normals[i] + m_tau * normalCorrection;
				m_mc_normals_smooth[i].normalize();

				m_mc_curv_smooth[i] =
					((static_cast<Real>(1.0) - m_tau) * m_mc_curv[i] + m_tau * correctionForCurvature) /
					((static_cast<Real>(1.0) - m_tau) + m_tau * correctionFactor);
			}
		}
	}


	// ################################################################################################
	// ## third pass, final blending and temporal smoothing
	// ################################################################################################
	
	m_CS_smooth_passes = std::max( 1, m_CS_smooth_passes );

	for (int si = 0; si < m_CS_smooth_passes; si++)
	{
		// smoothing pass 2 for sphericity
#pragma omp parallel default(shared)
		{
#pragma omp for schedule(static)  
			for (int i = 0; i < (int)numParticles; i++)
			{
				if (m_mc_normals[i] != Vector3r::Zero())
				{
					int count = 0;
					Real CsCorr = 0.0;

					const Vector3r& xi = m_model->getPosition( i );

					forall_fluid_neighbors_in_same_phase( 
						if (m_mc_normals[neighborIndex] != Vector3r::Zero())
						{
							CsCorr += m_pca_curv[neighborIndex];
							count++;
						} )
					

					if (count > 0)
						m_pca_curv_smooth[i] = static_cast<Real>(0.25) * m_pca_curv_smooth[i] + static_cast<Real>(0.75) * CsCorr / static_cast<Real>(count);
					else
						m_pca_curv_smooth[i] = m_pca_curv[i];

					m_pca_curv_smooth[i] /= supportRadius;

					m_pca_curv_smooth[i] *= 20.0;

					if (m_pca_curv_smooth[i] > 0.0)
						m_pca_curv_smooth[i] = std::min( 0.5f / supportRadius, m_pca_curv_smooth[i] );
					else
						m_pca_curv_smooth[i] = std::max( -0.5f / supportRadius, m_pca_curv_smooth[i] );


					Vector3r final_normal = Vector3r::Zero();
					Real     final_curvature = m_mc_curv_smooth[i];

					if (m_normal_mode == NormalMethod::MC)
					{
						final_normal = m_mc_normals_smooth[i];
						final_curvature = m_mc_curv_smooth[i];

					}
					else if (m_normal_mode == NormalMethod::PCA)
					{
						final_normal = m_pca_normals[i];
						final_curvature = m_pca_curv_smooth[i];
					}
					else if (m_normal_mode == NormalMethod::MIX) // <------------------------------
					{
						bool no_N = (m_pca_normals[i] == Vector3r::Zero());
						bool no_C = no_N || (m_pca_curv_smooth[i] != m_pca_curv_smooth[i]);

						if (no_N)
							final_normal = m_mc_normals_smooth[i];
						else
							final_normal = m_pca_N_mix * m_pca_normals[i] + (1.0 - m_pca_N_mix) * m_mc_normals_smooth[i];

						final_normal.normalize();


						if (no_C)
							m_pca_curv_smooth[i] = m_mc_curv_smooth[i];
						else
						{
							if (m_mc_curv_smooth[i] < 0.0)
								final_curvature = m_mc_curv_smooth[i];
							else
								final_curvature =
									m_pca_C_mix * m_pca_curv_smooth[i] +
									(static_cast<Real>(1.0) - m_pca_C_mix) * m_mc_curv_smooth[i];
						}
					}

					if (m_temporal_smoothing)
						m_final_curvatures[i] = static_cast<Real>(0.05) * final_curvature + static_cast<Real>(0.95) * m_final_curvatures_old[i];
					else
						m_final_curvatures[i] = final_curvature;

					Vector3r force = final_normal * k * m_final_curvatures[i];

					Vector3r& ai = m_model->getAcceleration( i );
					ai -= (1 / m_model->getMass( i )) * force;

					m_final_curvatures_old[i] = m_final_curvatures[i];
				}
				else // non surface particle blend 0.0 curvature
				{
					
					if (m_temporal_smoothing)
						m_final_curvatures[i] = static_cast<Real>(0.95) * m_final_curvatures_old[i];
					else
						m_final_curvatures[i] = 0.0;

					m_final_curvatures_old[i] = m_final_curvatures[i];
				}


			}
		}
	}

}


void SurfaceTension_ZorillaRitter2020::reset()
{
}

