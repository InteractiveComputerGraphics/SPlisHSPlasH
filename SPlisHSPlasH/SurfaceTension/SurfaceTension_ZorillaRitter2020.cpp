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

#include "SplisHSPlash/Simulation.h"
#include "SplisHSPlash/TimeManager.h"
//#include "SplisHSPlash/Viscosity/Viscosity_Standard.h"
//#include "Utilities/Timing.h"

#include "SurfaceTension_ZorillaRitter2020_haltonVec323.h"
#include "SurfaceTension_ZorillaRitter2020_EigenDecomposition.h"
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


vector<string> SurfaceTension_ZorillaRitter2020::split(string txt, char delimiter)
{
	string tmp;
	stringstream ss(txt);
	vector<string> tokens;

	while (getline(ss, tmp, delimiter))
		tokens.push_back(tmp);

	return tokens;
}


void SurfaceTension_ZorillaRitter2020::resizeV19States( size_t N )
{
	m_mc_normals.resize( N, Vector3r::Zero() );
	m_final_curvatures.resize( N, 0.0 );
}

void SPH::SurfaceTension_ZorillaRitter2020::setupGUIEnum(int& PARAMID, std::string name, std::string group, std::string description, 
	std::vector<pair<string, int*>> enum_ids, 
	const GenParam::ParameterBase::GetFunc<int>& getter,
	const GenParam::ParameterBase::SetFunc<int>& setter )
{
	std::string tmp = split(name, ' ')[0];
	PARAMID = createEnumParameter("surfTZR" + tmp, name, getter, setter);
	setGroup(PARAMID, group);
	setDescription(PARAMID, description);
	GenParam::EnumParameter* enumParam = static_cast<GenParam::EnumParameter*>(getParameter(PARAMID));
	for (auto& s : enum_ids)
		enumParam->addEnumValue(s.first, *s.second);
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

#ifdef RICH_OUTPUT
	m_classifier_input2.resize( N, 0.0 );

	m_nr_all_samples.resize(N, 0);
	m_nr_area_samples.resize(N, 0);
#endif
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
	, m_stateChange( -1.0 )
	, m_simulationStop( 5.0 )
	, m_CsdPostState( 200000 )
	, m_CsdFixPostState( 240 )
	, m_sTPostState( 0.02 )
	, m_start_export( -1.0 )
{
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();
	for( unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++ )
	{
		FluidModel* model = sim->getFluidModel( fluidModelIndex );

		model->addField( { "Curv.Final (SurfZR20)", FieldType::Scalar, [this, fluidModelIndex]( const unsigned int i ) -> Real* { return &this->getFinalCurvature( fluidModelIndex, i ); } } );
		model->addField( { "Surf.Class (SurfZR20)", FieldType::Scalar, [this, fluidModelIndex](const unsigned int i) -> Real* { return &this->getClassifierOutput(fluidModelIndex, i); } });

#ifdef RICH_OUTPUT
		model->addField({ "Curv.MC (SurfZR20)", FieldType::Scalar, [this, fluidModelIndex](const unsigned int i) -> Real* { return &this->getMcCurvature(fluidModelIndex, i); } });
		model->addField({ "Curv.MCSmooth (SurfZR20)", FieldType::Scalar, [this, fluidModelIndex](const unsigned int i) -> Real* { return &this->getMcCurvatureSmooth(fluidModelIndex, i); } });
		model->addField({ "Curv.PCA (SurfZR20)", FieldType::Scalar, [this, fluidModelIndex](const unsigned int i) -> Real* { return &this->getCs(fluidModelIndex, i); } });
		model->addField({ "Curv.PCASmooth (SurfZR20)", FieldType::Scalar, [this, fluidModelIndex](const unsigned int i) -> Real* { return &this->getCsCorr(fluidModelIndex, i); } });
		model->addField({ "Surf.Classif.In (SurfZR20)", FieldType::Scalar, [this, fluidModelIndex]( const unsigned int i ) -> Real* { return &this->getClassifierInput( fluidModelIndex, i ); } } );
#endif
	}


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

#ifdef MORE_CONTROL_VARS
	setupGUIParam<int>(NEIGH_LIMIT,  "NeighLimit (20)", grp20, "Limit nr of neighbors.", &m_neighs_limit);
	setupGUIParam<int>(SMOOTH_PASSES,"smoothPasses (20)", grp20, "Nr of smoothing passes.", &m_CS_smooth_passes);
	setupGUIParam<Real>(CLASS_D_OFF, "d-offset (20)",  grp20, "d-offset for pca.", &m_class_d_off);
#endif

	TEMPORAL_SMOOTH = createBoolParameter("surfTZRtemporalSmooth" , "temporalSmoothing (20)", &m_temporal_smoothing);
	setGroup(TEMPORAL_SMOOTH, grp20);
	setDescription(TEMPORAL_SMOOTH, "Enable temporal smoothing");
}



SurfaceTension_ZorillaRitter2020::~SurfaceTension_ZorillaRitter2020( void )
{
	Simulation* sim = Simulation::getCurrent();
	const unsigned int nModels = sim->numberOfFluidModels();
	for( unsigned int fluidModelIndex = 0; fluidModelIndex < nModels; fluidModelIndex++ )
	{
		FluidModel* model = sim->getFluidModel( fluidModelIndex );
		model->removeFieldByName( "SurfZor20 curvature smoothed" );
		model->removeFieldByName( "SurfZor20 curvature corr" );
		model->removeFieldByName( "SurfZor20 curvature" );
		model->removeFieldByName( "SurfZor20 Cs" );
		model->removeFieldByName( "SurfZor20 CsCorr" );
		model->removeFieldByName( "SurfZor20 CsFinal" );
	}
}


#ifdef _WIN32
	#undef max
#endif

/** Neural network classifier. Divides into surface or non-surface particle.
    Weights are pre-trained.
 * \param com       normalized center of mass / number of neighbors
 * \param non       number of neighbors
 * \param d_offset  constant parameter d
 * \return          true if surface, false otherwise
 **/
bool EvaluateNetwork( double com, int non )
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


bool SurfaceTension_ZorillaRitter2020::classifyParticle( double com, int non, double d_offset )
{
	//	double neighborsOnTheLine = ((360.0 / 241.0) * (com * 500.0) + 120) / 10;
	//	double neighborsOnTheLine = 74.688796680497925 * com + 12; // pre-multiplied

	double neighborsOnTheLine = m_class_k * com + m_class_d + d_offset; // pre-multiplied

	if (non <= neighborsOnTheLine) 
		return true;
	else 
		return false;
}

void SurfaceTension_ZorillaRitter2020::step()
{
	// -- switch version
	if( m_step_version == StepVersion::V2019 )
		stepZorilla();
	else
		stepRitter();
}



Vector3r anglesToVec( double theta, double phi, double radius )
{
	const double r_sin_phi = radius * sin( phi );

	return Vector3r(
		r_sin_phi * cos( theta ),
		r_sin_phi * sin( theta ),
		radius * cos( phi ) );
}

std::vector<Vector3r> SurfaceTension_ZorillaRitter2020::getSphereSamplesRnd( int N, Real supportRadius )
{
	std::vector<Vector3r> points;

	for (int i = 0; i < N; i++)
	{
		double theta = 2.0 * M_PI * rand() / double( RAND_MAX );
		double phi = acos( 1.0 - 2.0 * rand() / double( RAND_MAX ) );

		points.push_back( anglesToVec( theta, phi, supportRadius ) );
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
	double r1 = supportRadius;
	double r2 = m_r2mult * r1;        // m_r2mult=0.8 best results empirically
	const unsigned int NumberOfPoints = int( m_Csd * timeStep ); // was m_Csd=100000

	const int fluidModelIndex = 0; // hardcode for now

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


			if (sim->numberOfNeighbors(fluidModelIndex, fluidModelIndex, i ) == 0) { goto nextParticle; }
			const Vector3r& xi = m_model->getPosition( i );

			forall_fluid_neighbors_in_same_phase(
				Vector3r xjxi = (xj - xi);
				centerofMasses += xjxi / supportRadius; )

			if (classifyParticle( centerofMasses.norm() / double( numberOfNeighbours ), numberOfNeighbours )) //EvaluateNetwork also possible
			{
				std::vector<Vector3r> points = getSphereSamplesRnd( NumberOfPoints, supportRadius );

				forall_fluid_neighbors_in_same_phase(
					Vector3r xjxi = (xj - xi);
					for (int p = points.size() - 1; p >= 0; --p)
					{
						Vector3r vec = (points[p] - xjxi);
						Real dist = vec.squaredNorm();
						if (dist <= pow((r2 / r1), 2) * supportRadius * supportRadius)
						{
							points.erase(points.begin() + p);
						}
					} )

				for (int p = points.size() - 1; p >= 0; --p)
				{
					m_mc_normals[i] += points[p];
				}

				if (points.size() > 0)
				{
					m_mc_normals[i].normalize();

					m_mc_curv[i] = (1 / supportRadius) * (-2.0) * pow( (1 - (r2 * r2 / (r1 * r1))), -0.5 ) *
						cos( acos( 1 - 2 * (points.size() / double( NumberOfPoints )) ) + asin( r2 / r1 ) );

					m_classifier_output[i] = 1;
				}
				else
				{
					m_mc_curv[i] = 0.0;
					m_mc_normals[i] = Vector3r::Zero();
					m_classifier_output[i] = 0.5;
				}
			}

		nextParticle:;

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
				double curvatureCorrection = 0;
				double correctionFactor = 0.0;

				forall_fluid_neighbors_in_same_phase(
					if (m_mc_normals[neighborIndex] != Vector3r::Zero())
					{
						Vector3r xjxi = (xj - xi);
						double distanceji = xjxi.norm();

						normalCorrection += m_mc_normals[neighborIndex] * (1 - distanceji / supportRadius);
						curvatureCorrection += m_mc_curv[neighborIndex] * (1 - distanceji / supportRadius);
						correctionFactor += (1 - distanceji / supportRadius);
					})

					normalCorrection.normalize();
					normalCorrection = (1 - m_tau) * m_mc_normals[i] + m_tau * normalCorrection;
					normalCorrection.normalize();


					m_final_curvatures[i] =
						((1 - m_tau) * m_mc_curv[i] + m_tau * curvatureCorrection)
						/ ((1 - m_tau) * 1 + m_tau * correctionFactor);

					ai -= (1 / m_model->getMass(i)) * normalCorrection * k * m_final_curvatures[i];

			}
			else
				m_final_curvatures[i] = 0.0;
		}

	}

}


// -- Helper functions for extended step function

/**
 * Dyadic or tensor product of two vectors.
 * \return a symmetric second order tensor.
*/
Eigen::Matrix<double, 3, 3> dyadic(
	const Vector3r& a,
	const Vector3r& b )
{
	Eigen::Matrix<double, 3, 3> tmp;

	tmp( 0, 0 ) = a( 0 ) * b( 0 );
	tmp( 1, 1 ) = a( 1 ) * b( 1 );
	tmp( 2, 2 ) = a( 2 ) * b( 2 );

	tmp( 0, 1 ) = a( 0 ) * b( 1 );
	tmp( 1, 0 ) = tmp( 0, 1 );
	tmp( 0, 2 ) = a( 0 ) * b( 2 );
	tmp( 2, 0 ) = tmp( 0, 2 );
	tmp( 1, 2 ) = a( 1 ) * b( 2 );
	tmp( 2, 1 ) = tmp( 1, 2 );

	return tmp;
}

/** 
 * Snap very small values to zero.
 * \param x    a value
 * \param eps  the epsilon
 * \return     snapped value 
*/
Real zeroSnap( Real x, Real eps )
{
	if ((x < eps) && (x > -eps))
		return 0.0;
	return x;
}


/**
 * Eigen decomposition of a 3x3 matrix.
 * \param eigen_vectors    eigenvectors sorted by values, return by reference
 * \param eigen_values     eigenvalues sorted, return by reference
 * \param t                matrix to be decomposed
 * \param abs_sort         true sort by absolute, false don't sort
*/
void eigenDecomposition3(
	vector<Vector3r>& eigen_vectors,
	vector<Real>& eigen_values,
	const Eigen::Matrix<double, 3, 3>& t,
	bool abs_sort )
{
	eigen_vectors.resize( 3 );
	eigen_values.resize( 3 );

	double A[3][3] = {
		{ t( 0,0 ), t( 1,0 ), t( 2,0 ) },   // row 0
		{ t( 0,1 ), t( 1,1 ), t( 2,1 ) },   // row 1
		{ t( 0,2 ), t( 1,2 ), t( 2,2 ) } }; // row 2

	double V[3][3] = { {0.,0.,0.}, {0.,0.,0.}, {0.,0.,0.} };
	double D[3] = { 0.,0.,0. };

	eigen_decomposition<3>( A, V, D );

	Vector3r unsorted_eva( D[0], D[1], D[2] );

	vector<Vector3r> unsorted_eve = {
		Vector3r( V[0][0], V[1][0], V[2][0] ),
		Vector3r( V[0][1], V[1][1], V[2][1] ),
		Vector3r( V[0][2], V[1][2], V[2][2] ) };

	unsorted_eva[0] = zeroSnap( unsorted_eva[0], 1e-6 );
	unsorted_eva[1] = zeroSnap( unsorted_eva[1], 1e-6 );
	unsorted_eva[2] = zeroSnap( unsorted_eva[2], 1e-6 );

	vector<pair<Real, int>> per = {
		{unsorted_eva[0], 0},
		{unsorted_eva[1], 1},
		{unsorted_eva[2], 2} };

	if (abs_sort)
		sort( per.begin(), per.end(),
			[]( pair<Real, int>x, pair<Real, int>y )
			{ return abs( x.first ) <= abs( y.first ); } );

	eigen_values = {
		unsorted_eva[per[0].second],
		unsorted_eva[per[1].second],
		unsorted_eva[per[2].second] };

	eigen_vectors = {
		unsorted_eve[per[0].second],
		unsorted_eve[per[1].second],
		unsorted_eve[per[2].second] };
}


vector<Vector3r> SurfaceTension_ZorillaRitter2020::getSphereSamplesLookUp(
	int N, Real supportRadius, int start, 
	const vector<float>& vec3, int mod )
{
	vector<Vector3r> points;
	int s = (start / 3) * 3; // ensure to be dividable by 3
	for (int i = 0; i < N; i++)
	{
		int i3 = s + 3 * i;
		points.push_back( supportRadius * Vector3r( vec3[i3 % mod], vec3[(i3 + 1) % mod], vec3[(i3 + 2) % mod] ) );
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
	size_t step = floor( tm.getTime() / timeStep );

	m_r2 = m_r1 * m_r2mult;

#ifdef ENABLE_STATE_CHANGE
	if( m_stateChange > 0.0 )
		changeState(m_stateChange, m_simulationStop );
#endif

	Simulation* sim = Simulation::getCurrent();

	const Real supportRadius = sim->getSupportRadius();
	const unsigned int numParticles = m_model->numActiveParticles();
	const Real k = m_surfaceTension;

	unsigned int NrOfSamples;

	const int fluidModelIndex = 0; // hardcode for now
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

			if (sim->numberOfNeighbors( 0, 0, i ) == 0)
			{
				m_mc_curv[i] = 1.0 / supportRadius;
				goto nextParticle;
			}

			

			const Vector3r& xi = m_model->getPosition( i );

			forall_fluid_neighbors_in_same_phase(
				Vector3r xjxi = (xj - xi);
				centerofMasses += xjxi; )
			
			centerofMasses /= supportRadius;

			// cache classifier input, could also be recomputed later to avoid caching
			m_classifier_input[i] = centerofMasses.norm() / double( numberOfNeighbours );
			
#ifdef RICH_OUTPUT			
			m_classifier_input2[i] = double( numberOfNeighbours );
#endif

			// -- if it is a surface classified particle
			if (classifyParticle( m_classifier_input[i], numberOfNeighbours )) //EvaluateNetwork also possible
			{

				// -- create monte carlo samples on particle
				std::vector<Vector3r> points;

				if( m_halton_sampling == RandomMethod::HALTON )
					points = getSphereSamplesLookUp(
						NrOfSamples, supportRadius, i * NrOfSamples, haltonVec323, haltonVec323.size() ); // 8.5 // 15.0(double) // 9.0(float)
				else
					points = getSphereSamplesRnd( NrOfSamples, supportRadius );


				//  -- remove samples covered by neighbor spheres
				forall_fluid_neighbors_in_same_phase(
					Vector3r xjxi = (xj - xi);
					for (int p = points.size() - 1; p >= 0; --p)
					{
						Vector3r vec = (points[p] - xjxi);
						Real dist = vec.squaredNorm();

						if (dist <= pow( (m_r2 / m_r1), 2 ) * supportRadius * supportRadius)
							points.erase( points.begin() + p );
					} )

				// -- estimate normal by left over sample directions
				for (int p = points.size() - 1; p >= 0; --p)
					m_mc_normals[i] += points[p];


				// -- if surface classified and non-overlapping neighborhood spheres
				if (points.size() > 0)
				{
					m_mc_normals[i].normalize();

					// -- estimate curvature by sample ratio and particle radii
					m_mc_curv[i] = (1 / supportRadius) * (-2.0) * pow( (1 - (m_r2 * m_r2 / (m_r1 * m_r1))), -0.5 ) *
						cos( acos( 1.0 - 2.0 * (double( points.size() ) / double( NrOfSamples )) ) + asin( m_r2 / m_r1 ) );

					m_classifier_output[i] = 1.0; // -- used to visualize surface points (blue in the paper)
				}
				else
				{
					// -- correct false positives to inner points
					m_mc_normals[i] = Vector3r::Zero();
					m_mc_curv[i] = 0.0;
					m_classifier_output[i] = 0.5; // -- used for visualize post-correction points (white in the paper)
				}

				// -- export vars
#ifdef RICH_OUTPUT
				m_nr_all_samples[i] = NrOfSamples;
				m_nr_area_samples[i] = points.size();
#endif
			}
			else
			{
				// -- used to visualize inner points (green in the paper)
				m_classifier_output[i] = 0.0;
			}

		nextParticle:;

		}
	}

	// ################################################################################################
	// ## second pass, compute normals and curvature and compute PCA normal 
	// ################################################################################################

	CubicKernel kernel;
	kernel.setRadius( supportRadius );

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


				double correctionForCurvature = 0;
				double correctionFactor = 0.0;

				Vector3r centroid = xi;
				Vector3r surfCentDir = Vector3r::Zero();

				// collect neighbors
				multimap<double, size_t> neighs;

				Eigen::Matrix<double, 3, 3> t = Eigen::Matrix<double, 3, 3>::Zero();
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
						&& classifyParticle( m_classifier_input[neighborIndex], nrNeighhbors, m_class_d_off ))
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
				//else if (m_normal_mode == NormalMethod::Conditional)
				//{
				//	if (m_nr_area_samples[i] < 7 &&
				//		float( m_nr_area_samples[i] ) / float ( m_nr_all_samples[i] ) < 0.33 &&
				//		neighs.size() >= 6)
				//		do_pca = true;
				//}
				else if (m_normal_mode == NormalMethod::MIX)
				{
					if (m_pca_N_mix >= 0.0 && m_pca_C_mix >= 0.0)
						do_pca = true;
				}

				if (do_pca)
				{
					centroid /= (t_count + 1);
					surfCentDir.normalize();

					const size_t neigh_limit = m_neighs_limit;
					size_t neigh_count = 0;

					for (pair<Real, size_t> nn : neighs)
					{
						size_t j = nn.second;
						Vector3r& xj = m_model->getPosition( j );
						Vector3r xjxi = (xj - centroid);

						double distanceji = xjxi.norm();

						double w = (kernel.W( distanceji ) / kernel.W_zero());
						neighCent += w * xj;

						Eigen::Matrix<double, 3, 3> dy = dyadic( xjxi, xjxi );

						t = t + w * dy;
						t_count++;

						if (neigh_count >= neigh_limit - 1)
							break;

						neigh_count++;
					}

					if (neighs.size() >= 4)
					{
						vector<Vector3r> pdt_ec( 3 );
						vector<Real> pdt_ev( 3 );

						eigenDecomposition3( pdt_ec, pdt_ev, t, true );

						if (pdt_ec[0].dot( m_mc_normals[i] ) < 0.0)
							m_pca_normals[i] = -1.0 * pdt_ec[0];
						else
							m_pca_normals[i] = pdt_ec[0];

						m_pca_curv[i] = 3.0 * pdt_ev[0] / (pdt_ev[0] + pdt_ev[1] + pdt_ev[2]);

						if (surfCentDir.dot( m_pca_normals[i] ) > 0.0)
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
					((1 - m_tau) * m_mc_curv[i] + m_tau * correctionForCurvature) /
					((1 - m_tau) * 1 + m_tau * correctionFactor);
					//-0.075; // no offset found on the saddle example, but not used generally.

				//Vector3r force = m_normals_smoothed[i] * k * m_curvatures_smoothed[i];

//				if( i % 200 == 0 )
//				{
//					printf( "m_curvatures          [%d]: %.3f\n", i, m_curvatures[i] );
//					printf( "m_curvatures_corrected[%d]: %.3f\n", i, m_curvatures_corrected[i] );
//					printf( "m_curvatures_smoothed [%d]: %.3f\n", i, m_curvatures_smoothed[i] );
//				} 

				//ai -= (1 / m_model->getMass( i )) * force;

				//ai -= (1 / m_model->getMass( i )) * final_normal * k * final_curvature;

/*
				if( do_pca || m_ai_mode > 0 )
				{

					if( m_ai_mode == 1 ) /// pca normal, mc curvature
					{
						ai -= (1 / m_model->getMass( i )) * m_normals_pca[i] * k *
							((1 - m_tau) * m_curvatures[i] + m_tau * correctionForCurvature)
							/ ((1 - m_tau) * 1 + m_tau * correctionFactor);
					}
					else if( m_ai_mode == 2 ) // pca normal, sph curvature
					{
						ai -= (1 / m_model->getMass( i )) * m_normals_pca[i] * k * m_curv_by_sphericity[i];
					}
					else //if( m_ai_mode == 3 ) // pca normal, mixed sph curvature? // to be done
					{
						ai -= (1 / m_model->getMass( i )) * m_normals_pca[i] * k * m_curv_by_sphericity[i];
					}
				}
				else // full mc
				{
					ai -= (1 / m_model->getMass( i )) * normalCorrection * k *
						((1 - m_tau) * m_curvatures[i] + m_tau * correctionForCurvature)
						/ ((1 - m_tau) * 1 + m_tau * correctionFactor);
				}
				*/
			}
		}
	}


	// ################################################################################################
	// ## third pass, final blending and temporal smoothing
	// ################################################################################################
	
	m_CS_smooth_passes = max( 1, m_CS_smooth_passes );

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
					double CsCorr = 0.0;

					const Vector3r& xi = m_model->getPosition( i );

					forall_fluid_neighbors_in_same_phase( 
						if (m_mc_normals[neighborIndex] != Vector3r::Zero())
						{
							CsCorr += m_pca_curv[neighborIndex];
							count++;
						} )
					

					if (count > 0)
						m_pca_curv_smooth[i] = 0.25 * m_pca_curv_smooth[i] + 0.75 * CsCorr / count;
					else
						m_pca_curv_smooth[i] = m_pca_curv[i];

					m_pca_curv_smooth[i] /= supportRadius;

					m_pca_curv_smooth[i] *= 20.0;

					if (m_pca_curv_smooth[i] > 0.0)
						m_pca_curv_smooth[i] = min( 0.5f / supportRadius, m_pca_curv_smooth[i] );
					else
						m_pca_curv_smooth[i] = max( -0.5f / supportRadius, m_pca_curv_smooth[i] );


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
						//final_curvature = m_curvatures_smoothed[i];
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
									(1.0 - m_pca_C_mix) * m_mc_curv_smooth[i];
						}
						//						if( i % 200 == 0 )
						//							printf( "c_mc_corr: %.3f, c_pca: %.3f, c_pca_smo: %.3f, mixed: %.3f, smr: %.3f\n",
						//								m_curvatures_smoothed[i], m_curv_by_sphericity[i], m_curv_by_sphericity_corr[i], final_curvature, supportRadius );
					}
					/*
					else if (m_normal_mode == NormalMethod::Conditional)
					{
						if ( //m_nr_area_samples[i] < 7 &&
							float( m_nr_area_samples[i] ) / float ( m_nr_all_samples[i] ) < 0.33 &&
							m_normals_pca[i] != Vector3r::Zero())
						{
							if (i % 100 == 0) puts( "Switched normal PCA" );
							final_normal = m_normals_pca[i];
							final_curvature = m_curv_by_sphericity_corr[i];
						}
						else
						{
							if (i % 100 == 0) puts( "Switched normal MC" );
							final_normal = m_normals_smoothed[i];
							final_curvature = m_curvatures_smoothed[i];
						}
					}
					*/

					//m_normals_smoothed[i] = final_normal;

					if (m_temporal_smoothing)
						m_final_curvatures[i] = 0.05 * final_curvature + 0.95 * m_final_curvatures_old[i];
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
						m_final_curvatures[i] = 0.95 * m_final_curvatures_old[i];
					else
						m_final_curvatures[i] = 0.0;

					m_final_curvatures_old[i] = m_final_curvatures[i];
				}


			}
		}
	}


	// ################################################################################################
	// ## message printing and direct ascii export
	// ################################################################################################

	static int count = 0;

	/* // debug message
	if( count != 0 && count % 100 == 0 )
		printf( "c_avg_mc: %.3f, c_avg_pca: %.3f, Nr particels: %d, nr samples: %d, nr surf part: %d, sum forces: %f %f %f\n",
			c_mc_smoothed, c_pca_smoothed, (int) numParticles, NrOfSamples, m_nr_surface_particles_reduce,
			m_sum_surface_forces_reduce.x(), m_sum_surface_forces_reduce.y(), m_sum_surface_forces_reduce.z() );
	*/

#ifdef RICH_OUTPUT
	if (count != 0 && count % 100 == 0 && m_start_export > 0.0 && tm.getTime() > m_start_export)
		storeASCIIData( step );
#endif

	count++;
}


// -- Helper functions for evaluation. Export data, state change, ...
void SurfaceTension_ZorillaRitter2020::storeASCIIData( size_t step )
{
	Simulation* sim = Simulation::getCurrent();
	const Real R = sim->getParticleRadius();

	string path = "C:/Source/2019/surfacetensionextension/Measure/";

	string r_txt = to_string( R ).substr(0, 5);
	string t_txt = to_string( m_tau ).substr( 0, 3 );

	std::replace( r_txt.begin(), r_txt.end(), '.', '_' );
	std::replace( t_txt.begin(), t_txt.end(), '.', '_' );

	string file_name = 
		"saddle-" + r_txt + "-" + to_string( m_CsdFix ) + 
		"-d" + to_string( int( m_class_d ) ) + 
		"-t" + t_txt +
		"-" + to_string( step ) + ".csv";

	string file_name_class =
		"saddle-" + r_txt + "-" + to_string( m_CsdFix ) +
		"-d" + to_string( int( m_class_d ) ) +
		"-t" + t_txt +
		"-" + to_string( step ) + "_class.csv";

	FILE* file = fopen( (path + file_name).c_str(), "w" );
	FILE* file_class = fopen( (path + file_name_class).c_str(), "w" );


	fprintf( file, "X, Y, Z, NX, NY, NZ, NSX, NSY, NSZ, NPX, NPY, NPZ, C, C_Corr, C_Smooted, CN, CN_Smoothed, Cs, Cs_Smoothed, rAreaSample, NrSamples, ClassIn1, ClassIn2\n" );
	
	fprintf( file_class, "X, Y, Z, ClassIn1, ClassIn2\n" );

	for( size_t i = 0; i < m_mc_normals.size(); i++ )
	{
		const Vector3r& ni  = m_mc_normals[i];
		const Vector3r& nsi = m_mc_normals_smooth[i];
		const Vector3r& npi = m_pca_normals[i];
		const Vector3r& xi = m_model->getPosition( i );

#ifdef RICH_OUTPUT
		fprintf(
			file_class, "%.6f, %.6f, %.6f, %.6f, %.6f\n",
			xi.x(), xi.y(), xi.z(),
			m_classifier_input[i],
			m_classifier_input2[i] );

		if( m_mc_normals[i] != Vector3r::Zero() ) 
		{
			if( abs( xi.x() ) <= 0.8 && 
				abs( xi.z() ) <= 0.8 && 
				xi.y() < 1.0 )
				fprintf( 
					file, "%.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %d, %d\n",
					xi.x(), xi.y(), xi.z(), 
					ni.x(), ni.y(), ni.z(),
					nsi.x(), nsi.y(), nsi.z(),
					npi.x(), npi.y(), npi.z(),
					m_mc_curv[i], 
					m_mc_curv_smooth[i],
					m_pca_curv[i],
					m_pca_curv_smooth[i],
					m_nr_area_samples[i], 
					m_nr_all_samples[i] );
		}
#else
		puts( "!!!!  RICH_OUTPUT has to be enabled (#define, or compile flag -DRICH_OUTPUT)" );
		cin.ignore();
#endif
	}

	fclose( file );
	fclose( file_class );
	printf( "stored: %s\n", (path + file_name).c_str() );
}


#ifdef ENABLE_STATE_CHANGE

/** Change state and quit simulation at given time steps. 
 *	For faster convergence in the saddle example, viscocity and sampling configuration
 *  was changed at t=1.2. Higher viscocity and lower sampling at the start. 
 */
void SurfaceTension_ZorillaRitter2020::changeState( double time_change, double time_stop )
{
	auto& tm = *SPH::TimeManager::getCurrent();

	m_Csd = m_CsdPostState;

	static bool visc_consolidation = true;
	if( tm.getTime() > time_change&& visc_consolidation )
	{
		visc_consolidation = false;
		Simulation* sim = Simulation::getCurrent();
		FluidModel* model = sim->getFluidModel( 0 );

		Viscosity_Standard* v = dynamic_cast<Viscosity_Standard*>(model->getViscosityBase());
		if( v )
			v->setViscosity( 0.01 );

		m_surfaceTension = m_sTPostState;

		m_Csd = m_CsdPostState;
		m_CsdFix = m_CsdFixPostState;
	}

	if( time_stop > 0.0 && tm.getTime() > time_stop )
		exit( 0 );
}
#endif


/*
// surface point creation variants.
std::vector<Vector3r> SurfaceTension_ZorillaRitter2020::GetUniformSurfacePoints( int N, Real supportRadius )
{
	std::vector<Vector3r> points;

	int n = round( sqrt( float( N ) ) );

	for( int i = 0; i < n; i++ )
	{
		for( int j = 0; j < n; j++ )
		{
			double theta = 2.0 * M_PI * (float(i) / n);
			double phi   = acos( 1.0 - 2.0 * (float( j ) / n) );

			points.push_back( anglesToVec( theta, phi, supportRadius ) );
		}
	}

	return points;
}

std::vector<Vector3r> SurfaceTension_ZorillaRitter2020::GetLookupTableSurfacePoints(
	int N, 
	Real supportRadius, 
	int start, 
	const std::vector<double>& hal,
	int mod)
{
	std::vector<Vector3r> points;

	int s = (start / 2) * 2; // ensure to be even

	for( int i = 0; i < N; i++ )
	{
		double theta = 2.0 * M_PI * hal[(s + (2 * i)) % mod];
		double phi   = acos( 1.0 - 2.0 * hal[(s + (2 * i + 1)) % mod] );

		points.push_back( anglesToVec( theta, phi, supportRadius ) );
	}

	return points;
}

std::vector<Vector3r> SurfaceTension_ZorillaRitter2020::GetLookupTableSurfacePoints(
	int N,
	Real supportRadius,
	int start,
	const std::vector<float>& hal,
	int mod )
{
	std::vector<Vector3r> points;

	int s = (start / 2) * 2; // essure to be even

	for( int i = 0; i < N; i++ )
	{
		double theta = 2.0 * M_PI * hal[(s + (2 * i)) % mod];
		double phi = acos( 1.0 - 2.0 * hal[(s + (2 * i + 1)) % mod] );

		points.push_back( anglesToVec( theta, phi, supportRadius ) );
	}

	return points;
}

std::vector<Vector3r> SurfaceTension_ZorillaRitter2020::getSphereSamplesLookUp(
	int N, Real supportRadius, int start, const std::vector<double>& vec3, int mod )
{
	std::vector<Vector3r> points;
	int s = (start / 3) * 3; // ensure to be dividable by 3
	for( int i = 0; i < N; i++ )
	{
		int i3 = s + 3 * i;
		points.push_back( supportRadius * Vector3r( vec3[i3%mod], vec3[(i3 + 1)%mod], vec3[(i3+2)%mod] ) );
	}
	return points;
}
*/


/*
// some time measure
void SurfaceTension_ZorillaRitter2020::classificationTiming()
{


	const Real supportRadius = m_model->getSupportRadius();
	const unsigned int numParticles = m_model->numActiveParticles();
	const Real k = m_surfaceTension;
	Real timeStep = SPH::TimeManager::getCurrent()->getTimeStepSize();
	double r1 = supportRadius;
	double r2 = 0.8 * r1;// 0.8 best results empirically according to Fernando
	const unsigned int NrOfSamples = int(100000 * timeStep);

	double smoothingFactor = 0.5;


#pragma omp parallel default(shared)
	{
#pragma omp for schedule(static)
		for (int i = 0; i < (int)numParticles; i++)
		{

			m_normals[i] = Vector3r::Zero();
			m_model->setCurvature(i, 0);

			Vector3r centerofMasses = Vector3r::Zero();
			int numberOfNeighbours = m_model->numberOfNeighbors(0, i);


			if (m_model->numberOfNeighbors(0, i) == 0) { goto nextParticle; }
			const Vector3r& xi = m_model->getPosition(0, i);

			for (unsigned int j = 0; j < m_model->numberOfNeighbors(0, i); j++)
			{
				unsigned int neighborIndex = m_model->getNeighbor(0, i, j);
				Vector3r& xj = m_model->getPosition(0, neighborIndex);
				Vector3r xjxi = (xj - xi);

				centerofMasses += xjxi / supportRadius;
			}
			if (classifyParticle(centerofMasses.norm() / double(numberOfNeighbours), numberOfNeighbours)) //EvaluateNetwork also possible
			{

			}

		nextParticle:;

		}
	}

}
*/


/*
// this was used for training data generation
void SurfaceTension_ZorillaRitter2020::stepData()
{
	Simulation* sim = Simulation::getCurrent();

	const Real supportRadius = sim->getSupportRadius();

	const unsigned int numParticles = m_model->numActiveParticles();
	const Real k = m_surfaceTension;
	Real timeStep = SPH::TimeManager::getCurrent()->getTimeStepSize();

	const unsigned int NrOfSamples = int( 1000000 * timeStep );

#pragma omp parallel default(shared)
	{
#pragma omp for schedule(static)  
		for (int i = 0; i < (int)numParticles; i++)
		{
			std::vector<Vector3r> points = getSphereSamplesRnd( NrOfSamples, supportRadius );

			m_normals[i] = Vector3r::Zero();

			Vector3r centerofMasses = Vector3r::Zero();
			int numberOfNeighbours = sim->numberOfNeighbors( 0, i );
			int SurfaceFlag = 0;

			if (sim->numberOfNeighbors( 0, i ) == 0) { goto nextParticle; }
			const Vector3r& xi = m_model->getPosition( i );

			for (unsigned int j = 0; j < sim->numberOfNeighbors( 0, i ); j++)
			{
				unsigned int neighborIndex = sim->getNeighbor( 0, i, j );
				Vector3r& xj = m_model->getPosition( neighborIndex );
				Vector3r xjxi = (xj - xi);

				centerofMasses += xjxi / supportRadius;

				for (int p = points.size() - 1; p >= 0; --p)
				{
					Vector3r vec = (points[p] - xjxi);
					Real dist = vec.squaredNorm();
					if (dist <= 0.16 * supportRadius * supportRadius)
					{
						points.erase( points.begin() + p );
					}
				}
			}

			for (int p = points.size() - 1; p >= 0; --p)
			{
				m_normals[i] += points[p];
			}
			if (points.size() > 0)
			{
				SurfaceFlag = 1;
				m_normals[i].normalize();
				//m_curvatures[i] = 0.917*(points.size() / double(NrOfSamples)-0.13);
				//m_curvatures[i] = 0.8663*(points.size() / double(NrOfSamples) - 0.1) - 0.0115;
				//m_curvatures[i] = 2.26*(points.size() / double(NrOfSamples)-0.082);
				m_curvatures[i] =
					-1.5287 * pow( (points.size() / double( NrOfSamples )) - 0.1, 2 )
					+ 1.0134 * pow( (points.size() / double( NrOfSamples )) - 0.1, 1 )
					- 0.0044;
				//m_curvatures[i] =
				//	12.924* pow((points.size() / double(NrOfSamples)), 3)
				//	- 9.6294* pow((points.size() / double(NrOfSamples)), 2)
				//	+ 3.8219* pow((points.size() / double(NrOfSamples)), 1)
				//	-0.347;
			}

			// some output to be enabled
			//fprintf(SurfaceParticlesFile, "%f;%d;%d\n", centerofMasses.norm()/double(numberOfNeighbours), numberOfNeighbours, SurfaceFlag);

		nextParticle:;
		}
	}
}
*/

void SurfaceTension_ZorillaRitter2020::reset()
{
}

