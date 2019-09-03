#include "PartioReaderWriter.h"
#include "extern/partio/src/lib/Partio.h"
#include "FileSystem.h"

using namespace Utilities;


bool PartioReaderWriter::readParticles(const std::string &fileName, const Vector3r &translation, const Matrix3r &rotation, const Real scale,
	std::vector<Vector3r> &positions, std::vector<Vector3r> &velocities)
{
	if (!FileSystem::fileExists(fileName))
		return false;

	Partio::ParticlesDataMutable* data = Partio::read(fileName.c_str());
	if (!data)
		return false;

	unsigned int posIndex = 0xffffffff;
	unsigned int velIndex = 0xffffffff;

	for (int i = 0; i < data->numAttributes(); i++)
	{
		Partio::ParticleAttribute attr;
		data->attributeInfo(i, attr);
		if (attr.name == "position")
			posIndex = i;
		else if (attr.name == "velocity")
			velIndex = i;
	}

	Partio::ParticleAttribute attr;

	if (posIndex != 0xffffffff)
	{
		unsigned int fSize = (unsigned int) positions.size();
		positions.resize(fSize + data->numParticles());
		data->attributeInfo(posIndex, attr);
		for (int i = 0; i < data->numParticles(); i++)
		{
			const float *pos = data->data<float>(attr, i);
			Vector3r x(pos[0], pos[1], pos[2]);
			x = rotation * (x*scale) + translation;
			positions[i + fSize] = x;
		}
	}

	if (velIndex != 0xffffffff)
	{
		unsigned int fSize = (unsigned int) velocities.size();
		velocities.resize(fSize + data->numParticles());
		data->attributeInfo(velIndex, attr);
		for (int i = 0; i < data->numParticles(); i++)
		{
			const float *vel = data->data<float>(attr, i);
			Vector3r v(vel[0], vel[1], vel[2]);
			velocities[i + fSize] = v;
		}
	}
	else
	{
		unsigned int fSize = (unsigned int) velocities.size();
		velocities.resize(fSize + data->numParticles());
		for (int i = 0; i < data->numParticles(); i++)
			velocities[i + fSize].setZero();
	}

	data->release();
	return true;
}

bool PartioReaderWriter::readParticles(const std::string &fileName, const Vector3r &translation, const Matrix3r &rotation, const Real scale,
	std::vector<Vector3r> &positions, std::vector<Vector3r> &velocities, Real &particleRadius)
{
	if (!FileSystem::fileExists(fileName))
		return false;

	Partio::ParticlesDataMutable* data = Partio::read(fileName.c_str());

	if (!data)
		return false;

	unsigned int posIndex = 0xffffffff;
	unsigned int velIndex = 0xffffffff;
	unsigned int radiusIndex = 0xffffffff;

	for (int i = 0; i < data->numAttributes(); i++)
	{
		Partio::ParticleAttribute attr;
		data->attributeInfo(i, attr);
		if (attr.name == "position")
			posIndex = i;
		else if (attr.name == "velocity")
			velIndex = i;
		else if (attr.name == "pscale")
			radiusIndex = i;
	}

	Partio::ParticleAttribute attr;

	if (posIndex != 0xffffffff)
	{
		unsigned int fSize = (unsigned int)positions.size();
		positions.resize(fSize + data->numParticles());
		data->attributeInfo(posIndex, attr);
		for (int i = 0; i < data->numParticles(); i++)
		{
			const float *pos = data->data<float>(attr, i);
			Vector3r x(pos[0], pos[1], pos[2]);
			x = rotation * (x*scale) + translation;
			positions[i + fSize] = x;
		}
	}

	if (velIndex != 0xffffffff)
	{
		unsigned int fSize = (unsigned int)velocities.size();
		velocities.resize(fSize + data->numParticles());
		data->attributeInfo(velIndex, attr);
		for (int i = 0; i < data->numParticles(); i++)
		{
			const float *vel = data->data<float>(attr, i);
			Vector3r v(vel[0], vel[1], vel[2]);
			velocities[i + fSize] = v;
		}
	}
	else
	{
		unsigned int fSize = (unsigned int)velocities.size();
		velocities.resize(fSize + data->numParticles());
		for (int i = 0; i < data->numParticles(); i++)
			velocities[i + fSize].setZero();
	}

	if (radiusIndex != 0xffffffff)
	{
		data->attributeInfo(radiusIndex, attr);
		const float *radius = data->data<float>(attr, 0);
		particleRadius = radius[0];
	}

	data->release();

	return true;
}

bool PartioReaderWriter::readParticles(const std::string &fileName, const Vector3r &translation, const Matrix3r &rotation, const Real scale,
	std::vector<Vector3r> &positions)
{
	if (!FileSystem::fileExists(fileName))
		return false;

	Partio::ParticlesDataMutable* data = Partio::read(fileName.c_str());

	if (!data)
		return false;

	unsigned int posIndex = 0xffffffff;

	for (int i = 0; i < data->numAttributes(); i++)
	{
		Partio::ParticleAttribute attr;
		data->attributeInfo(i, attr);
		if (attr.name == "position")
		{
			posIndex = i;
			break;
		}
	}

	Partio::ParticleAttribute attr;
	if (posIndex != 0xffffffff)
	{
		unsigned int fSize = (unsigned int)positions.size();
		positions.resize(fSize + data->numParticles());
		data->attributeInfo(posIndex, attr);
		for (int i = 0; i < data->numParticles(); i++)
		{
			const float *pos = data->data<float>(attr, i);
			Vector3r x(pos[0], pos[1], pos[2]);
			x = rotation * (x*scale) + translation;
			positions[i + fSize] = x;
		}
	}

	data->release();

	return true;
}


void PartioReaderWriter::writeParticles(const std::string &fileName, const unsigned int numParticles, const Vector3r *particlePositions,
	const Vector3r *particleVelocities, const Real particleRadius)
{
	if (numParticles == 0)
		return;

	Partio::ParticlesDataMutable& particleData = *Partio::create();
	Partio::ParticleAttribute posAttr = particleData.addAttribute("position", Partio::VECTOR, 3);
	Partio::ParticleAttribute velAttr;
	if (particleVelocities != NULL)
		velAttr = particleData.addAttribute("velocity", Partio::VECTOR, 3);
	Partio::ParticleAttribute scaleAttr;
	if (particleRadius != 0.0)
		scaleAttr = particleData.addAttribute("pscale", Partio::FLOAT, 1);
	Partio::ParticleAttribute idAttr = particleData.addAttribute("id", Partio::INT, 1);

	for (unsigned int i = 0; i < numParticles; i++)
	{
		Partio::ParticleIndex index = particleData.addParticle();
		float* pos = particleData.dataWrite<float>(posAttr, index);
		int* id = particleData.dataWrite<int>(idAttr, index);

		const Vector3r &x = particlePositions[i];
		pos[0] = (float)x[0];
		pos[1] = (float)x[1];
		pos[2] = (float)x[2];

		if (particleVelocities != NULL)
		{
			float* vel = particleData.dataWrite<float>(velAttr, index);
			const Vector3r &v = particleVelocities[i];
			vel[0] = (float)v[0];
			vel[1] = (float)v[1];
			vel[2] = (float)v[2];
		}

		if (particleRadius != 0.0)
		{
			float* scale = particleData.dataWrite<float>(scaleAttr, index);
			scale[0] = (float)particleRadius;
		}
		id[0] = i;
	}

	Partio::write(fileName.c_str(), particleData, true);
	particleData.release();
}
