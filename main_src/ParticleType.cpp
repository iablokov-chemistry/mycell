//
// Created by root on 11.02.2019.
//
using namespace std;

#include <iostream>
#include "ParticleType.h"
#include "Reaction.h"
#include "System.h"


ParticleType::ParticleType(string _type)
{
	this->type = _type;
}

void ParticleType::addInteraction(ParticleType* _particleType, Interaction* _interaction)
{
	this->interactions[_particleType] = _interaction;
}

void ParticleType::addCollision(ParticleType* _particleType, Reaction* _reaction)
{
	this->collisions[_particleType].emplace_back(_reaction);
}

void ParticleType::addDecay(Reaction* _reaction)
{
	this->decays.emplace_back(_reaction);
}

void ParticleType::addExistingParticle(Particle* _particle)
{
	this->existingParticles.emplace_back(_particle);
}

void ParticleType::removeExistingParticle(Particle* _particle)
{
	this->existingParticles.remove(_particle);
}

void ParticleType::calcParams()
{
	System* system = System::getInstance();

	this->xi_medium = 6 * system->PI * system->viscosity * this->r * 1E-10; // kg /s * m / m
	this->xi_micelle = 6 * system->PI * system->micelle_viscosity * this->r * 1E-10; // kg /s * m / m
}