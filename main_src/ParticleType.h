//
// Created by root on 11.02.2019.
//

#ifndef MOLECULARDYNAMICS_PARTICLETYPE_H
#define MOLECULARDYNAMICS_PARTICLETYPE_H


#include <string>
#include <list>
#include <map>
#include "Reaction.h"
#include "Interaction.h"

class Particle;

using namespace std;

class ParticleType
{
public:
	ParticleType(string _type);
	ParticleType() = default; // default constructor

	void addCollision(ParticleType* _particleType, Reaction* _reaction);
	void addDecay(Reaction* _reaction);
	void addInteraction(ParticleType* _particleType, Interaction* _interaction);

	void addExistingParticle(Particle* _particle);
	void removeExistingParticle(Particle* _particle);

	void calcParams();

	string type;
	double xi_medium = 0.0, xi_micelle = 0.0;
	double r = 0.0;

	bool medium_bound = false, micelle_bound = false, fixed_concentration = false, final_product = false;
	int fixed_concentration_count = 0, fixed_added = 0;

	map<ParticleType*, list<Reaction*>> collisions;
	list<Reaction*> decays;

	map<ParticleType*, Interaction*> interactions;
	list<Particle*> existingParticles;

};


#endif //MOLECULARDYNAMICS_PARTICLETYPE_H
