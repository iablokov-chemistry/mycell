//
// Created by root on 15.07.2019.
//

#ifndef MOLECULARDYNAMICS_INTERACTION_H
#define MOLECULARDYNAMICS_INTERACTION_H

#include <string>
#include <map>
#include <vector>
#include <functional>

//#include "Particle.h"
//#include "ParticleType.h"
class System;
class Particle;
class ParticleType;

using namespace std;

class Interaction
{
public:
	Interaction(string _name);
	Interaction() {}; // for iterations in containers

	void setParticleTypes(string _particleType1, string _particleType2);
	void setInteractionType(string _type);
	void addParam(double _paramValue);

	string name;
	string type = "";
	double cutoff = 0.0;
	double params[5];
	int __paramCount = 0;

	ParticleType* particleType1, * particleType2;

	function<void(Particle*, Particle*, Interaction*, double&, double&, double&, double&)> interactionFunction = nullptr;

};


#endif //MOLECULARDYNAMICS_INTERACTION_H
