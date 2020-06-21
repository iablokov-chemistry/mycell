//
// Created by root on 11.02.2019.
//

#ifndef MOLECULARDYNAMICS_REACTION_H
#define MOLECULARDYNAMICS_REACTION_H


#include <string>
#include <map>
#include <vector>
#include <functional>

//#include "Particle.h"
class Particle;
class ParticleType;

using namespace std;

class Reaction
{
public:
	Reaction(string _name);
	Reaction() {}; // for iterations in containers
	void addReactant(string _reactantName);
	void addProduct(string _productName);

	string name;
	double probability = 0.0;
	double molecularity = 0;
	double energy = 0.0;
	double cutoff = 0.0;

	int collision_count = 0;
	int reaction_count = 0;

	vector<ParticleType*> reactants;
	vector<ParticleType*> products;
};



#endif //MOLECULARDYNAMICS_REACTION_H
