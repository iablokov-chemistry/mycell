//
// Created by root on 15.07.2019.
//

#include <iostream>
#include <cmath>
#include <algorithm>

#include "System.h"
#include "Reaction.h"
#include "Interaction.h"

using namespace std;

static void interactionFunctionLJ(Particle* _particle1, Particle* _particle2, Interaction* _interaction, double& Fx, double& Fy, double& Fz, double& E)
{
	System* system = System::getInstance();

	double eps = _interaction->params[0];
	double sigma = _interaction->params[1];

	E = 0.0;
	Fx = 0.0;
	Fy = 0.0;
	Fz = 0.0;

	if (system->__dr[_particle2->idx] <= _interaction->cutoff)
	{
		// dr can't exceed 2
		//double dr = max(system->__dr[_particle2->idx], (_particle1->particleType->r + _particle2->particleType->r) / 2.0);
		double dr = max(system->__dr[_particle2->idx], 2.0);

		double p = sigma / dr;
		double p2 = p * p;
		double p6 = p2 * p2 * p2;
		double p12 = p6 * p6;

		E = 4 * eps * (p12 - p6);
		//cout << _particle1->particleId << " " << _particle2->particleId << endl;
		//cout << dx * 1E10 << " " << dy * 1E10 << " " << dz * 1E10 << endl;
		//cout << "> p=" << p << "; sigma=" << sigma * 1E10 << endl;
		//cout << "> r=" << r * 1E10 << "; cut=" << _reaction->interactionCutoff * 1E10 << "; E=" << E << endl;

		auto eedr2 = 24 * eps * (2 * p12 - p6) / (dr * dr);
		Fx = eedr2 * system->__dx[_particle2->idx]; // eV / A
		Fy = eedr2 * system->__dy[_particle2->idx]; // eV / A
		Fz = eedr2 * system->__dz[_particle2->idx]; // eV / A


		// old
		// Fx = E * system->__dx[_particle2->idx] / dr;
		// Fy = E * system->__dy[_particle2->idx] / dr;
		// Fz = E * system->__dz[_particle2->idx] / dr;
	}
}

Interaction::Interaction(string _name)
{
	this->name = _name;
};

void Interaction::setParticleTypes(string _particleType1, string _particleType2)
{
	System* system = System::getInstance();
	this->particleType1 = system->particleTypes[_particleType1];
	this->particleType2 = system->particleTypes[_particleType2];
}

void Interaction::setInteractionType(string _type)
{
	this->type = _type;
	if (_type == "LJ")
		this->interactionFunction = interactionFunctionLJ;

}

void Interaction::addParam(double _paramValue)
{
	this->params[this->__paramCount++] = _paramValue;
}