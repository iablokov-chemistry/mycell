//
// Created by root on 11.02.2019.
//

#ifndef MOLECULARDYNAMICS_PARTICLE_H
#define MOLECULARDYNAMICS_PARTICLE_H


#include <string>
#include <list>
#include "Reaction.h"

using namespace std;

// ...
class ParticleType;

class Particle
{
public:
	Particle(ParticleType* _particleType);

	ParticleType* particleType;
	bool reacted = false;

	int nearby = 0;

	long int id = 0;
	long int idx = 0;

	list<Reaction> reactions; // list of reactions this particle participates in

	double x();
	double y();
	double z();

	unsigned int subBoxId();
	unsigned int subBoxW();
	unsigned int subBoxD();
	unsigned int subBoxH();

	void setSubBox(unsigned int _W, unsigned int _D, unsigned int _H);

	void setX(double _x);
	void setY(double _y);
	void setZ(double _z);

	void updateX(double _x);
	void updateY(double _y);
	void updateZ(double _z);

	bool operator == (const Particle& _particle) const { return this->id == _particle.id; }
	bool operator != (const Particle& _particle) const { return this->id != _particle.id; }
private:
	unsigned int __subBoxId, __subBoxW, __subBoxD, __subBoxH;

};


#endif //MOLECULARDYNAMICS_PARTICLE_H
