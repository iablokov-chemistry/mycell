//
// Created by root on 11.02.2019.
//

#include <iostream>
#include <string>
#include <map>
#include <cmath>
#include "ParticleType.h"
#include "Particle.h"
#include "System.h"

using namespace std;

Particle::Particle(ParticleType* _particleType)
{
	this->particleType = _particleType;
};

double Particle::x()
{
	static System* system = System::getInstance();
	return system->rx[this->idx];
}

double Particle::y()
{
	static System* system = System::getInstance();
	return system->ry[this->idx];
}

double Particle::z()
{
	static System* system = System::getInstance();
	return system->rz[this->idx];
}


void Particle::setX(double _x)
{
	static System* system = System::getInstance();

	system->rx[this->idx] = _x;
	// fix for boundaries
	if (system->bcX)
		system->rx[this->idx] -= round(this->x() / system->boxWidth - 0.5) * system->boxWidth;
}

void Particle::setY(double _y)
{
	static System* system = System::getInstance();

	system->ry[this->idx] = _y;
	// fix for boundaries
	if (system->bcY)
		system->ry[this->idx] -= round(this->y() / system->boxDepth - 0.5) * system->boxDepth;
}

void Particle::setZ(double _z)
{
	static System* system = System::getInstance();

	system->rz[this->idx] = _z;
	// fix for boundaries
	if (system->bcZ)
		system->rz[this->idx] -= round(this->z() / system->boxHeight - 0.5) * system->boxHeight;
}

void Particle::updateX(double _dx)
{
	this->setX(this->x() + _dx);
}

void Particle::updateY(double _dy)
{
	this->setY(this->y() + _dy);
}

void Particle::updateZ(double _dz)
{
	this->setZ(this->z() + _dz);
}

unsigned int Particle::subBoxId() { return this->__subBoxId; }
unsigned int Particle::subBoxW() { return this->__subBoxW; }
unsigned int Particle::subBoxD() { return this->__subBoxD; }
unsigned int Particle::subBoxH() { return this->__subBoxH; }

void Particle::setSubBox(unsigned int _W, unsigned int _D, unsigned int _H)
{
	static System* system = System::getInstance();

	this->__subBoxW = _W;
	this->__subBoxD = _D;
	this->__subBoxH = _H;

	this->__subBoxId = system->subBoxIndex(_W, _D, _H);
}
