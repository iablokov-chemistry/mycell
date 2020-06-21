//
// Created by root on 11.02.2019.
//
using namespace std;

#include "System.h"
#include "Reaction.h"


Reaction::Reaction(string _name)
{
	this->name = _name;
};

void Reaction::addReactant(string _reactantName)
{
	System* system = System::getInstance();
	this->reactants.push_back(system->particleTypes[_reactantName]);
	this->molecularity++;
};

void Reaction::addProduct(string _productName)
{
	System* system = System::getInstance();
	this->products.push_back(system->particleTypes[_productName]);
};


