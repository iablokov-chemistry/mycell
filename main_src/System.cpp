//
// Created by root on 11.02.2019.
//

#include <iostream>
#include <ctime>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <tuple>
#include <set>
#include <map>
#include <algorithm>
#include <random>
#include <cstdlib>
#include <algorithm>

#include "System.h"
#include "Particle.h"
#include "Interaction.h"
#include "Reaction.h"
#include "ParticleType.h"

#include "json.h"

System* System::__instance = nullptr;

System::System()
{
	std::random_device rd;
	srand(rd());
}

System* System::getInstance()
{
	if (__instance == nullptr)
	{
		__instance = new System();
	}

	return __instance;
}

void System::setThreadId(string _threadId)
{
	this->threadId = _threadId;
}

void System::__init_micelle()
{
	if (this->micelle_position == "center")
	{
		this->micelle_geometry_x = this->boxWidth / 2.0;
		this->micelle_geometry_y = this->boxDepth / 2.0;
		this->micelle_geometry_z = this->boxHeight / 2.0;
	}
}

bool System::__is_in_micelle(double _x, double _y, double _z)
{
	double dx = this->micelle_geometry_x - _x;
	double dy = this->micelle_geometry_y - _y;
	double dz = this->micelle_geometry_z - _z;
	double r = sqrt(dx * dx + dy * dy + dz * dz);
	
	if (r < this->micelle_geometry_r)
		return true;
	else
		return false;
}

void System::load_config()
{
	Json::Value config;

	ifstream fConfig(this->working_path + this->cfg_file, ifstream::binary);
	fConfig >> config;
	fConfig.close();

	// loading run settings

	this->dt = config["run"]["dt"].asDouble();
	this->steps = config["run"]["steps"].asUInt();
	this->tsteps = config["run"]["tsteps"].asUInt();
	this->maxParticles = config["run"]["max-particles"].asUInt();
	this->logFrequency = config["run"]["log-frequency"].asUInt();
	this->log_coordinates = config["run"]["log-coordinates"].asBool();
	this->log_concentrations = config["run"]["log-concentrations"].asBool();
	this->log_collisions = config["run"]["log-collisions"].asBool();
	this->log_reactions = config["run"]["log-reactions"].asBool();

	// probability mode
	this->probability_mode_enabled = config["run"]["probability-mode"]["enabled"].asBool();
	this->probability_mode_product = config["run"]["probability-mode"]["product"].asString();
	this->probability_mode_limit = config["run"]["probability-mode"]["limit"].asUInt();

	// loading medium settings

	this->viscosity = config["medium"]["viscosity"].asDouble();
	this->density = config["medium"]["density"].asDouble();
	this->temperature = config["medium"]["temperature"].asDouble();

	// loading geometry

	this->boxWidth = config["geometry"]["width"]["size"].asDouble();
	this->boxDepth = config["geometry"]["depth"]["size"].asDouble();
	this->boxHeight = config["geometry"]["height"]["size"].asDouble();

	this->boxWidthDiv = config["geometry"]["width"]["div"].asUInt();
	this->boxDepthDiv = config["geometry"]["depth"]["div"].asUInt();
	this->boxHeightDiv = config["geometry"]["height"]["div"].asUInt();

	this->subBoxWidth = this->boxWidth / (double)this->boxWidthDiv;
	this->subBoxDepth = this->boxDepth / (double)this->boxDepthDiv;
	this->subBoxHeight = this->boxHeight / (double)this->boxHeightDiv;

	this->bcX = config["geometry"]["width"]["periodic"].asBool();
	this->bcY = config["geometry"]["depth"]["periodic"].asBool();
	this->bcZ = config["geometry"]["height"]["periodic"].asBool();

	// loading micelle settings
	// if "micelle" in ... {}
	this->use_micelle = true;

	this->micelle_density = config["micelle"]["density"].asDouble();
	this->micelle_viscosity = config["micelle"]["viscosity"].asDouble();
	this->micelle_position = config["micelle"]["position"].asString();
	this->micelle_geometry_r = config["micelle"]["geometry"]["r"].asDouble();

	this->__init_micelle();

	// loading particle types
	for (unsigned int iPt = 0; iPt < config["particle-types"].size(); iPt++)
	{
		string particle_type_name = config["particle-types"][iPt]["type"].asString();
		cout << "Particle type: " << particle_type_name << endl;
		ParticleType* particleType = new ParticleType(config["particle-types"][iPt]["type"].asString());
		particleType->r = config["particle-types"][iPt]["r"].asDouble();

		if (this->use_micelle)
		{
			particleType->medium_bound = config["particle-types"][iPt]["medium-bound"].asBool();
			particleType->micelle_bound = config["particle-types"][iPt]["micelle-bound"].asBool();
			particleType->fixed_concentration = config["particle-types"][iPt]["fixed-concentration"].asBool();
			particleType->final_product = config["particle-types"][iPt]["final-product"].asBool();
		}
		//cout << config["particle-types"][iPt]["type"].asString() << endl;
		particleType->calcParams();
		this->addParticleType(particleType); 
	}

	// loading interactions
	for (unsigned int iInt = 0; iInt < config["interactions"].size(); iInt++)
	{
		Interaction* interaction = new Interaction(config["interactions"][iInt]["name"].asString());

		Json::Value particles = config["interactions"][iInt]["particles"];
		interaction->setParticleTypes(particles[0].asString(), particles[1].asString());
		interaction->setInteractionType(config["interactions"][iInt]["type"].asString());
		interaction->cutoff = config["interactions"][iInt]["cutoff"].asDouble();

		Json::Value params = config["interactions"][iInt]["params"];
		for (unsigned int iPrm = 0; iPrm < params.size(); iPrm++)
			interaction->addParam(params[iPrm].asDouble());

		this->addInteraction(interaction);
	}


	// loading reactions
	for (unsigned int iRxn = 0; iRxn < config["reactions"].size(); iRxn++)
	{
		Reaction* reaction = new Reaction(config["reactions"][iRxn]["name"].asString());

		reaction->probability = config["reactions"][iRxn]["probability"].asDouble();
		reaction->energy = config["reactions"][iRxn]["energy"].asDouble();
		reaction->cutoff = config["reactions"][iRxn]["cutoff"].asDouble();

		Json::Value reactants = config["reactions"][iRxn]["reactants"];
		for (unsigned int iRct = 0; iRct < reactants.size(); iRct++)
			reaction->addReactant(reactants[iRct].asString());

		Json::Value products = config["reactions"][iRxn]["products"];
		for (unsigned int iPrd = 0; iPrd < products.size(); iPrd++)
			reaction->addProduct(products[iPrd].asString());

		cout << "Reaction: " << config["reactions"][iRxn]["name"].asString() << endl;
		this->addReaction(reaction);
	}

	// allocate memory for particles coordinates
	this->rx = new double[this->maxParticles];
	this->ry = new double[this->maxParticles];
	this->rz = new double[this->maxParticles];
	this->particles = new Particle * [this->maxParticles];
	this->lastParticleIdx = -1;
	this->lastParticleId = -1;

	// total number of subboxes
	int subBoxCount = this->boxWidthDiv * this->boxHeightDiv * this->boxDepthDiv;

	// storage for indices of particles in each subbox
	for (int di = 0; di < subBoxCount; di++)
	{
		this->subBox.emplace_back(list<unsigned int>());
	}

	// adding particles
	for (unsigned int iPt = 0; iPt < config["particles"].size(); iPt++)
	{
		string type = config["particles"][iPt]["type"].asString();
		string position = config["particles"][iPt]["position"].asString();
		unsigned int amount = config["particles"][iPt]["amount"].asUInt();
		ParticleType* particleType = this->particleTypes[type];

		if (particleType->fixed_concentration)
			particleType->fixed_concentration_count = amount;

		for (unsigned int iA = 0; iA < amount; iA++)
		{
			Particle* particle = new Particle(particleType);

			if (position == "random" || position == "random_micelle" || position == "random_medium")
				this->addParticle(particle, position);   // true for random position
			else if (position == "center")
				this->addParticle(particle, "fixed", this->boxWidth / 2.0, this->boxHeight / 2.0, this->boxDepth / 2.0);
		}
	}

	cout << "Added total " << this->lastParticleIdx + 1 << " particles" << endl;
	cout << "configuration successfully loaded" << endl;
}

void System::test_system()
{
	cout << "========================================" << endl;
	
	for (pair<string, ParticleType*> el : this->particleTypes)
	{
		string particle_type_name = el.first;
		ParticleType* particle_type = el.second;
		map<ParticleType*, list<Reaction*>> collisions = particle_type->collisions;

		cout << "Particle type: " << particle_type_name << endl;
		cout << "  collisions: " << endl;
		
		for (pair<ParticleType*, list<Reaction*>> col : collisions)
			for (Reaction* reaction : col.second)
				cout << "    " << reaction->name << endl;
	}
	
	cout << "========================================" << endl;
}

unsigned int System::subBoxIndex(int _iW, int _iD, int _iH)
{
	unsigned int i = 0;
	i += _iW * this->boxDepthDiv * this->boxHeightDiv;
	i += _iD * this->boxHeightDiv;
	i += _iH;
	return i;
}

// subbox methods

void System::__addParticleToSubBox(Particle* _particle)
{
	// this->subBoxWidth = this->boxWidth / (double)this->boxWidthDiv;

	auto iW = (unsigned int)(_particle->x() / this->subBoxWidth) % this->boxWidthDiv;
	auto iD = (unsigned int)(_particle->y() / this->subBoxDepth) % this->boxDepthDiv;
	auto iH = (unsigned int)(_particle->z() / this->subBoxHeight) % this->boxHeightDiv;

	_particle->setSubBox(iW, iD, iH);

	int sub_box_id = _particle->subBoxId();

	this->subBox[sub_box_id].emplace_back(_particle->idx);
}

void System::__removeParticleFromSubBox(Particle* _particle)
{
	//cout << _particle->idx << " ";
	//cout << this->subBox[_particle->subBoxId()].size();
	this->subBox[_particle->subBoxId()].remove(_particle->idx);
	//cout << " " << this->subBox[_particle->subBoxId()].size() << endl;
}

void System::__updateParticleSubBox(Particle* _particle)
{
	auto iW = (unsigned int)(_particle->x() / this->subBoxWidth) % this->boxWidthDiv;
	auto iD = (unsigned int)(_particle->y() / this->subBoxDepth) % this->boxWidthDiv;
	auto iH = (unsigned int)(_particle->z() / this->subBoxHeight) % this->boxWidthDiv;

	auto oldId = _particle->subBoxId();
	auto newId = this->subBoxIndex(iW, iD, iH);

	if (oldId != newId)
	{
		this->__removeParticleFromSubBox(_particle);
		this->__addParticleToSubBox(_particle);
	}
}



void System::addParticle(Particle* _particle, string _position, double _x, double _y, double _z)
{
	// set particle id
	_particle->id = ++this->lastParticleId;

	// set particle index in the list
	_particle->idx = ++this->lastParticleIdx;
	this->particles[this->lastParticleIdx] = _particle;

	_particle->particleType->addExistingParticle(_particle);

	// set particle's coordinates

	//cout << "addParticle [1]: " << _particle->x << " " << _particle->y << " " << _particle->z << endl;
	if (_position == "random" || _position == "random_micelle" || _position == "random_medium")
	{
		double rnd;
		double x, y, z;
		while (true)
		{
			rnd = this->__random();
			x = this->boxWidth * rnd;
			rnd = this->__random();
			y = this->boxDepth * rnd;
			rnd = this->__random();
			z = this->boxHeight * rnd;

			// place particles appropriatly
			if (_position == "random" || 
				_position == "random_micelle" &&  this->__is_in_micelle(x, y, z) || 
				_position == "random_medium"  && !this->__is_in_micelle(x, y, z))
				break;
		}

		//if (_particle->particleType->type == "rO2*")
		//	cout << _particle->particleType->type << " " << _position << " " << x << " " << y << " " << z << " " << this->__is_in_micelle(x, y, z) << endl;

		_particle->setX(x);
		_particle->setY(y);
		_particle->setZ(z);
	}
	else if (_position == "fixed")
	{
		_particle->setX(_x);
		_particle->setY(_y);
		_particle->setZ(_z);
	}
	else
	{
		//???
	}

	this->__addParticleToSubBox(_particle);

	//cout << "addParticle: " << _particle->x() << " " << _particle->y() << " " << _particle->z() << endl;

	//_particle->particleType->addExistingParticle(_particle);  // add to ParticleType list
};

void System::removeParticle(Particle* _particle)
{
	//cout << "remove idx: " << _particle->idx << endl;
	_particle->particleType->removeExistingParticle(_particle);
	this->__removeParticleFromSubBox(_particle);

	moveParticle(this->lastParticleIdx, _particle->idx);
	//cout << "last idx: " << lastParticleIdx << endl;
	this->lastParticleIdx--;
	

	// destroy _particle?
	//_particle->~Particle();
	//delete _particle;

	// do we really need to delete/destroy it?
}

void System::moveParticle(long int _idxFrom, long int _idxTo)
{
	// when moving particles inside the list their idx is changing
	// therefore, one needs to put updated idx into particle's subbox

	int sbid = this->particles[_idxTo]->subBoxId();
	this->subBox[sbid].remove(_idxFrom);
	if (_idxFrom != _idxTo)
	{
		this->particles[_idxTo] = this->particles[_idxFrom];
		this->particles[_idxTo]->idx = _idxTo;
		this->rx[_idxTo] = this->rx[_idxFrom];
		this->ry[_idxTo] = this->ry[_idxFrom];
		this->rz[_idxTo] = this->rz[_idxFrom];
		
		this->subBox[sbid].emplace_back(_idxTo);
	}
}

// add Interaction
void System::addInteraction(Interaction* _interaction)
{
	// add to interactions list
	this->interactions[_interaction->name] = _interaction;

	// cross reference of interactions between particles 1 & 2
	_interaction->particleType1->addInteraction(_interaction->particleType2, _interaction);
	_interaction->particleType2->addInteraction(_interaction->particleType1, _interaction);

}

// add Reaction
void System::addReaction(Reaction* _reaction)
{
	// add to reactions list
	this->reactions[_reaction->name] = _reaction;

	// assign reaction to particle types
	if (_reaction->reactants.size() == 1)
	{
		_reaction->reactants[0]->addDecay(_reaction);
	}
	else if (_reaction->reactants.size() == 2)
	{
		_reaction->reactants[0]->addCollision(_reaction->reactants[1], _reaction);
		_reaction->reactants[1]->addCollision(_reaction->reactants[0], _reaction);
	}
	else
	{
		cout << "Error adding reaction " << _reaction->name << endl;
	}
};

// add ParticleType 
void System::addParticleType(ParticleType* _particleType)
{
	this->particleTypes[_particleType->type] = _particleType;
};



///////////////////////////////////////////////////////////////

double System::__random()
{
	static std::random_device rd;
	static std::mt19937 e2(rd());
	
	static unsigned long rnd_int;
	static double rnd;

	rnd_int = e2();
	rnd = (double)rnd_int / e2.max();

	//cout.precision(12);
	//cout << std::fixed << rnd << endl;
	return rnd;
}

double System::__randomBD(double _sigma2)
{
	return pow(abs(-2.0 * _sigma2 * log(this->__random() + 0.000000001)), 0.5) * cos(2.0 * this->PI * this->__random());
}

void System::__save_coordinates(string _fname)
{
	ofstream fOut(_fname, ofstream::out);
	string sep = ",";
	fOut << "particle-id" << sep 
		 << "particle-type" << sep 
		 << "x" << sep << "y" << sep << "z" << sep
		 << "box_w" << sep << "box_d" << sep << "box_h" << sep << "box_id" << sep
		 << "in-micelle" << sep
	     << "reacted" << endl;
	for (int i = 0; i <= this->lastParticleIdx; i++)
	{
		Particle* particle = this->particles[i];
		fOut << particle->id << sep << particle->particleType->type << sep 
			 << particle->x() << sep << particle->y() << sep << particle->z() << sep
  			 << particle->subBoxW() << sep << particle->subBoxD() << sep << particle->subBoxH() << sep << particle->subBoxId() << sep
			 << int(this->__is_in_micelle(particle->x(), particle->y(), particle->z())) << sep
 			 << int(particle->reacted)
			 << endl;
	}

	fOut.close();
}

void System::__save_concentrations(string _thread_id, double _t)
{
	static bool firstLine = true;
	if (firstLine)
	{
		ofstream f_out_all(this->working_path + "log/" + _thread_id + "/conc.csv", ofstream::out);
		ofstream f_out_fixed(this->working_path + "log/" + _thread_id + "/conc_fixed.csv", ofstream::out);
		ofstream f_out_micelle(this->working_path + "log/" + _thread_id + "/conc_micelle.csv", ofstream::out);
		ofstream f_out_medium(this->working_path + "log/" + _thread_id + "/conc_medium.csv", ofstream::out);
		firstLine = false;
	
		// print header
		f_out_all << "Time";
		f_out_fixed << "Time";
		f_out_micelle << "Time";
		f_out_medium << "Time";
		
		for (auto particleTypeItem : this->particleTypes)
		{
			f_out_all << "," << particleTypeItem.first;
			f_out_fixed << "," << particleTypeItem.first;
			f_out_micelle << "," << particleTypeItem.first;
			f_out_medium << "," << particleTypeItem.first;
		}
		
		f_out_all << endl;
		f_out_fixed << endl;
		f_out_micelle << endl;
		f_out_medium << endl;
		
		f_out_all.close();
		f_out_fixed.close();
		f_out_micelle.close();
		f_out_medium.close();
	}

	ofstream f_out_all(this->working_path + "log/" + _thread_id + "/conc.csv", ofstream::app);
	ofstream f_out_fixed(this->working_path + "log/" + _thread_id + "/conc_fixed.csv", ofstream::app);
	ofstream f_out_micelle(this->working_path + "log/" + _thread_id + "/conc_micelle.csv", ofstream::app);
	ofstream f_out_medium(this->working_path + "log/" + _thread_id + "/conc_medium.csv", ofstream::app);

	f_out_all << _t;
	f_out_fixed << _t;
	f_out_micelle << _t;
	f_out_medium << _t;

	for (auto particleTypeItem : this->particleTypes)
	{
		auto exps = particleTypeItem.second->existingParticles;
		int particles_all = exps.size();
		int particles_fixed = particleTypeItem.second->fixed_added;
		int particles_micelle = 0;
		int particles_medium = 0;

		for (auto p : exps)
		{
			if (this->__is_in_micelle(p->x(), p->y(), p->z())) 
				particles_micelle++;
			else 
				particles_medium++;
		}
		f_out_all << "," << particles_all;
		f_out_fixed << "," << particles_fixed;
		f_out_micelle << "," << particles_micelle;
		f_out_medium << "," << particles_medium;
	}

	f_out_all << endl;
	f_out_fixed << endl;
	f_out_micelle << endl;
	f_out_medium << endl;

	f_out_all.close();
	f_out_fixed.close();
	f_out_micelle.close();
	f_out_medium.close();
}

void System::__save_collisions(string _thread_id, double _t)
{
	static bool firstLine = true;
	if (firstLine)
	{
		ofstream f_out_collisions(this->working_path + "log/" + _thread_id + "/collisions.csv", ofstream::out);
		firstLine = false;

		// print header
		f_out_collisions << "Time";

		for (auto reactionItem : this->reactions)
		{
			f_out_collisions << "," << reactionItem.first;
		}

		f_out_collisions << endl;
		f_out_collisions.close();
	}

	ofstream f_out_collisions(this->working_path + "log/" + _thread_id + "/collisions.csv", ofstream::app);

	f_out_collisions << _t;

	for (auto reactionItem : this->reactions)
	{
		f_out_collisions << "," << reactionItem.second->collision_count;
	}

	f_out_collisions << endl;
	f_out_collisions.close();
}

void System::__save_reactions(string _thread_id, double _t)
{
	static bool firstLine = true;
	if (firstLine)
	{
		ofstream f_out_reactions(this->working_path + "log/" + _thread_id + "/reactions.csv", ofstream::out);
		firstLine = false;

		// print header
		f_out_reactions << "Time";

		for (auto reactionItem : this->reactions)
		{
			f_out_reactions << "," << reactionItem.first;
		}

		f_out_reactions << endl;
		f_out_reactions.close();
	}

	ofstream f_out_reactions(this->working_path + "log/" + _thread_id + "/reactions.csv", ofstream::app);

	f_out_reactions << _t;

	for (auto reactionItem : this->reactions)
	{
		f_out_reactions << "," << reactionItem.second->reaction_count;
	}

	f_out_reactions << endl;
	f_out_reactions.close();
}


// not used, for debugging purposes...
void System::__printDistance(Particle* _particle1, Particle* _particle2)
{
	double dx = (_particle1->x() - _particle2->x());
	double dy = (_particle1->y() - _particle2->y());
	double dz = (_particle1->z() - _particle2->z());

	double r = pow(dx * dx + dy * dy + dz * dz, 0.5);
	//cout << _particle1->id << " " << _particle2->id << endl;
	cout << "r=" << r << "; dx=" << dx << "; dy=" << dy << "; dz=" << dz << endl;
}

// 
//double System::__

// ==================================
// main routine
// ==================================

void System::__debug_box()
{
	int cnt = 0, cnt_final = 0;
	for (auto u : this->subBox[13])
	{
		auto pt = this->particles[u];
		std::cout << u << " " << pt->id << " " << pt->idx << " " << pt->particleType->type << " " << pt->particleType->final_product << endl;

		if (!this->particles[u]->particleType->final_product)
			cnt++;
		else
			cnt_final++;
	}
	cout << cnt << " " << cnt_final << endl;
}

void System::run()
{
	// temporary storage for coordinates/forces/energies/etc...
	this->__dx = new double[this->maxParticles];
	this->__dy = new double[this->maxParticles];
	this->__dz = new double[this->maxParticles];
	this->__dr = new double[this->maxParticles];
	this->__dr2 = new double[this->maxParticles];

	double dFx, dFy, dFz, E;

	double* Fx = new double[this->maxParticles];
	double* Fy = new double[this->maxParticles];
	double* Fz = new double[this->maxParticles];


	// timestep cycle
	for (unsigned long int step = 0; step < this->steps; step++)
	{
		// stop simulation when reached product's limit in probability mode
		if (this->probability_mode_enabled && this->particleTypes[this->probability_mode_product]->existingParticles.size() >= this->probability_mode_limit)
			break;

		// logging
		if (step % this->logFrequency == 0)
		{
			auto padded = std::to_string(step + 1);
			padded.insert(0, 14U - std::min(std::string::size_type(14), padded.length()), '0');

			cout << "Step: " << step << " / " << this->steps << endl;
			if (this->log_coordinates)
				this->__save_coordinates(this->working_path + "log/" + this->threadId + "/coords/" + padded + ".csv");

			if (this->log_concentrations)
				this->__save_concentrations(this->threadId, step * this->dt);

			if (this->log_collisions)
				this->__save_collisions(this->threadId, step * this->dt);

			if (this->log_reactions)
				this->__save_reactions(this->threadId, step * this->dt);
		}

		// init force values
		for (unsigned int f = 0; f < this->maxParticles; f++)
		{
			Fx[f] = 0.0; Fy[f] = 0.0; Fz[f] = 0.0;
		}

		vector<tuple<Particle*, Particle*, Reaction*>> possibleReactions;    // possible reactions during one cycle
		list<Particle*> reactedParticles;                                    // reacted particles, to remove
		list<tuple<Particle*, double, double, double>> createdParticles;     // created particles, to add

		// reaction cycle

		//cout << "1) begin cycle" << endl;
		//this->__debug_box();

		// iterate through all subboxes
		for (unsigned int iW = 0; iW < this->boxWidthDiv; iW++)
			for (unsigned int iD = 0; iD < this->boxDepthDiv; iD++)
				for (unsigned int iH = 0; iH < this->boxHeightDiv; iH++)
				{
					auto sbi1 = this->subBoxIndex(iW, iD, iH);

					unsigned int dbW_ = 1, dbD_ = 1, dbH_ = 1;

					// accounting for all particles in the proximity
					vector<pair<unsigned int, bool>> nearbyParticles;
					for (unsigned int jW_ = iW; jW_ <= iW + dbW_; jW_++)
						for (unsigned int jD_ = iD; jD_ <= iD + dbD_; jD_++)
							for (unsigned int jH_ = iH; jH_ <= iH + dbH_; jH_++)
							{
								auto jW = jW_ % this->boxWidthDiv;
								auto jD = jD_ % this->boxDepthDiv;
								auto jH = jH_ % this->boxHeightDiv;

								auto sbi2 = this->subBoxIndex(jW, jD, jH);

								// particles from this box (sbi1) should not be used twice
								auto isThisBox = (jW == iW && jD == iD && jH == iH);

								for (auto j : this->subBox[sbi2])
									if (!this->particles[j]->particleType->final_product)
										nearbyParticles.emplace_back(make_pair(j, isThisBox));
							}

					

					for (auto i : this->subBox[sbi1])
					{
						if (!this->particles[i]->particleType->final_product)
						{
							Particle* particle1 = this->particles[i];
							auto interactions = particle1->particleType->interactions;

							// !!! this step could be automated with CUDA
							// calculating dx/dy/dz/dr2/dr

							// deacys due to interaction with medium
							if (step >= this->tsteps)
							{
								auto decays = particle1->particleType->decays;
								// each reaction is allowed to have multiple channels (A->C, A->E+F, ...)
								for (Reaction* reaction : decays)
								{
									if (this->__random() < reaction->probability)
										possibleReactions.emplace_back(make_tuple(particle1, nullptr, reaction));
								}
							}


							if (interactions.size() > 0)
							{
								// fill position data: dr, dr^2, ...
								// CUDA?
								for (auto particleItem : nearbyParticles)
								{
									unsigned int j = particleItem.first;
									bool isThisBox = particleItem.second;

									if (!isThisBox || j > i)
									{
										this->__dx[j] = this->rx[i] - this->rx[j];
										this->__dy[j] = this->ry[i] - this->ry[j];
										this->__dz[j] = this->rz[i] - this->rz[j];

										// accounting for periodic boundary conditions along each axis
										if (this->bcX)
											this->__dx[j] -= rint(this->__dx[j] / this->boxWidth) * this->boxWidth;
										if (this->bcY)
											this->__dy[j] -= rint(this->__dy[j] / this->boxDepth) * this->boxDepth;
										if (this->bcZ)
											this->__dz[j] -= rint(this->__dz[j] / this->boxHeight) * this->boxHeight;

										this->__dr2[j] = this->__dx[j] * this->__dx[j] + this->__dy[j] * this->__dy[j] +
											this->__dz[j] * this->__dz[j];
										this->__dr[j] = sqrt(this->__dr2[j]);
									}
								}

								// interactions with medium
								// ...

						

								// interactions & reactions (cycle over neraby boxes)
								for (auto particleItem : nearbyParticles)
								{
									unsigned int j = particleItem.first;
									bool isThisBox = particleItem.second;

									if (!isThisBox || j > i)
									{
										Particle* particle2 = this->particles[j];

										// interaction
										if (interactions.find(particle2->particleType) != interactions.end())
										{
											auto interaction = particle1->particleType->interactions[particle2->particleType];

											auto interactionFunction = interaction->interactionFunction;
											interactionFunction(particle1, particle2, interaction, dFx, dFy, dFz, E);

											// update total force acting on i-th and j-th particles due to interaction
											Fx[i] += dFx;
											Fy[i] += dFy;
											Fz[i] += dFz;
											Fx[j] -= dFx;
											Fy[j] -= dFy;
											Fz[j] -= dFz;

											// reactions are possible only after initial thermalization
											if (step >= this->tsteps)
											{
												// reactions: collisions
												auto collisions = particle1->particleType->collisions;
												if (collisions.find(particle2->particleType) != collisions.end())
													// each reaction is allowed to have multiple channels (A+B->C+D, A+B->E, ...)
													for (Reaction* reaction : collisions[particle2->particleType])
													{
														//if (reaction->name == "rO2* + RCH2CH=RH -> RCH2CH*=RHOOr" && E > 0.0001)
														//	cout << "E=" << E << endl;
														if ( 
															(reaction->energy == 0.0 && this->__dr[particle2->idx] < reaction->cutoff) ||
															(reaction->energy > 0 && E > reaction->energy) ||
															(reaction->energy < 0 && E < reaction->energy)
														   )
															if (this->__is_in_micelle(particle1->x(), particle1->y(), particle1->z()) == this->__is_in_micelle(particle2->x(), particle2->y(), particle2->z()))
															{
																reaction->collision_count++;
																if (this->__random() < reaction->probability)
																	possibleReactions.emplace_back(make_tuple(particle1, particle2, reaction));
															}
													}

												// reactions: decays due to medium
												// [?] what happens with decays due to interaction with medium?
											}
										}
									}
								}
							}
						}
					}
				}

		// all interactions/reactions are done...

		// shuffle reactions so that there is no bias towards particles with lower ids
		auto rng = std::default_random_engine{};
		std::shuffle(std::begin(possibleReactions), std::end(possibleReactions), rng);

		//cout << "Possible reactions: " << possibleReactions.size() << endl;

		// perform reactions
		for (auto reactionItem : possibleReactions)
		{
			Particle* particle1 = get<0>(reactionItem);
			Particle* particle2 = get<1>(reactionItem);
			Reaction* reaction = get<2>(reactionItem);

			bool rc;
			if (particle2 == nullptr)
				rc = !particle1->reacted;
			else
				rc = !particle1->reacted && !particle2->reacted;

			if (rc)
			{
				// average coordinates where reaction took place
				double  ax = this->rx[particle1->idx],
						ay = this->ry[particle1->idx],
						az = this->rz[particle1->idx];

				particle1->reacted = true;
				reactedParticles.emplace_back(particle1);

				reaction->reaction_count++;

				if (particle2 != nullptr)
				{
					particle2->reacted = true;
					reactedParticles.emplace_back(particle2);

					ax = (ax + this->rx[particle2->idx]) / 2 - floor(abs(ax - this->rx[particle2->idx]) / (this->boxWidth / 2)) * this->boxWidth / 2;
					ay = (ay + this->ry[particle2->idx]) / 2 - floor(abs(ay - this->ry[particle2->idx]) / (this->boxDepth / 2)) * this->boxDepth / 2;
					az = (az + this->rz[particle2->idx]) / 2 - floor(abs(az - this->rz[particle2->idx]) / (this->boxHeight / 2)) * this->boxHeight / 2;
				}

				// adding products
				double dr = 2.0, dx, dy, dz;
				if (particle2 != nullptr)
				{
					dx = particle1->x() - particle2->x();
					dy = particle1->y() - particle2->y();
					dz = particle1->z() - particle2->z();

					// accounting for periodic boundary conditions
					if (this->bcX)
						dx -= rint(dx / this->boxWidth) * this->boxWidth;
					if (this->bcY)
						dy -= rint(dy / this->boxDepth) * this->boxDepth;
					if (this->bcZ)
						dz -= rint(dz / this->boxHeight) * this->boxHeight;

					dr = sqrt(dx*dx + dy*dy + dz*dz);
				}

				for (ParticleType* productParticleType : reaction->products)
				{
					if (productParticleType->type == "<void>") continue;
					Particle* product = new Particle(productParticleType);

					// creating new particles nearby the average point
					// in the range [2, 8] * radius for each coordinate
					// [?] rethink...

					//here we create both particles either in micelle or outside depending on their initial positions

					double pax, pay, paz;
					int try_count = 0;
					bool pim_1, pim_2, pim_product;

					while (true)
					{
						pax = ax + (0.5 - this->__random()) * dr;
						pay = ay + (0.5 - this->__random()) * dr;
						paz = az + (0.5 - this->__random()) * dr;

						pax -= floor(pax / this->boxWidth) * this->boxWidth;
						pay -= floor(pay / this->boxDepth) * this->boxDepth;
						paz -= floor(paz / this->boxHeight) * this->boxHeight;

						pim_1 = this->__is_in_micelle(particle1->x(), particle1->y(), particle1->z());
						
						if (particle2 != nullptr)
							pim_2 = this->__is_in_micelle(particle2->x(), particle2->y(), particle2->z());

						pim_product = this->__is_in_micelle(pax, pay, paz);
						
						if (pim_1 == pim_product) break;

						try_count++;
						if (try_count > 100)
						{
							cout << "e1    " << pim_1 << " " << pim_product << " " << endl;
							if (particle2 != nullptr)
								cout << "e2    " << pim_2 << " " << pim_product << " " << endl;
							cout << "pa    " << pax << " " << pay << " " << paz << endl;
							cout << "a     " <<  ax << " " <<  ay << " " <<  az << endl;
							cout << "1: " << particle1->particleType->type << " " << particle1->x() << " " << particle1->y() << " " << particle1->z() << endl;
							if (particle2 != nullptr)
								cout << "2: " << particle2->particleType->type << " " << particle2->x() << " " << particle2->y() << " " << particle2->z() << endl;
						}
					}
					if (try_count > 100) cout << "bad_rand " << reaction->name << endl;

					createdParticles.emplace_back(make_tuple(product, pax, pay, paz));
				}
			}
		}

		//cout << "2) before brownian" << endl;
		//this->__debug_box();

		///////////////////////////////////////
		// do Brownian & force motion
		//////////////////////////////////////////

		double upd_x, upd_y, upd_z;

		for (int i = 0; i <= this->lastParticleIdx; i++)
		{
			Particle* particle = this->particles[i];
			//if (step > 900)
			//   cout << ">>> " << particle->id << " " << particle->idx << " " << particle->particleType->type << endl;
			if (!particle->reacted && !this->particles[i]->particleType->final_product)
			{
				double xi;
				bool pim = this->__is_in_micelle(particle->x(), particle->y(), particle->z());

				// use different viscosity values for micelle/medium
				if (pim)
					xi = particle->particleType->xi_micelle;
				else
					xi = particle->particleType->xi_medium;

				double sigma2 = (1E+20) * 2 * this->KB * this->temperature * (this->dt * 1E-12) / xi; // A^2
				// !!! redundant calculations
				// should be calculated once and saved as particle->particleType->sigma
				// 1E-12 - ps
				// 1E+10 - A

				int try_count = 0;
				int try_count_max = 10;
				while (true)
				{
					double ux = this->dt * 1E-12 * Fx[i] * (1.6E-19 * 1E+20) / xi;
					double uy = this->dt * 1E-12 * Fy[i] * (1.6E-19 * 1E+20) / xi;
					double uz = this->dt * 1E-12 * Fz[i] * (1.6E-19 * 1E+20) / xi;

					//if (!particle->particleType->micelle_bound)
					//cout << "0) " << xi << " " << sigma2 << " " << sqrt(sigma2) << endl;
					//if (!particle->particleType->micelle_bound)
					//cout << "1) " << particle->id << " " << ux << " " << uy << " " << uz << endl;

					// setting threshold for displacement in any direction
					ux = copysign(min(abs(ux), 3.0), ux);
					uy = copysign(min(abs(uy), 3.0), uy);
					uz = copysign(min(abs(uz), 3.0), uz);

					ux += this->__randomBD(sigma2);
					uy += this->__randomBD(sigma2);
					uz += this->__randomBD(sigma2);

					//if (!particle->particleType->micelle_bound)
					//cout << "2) " << particle->id << " " << ux << " " << uy << " " << uz << endl;

					// diffusion inside micelle 
					if (pim && particle->particleType->micelle_bound)
					{
						ux = fmod(ux, this->micelle_geometry_r / 2);
						uy = fmod(uy, this->micelle_geometry_r / 2);
						uz = fmod(uz, this->micelle_geometry_r / 2);
					}
					//if (!particle->particleType->micelle_bound)
					//	cout << "3) " << particle->id << " " << ux << " " << uy << " " << uz << endl;

					upd_x = particle->x() + ux;
					upd_y = particle->y() + uy;
					upd_z = particle->z() + uz;

					//cout << "2) " << particle->id << " " << ux << " " << uy << " " << uz << endl;

					bool pim_upd = this->__is_in_micelle(upd_x, upd_y, upd_z);

					if (!particle->particleType->micelle_bound) 
						// non-bound particles are not affected
						break;
					else
						if (!(pim && !pim_upd)) 
							// bound particles can't leave micelle, other displacements are allowed
							break;
					try_count++;
					if (try_count > try_count_max) break;
				}

				if (try_count <= try_count_max)
				{
					particle->setX(upd_x);
					particle->setY(upd_y);
					particle->setZ(upd_z);
				}
				// should also be optimized?

				this->__updateParticleSubBox(particle);
			}
		}

		//cout << "3) after brownian" << endl;
		//this->__debug_box();


		// remove all reacted particles
		for (Particle* particle : reactedParticles)
		{
			//cout << "Removed: " << particle->particleType->type << " " << particle->id <<  " " << particle->idx << endl;
			this->removeParticle(particle);
		}

		//cout << "4) after remove" << endl;
		//this->__debug_box();

		// add all created particles
		for (tuple<Particle*, double, double, double>& particleItem : createdParticles)
		{
			//cout << "Added: " << get<0>(particleItem)->particleType->type << " " << get<0>(particleItem)->id << endl;
			this->addParticle(get<0>(particleItem), "fixed", get<1>(particleItem), get<2>(particleItem), get<3>(particleItem));
		}

		//cout << "5) after add" << endl;
		//this->__debug_box();


		// regenerate fixed concentration

		for (pair<string, ParticleType*> element : this->particleTypes)
		{
			ParticleType* particleType = element.second;
			if (particleType->fixed_concentration)
			{
				int ex = particleType->existingParticles.size();
				int fc = particleType->fixed_concentration_count;
				if (ex < fc)
				{
					for (int k = 0; k < fc - ex; k++)
					{
						Particle* particle = new Particle(particleType);
						if (particleType->micelle_bound)
							this->addParticle(particle, "random_micelle"); 
						else if (particleType->medium_bound)
							this->addParticle(particle, "random_medium");
						else 
							this->addParticle(particle, "random");
						
						particleType->fixed_added++;   // counting how many particles were added to maintain constant concentrations
					}
				}
			}

		}

		//cout << "6) after reg" << endl;
		//this->__debug_box();


	}
}

