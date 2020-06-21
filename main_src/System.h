//
// Created by root on 11.02.2019.
//

#ifndef MOLECULARDYNAMICS_SYSTEM_H
#define MOLECULARDYNAMICS_SYSTEM_H


#include <list>
#include <string>
#include <vector>
#include <map>
#include "Particle.h"
#include "ParticleType.h"
#include "Reaction.h"
#include "Interaction.h"


using namespace std;

class System
{
public:
	// singleton definitions
	static System* getInstance();

	void addParticle(Particle* _particle, string _position, double _x = 0.0, double _y = 0.0, double _z = 0.0);
	void removeParticle(Particle* _particle);
	void moveParticle(long int _idxFrom, long int _idxTo);
	void addParticleType(ParticleType* _particleType);

	void addReaction(Reaction* _reaction);
	void addInteraction(Interaction* _interaction);

	void run();
	void load_config();
	void test_system();

	void setThreadId(string _threadId);

	map<string, Reaction*> reactions = {};
	map<string, Interaction*> interactions = {};
	map<string, ParticleType*> particleTypes = {};

	// Particle storage
	Particle** particles;
	double* rx, * ry, * rz;
	vector<list<unsigned int>> subBox;
	long int lastParticleIdx = -1;
	long int lastParticleId = -1;

	// probability mode
	bool probability_mode_enabled = false;
	string probability_mode_product = "";
	int probability_mode_limit = 0;

	// intermediate param storage
	double* __dx, * __dy, * __dz, * __dr, * __dr2;

	double boxWidth = 10.0, boxDepth = 10.0, boxHeight = 10.0;        // x, y, z
	double subBoxWidth = 10.0, subBoxDepth = 10.0, subBoxHeight = 10.0;
	unsigned int boxWidthDiv = 1, boxHeightDiv = 1, boxDepthDiv = 1;
	double temperature = 300.0, density = 1.0, viscosity = 8e-4;
	double dt = 1;

	// micelle params
	bool use_micelle = false;
	double micelle_density = 1.0;
	double micelle_viscosity = 8E-4;
	string micelle_position = "center";
	double micelle_geometry_r = 4.0;
	double micelle_geometry_x = 0.0;
	double micelle_geometry_y = 0.0;
	double micelle_geometry_z = 0.0;


	// run settings
	unsigned long tsteps = 100, steps = 1000, logFrequency = 10000;
	bool log_coordinates = false, log_concentrations = true, log_collisions = false, log_reactions = false;
	unsigned int maxParticles = 10000;
	
	string working_path = "";
	string threadId = "";
	string cfg_file = "config.json";


	bool bcX = true, bcY = true, bcZ = true;

	// constants
	double PI = 3.1415926, KB = 1.38064852e-23, NA = 6.02214086e-23;

	unsigned int subBoxIndex(int _iW, int _iD, int _iH);

private:

	// singleton definitions
	static System* __instance;
	System();

	void __debug_box();

	// methods
	double __random();
	double __randomBD(double _sigma);

	void __save_coordinates(string _fname);
	void __save_concentrations(string _thread_id, double _t);
	void __save_collisions(string _thread_id, double _t);
	void __save_reactions(string _thread_id, double _t);

	void __addParticleToSubBox(Particle* _particle);
	void __removeParticleFromSubBox(Particle* _particle);
	void __updateParticleSubBox(Particle* _particle);

	void __printDistance(Particle* _particle1, Particle* _particle2);

	double __interactionEnergy(Particle* _particle1, Particle* _particle2);

	// micelle
	void __init_micelle();
	bool __is_in_micelle(double _x, double _y, double _z);

	unsigned long __randomSequenceLength, __randomSequenceIndex;
	vector<double> __randomSequence;

};



#endif //MOLECULARDYNAMICS_SYSTEM_H
