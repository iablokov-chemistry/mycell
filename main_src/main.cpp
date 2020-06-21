#include <iostream>
#include <vector>
#include <fstream>
#include <string>

#include "Particle.h"
#include "System.h"


using namespace std;


int main(int argc, char* argv[])
{
	System* system = System::getInstance();
	if (argc > 1)
	{
		system->setThreadId((string)argv[1]);
		system->cfg_file = (string)argv[2];
	}
	else
	{
		system->setThreadId("th5");
		system->cfg_file = "config.json";
	}
	
	if (false)
		system->working_path = "/study/univer/chemistry/masters/projects/MolecularDynamics.VS/MD/x64/Debug/";
	else
	{
		system->working_path = "./";
	}

	system->load_config();
	system->test_system();
	system->run();
	return 0;
}