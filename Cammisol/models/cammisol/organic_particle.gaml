/**
* Name: organic_particle
* Based on the internal empty template. 
* Author: pbreugno
* Tags: 
*/


model organic_particle

import "local_time.gaml"

/* Insert your model definition here */

global {
	
}

species OrganicParticle {
	int grid_x;
	int grid_y;
	float C_labile;
	float N_labile;
	float P_labile;
	float C_recalcitrant;
	float N_recalcitrant;
	float P_recalcitrant;
	
	bool in_pore <- false;
	
	list<OrganicParticle> organic_neighbors;
}
