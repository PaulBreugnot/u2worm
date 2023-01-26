/**
* Name: organic_particle
* Based on the internal empty template. 
* Author: pbreugno
* Tags: 
*/


model organic_particle

/* Insert your model definition here */

global {
	
}

species OrganicParticle {
	int grid_x;
	int grid_y;
	float C_labile;
	float C_recalcitrant;
	float N;
	float P;
	bool in_pore <- false;
	
	list<OrganicParticle> organic_neighbors;
}
