/**
* Name: grid
* Based on the internal skeleton template. 
* Author: pbreugno
* Tags: 
*/

model cammisol_grid

import "nematode/nematode.gaml"
import "mineral_particle.gaml"

global {
	string ORGANIC <- "organic";
	string MINERAL <- "mineral";
	string PORE <- "pore";

	float organic_rate <- 1/3;
	float mineral_rate <- 1/3;
	float pore_rate <- 1 - organic_rate - mineral_rate;

	int grid_size <- 3;
	
	float soil_surface <- 1#cm;
	
	geometry shape <- square(soil_surface);
	float model_area <- shape.area;
	float cell_size <- soil_surface/grid_size;
	float cell_area <- cell_size^2;
	
	action init_grid {
		ask Particle {
			type <- rnd_choice([ORGANIC::organic_rate, MINERAL::mineral_rate, PORE::pore_rate]);
			do init_particle;
		}
		// Counts the number of PORES after the initialization of the grid
		int pores_count <- length(PoreParticle);
		ask PoreParticle {
			add self to: pore_neighbors;
			ask Particle[grid_x, grid_y] {
				loop neighbor_cell over: neighbors {
					if(neighbor_cell.type = PORE) {
						add PoreParticle(neighbor_cell.particle) to: myself.pore_neighbors;
					}
				}
			}
			
			add self.organic_particle to: accessible_organics;
			ask Particle[grid_x, grid_y] {
				loop neighbor_cell over: neighbors {
					if(neighbor_cell.type = ORGANIC) { // TODO: and not from neighbor pores?
						add OrganicParticle(neighbor_cell.particle) to: myself.accessible_organics;
					}
				}
			}
		}
	}
}

grid Particle width: grid_size height: grid_size neighbors: 4 {
	string type;
	agent particle;
	
	action init_particle {
		switch type {
			match MINERAL 
			{
				color <- #yellow;
				create MineralParticle with: (
					grid_x: self.grid_x,
					grid_y: self.grid_y,
					location: self.location,
					shape: self.shape
				) {
					myself.particle <- self;
				}
			}
			match ORGANIC {
				color <- #green;
				create OrganicParticle with: (
					grid_x: self.grid_x,
					grid_y: self.grid_y,
					location: self.location,
					shape: self.shape
				) {
					myself.particle <- self;
				}
			}
			match PORE { 
				color <- #black;
				// TODO: autant de nutriment dans la partie organique du pore que dans une vraie particule organique?
				OrganicParticle pore_organic_particle;
				Dam pore_dam;

				create OrganicParticle with: (
					C_labile: 0,
					N_labile: 0,
					P_labile: 0,
					C_recalcitrant: 0,
					N_recalcitrant: 0,
					P_recalcitrant: 0,
					location: self.location,
					grid_x: self.grid_x,
					grid_y: self.grid_y,
					in_pore: true,
					location: self.location,
					shape: self.shape
				) {
					pore_organic_particle <- self;
				}
				
				create Dam {
					pore_dam <- self;
				}
				create PoreParticle with: (
					organic_particle: pore_organic_particle,
					dam: pore_dam,
					grid_x: self.grid_x,
					grid_y: self.grid_y,
					location: self.location,
					shape: self.shape
				) {
					myself.particle <- self;
				}
			}
			
		}
	}
}

experiment grid type: gui {
	/** Insert here the definition of the input and output of the model */
	output {
	}
}
