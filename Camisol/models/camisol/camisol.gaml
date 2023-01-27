/**
* Name: camisol
* Based on the internal empty template. 
* Author: pbreugno
* Tags: 
*/


model camisol

import "nematode.gaml"
import "mineral_particle.gaml"

global {
	string ORGANIC <- "organic";
	string MINERAL <- "mineral";
	string PORE <- "pore";

	float organic_rate <- 1/3;
	float mineral_rate <- 1/3;
	float pore_rate <- 1 - organic_rate - mineral_rate;

	int grid_width <- 30;
	int grid_height <- 30;
	
	int nematodes_count <- 30;
	
	geometry shape <- square(1#cm);
	float model_area <- shape.area;

	/**
	 * Solid matter weight by volume of soil.
	 */
	float bulk_density <- 1.17#gram/(#cm*#cm*#cm);
	/**
	 * The total weight of matter in the model.
	 */
	float total_model_weight <- bulk_density * world.shape.area * 1#cm;
	
	float carbone_concentration_in_dam <- (729.0#gram * 10^-6)/#gram;
	float azote_concentration_in_dam <- (60.0#gram * 10^-6)/#gram;
	float azote_concentration_in_dim <- (4.74#gram * 10^-6)/#gram;
	float phosphore_concentration_in_dam <- (400.0#gram * 10^-6)/#gram;
	float phosphore_concentration_in_dim <- (1.43#gram * 10^-6)/#gram;
	
	map<string,float> soil_characteristics <-   map(["C"::0.02157#gram/(#cm*#cm),"N"::0.00132#gram/(#cm*#cm),"P"::0.00077#gram/(#cm*#cm)]);
	
	init {
		step <- 1#h;
		// Counts the number of PORES after the initialization of the grid
		int pores_count <- length(PoreParticle);
		ask PoreParticle {
			// The carrying capacity of each pore is equal to 10 times the initial bacteria population
			carrying_capacity <- 10 * total_initial_bacteria_weight / pores_count;
			
			ask dam {
				/*Jean*/
				float carbon_in_all_pore <- total_model_weight * carbone_concentration_in_dam; 
				float carbon_in_pore <- carbon_in_all_pore / pores_count;
				
				float azote_in_all_pore <- total_model_weight * azote_concentration_in_dam;
				float azote_in_pore <- azote_in_all_pore / pores_count;
				
				float phosphore_in_all_pore <- total_model_weight * phosphore_concentration_in_dam;
				float phosphore_in_pore <- phosphore_in_all_pore / pores_count;
				/*Jean*/
		
				/*Bernard Et Al 2021 -> Bouzac 2015*/
				float phosphore_in_dim <- phosphore_concentration_in_dim * total_model_weight;
				float phosphore_in_dim_in_pore <- phosphore_in_dim / pores_count;
				
				float azote_in_dim <- azote_concentration_in_dim * total_model_weight;
				float azote_in_dim_in_pore <- phosphore_in_dim / pores_count;
				/*Bernard Et Al 2021 -> Bouzac 2015*/
				
				dom <- [azote_in_pore, phosphore_in_pore, carbon_in_pore];
				dim <- [azote_in_dim_in_pore, phosphore_in_dim_in_pore];
			}
			
			create Copiotrophe_R with: (dam: self.dam) {
				add self to:myself.populations;
			}
			create Copiotrophe_K with: (dam: self.dam) {
				add self to:myself.populations;
			}
			create Oligotrophe_K with: (dam: self.dam) {
				add self to:myself.populations;
			}
			
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
		create Nematode number: nematodes_count
		{
			current_pore <- one_of(PoreParticle); 
			location <- current_pore.location;
		}
	}
}

grid Particle width: grid_width height: grid_height neighbors: 4 {
	string type <- rnd_choice([ORGANIC::organic_rate, MINERAL::mineral_rate, PORE::pore_rate]);
	agent particle;
	
	init {
		float cell_area <- self.shape.area;
		switch type {
			match MINERAL 
			{
				color <- #yellow;
				create MineralParticle with: (
					N: (soil_characteristics["N"])*cell_area,
					P: (soil_characteristics["P"])*cell_area,
					grid_x: self.grid_x,
					grid_y: self.grid_y
				) {
					myself.particle <- self;
				}
			}
			match ORGANIC {
				color <- #green;
				create OrganicParticle with: (
					N: (soil_characteristics["N"])*cell_area,
					P: (soil_characteristics["P"])*cell_area,
					C_labile: (soil_characteristics["C"]/2)*cell_area, // TODO: variable factor
					C_recalcitrant: (soil_characteristics["C"]/2)*cell_area,
					grid_x: self.grid_x,
					grid_y: self.grid_y
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
					N: 0,
					P: 0,
					C_labile: 0,
					C_recalcitrant: 0,
					location: self.location,
					grid_x: self.grid_x,
					grid_y: self.grid_y,
					in_pore: true
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
					grid_y: self.grid_y
				) {
					myself.particle <- self;
				}
			}
			
		}
	}
}

experiment display_grid {
	output {
		display grid {
			grid Particle;
		}
	}
}

experiment camisol_no_output {
	reflex state when: cycle mod 100 = 0 {
				ask simulation {
			write "Time: " + time/#day;
			write "Dom: N=" + sum(Dam collect each.dom[0]) + ", P=" + sum(Dam collect each.dom[1]) + ", C=" + sum(Dam collect each.dom[2]);
			write "Dim: N=" + sum(Dam collect each.dim[0]) + ", P=" + sum(Dam collect each.dim[1]);
			write "Active nematodes: " + length(Nematode where each.awake);
			write "";
		}
	}
	
	reflex when: (time>6#month) {
		ask simulation {
			do pause;
		}
	}
}

experiment camisol_output {
	// Test to reactivate bacterias
	reflex add_N_P when: cycle mod 300 = 0 {
		ask simulation {
			ask PoreParticle {
				ask dam {
					int pores_count <- length(PoreParticle);
					float carbon_in_all_pore <- total_model_weight * carbone_concentration_in_dam; 
					float carbon_in_pore <- carbon_in_all_pore / pores_count;
					
					float azote_in_all_pore <- total_model_weight * azote_concentration_in_dam;
					float azote_in_pore <- azote_in_all_pore / pores_count;
					
					float phosphore_in_all_pore <- total_model_weight * phosphore_concentration_in_dam;
					float phosphore_in_pore <- phosphore_in_all_pore / pores_count;
					dom[0] <- dom[0] + azote_in_pore;
					dom[1] <- dom[1] + phosphore_in_pore;
					dom[2] <- dom[2] + carbon_in_pore;
				}
			}
		}
	}
	
	output {
//		display grid {
//			grid Particle;
//		}
		
		display "awake population" type: java2D {
			chart "awake population" type: series {
				data "CopioR awake %" value: (sum(Copiotrophe_R collect (each.awake_population)) / length(Copiotrophe_R)) * 100 style:spline;
				data "CopioK awake %" value: (sum(Copiotrophe_K collect (each.awake_population)) / length(Copiotrophe_K)) * 100 style:spline;
				data "Oligo awake %" value: (sum(Oligotrophe_K collect (each.awake_population)) / length(Oligotrophe_K)) * 100 style:spline;
				data "Nematode awake %" value: (sum(Nematode collect (each.awake as int)) / length(Nematode)) * 100 style:spline;
			}
		}
		
		display "dam" type: java2D {
			chart "dam" type:series {
				data "N (dom)" value: (sum(Dam collect each.dom[0])) style:spline;
				data "P (dom)" value: (sum(Dam collect each.dom[1])) style:spline;
				data "C (dom)" value: (sum(Dam collect each.dom[2])) style:spline;
				data "N (dim)" value: (sum(Dam collect each.dim[0])) style:spline;
				data "P (dim)" value: (sum(Dam collect each.dim[1])) style:spline;
			}
		}
		
//		display "CO2" type: java2D {
//			chart "Emitted CO2" type:series {
//				data "CO2" value: total_CO2_produced;
//			}
//		}
		
		display "organics" type:java2D {
			chart "Organics composition" type:series {
				data "C_labile" value: (sum(OrganicParticle collect each.C_labile)) style:spline;
				data "C_recalcitrant" value: (sum(OrganicParticle collect each.C_recalcitrant)) style:spline;
				data "N" value: (sum(OrganicParticle collect each.N)) style:spline;
				data "P" value: (sum(OrganicParticle collect each.P)) style:spline;
			}
		}
		
		display "populations" type:java2D {
			chart "Bacteria populations" type:series {
				data "Copiotrophe K" value: (sum(Copiotrophe_K collect each.C)) style:spline;
				data "Copiotrophe R" value: (sum(Copiotrophe_R collect each.C)) style:spline;
				data "Oligotrophe K" value: (sum(Oligotrophe_K collect each.C)) style:spline;
			}
		}
	}
}

