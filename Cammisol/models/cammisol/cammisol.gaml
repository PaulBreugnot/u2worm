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

	int grid_size <- 30;
	float env_size <- 1#cm;
	
	int nematodes_count <- 50;
	
	geometry shape <- square(env_size);
	float model_area <- shape.area;
	float cell_size <- env_size/grid_size;
	float cell_area <- cell_size^2;

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
		
	// 5E8 -> 5E9 bacterie / gramme de sol
	/*
	 * A modifier................................................ 
	 */
	 // TODO: poids total de bactérie par gramme de sol * surface du modèle
	/**
	 * Total bacteria weight in the model.
	 */
	float total_initial_bacteria_weight <-  0.05*1.5#gram/(#cm*#cm)*world.shape.area;
	// TODO: what rates?
	float copiotrophe_R_rate <- 0.1;
	float copiotrophe_K_rate <- 0.2;
	float oligotrophe_K_rate <- 0.7;
	
	float rain_diffusion_rate <- 0.1;
	float rain_period <- 7#days;
	
	init {
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
			
			create Y_Strategist with: [C::copiotrophe_R_rate * total_initial_bacteria_weight / length(PoreParticle)] {
				add self to:myself.populations;
			}
			create A_Strategist with: [C::copiotrophe_K_rate * total_initial_bacteria_weight / length(PoreParticle)]{
				add self to:myself.populations;
			}
			create S_Strategist with: [C::oligotrophe_K_rate * total_initial_bacteria_weight / length(PoreParticle)]{
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
			location <- any_location_in(current_pore);
		}
	 	// Default soil initialization. Can be safely overriden with subsequent init_soil() calls.
		do init_soil(
			C: 0.02157#gram/(#cm*#cm),
			N: 0.00132#gram/(#cm*#cm),
			P: 0.00077#gram/(#cm*#cm)
		);
	}
	
	action init_soil(float C, float N, float P) {
		write "Soil C: " + C + "kg/m2";
		write "Soil N: " + N + "kg/m2";
		write "Soil P: " + P + "kg/m2";
		ask MineralParticle {
			N <- N*cell_area;
			P <- P*cell_area;
		}
		ask OrganicParticle {
			if(!in_pore) {
				N <- N*cell_area;
				P <- P*cell_area;
				C_labile <- (C/2)*cell_area; // TODO: variable factor
				C_recalcitrant <- (C/2)*cell_area;
			}
		}
	}
		
	reflex rain when: rnd(1.0) < (1/rain_period)*local_step {
//		write "It's raining today!";
		map<species<MicrobePopulation>, MicrobePopulation> extracted_populations;
		loop bacteria_type over: bacteria_types {
			create MicrobePopulation with: (
				C: 0.0, cytosol_C: 0.0,
				N: 0.0, cytosol_N: 0.0,
				P: 0.0, cytosol_P: 0.0
			) {
				extracted_populations[bacteria_type] <- self;
			}
		}
		
		ask PoreParticle {
			ask populations {
				ask extracted_populations[species(self)] {
					float extracted_C <- myself.C * rain_diffusion_rate;
					self.C <- self.C + extracted_C;
					myself.C <- myself.C - extracted_C;
					
					float extracted_P <- myself.P * rain_diffusion_rate;
					self.P <- self.P + extracted_P;
					myself.P <- myself.P - extracted_P;
					
					float extracted_N <- myself.N * rain_diffusion_rate;
					self.N <- self.N + extracted_N;
					myself.N <- myself.N - extracted_N;
					
					float extracted_cytosol_C <- myself.cytosol_C * rain_diffusion_rate;
					self.cytosol_C <- self.cytosol_C + extracted_cytosol_C;
					myself.cytosol_C <- myself.cytosol_C - extracted_cytosol_C;
					
					float extracted_cytosol_P <- myself.cytosol_P * rain_diffusion_rate;
					self.cytosol_P <- self.cytosol_P + extracted_cytosol_P;
					myself.cytosol_P <- myself.cytosol_P - extracted_cytosol_P;
					
					float extracted_cytosol_N <- myself.cytosol_N * rain_diffusion_rate;
					self.cytosol_N <- self.cytosol_N + extracted_cytosol_N;
					myself.cytosol_N <- myself.cytosol_N - extracted_cytosol_N;
				}
			}
		}
		
		ask PoreParticle {
			ask populations {
				MicrobePopulation extracted_pop <- extracted_populations[species(self)];
				self.C <- self.C + extracted_pop.C/length(PoreParticle);
				self.N <- self.N + extracted_pop.N/length(PoreParticle);
				self.P <- self.P + extracted_pop.P/length(PoreParticle);
				self.cytosol_C <- self.cytosol_C + extracted_pop.cytosol_C/length(PoreParticle);
				self.cytosol_N <- self.cytosol_N + extracted_pop.cytosol_N/length(PoreParticle);
				self.cytosol_P <- self.cytosol_P + extracted_pop.cytosol_P/length(PoreParticle);
			}
		}
	
		loop bacteria_type over: bacteria_types {
			ask extracted_populations[bacteria_type] {
				do die;
			}
		}
	}
}

grid Particle width: grid_size height: grid_size neighbors: 4 {
	string type <- rnd_choice([ORGANIC::organic_rate, MINERAL::mineral_rate, PORE::pore_rate]);
	agent particle;
	
	init {
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

experiment display_grid {
	bool show_nematodes;
	parameter "Show nematodes" var: show_nematodes <- false;
	
	output {
		display grid type:2d axes:false {
			grid Particle;
			graphics legend {
				draw square(cell_size) at: {env_size + cell_size, 1.5*cell_size, 0} color: #yellow;
				draw square(cell_size) at: {env_size + cell_size, 2.8*cell_size, 0} color: #green;
				draw square(cell_size) at: {env_size + cell_size, 4.1*cell_size, 0} color: #black;
				draw "Mineral particle" font:font("Helvetica", 20 , #bold) at: {env_size + 2.1*cell_size, 1.7*cell_size, 0} color: #black;
				draw "Organic particle" font:font("Helvetica", 20 , #bold) at: {env_size + 2.1*cell_size, 3*cell_size, 0} color: #black;
				draw "Pore particle" font:font("Helvetica", 20 , #bold) at: {env_size + 2.1*cell_size, 4.3*cell_size, 0} color: #black;
				if show_nematodes {
					draw circle(cell_size/4) at: {env_size + cell_size, 5.4*cell_size, 0} color: #red;
					draw "Nematode" font:font("Helvetica", 20 , #bold) at: {env_size + 2.1*cell_size, 5.5*cell_size, 0} color: #black;
				}
			}
			species Nematode aspect: red_dot visible: show_nematodes;
		}
	}
}

experiment camisol_no_output {

	reflex state when: local_cycle mod 100 = 0 {
		ask simulation {
			write "Time: " + time/#day;
			write "Dom: N=" + sum(Dam collect each.dom[0]) + ", P=" + sum(Dam collect each.dom[1]) + ", C=" + sum(Dam collect each.dom[2]);
			write "Dim: N=" + sum(Dam collect each.dim[0]) + ", P=" + sum(Dam collect each.dim[1]);
			write "Active nematodes: " + length(Nematode where each.awake);
			write "";
		}
	}
	
	reflex when: (local_time>6#month) {
		ask simulation {
			do pause;
		}
	}
}

experiment camisol {
	// Test to reactivate bacterias
	reflex add_N_P when: local_cycle mod 300 = 0 {
		ask simulation {
			// Redistributes the initial quantity of C/N/P in all pores.
			// Notice that this has no scientific meaning, and is used only for test purpose.
			// This might however represent a kind of fertilizer added to the soil.
			int pores_count <- length(PoreParticle);
			float carbon_in_all_pore <- total_model_weight * carbone_concentration_in_dam; 
			float azote_in_all_pore <- total_model_weight * azote_concentration_in_dam;
			float phosphore_in_all_pore <- total_model_weight * phosphore_concentration_in_dam;
			
			ask PoreParticle {
				ask dam {
					float carbon_in_pore <- carbon_in_all_pore / pores_count;
					float azote_in_pore <- azote_in_all_pore / pores_count;
					float phosphore_in_pore <- phosphore_in_all_pore / pores_count;
					
					dom[0] <- dom[0] + azote_in_pore;
					dom[1] <- dom[1] + phosphore_in_pore;
					dom[2] <- dom[2] + carbon_in_pore;
				}
			}
		}
	}
	
	reflex update_particle_color {
		float max_population <- 0.0;
		ask PoreParticle {
			float population <- sum(populations collect (each.C + each.cytosol_C));
			if(population > max_population) {
				max_population <- population;
			}
		}
		ask Particle {
			if(type = PORE) {
				color <- rgb(0, 0, 255 * sum(PoreParticle(particle).populations collect (each.C + each.cytosol_C))/max_population);
			}
		}
	}
	
	output {
		display grid {
			grid Particle;
			species Nematode aspect: red_dot;
		}
		
		display "Awoken population" type: java2D {
			chart "Awoken population" type: series {
				data "CopioR awake %" value: (sum(Y_Strategist collect (each.awake_population)) / length(Y_Strategist)) * 100 style:spline color: #red marker:false thickness:3;
				data "CopioK awake %" value: (sum(A_Strategist collect (each.awake_population)) / length(A_Strategist)) * 100 style:spline color: #green marker:false thickness:3;
				data "Oligo awake %" value: (sum(S_Strategist collect (each.awake_population)) / length(S_Strategist)) * 100 style:spline color: #blue marker:false thickness:3;
				data "Nematode awake %" value: (sum(Nematode collect (each.awake as int)) / length(Nematode)) * 100 style:spline color: #yellow marker:false thickness:3;
			}
		}
			
		display "dam" type: java2D {
			chart "dam" type:series {
				data "N (dom)" value: (sum(Dam collect each.dom[0])) style:spline marker:false thickness:3;
				data "P (dom)" value: (sum(Dam collect each.dom[1])) style:spline marker:false thickness:3;
				data "C (dom)" value: (sum(Dam collect each.dom[2])) style:spline marker:false thickness:3;
				data "N (dim)" value: (sum(Dam collect each.dim[0])) style:spline marker:false thickness:3;
				data "P (dim)" value: (sum(Dam collect each.dim[1])) style:spline marker:false thickness:3;
			}
		}
		
		display "organics" type:java2D {
			chart "Organics composition" type:series {
				data "C_labile" value: (sum(OrganicParticle collect each.C_labile)) style:spline marker:false thickness:3;
				data "N_labile" value: (sum(OrganicParticle collect each.N_labile)) style:spline marker:false thickness:3;
				data "P_labile" value: (sum(OrganicParticle collect each.P_labile)) style:spline marker:false thickness:3;
				data "C_recalcitrant" value: (sum(OrganicParticle collect each.C_recalcitrant)) style:spline marker:false thickness:3;
				data "N_recalcitrant" value: (sum(OrganicParticle collect each.N_recalcitrant)) style:spline marker:false thickness:3;
				data "P_recalcitrant" value: (sum(OrganicParticle collect each.P_recalcitrant)) style:spline marker:false thickness:3;
			}
		}
		
		display "populations" type:java2D {
			chart "Bacteria populations" type:series {
				data "Copiotrophe R" value: (sum(Y_Strategist collect each.C)) style:spline color: #red marker:false thickness:3;
				data "Copiotrophe K" value: (sum(A_Strategist collect each.C)) style:spline color: #green marker:false thickness:3;
				data "Oligotrophe K" value: (sum(S_Strategist collect each.C)) style:spline color: #blue marker:false thickness:3;
			}
		}
	}
}

