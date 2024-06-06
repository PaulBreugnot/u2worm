/**
* Name: camisol
* Based on the internal empty template. 
* Author: pbreugno
* Tags: 
*/


model cammisol

import "environment/grid.gaml"
import "nematode/nematode.gaml"

global {
	int nematodes_count <- 50;


		
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
	float init_Y_rate <- 0.1;
	float init_A_rate <- 0.2;
	float init_S_rate <- 0.7;
	
	float rain_diffusion_rate <- 0.1;
	float rain_period <- 7#days;
	
	init {
		do init_grid;
		do init_enzymatic_optimisation;
		do init_enzymes;
		// Counts the number of PORES after the initialization of the grid
		int pores_count <- length(PoreParticle);
		ask PoreParticle {
			// The carrying capacity of each pore is equal to 10 times the initial bacteria population
			carrying_capacity <- 10 * total_initial_bacteria_weight / pores_count;
			
			create Y_Strategist with: [C::init_Y_rate * total_initial_bacteria_weight / length(PoreParticle)] {
				add self to:myself.populations;
			}
			create A_Strategist with: [C::init_A_rate * total_initial_bacteria_weight / length(PoreParticle)]{
				add self to:myself.populations;
			}
			create S_Strategist with: [C::init_S_rate * total_initial_bacteria_weight / length(PoreParticle)]{
				add self to:myself.populations;
			}
		}
		
		create Nematode number: nematodes_count
		{
			current_pore <- one_of(PoreParticle); 
			location <- any_location_in(current_pore);
		}
		
		ask PoreParticle {
			ask populations {
				do update;
				write "Initial enzymes optimization for " + self;
				do optimize_enzymes(myself.dam, myself.accessible_organics);
			}
		}
	}
	
//	reflex rain when: rnd(1.0) < (1/rain_period)*local_step {
////		write "It's raining today!";
//		map<species<MicrobePopulation>, MicrobePopulation> extracted_populations;
//		loop bacteria_type over: bacteria_types {
//			create MicrobePopulation with: (
//				C: 0.0, cytosol_C: 0.0,
//				N: 0.0, cytosol_N: 0.0,
//				P: 0.0, cytosol_P: 0.0
//			) {
//				extracted_populations[bacteria_type] <- self;
//			}
//		}
//		
//		ask PoreParticle {
//			ask populations {
//				ask extracted_populations[species(self)] {
//					float extracted_C <- myself.C * rain_diffusion_rate;
//					self.C <- self.C + extracted_C;
//					myself.C <- myself.C - extracted_C;
//					
//					float extracted_P <- myself.P * rain_diffusion_rate;
//					self.P <- self.P + extracted_P;
//					myself.P <- myself.P - extracted_P;
//					
//					float extracted_N <- myself.N * rain_diffusion_rate;
//					self.N <- self.N + extracted_N;
//					myself.N <- myself.N - extracted_N;
//					
//					float extracted_cytosol_C <- myself.cytosol_C * rain_diffusion_rate;
//					self.cytosol_C <- self.cytosol_C + extracted_cytosol_C;
//					myself.cytosol_C <- myself.cytosol_C - extracted_cytosol_C;
//					
//					float extracted_cytosol_P <- myself.cytosol_P * rain_diffusion_rate;
//					self.cytosol_P <- self.cytosol_P + extracted_cytosol_P;
//					myself.cytosol_P <- myself.cytosol_P - extracted_cytosol_P;
//					
//					float extracted_cytosol_N <- myself.cytosol_N * rain_diffusion_rate;
//					self.cytosol_N <- self.cytosol_N + extracted_cytosol_N;
//					myself.cytosol_N <- myself.cytosol_N - extracted_cytosol_N;
//				}
//			}
//		}
//		
//		ask PoreParticle {
//			ask populations {
//				MicrobePopulation extracted_pop <- extracted_populations[species(self)];
//				self.C <- self.C + extracted_pop.C/length(PoreParticle);
//				self.N <- self.N + extracted_pop.N/length(PoreParticle);
//				self.P <- self.P + extracted_pop.P/length(PoreParticle);
//				self.cytosol_C <- self.cytosol_C + extracted_pop.cytosol_C/length(PoreParticle);
//				self.cytosol_N <- self.cytosol_N + extracted_pop.cytosol_N/length(PoreParticle);
//				self.cytosol_P <- self.cytosol_P + extracted_pop.cytosol_P/length(PoreParticle);
//			}
//		}
//	
//		loop bacteria_type over: bacteria_types {
//			ask extracted_populations[bacteria_type] {
//				do die;
//			}
//		}
//	}

	reflex {
		ask shuffle(Nematode) {
			do life;
		}
		ask shuffle(PoreParticle) {
			ask populations {
				if flip(local_step / enzymes_optimization_period) {
					do update;
					do optimize_enzymes(myself.dam, myself.accessible_organics);
				}
			}
			do decompose;
			do microbe_life;
		}
	}
}

experiment base_cammisol_output {
	parameter "Grid size" category: "Environment" var:grid_size;
	parameter "Nematodes count" category: "Environment" var:nematodes_count;
	parameter "Organic particle rate" category: "Environment" var:organic_rate;
	parameter "Mineral particle rate" category: "Environment" var:mineral_rate;
	parameter "Init C (FOM)" category: "Environment" var:C_concentration_in_pom;
	parameter "Init N (FOM)" category: "Environment" var:N_concentration_in_pom;
	parameter "Init P (FOM)" category: "Environment" var:P_concentration_in_pom;
	parameter "Init labile rate (FOM)" category: "Environment" var:labile_rate_pom;
	
	parameter "Nematodes predation rate" category: "Nematode" var:nematode_predation_rate;
	parameter "Nematodes CUE" category: "Nematode" var:nematode_CUE;
	parameter "Nematodes C/N" category: "Nematode" var:nematode_C_N;
	parameter "Nematodes C/P" category: "Nematode" var:nematode_C_P;
	
	parameter "Dividing time (Y)" category: "Y Strategists" var:dividing_time_Y;
	parameter "CUE (Y)" category: "Y Strategists" var:carbon_use_efficiency_Y;
	parameter "Minimum awake rate (Y)" category: "Y Strategists" var:minimum_awake_rate_Y;
	
	parameter "Dividing time (A)" category: "A Strategists" var:dividing_time_A;
	parameter "CUE (A)" category: "A Strategists" var:carbon_use_efficiency_A;
	parameter "Minimum awake rate (A)" category: "A Strategists" var:minimum_awake_rate_A;
	
	parameter "Dividing time (S)" category: "S Strategists" var:dividing_time_S;
	parameter "CUE (S)" category: "S Strategists" var:carbon_use_efficiency_S;
	parameter "Minimum awake rate (S)" category: "S Strategists" var:minimum_awake_rate_S;
	
	reflex update_particle_color {
		float max_population <- 0.0;
		ask PoreParticle {
			float population <- sum(populations collect (each.C + each.cytosol_C));
			if(population > max_population) {
				max_population <- population;
			}
		}
		if(max_population > 0.0) {
			ask Particle {
				if(type = PORE) {
					color <- rgb(0, 0, 255 * sum(PoreParticle(particle).populations collect (each.C + each.cytosol_C))/max_population);
				}
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
				if (length(PoreParticle) > 0) {
					data "Y awake (%)" value: sum(Y_Strategist collect (each.awake_population))/length(Y_Strategist) * 100 style:spline color: #red marker:false thickness:3;
					data "A awake (%)" value: sum(A_Strategist collect (each.awake_population))/length(A_Strategist) * 100 style:spline color: #green marker:false thickness:3;
					data "S awake (%)" value: sum(S_Strategist collect (each.awake_population))/length(S_Strategist) * 100 style:spline color: #blue marker:false thickness:3;
				}
				if(nematodes_count > 0) {
					data "Nematode awake (%)" value: (sum(Nematode collect (each.awake as int)) / length(Nematode)) * 100 style:spline color: #yellow marker:false thickness:3;				
				}
			}
		}
			
		display "dam" type: java2D {
			chart "dam" type:series {
				data "N dom (g)" value: sum(Dam collect each.dom[0])/#gram style:spline marker:false thickness:3;
				data "P dom (g)" value: sum(Dam collect each.dom[1])/#gram style:spline marker:false thickness:3;
				data "C dom (g)" value: sum(Dam collect each.dom[2])/#gram style:spline marker:false thickness:3;
				data "N dim (g)" value: sum(Dam collect each.dim[0])/#gram style:spline marker:false thickness:3;
				data "P dim (g)" value: sum(Dam collect each.dim[1])/#gram style:spline marker:false thickness:3;
			}
		}
		
		display "organics" type:java2D {
			chart "Organics composition" type:series {
				data "C labile (g)" value: sum(OrganicParticle collect each.C_labile)/#gram style:spline marker:false thickness:3;
				data "N labile (g)" value: sum(OrganicParticle collect each.N_labile)/#gram style:spline marker:false thickness:3;
				data "P labile (g)" value: sum(OrganicParticle collect each.P_labile)/#gram style:spline marker:false thickness:3;
				data "C recalcitrant (g)" value: sum(OrganicParticle collect each.C_recalcitrant)/#gram style:spline marker:false thickness:3;
				data "N recalcitrant (g)" value: sum(OrganicParticle collect each.N_recalcitrant)/#gram style:spline marker:false thickness:3;
				data "P recalcitrant (g)" value: sum(OrganicParticle collect each.P_recalcitrant)/#gram style:spline marker:false thickness:3;
			}
		}
		
		display "populations" type:java2D {
			chart "Bacteria populations" type:series {
				data "Y (g)" value: sum(Y_Strategist collect each.C)/#gram style:spline color: #red marker:false thickness:3;
				data "A (g)" value: sum(A_Strategist collect each.C)/#gram style:spline color: #green marker:false thickness:3;
				data "S (g)" value: sum(S_Strategist collect each.C)/#gram style:spline color: #blue marker:false thickness:3;
			}
		}
	}
}

experiment test_microbes parent:base_cammisol_output {
	float C_soil <- 1#gram/(#cm#cm);
	float N_soil <- 0.1#gram/(#cm#cm);
	float P_soil <- 0.05#gram/(#cm#cm);
	
	parameter "Grid size" var:grid_size;
	parameter "Nematodes" var:nematodes_count;
	parameter "C soil" var:C_soil;
	parameter "N soil" var:N_soil;
	parameter "P soil" var:P_soil;
	parameter "Initial Y population rate" var:init_Y_rate;
	parameter "Initial A population rate" var:init_A_rate;
	parameter "Initial S population rate" var:init_S_rate;
	
	parameter "Enzymes optimization period" var:enzymes_optimization_period init:local_step;
	
	init {
		ask simulation {
			do init_POM(myself.C_soil, myself.N_soil, myself.P_soil);		
		}
	}
}

experiment display_grid {
	bool show_nematodes;
	parameter "Show nematodes" var: show_nematodes <- false;
	parameter "Grid size" category: "Environment" var:grid_size;
	
	output {
		display grid type:2d axes:false {
			grid Particle;
			graphics legend {
				draw square(cell_size) at: {soil_size + cell_size, 1.5*cell_size, 0} color: #yellow;
				draw square(cell_size) at: {soil_size + cell_size, 2.8*cell_size, 0} color: #green;
				draw square(cell_size) at: {soil_size + cell_size, 4.1*cell_size, 0} color: #black;
				draw "Mineral particle" font:font("Helvetica", 20 , #bold) at: {soil_size + 2.1*cell_size, 1.7*cell_size, 0} color: #black;
				draw "Organic particle" font:font("Helvetica", 20 , #bold) at: {soil_size + 2.1*cell_size, 3*cell_size, 0} color: #black;
				draw "Pore particle" font:font("Helvetica", 20 , #bold) at: {soil_size + 2.1*cell_size, 4.3*cell_size, 0} color: #black;
				if show_nematodes {
					draw circle(cell_size/4) at: {soil_size + cell_size, 5.4*cell_size, 0} color: #red;
					draw "Nematode" font:font("Helvetica", 20 , #bold) at: {soil_size + 2.1*cell_size, 5.5*cell_size, 0} color: #black;
				}
			}
			species Nematode aspect: red_dot visible: show_nematodes;
		}
	}
}

experiment cammisol_no_output {

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

experiment cammisol parent:base_cammisol_output {
	// Test to reactivate bacterias
	reflex add_N_P when: local_cycle mod 300 = 0 {
		ask simulation {
			// Redistributes the initial quantity of C/N/P in all pores.
			// Notice that this has no scientific meaning, and is used only for test purpose.
			// This might however represent a kind of fertilizer added to the soil.
			int pores_count <- length(PoreParticle);
			float carbon_in_all_pore <- soil_weight * C_concentration_in_dom; 
			float azote_in_all_pore <- soil_weight * N_concentration_in_dom;
			float phosphore_in_all_pore <- soil_weight * P_concentration_in_dom;
			
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
}

