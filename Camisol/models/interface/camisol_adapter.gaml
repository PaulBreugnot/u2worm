/**
* Name: camisoladapter
* Based on the internal empty template. 
* Author: pbreugno
* Tags: 
*/


model camisoladapter

import "../camisol/camisol.gaml"

global {
	float N_P_dim_available_to_plant <- 0.08;
	
	float production(float N_seed, float P_seed, float N_plant, float P_plant, float harvest_index, float N_from_soil, float plot_surface) {
		// N_from_soil = 0 => an infinite quantity of N can be retrieved from air
		write "N_seed: " + N_seed;
		write "P_seed: " + P_seed;
		write "N_plant: " + N_plant;
		write "P_plant: " + P_plant;
		write "Harvest index: " + harvest_index;
		write "Plot surface: " + plot_surface;
		float required_N <- (harvest_index * N_seed + (1-harvest_index) * N_plant) * N_from_soil;
		float required_P <- harvest_index * P_seed + (1-harvest_index) * P_plant;
		write "Required N: " + required_N;
		write "Required P: " + required_P;
		float total_biomass_production <- 0.0;
		ask Dam {
			// The rest of the dim is "lost" in the nature over the season (for example flushed by water)
			// TODO: remove this hack once the model is complete
			dim[0] <- dim[0] * N_P_dim_available_to_plant;
			dim[1] <- dim[1] * N_P_dim_available_to_plant;
			
			float available_N <- dim[0];
			float available_P <- dim[1];
			float local_biomass_production <- min([
				required_N > 0 ? available_N / required_N : #infinity,
				required_P > 0 ? available_P / required_P : #infinity
			]);
			total_biomass_production <- total_biomass_production + local_biomass_production;
			dim[0] <- max([0.0, dim[0] - total_biomass_production * required_N]);
			dim[1] <- max([0.0, dim[1] - total_biomass_production * required_P]);
		}
		
		write "Total biomass production: " + total_biomass_production;
		// Returns the total mass of seed production
		return harvest_index * total_biomass_production * plot_surface / model_area;
	}

	action fertilize(
		float C_rate,
		float N_dim,
		float P_dim,
		float solubles,
		float hemicellulose,
		float cellulose,
		float lignine,
		float C_N,
		float C_P,
		float sample_dose
	) {
		write "C%: " + C_rate;
		write "N_dim: " + N_dim;
		write "P_dim: " + P_dim;
		write "Solubles: " + solubles;
		write "Hemicellulose: " + hemicellulose;
		write "Cellulose: " + cellulose;
		write "Lignine: " + lignine;
		write "C_N: " + C_N;
		write "C_P: " + C_P;
		write "Sample dose: " + sample_dose;
		// 80% soluble -> C_labile 
		// 20% soluble -> DOM
		
		// hemi + celluose -> C_labile_organics
		// lignine -> C_recalcitrant_organics
		
		float soluble_labile_rate <- 0.8;
		
		//10% de carbone
		write "Sample dose: " + sample_dose;
		write "Model area: " + model_area;
		float quantity_for_model <- sample_dose * model_area * C_rate;
		write "Quantity of C for model: " + quantity_for_model;
		
		// write("Fertilizer quantity for model: "+ qty_for_model);
		
		if(C_N > 0 and C_P > 0) {
			float soluble_Cl <- ((quantity_for_model * solubles / 100)) * soluble_labile_rate;
			float N_pore <- soluble_Cl / C_N;
			float P_pore <- soluble_Cl / C_P;
			
			float soluble_dom <- ((quantity_for_model * solubles / 100)) * (1 - soluble_labile_rate);
			float N_dom <- soluble_dom / C_N;
			float P_dom <- soluble_dom / C_P;
			
			float labile_organic <- ((quantity_for_model * hemicellulose/ 100)) + ((quantity_for_model * cellulose/ 100));
			float recalcitrant_organic <- ((quantity_for_model * lignine/ 100));
			float N_organic <- (labile_organic+recalcitrant_organic) / C_N;
			float P_organic <- (labile_organic+recalcitrant_organic) / C_P;
			
			ask PoreParticle {
				ask organic_particle {
					self.C_labile <- self.C_labile + soluble_Cl/length(PoreParticle);
					self.N <- self.N + N_pore/length(PoreParticle);
					self.P <- self.P + P_pore/length(PoreParticle);
				}
				
				self.dam.dom <- [
					self.dam.dom[0] + (N_dom/length(PoreParticle)),
					self.dam.dom[1] + (P_dom/length(PoreParticle)),
					self.dam.dom[2] + (soluble_dom/length(PoreParticle))
				];
			}
			int num_organics <- length(OrganicParticle)-length(PoreParticle);
			ask OrganicParticle{
				if(!in_pore) {
					self.C_labile <- self.C_labile + labile_organic/num_organics;
					self.C_recalcitrant <- self.C_recalcitrant + recalcitrant_organic/num_organics;
					self.N <- self.N + N_organic/num_organics;
					self.P <- self.P + P_organic/num_organics;
				}
			}
		}
		ask Dam {
			self.dim <- [
				dim[0] + sample_dose * model_area * N_dim/100 / length(PoreParticle),
				dim[1] + sample_dose * model_area * P_dim/100 / length(PoreParticle)
			];
		}
	}
	
//	reflex init_camisol when: local_cycle=1 {
//		write "[Camisol] Time step duration: " + local_step/(1#h) + "h";
//	}
}

experiment TestCamisolWithFertilizer {
	action _init_ {
		create simulation;
		ask simulation {
			// Mada Compost
//			do fertilize
//				C_rate: 16.6
//				N_dim: 0.0
//				P_dim: 0.0
//				solubles: 26.66
//				hemicellulose: 4.72
//				cellulose: 55.19
//				lignine: 13.43
//				C_N: 12.25
//				C_P: 35.0
//				sample_dose: 3000.0#kg/(10000#m2);
			// NPK
			do fertilize
				C_rate: 0.0
				N_dim: 11.0
				P_dim: 22.0
				solubles: 0.0
				hemicellulose: 0.0
				cellulose: 0.0
				lignine: 0.0
				C_N: 0.0
				C_P: 0.0
				sample_dose: 100.0#kg/(10000#m2);
			// Guanomad
//			do fertilize
//				solubles: 82.2
//				hemicellulose: 9.6
//				cellulose: 7.8
//				lignine: 0.43
//				C_N: 3.55
//				C_P: 0.65
//				sample_dose: 500.0;
		}
	}
	
	reflex fert when: local_cycle > 1 and local_cycle mod 4464 = 0 {
		ask simulation {
//			// Rice
			write "Rice production: " + production(
				N_seed: 16.0#kg/#ton,
				P_seed: 3.5#kg/#ton,
				N_plant: 12.5#kg/#ton,
				P_plant: 1.5#kg/#ton,
				harvest_index: 0.44,
				N_from_soil: 1.0,
				plot_surface: 100#m * 100#m
			) + "kg";
			// Tomato
//			write "Tomato production: " + production(
//				N_seed: 1.5#kg/#ton,
//				P_seed: 0.35#kg/#ton,
//				N_plant: 1.5#kg/#ton,
//				P_plant: 0.35#kg/#ton,
//				harvest_index: 0.85#kg/#ton,
//				N_from_soil: 1.0,
//				plot_surface: 10000#m2
//			) + "kg";
			loop i from: 0 to: 0 {
				// Mada Compost
				do fertilize
					C_rate: 16.6
					N_dim: 0.0
					P_dim: 0.0
					solubles: 26.66
					hemicellulose: 4.72
					cellulose: 55.19
					lignine: 13.43
					C_N: 12.25
					C_P: 35.0
					sample_dose: 3000.0#kg/(10000#m2);
//				do fertilize
//					solubles: 82.2
//					hemicellulose: 9.6
//					cellulose: 7.8
//					lignine: 0.43
//					C_N: 3.55
//					C_P: 0.65
//					sample_dose: 500.0#kg/(10000#m2);
			}
			do pause;
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
		
		display "awake population" type: java2D {
			chart "awake population" type: series {
				data "CopioR awake %" value: (sum(Copiotrophe_R collect (each.awake_population)) / length(Copiotrophe_R)) * 100 style:spline color: #red marker:false thickness: 3.0;
				data "CopioK awake %" value: (sum(Copiotrophe_K collect (each.awake_population)) / length(Copiotrophe_K)) * 100 style:spline color: #green marker:false thickness: 3.0;
				data "Oligo awake %" value: (sum(Oligotrophe_K collect (each.awake_population)) / length(Oligotrophe_K)) * 100 style:spline color: #blue marker:false thickness: 3.0;
				data "Nematode awake %" value: (sum(Nematode collect (each.awake as int)) / length(Nematode)) * 100 style:line color: #yellow marker:false thickness: 3.0;
			}
		}
		
//		display "Nematode" type: java2D {
//			chart "Awake nematodes" type: series {
//				data "Nematode awake %" value: (sum(Nematode collect (each.awake as int)) / length(Nematode)) * 100 style:spline color: #yellow;
//			}
//		}
		
		display "dam" type: java2D {
			chart "dam" type:series {
				data "N (dom)" value: (sum(Dam collect each.dom[0])) style:spline marker:false thickness: 3.0;
				data "P (dom)" value: (sum(Dam collect each.dom[1])) style:spline marker:false thickness: 3.0;
				data "C (dom)" value: (sum(Dam collect each.dom[2])) style:spline marker:false thickness: 3.0;
				data "N (dim)" value: (sum(Dam collect each.dim[0])) style:spline marker:false thickness: 3.0;
				data "P (dim)" value: (sum(Dam collect each.dim[1])) style:spline marker:false thickness: 3.0;
			}
		}
		
//		display "CO2" type: java2D {
//			chart "Emitted CO2" type:series {
//				data "CO2" value: total_CO2_produced;
//			}
//		}
		
		display "organics" type:java2D {
			chart "Organics composition" type:series {
				data "C_labile" value: (sum(OrganicParticle collect each.C_labile)) style:spline marker:false thickness: 3.0;
				data "C_recalcitrant" value: (sum(OrganicParticle collect each.C_recalcitrant)) style:spline marker:false thickness: 3.0;
				data "N" value: (sum(OrganicParticle collect each.N)) style:spline marker:false thickness: 3.0;
				data "P" value: (sum(OrganicParticle collect each.P)) style:spline marker:false thickness: 3.0;
			}
		}
		
		display "populations" type:java2D {
			chart "Bacteria populations" type:series {
				data "Copiotrophe R (C)" value: (sum(Copiotrophe_R collect each.C)) style:spline color: #red marker:false thickness: 3.0;
				data "Copiotrophe K (C)" value: (sum(Copiotrophe_K collect each.C)) style:spline color: #green marker:false thickness: 3.0;
				data "Oligotrophe K (C)" value: (sum(Oligotrophe_K collect each.C)) style:spline color: #blue marker:false thickness: 3.0;
			}
		}
	}
}

experiment TestProduction type:gui {
	action _init_ {
		create simulation;
		ask simulation {
			loop i from: 0 to: 3 {
			// Mada Compost
			do fertilize
				C_rate: 16.6
				N_dim: 0.0
				P_dim: 0.0
				solubles: 26.66
				hemicellulose: 4.72
				cellulose: 55.19
				lignine: 13.43
				C_N: 12.25
				C_P: 35.0
				sample_dose: 3000.0#kg/(10000#m2);
			}
			// Guanomad
//			do fertilize
//				solubles: 82.2
//				hemicellulose: 9.6
//				cellulose: 7.8
//				lignine: 0.43
//				C_N: 3.55
//				C_P: 0.65
//				sample_dose: 500.0;
		}
	}
	
	reflex produce when: local_time > 6#month {
		ask simulation {
			do pause;
			// Rice
			write "Production: " + production(
				N_seed: 16.0#kg,
				P_seed: 3.5#kg,
				N_plant: 12.5#kg,
				P_plant: 1.5#kg,
				harvest_index: 0.44,
				N_from_soil: 1.0,
				plot_surface: 100#m * 100#m
			) + "kg";
		}
	}
//	output {
//		display Production {
//			chart "Production" type:series {
//				data "Production (kg)" value: estimate_production(
//					N_seed: 16.0,
//					P_seed: 3.5,
//					N_plant: 12.5,
//					P_plant: 1.5,
//					harvest_index: 0.44,
//					N_from_soil: 1.0,
//					plot_surface: 100#m * 100#m
//				);
//			}
//		}	
//	}
}

experiment Simple type:gui {

	output {
	}
}
