/**
* Name: camisoladapter
* Based on the internal empty template. 
* Author: pbreugno
* Tags: 
*/


model camisoladapter

import "../camisol/camisol.gaml"

global {
	
	float estimate_production(float N_seed, float P_seed, float N_plant, float P_plant, float harvest_index, float N_from_soil, float plot_surface) {
		// N_from_soil = 0 => an infinite quantity of N can be retrieved from air
		write "N_seed: " + N_seed;
		write "P_seed: " + P_seed;
		write "N_plant: " + N_plant;
		write "P_plant: " + P_plant;
		write "Harvest index: " + harvest_index;
		float required_N <- (harvest_index * N_seed + (1-harvest_index) * N_plant) * N_from_soil * #kg * #ton^-1;
		float required_P <- harvest_index * P_seed + (1-harvest_index) * P_seed * #kg * #ton^-1;
		write "Required N: " + required_N;
		write "Required P: " + required_P;
		float available_N <- sum(Dam collect each.dim[0]);
		float available_P <- sum(Dam collect each.dim[0]);
		write "Available N: " + available_N;
		write "Available P: " + available_P;
		float total_biomass_production <- min([
			required_N > 0 ? available_N / required_N : #infinity,
			required_P > 0 ? available_P / required_P : #infinity
		]);
		write "Total biomass production: " + total_biomass_production;

		write plot_surface / model_area;
		// Returns the total mass of seed production
		return harvest_index * total_biomass_production * plot_surface / model_area;
	}
	
	float production(float N_seed, float P_seed, float N_plant, float P_plant, float harvest_index, float N_from_soil, float plot_surface) {
		// N_from_soil = 0 => an infinite quantity of N can be retrieved from air
		write "N_seed: " + N_seed;
		write "P_seed: " + P_seed;
		write "N_plant: " + N_plant;
		write "P_plant: " + P_plant;
		write "Harvest index: " + harvest_index;
		float required_N <- (harvest_index * N_seed + (1-harvest_index) * N_plant) * N_from_soil * #kg * #ton^-1;
		float required_P <- harvest_index * P_seed + (1-harvest_index) * P_seed * #kg * #ton^-1;
		write "Required N: " + required_N;
		write "Required P: " + required_P;
		float available_N <- sum(Dam collect each.dim[0]);
		float available_P <- sum(Dam collect each.dim[1]);
		write "Available N: " + available_N;
		write "Available P: " + available_P;
		float total_biomass_production <- min([
			required_N > 0 ? available_N / required_N : #infinity,
			required_P > 0 ? available_P / required_P : #infinity
		]);
				
		ask Dam {
			dim <- dim - [total_biomass_production * required_N, total_biomass_production * required_P];
		}
		
		write "Total biomass production: " + total_biomass_production;
		// Returns the total mass of seed production
		return harvest_index * total_biomass_production * plot_surface / model_area;
	}

	action fertilize(
		float solubles,
		float hemicellulose,
		float cellulose,
		float lignine,
		float C_N,
		float C_P,
		float sample_dose
	) {
		// 80% soluble -> C_labile 
		// 20% soluble -> DOM
		
		// hemi + celluose -> C_labile_organics
		// lignine -> C_recalcitrant_organics
		
		float soluble_labile_rate <- 0.8;
		
		//10% de carbone
		
		float hectare <- 10000#m2;
		float model_size <- model_area; // todo selon l'echelle
		
		float model_ratio <- model_size / hectare;
		
		float quantity_for_model <- sample_dose * model_ratio;
		
		// write("Fertilizer quantity for model: "+ qty_for_model);
		
		float soluble_Cl <- ((quantity_for_model * solubles / 100)) * soluble_labile_rate;
		float N_pore <- soluble_Cl * C_N;
		float P_pore <- soluble_Cl * C_P;
		
		float soluble_dom <- ((quantity_for_model * solubles / 100)) * (1 - soluble_labile_rate);
		float N_dom <- soluble_dom * C_N;
		float P_dom <- soluble_dom * C_N;
		
		float labile_organic <- ((quantity_for_model * hemicellulose/ 100)) + ((quantity_for_model * cellulose/ 100));
		float recalcitrant_organic <- ((quantity_for_model * lignine/ 100));
		float N_organic <- (labile_organic+recalcitrant_organic) * C_N;
		float P_organic <- (labile_organic+recalcitrant_organic) * C_N;
		
		ask PoreParticle {
			ask organic_particle {
				self.C_labile <- self.C_labile + soluble_Cl;
				self.N <- self.N + N_pore;
				self.P <- self.P + P_pore;
			}
			
			self.dam.dom <- [
				self.dam.dom[0] + (N_dom/length(PoreParticle)),
				self.dam.dom[1] + (P_dom/length(PoreParticle)),
				self.dam.dom[2] + (soluble_dom/length(PoreParticle))
			];
		}
		ask OrganicParticle{
			if(!in_pore) {
				self.C_labile <- self.C_labile + labile_organic;
				self.C_recalcitrant <- self.C_recalcitrant + recalcitrant_organic;
				self.N <- self.N + N_organic;
				self.P <- self.P + P_organic;
			}
		}
	}
}

experiment TestProduction type:gui {
	action _init_ {
		create simulation;
		ask simulation {
			// Mada Compost
			do fertilize
				solubles: 26.66
				hemicellulose: 4.72
				cellulose: 55.19
				lignine: 13.43
				C_N: 12.25
				C_P: 35.0
				sample_dose: 3000.0;
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
	
	reflex produce when: time > 6#month {
		ask simulation {
		do pause;
		write "Production: " + production(
			N_seed: 16.0,
			P_seed: 3.5,
			N_plant: 12.5,
			P_plant: 1.5,
			harvest_index: 0.44,
			N_from_soil: 1.0,
			plot_surface: 10#m * 10#m
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