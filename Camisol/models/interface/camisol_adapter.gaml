/**
* Name: camisoladapter
* Based on the internal empty template. 
* Author: pbreugno
* Tags: 
*/


model camisoladapter

import "../camisol.gaml"

global {
	reflex print_step when: (local_time mod 100 = 0) {
		write "Step " + local_time;
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
		write "Masse totale N: " + masse_total_N;
		write "Masse totale P: " + masse_total_P;
		float total_biomass_production <- min([
			required_N > 0 ? masse_total_N / required_N : #infinity,
			required_P > 0 ? masse_total_P / required_P : #infinity
		]);
		write "Total biomass production: " + total_biomass_production;
		
		ask Dam {
			dim <- dim - [total_biomass_production * required_N, total_biomass_production * required_P];
		}

		write plot_surface / camisol_area();
		// Returns the total mass of seed production
		return harvest_index * total_biomass_production * plot_surface / camisol_area();
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
		float model_size <- camisol_area(); // todo selon l'echelle
		
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
		
		ask pores_particles{
			self.C_labile <- self.C_labile + soluble_Cl;
			self.N <- self.N + N_pore;
			self.P <- self.P + P_pore;
			
			self.dam.dom <- [self.dam.dom[0] + (N_dom/nb_pores),self.dam.dom[1] + (P_dom/nb_pores),self.dam.dom[2] + (soluble_dom/nb_pores)];
		}
		ask organics_particles{
			self.C_labile <- self.C_labile + labile_organic;
			self.C_recalcitrant <- self.C_recalcitrant + recalcitrant_organic;
			self.N <- self.N + N_organic;
			self.P <- self.P + P_organic;
		}
	}
}

experiment Simple type:gui {

	output {
	}
}