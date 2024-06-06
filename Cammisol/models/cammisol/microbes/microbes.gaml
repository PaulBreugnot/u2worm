/**
* Name: microbes
* Based on the internal empty template. 
* Author: pbreugno
* Tags: 
*/

model microbes

import "../environment/pore_particle.gaml"

global {
	list<species<MicrobePopulation>> bacteria_types <- [Y_Strategist, A_Strategist, S_Strategist];

	float dividing_time_Y <- 1#h;
	float dividing_time_A <- 24#h;
	float dividing_time_S <- 368#h;
	
	float carbon_use_efficiency_Y <- 0.7;
	float carbon_use_efficiency_A <- 0.3;
	float carbon_use_efficiency_S <- 0.5;
	
	float minimum_awake_rate_Y <- 0.002;
	float minimum_awake_rate_A <- 0.018;
	float minimum_awake_rate_S <- 1.0;
	
	float min_T_cellulolytic_Y <- 0.0;
	float min_T_amino_Y <- 0.0;
	float min_T_P_Y <- 0.0;
	float min_T_recal_Y <- 0.0;
	float max_T_cellulolytic_Y <- 1 #gram / #gram / #d;
	float max_T_amino_Y <- 0.2 #gram / #gram / #d;
	float max_T_P_Y <- 0.04 #gram / #gram / #d;
	float max_T_recal_Y <- 0.0;
	
	float min_T_cellulolytic_A <- 0.0;
	float min_T_amino_A <- 0.0;
	float min_T_P_A <- 0.0;
	float min_T_recal_A <- 0.0;
	float max_T_cellulolytic_A <- 5 #gram/ #gram / #d;
	float max_T_amino_A <- 1 #gram / #gram / #d;
	float max_T_P_A <- 0.2 #gram / #gram / #d;
	float max_T_recal_A <- 1 #gram / #gram / #d;
	
	float min_T_cellulolytic_S <- 0.1 #gram/ #gram / #d;
	float min_T_amino_S <- 0.02 #gram/ #gram / #d;
	float min_T_P_S <- 0.004 #gram/ #gram / #d;
	float min_T_recal_S <- 0.0;
	float max_T_cellulolytic_S <- 0.2 #gram/ #gram / #d;
	float max_T_amino_S <- 0.04 #gram / #gram / #d;
	float max_T_P_S <- 0.008 #gram / #gram / #d;
	float max_T_recal_S <- 0.0;
	
	Enzymes min_enzymes_Y;
	Enzymes max_enzymes_Y;
	Enzymes min_enzymes_A;
	Enzymes max_enzymes_A;
	Enzymes min_enzymes_S;
	Enzymes max_enzymes_S;
	
	action init_enzymes {
		create Enzymes with: [
			name::"Min enzymes (Y)",
			T_cellulolytic::min_T_cellulolytic_Y,
			T_amino::min_T_amino_Y,
			T_P::min_T_P_Y,
			T_recal::min_T_recal_Y
		] {
			myself.min_enzymes_Y <- self;
		}
		
		create Enzymes with: [
			name::"Max enzymes (Y)",
			T_cellulolytic::max_T_cellulolytic_Y,
			T_amino::max_T_amino_Y,
			T_P::max_T_P_Y,
			T_recal::max_T_recal_Y
		] {
			myself.max_enzymes_Y <- self;
		}
		
		create Enzymes with: [
			name::"Min enzymes (A)",
			T_cellulolytic::min_T_cellulolytic_A,
			T_amino::min_T_amino_A,
			T_P::min_T_P_A,
			T_recal::min_T_recal_A
		] {
			myself.min_enzymes_A <- self;
		}
		
		create Enzymes with: [
			name::"Max enzymes (A)",
			T_cellulolytic::max_T_cellulolytic_A,
			T_amino::max_T_amino_A,
			T_P::max_T_P_A,
			T_recal::max_T_recal_A
		] {
			myself.max_enzymes_A <- self;
		}	
		
		create Enzymes with: [
			name::"Min enzymes (S)",
			T_cellulolytic::min_T_cellulolytic_S,
			T_amino::min_T_amino_S,
			T_P::min_T_P_S,
			T_recal::min_T_recal_S
		] {
			myself.min_enzymes_S <- self;
		}
		
		create Enzymes with: [
			name::"Max enzymes (S)",
			T_cellulolytic::max_T_cellulolytic_S,
			T_amino::max_T_amino_S,
			T_P::max_T_P_S,
			T_recal::max_T_recal_S
		] {
			myself.max_enzymes_S <- self;
		}	
	}
}

species Y_Strategist parent:MicrobePopulation schedules:[] {
	init {
		dividing_time <- dividing_time_Y;
		carbon_use_efficiency <- carbon_use_efficiency_Y;
		minimum_awake_rate <- minimum_awake_rate_Y;
		
		cytosol_mineralization_rate <- 1.0;
		
		min_enzymes <- min_enzymes_Y;
		max_enzymes <- max_enzymes_Y;
	}
}
species A_Strategist parent:MicrobePopulation schedules:[] {
	init {
		dividing_time <- dividing_time_A;
		carbon_use_efficiency <- carbon_use_efficiency_A;
		minimum_awake_rate <- minimum_awake_rate_A;
		
		cytosol_mineralization_rate <- 0.0;
		
		min_enzymes <- min_enzymes_A;
		max_enzymes <- max_enzymes_A;
	}
}

species S_Strategist parent:MicrobePopulation schedules:[] {
	init {
		dividing_time <- dividing_time_S;
		carbon_use_efficiency <- carbon_use_efficiency_S;
		minimum_awake_rate <- minimum_awake_rate_S;

		cytosol_mineralization_rate <- 0.0;
		
		min_enzymes <- min_enzymes_S;
		max_enzymes <- max_enzymes_S;
	}
}
