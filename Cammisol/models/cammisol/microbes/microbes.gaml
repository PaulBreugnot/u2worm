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
	
	float carbon_use_efficiency_Y <- 7/10;
	float carbon_use_efficiency_A <- 3/10;
	float carbon_use_efficiency_S <- 5/10;
	
	float minimum_awake_rate_Y <- 0.002;
//	float minimum_awake_rate_A <- 0.018;
	float minimum_awake_rate_A <- 1.0;
	float minimum_awake_rate_S <- 1.0;
}

species Y_Strategist parent:MicrobePopulation schedules:[] {
	init {
		dividing_time <- dividing_time_Y;
		carbon_use_efficiency <- carbon_use_efficiency_Y;
		minimum_awake_rate <- minimum_awake_rate_Y;
		
		cytosol_mineralization_rate <- 1.0;
		
		Enzymes _min_enzymes;
		create Enzymes with: [
			T_cellulolytic::0.0,
			T_amino::0.0,
			T_P::0.0,
			T_recal::0.0
		] {
			_min_enzymes <- self;
		}
		
		Enzymes _max_enzymes;
		create Enzymes with: [
			T_cellulolytic::10 #gram/ #gram / #d,
			T_amino::1.0 #gram / #gram / #d,
			T_P::0.5 #gram / #gram / #d,
			T_recal::0.00001 #gram / #gram / #d
		] {
			_max_enzymes <- self;
		}
		do set_min_max_enzymes(_min_enzymes, _max_enzymes);
	}
}
species A_Strategist parent:MicrobePopulation schedules:[] {
	init {
		dividing_time <- dividing_time_A;
		carbon_use_efficiency <- carbon_use_efficiency_A;
		minimum_awake_rate <- minimum_awake_rate_A;
		
		cytosol_mineralization_rate <- 0.0;

		Enzymes _min_enzymes;
		create Enzymes with: [
			T_cellulolytic::0.0,
			T_amino::0.0,
			T_P::0.0,
			T_recal::0.0
		] {
			_min_enzymes <- self;
		}
		
		Enzymes _max_enzymes;
		create Enzymes with: [
			T_cellulolytic::40 #gram/ #gram / #d,
			T_amino::80 #gram / #gram / #d,
			T_P::20 #gram / #gram / #d,
			T_recal::80 #gram / #gram / #d
		] {
			_max_enzymes <- self;
		}
		do set_min_max_enzymes(_min_enzymes, _max_enzymes);
	}
}

species S_Strategist parent:MicrobePopulation schedules:[] {
	init {
		dividing_time <- dividing_time_S;
		carbon_use_efficiency <- carbon_use_efficiency_S;
		minimum_awake_rate <- minimum_awake_rate_S;

		cytosol_mineralization_rate <- 0.0;
		
		Enzymes _min_enzymes;
		create Enzymes with: [
			T_cellulolytic::0.0,
			T_amino::0.0,
			T_P::0.0,
			T_recal::0.0
		] {
			_min_enzymes <- self;
		}
		
		Enzymes _max_enzymes;
		create Enzymes with: [
			T_cellulolytic::2 #gram/ #gram / #d,
			T_amino::0.5 #gram / #gram / #d,
			T_P::0.1 #gram / #gram / #d,
			T_recal::0.2 #gram / #gram / #d
		] {
			_max_enzymes <- self;
		}
		do set_min_max_enzymes(_min_enzymes, _max_enzymes);
	}
}
