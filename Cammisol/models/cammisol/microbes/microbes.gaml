/**
* Name: microbes
* Based on the internal empty template. 
* Author: pbreugno
* Tags: 
*/

model microbes

import "../pore_particle.gaml"

global {
	list<species<MicrobePopulation>> bacteria_types <- [Y_Strategist, A_Strategist, S_Strategist];

	float dividing_time_Y <- 1#h;
	float dividing_time_A <- 24#h;
	float dividing_time_S <- 368#h;
	
	float carbon_use_efficiency_Y <- 7/10;
	float carbon_use_efficiency_A <- 3/10;
	float carbon_use_efficiency_S <- 5/10;
	
	float minimum_awake_rate_Y <- 0.002;
	float minimum_awake_rate_A <- 0.018;
	float minimum_awake_rate_S <- 1.0;
	
	float init_structural_cytosol_rate <- 1.0;
}

species Y_Strategist parent:MicrobePopulation {
	init {
		float total_C_init <- C;
		dividing_time <- dividing_time_Y;
		carbon_use_efficiency <- carbon_use_efficiency_Y;
		minimum_awake_rate <- minimum_awake_rate_Y;
		
		C <- init_structural_cytosol_rate * total_C_init; // TODO: /3 => car trois différentes espèces de bactéries, hypothèse toutes les bactéries ont le même poids
		N <- C / C_N;
		P <- C / C_P;
		
		cytosol_C <- total_C_init - C; // TODO: cytosol = C?
		cytosol_N <- cytosol_C / C_N;
		cytosol_P <- cytosol_C / C_P;
		
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
//			T_cellulolytic::0.5 #gram/ #gram / #d,
//			T_amino::0.1 #gram / #gram / #d,
//			T_P::0.05 #gram / #gram / #d,
//			T_recal::0.001 #gram / #gram / #d
			T_cellulolytic::1.0 #gram/ #gram / #d,
			T_amino::0.1 #gram / #gram / #d,
			T_P::0.05 #gram / #gram / #d,
			T_recal::0.00001 #gram / #gram / #d
		] {
			_max_enzymes <- self;
		}
		do set_min_max_enzymes(_min_enzymes, _max_enzymes);
	}
}
species A_Strategist parent:MicrobePopulation {
	init {
		float total_C_init <- C;
		dividing_time <- dividing_time_A;
		carbon_use_efficiency <- carbon_use_efficiency_A;
		minimum_awake_rate <- minimum_awake_rate_A;
		
		C <- init_structural_cytosol_rate * total_C_init; //*0.2;//0.000001#gram;
		N <- C/ C_N;
		P <- C/ C_P;
		
		cytosol_C <- total_C_init - C;
		cytosol_N <- cytosol_C / C_N;
		cytosol_P <- cytosol_C / C_P;
		
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
			T_cellulolytic::0.4 #gram/ #gram / #d,
			T_amino::0.08 #gram / #gram / #d,
//			T_P::0.08 #gram / #gram / #d,
			T_P::0.02 #gram / #gram / #d,
			T_recal::0.08 #gram / #gram / #d
		] {
			_max_enzymes <- self;
		}
		do set_min_max_enzymes(_min_enzymes, _max_enzymes);
	}
}

species S_Strategist parent:MicrobePopulation {
	init {
		float total_C_init <- C;
		dividing_time <- dividing_time_S;
		carbon_use_efficiency <- carbon_use_efficiency_S;
		minimum_awake_rate <- minimum_awake_rate_S;
		
		C <- init_structural_cytosol_rate * total_C_init;   //*0.7;//0.000001#gram;
		N <- C/ C_N;
		P <- C/ C_P;
		
		cytosol_C <- total_C_init - C;
		cytosol_N <- cytosol_C / C_N;
		cytosol_P <- cytosol_C / C_P;
		
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
			T_cellulolytic::0.2 #gram/ #gram / #d,
			T_amino::0.05 #gram / #gram / #d,
//			T_P::0.05 #gram / #gram / #d,
			T_P::0.01 #gram / #gram / #d,
			T_recal::0.02 #gram / #gram / #d
		] {
			_max_enzymes <- self;
		}
		do set_min_max_enzymes(_min_enzymes, _max_enzymes);
	}
}
