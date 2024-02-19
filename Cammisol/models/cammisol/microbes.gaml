/**
* Name: microbes
* Based on the internal empty template. 
* Author: pbreugno
* Tags: 
*/

model microbes

import "pore_particle.gaml"

global {
	string Y_STRATEGIST <- "Y Strategist";
	string A_STRATEGIST <- "A Strategist";
	string S_STRATEGIST <- "S Strategist";
	list<string> bacteria_types <- [Y_STRATEGIST, A_STRATEGIST, S_STRATEGIST];

	float dividing_time_copiotrophe_R <- 1#h;
	float dividing_time_copiotrophe_K <- 24#h;
	float dividing_time_oligotrophe <- 368#h;
	
	float init_structural_cytosol_rate <- 1.0;
}

species Y_Strategist parent:MicrobePopulation {
	init {
		float total_C_init <- C;
		bacteria_name <- Y_STRATEGIST;
		dividing_time <- dividing_time_copiotrophe_R;
		awake_population <- wakeup_factor;
		L_R_enzyme_rate <- 1.0;
		C <- init_structural_cytosol_rate * total_C_init; // TODO: /3 => car trois différentes espèces de bactéries, hypothèse toutes les bactéries ont le même poids
		N <- C / C_N;
		P <- C / C_P;
		
		cytosol_C <- total_C_init - C; // TODO: cytosol = C?
		cytosol_N <- cytosol_C / C_N;
		cytosol_P <- cytosol_C / C_P;
		
		wakeup_factor <- 0.002;

//		taux_respiration <- 0.2;
//		division_enzyme_rate <- 0.714;
		// 3/10 C used for respiration
		respiration_rate <- 3/10;
		// Among the 7 C left, 5 are used for division (and the 2 remaining are used for enzyme production)
		division_enzyme_rate <- 5/7;
		
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
			T_cellulolytic::0.5 #gram/ #gram / #d,
			T_amino::0.1 #gram / #gram / #d,
//			T_P::0.1 #gram / #gram / #d,
			T_P::0.05 #gram / #gram / #d,
			T_recal::0.001 #gram / #gram / #d
		] {
			_max_enzymes <- self;
		}
		do set_min_max_enzymes(_min_enzymes, _max_enzymes);
	}
}
species A_Strategist parent:MicrobePopulation {
	init {
		float total_C_init <- C;
		bacteria_name <- A_STRATEGIST;
		dividing_time <- dividing_time_copiotrophe_K;
		awake_population <- wakeup_factor;
		L_R_enzyme_rate <- 0.2;
		
		C <- init_structural_cytosol_rate * total_C_init; //*0.2;//0.000001#gram;
		N <- C/ C_N;
		P <- C/ C_P;
		
		cytosol_C <- total_C_init - C;
		cytosol_N <- cytosol_C / C_N;
		cytosol_P <- cytosol_C / C_P;
		
		wakeup_factor <- 0.018;

//		taux_respiration <- 0.4;
//		division_enzyme_rate <- 0.5;	
		// 7/10 C used for respiration
		respiration_rate <- 7/10;
		// Among the 3 C left, 1.5 are used for division (and the 1.5 remaining are used for enzyme production)
		division_enzyme_rate <- 1.5/3;

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
		bacteria_name <- S_STRATEGIST;
		dividing_time <- dividing_time_oligotrophe;
		awake_population <- wakeup_factor;
		L_R_enzyme_rate <- 0.8;
		C <- init_structural_cytosol_rate * total_C_init;   //*0.7;//0.000001#gram;
		N <- C/ C_N;
		P <- C/ C_P;
		
		cytosol_C <- total_C_init - C;
		cytosol_N <- cytosol_C / C_N;
		cytosol_P <- cytosol_C / C_P;
		
		wakeup_factor <- 1.0;

//		taux_respiration <- 0.7;
//		division_enzyme_rate <- 0.5;
		// 5/10 C used for respiration	
		respiration_rate <- 5/10;
		// Among the 5 C left, 2.5 are used for division (and the 2.5 remaining are used for enzyme production)
		division_enzyme_rate <- 2.5/5;
		
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
