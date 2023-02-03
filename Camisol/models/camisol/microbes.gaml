/**
* Name: microbes
* Based on the internal empty template. 
* Author: pbreugno
* Tags: 
*/

model microbes

import "pore_particle.gaml"

global {
	string COPIOTROPHE_R <- "Copiotrophe R";
	string COPIOTROPHE_K <- "Copiotrophe K";
	string OLIGOTROPHE_K <- "Oligotrophe K";
	list<string> bacteria_types <- [COPIOTROPHE_R, COPIOTROPHE_K, OLIGOTROPHE_K];
	
	// TODO: what rates?
	float copiotrophe_R_rate <- 0.1;
	float copiotrophe_K_rate <- 0.2;
	float oligotrophe_K_rate <- 0.7;

	float dividing_time_copiotrophe_R <- 1#h;
	float dividing_time_copiotrophe_K <- 24#h;
	float dividing_time_oligotrophe <- 368#h;
	
	float init_structural_cytosol_rate <- 1.0;
}

species Copiotrophe_R parent:MicrobePopulation {
	init {
		float total_C_init <- copiotrophe_R_rate * total_initial_bacteria_weight / length(PoreParticle);
		bacteria_name <- COPIOTROPHE_R;
		dividing_time <- dividing_time_copiotrophe_R;
		awake_population <- 0.0;
		L_R_rate <- 1.0;
		C_N <- 5.0;
		C_P <- 10.0;
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
	}
}
species Copiotrophe_K parent:MicrobePopulation {
	init {
		float total_C_init <- copiotrophe_K_rate * total_initial_bacteria_weight / length(PoreParticle);
		bacteria_name <- COPIOTROPHE_K;
		dividing_time <- dividing_time_copiotrophe_K;
		awake_population <- 0.0;
		L_R_rate <- 0.2;
		C_N <- 5.0;
		C_P <- 10.0;
		
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
	}
}

species Oligotrophe_K parent:MicrobePopulation {
	init {
		float total_C_init <- oligotrophe_K_rate * total_initial_bacteria_weight / length(PoreParticle);
		bacteria_name <- OLIGOTROPHE_K;
		dividing_time <- dividing_time_oligotrophe;
		awake_population <- 1.0;
		L_R_rate <- 0.8;
		C_N <- 5.0;
		C_P <- 10.0;
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
	}
}

