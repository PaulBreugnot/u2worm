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
	
	// TODO: what rates?
	float copiotrophe_R_rate <- 0.3;
	float copiotrophe_K_rate <- 0.3;
	float oligotrophe_K_rate <- 0.3;

	float dividing_time_copiotrophe_R <- 1#h;
	float dividing_time_copiotrophe_K <- 24#h;
	float dividing_time_oligotrophe <- 368#h;
}

species Copiotrophe_R parent:MicrobePopulation {
	init {
		bacteria_name <- COPIOTROPHE_R;
		dividing_time <- dividing_time_copiotrophe_R;
		awake_population <- 0.0;
		L_R_rate <- 1.0;
		C_N <- 5.0;
		C_P <- 10.0;
		C <- copiotrophe_R_rate * total_initial_bacteria_weight / length(PoreParticle); // TODO: /3 => car trois différentes espèces de bactéries, hypothèse toutes les bactéries ont le même poids
		N <- C / C_N;
		P <- C / C_P;
		
		cytosol_C <- C; // TODO: cytosol = C?
		cytosol_N <- C / C_N;
		cytosol_P <- C / C_P;
		
		wakeup_factor <- 0.002;

		taux_respiration <- 0.2;
		division_enzyme_rate <- 0.714;	
	}
}
species Copiotrophe_K parent:MicrobePopulation {
	init {
		bacteria_name <- COPIOTROPHE_K;
		dividing_time <- dividing_time_copiotrophe_K;
		awake_population <- 0.0;
		L_R_rate <- 0.2;
		C_N <- 5.0;
		C_P <- 10.0;
		
		C <- copiotrophe_K_rate * total_initial_bacteria_weight / length(PoreParticle); //*0.2;//0.000001#gram;
		N <- C/ C_N;
		P <- C/ C_P;
		
		cytosol_C <- C;
		cytosol_N <- C / C_N;
		cytosol_P <- C / C_P;
		
		wakeup_factor <- 0.018;

		taux_respiration <- 0.4;
		division_enzyme_rate <- 0.5;		
	}
}

species Oligotrophe_K parent:MicrobePopulation {
	init {
		bacteria_name <- OLIGOTROPHE_K;
		dividing_time <- dividing_time_oligotrophe;
		awake_population <- 1.0;
		L_R_rate <- 0.8;
		C_N <- 5.0;
		C_P <- 10.0;
		C <- oligotrophe_K_rate * total_initial_bacteria_weight / length(PoreParticle);   //*0.7;//0.000001#gram;
		N <- C/ C_N;
		P <- C/ C_P;
		
		cytosol_C <- C;
		cytosol_N <- C / C_N;
		cytosol_P <- C / C_P;
		
		wakeup_factor <- 1.0;

		taux_respiration <- 0.7;
		division_enzyme_rate <- 0.5;	
	}
}

