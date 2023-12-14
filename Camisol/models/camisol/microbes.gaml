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
		awake_population <- wakeup_factor;
		L_R_enzyme_rate <- 1.0;
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
		awake_population <- wakeup_factor;
		L_R_enzyme_rate <- 0.2;
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
		awake_population <- wakeup_factor;
		L_R_enzyme_rate <- 0.8;
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

experiment TestEnzymes type: gui {
	init {
		float carbone_concentration_in_dam <- (729.0#gram * 10^-6)/#gram;
		float azote_concentration_in_dam <- (60.0#gram * 10^-6)/#gram;
		float azote_concentration_in_dim <- (4.74#gram * 10^-6)/#gram;
		float phosphore_concentration_in_dam <- (400.0#gram * 10^-6)/#gram;
		float phosphore_concentration_in_dim <- (1.43#gram * 10^-6)/#gram;
		float model_surface <- 1#cm*1#cm;
		float model_weight <- 1.17#gram * model_surface;
		float C <- 0.02157#gram/(#cm*#cm);
		float N <- 0.00132#gram/(#cm*#cm);
		float P <- 0.00077#gram/(#cm*#cm);
		
		create Dam with: [
				dom: [
					azote_concentration_in_dam * model_weight,
					phosphore_concentration_in_dam * model_weight,
					carbone_concentration_in_dam * model_weight
				],
				dim: [
					azote_concentration_in_dim * model_weight,
					phosphore_concentration_in_dim * model_weight
				]
			] {
			}
		create PoreParticle with: [carrying_capacity::10*total_initial_bacteria_weight] {
			ask Dam {
				myself.dam <- self;
			}
		}
		
		create Copiotrophe_R {
			ask PoreParticle {
				add myself to: self.populations;
			}
		}
		create OrganicParticle with: [
			C_labile: (C/2)*model_surface,
			N_labile: (N/2)*model_surface,
			P_labile: (P/2)*model_surface,
			C_recalcitrant: (C/2)*model_surface,
			N_recalcitrant: (N/2)*model_surface,
			P_recalcitrant: (P/2)*model_surface
		] {
			ask PoreParticle {
				organic_particle <- myself;
				add myself to: accessible_organics;
			}
		}
	}
	
	reflex {
		ask Dam {
			write dom;
			write dim;
		}
	}
	
	/** Insert here the definition of the input and output of the model */
	output {
		display "dam" type: java2D {
			chart "dam" type:series {
				data "N (dom)" value: (sum(Dam collect each.dom[0])) style:spline marker:false thickness:3;
				data "P (dom)" value: (sum(Dam collect each.dom[1])) style:spline marker:false thickness:3;
				data "C (dom)" value: (sum(Dam collect each.dom[2])) style:spline marker:false thickness:3;
				data "N (dim)" value: (sum(Dam collect each.dim[0])) style:spline marker:false thickness:3;
				data "P (dim)" value: (sum(Dam collect each.dim[1])) style:spline marker:false thickness:3;
			}
		}
		
		display "organics" type:java2D {
			chart "Organics composition" type:series {
				data "C_labile" value: (sum(OrganicParticle collect each.C_labile)) style:spline marker:false thickness:3;
				data "N_labile" value: (sum(OrganicParticle collect each.N_labile)) style:spline marker:false thickness:3;
				data "P_labile" value: (sum(OrganicParticle collect each.P_labile)) style:spline marker:false thickness:3;
				data "C_recalcitrant" value: (sum(OrganicParticle collect each.C_recalcitrant)) style:spline marker:false thickness:3;
				data "N_recalcitrant" value: (sum(OrganicParticle collect each.N_recalcitrant)) style:spline marker:false thickness:3;
				data "P_recalcitrant" value: (sum(OrganicParticle collect each.P_recalcitrant)) style:spline marker:false thickness:3;
			}
		}
		
		display "populations" type:java2D {
			chart "Bacteria populations" type:series {
				data "Copiotrophe R" value: (sum(Copiotrophe_R collect each.C)) style:spline color: #red marker:false thickness:3;
				data "Copiotrophe K" value: (sum(Copiotrophe_K collect each.C)) style:spline color: #green marker:false thickness:3;
				data "Oligotrophe K" value: (sum(Oligotrophe_K collect each.C)) style:spline color: #blue marker:false thickness:3;
			}
		}
	}
}
