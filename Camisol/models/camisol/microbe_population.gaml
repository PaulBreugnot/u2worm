/**
* Name: microbe_population
* Based on the internal skeleton template. 
* Author: Nicolas Marilleau, Lucas Grosjean, Paul Breugnot
* Tags: 
*/

model microbe_population

import "dam.gaml"

global {
	// 5E8 -> 5E9 bacterie / gramme de sol
	/*
	 * A modifier................................................ 
	 */
	 // TODO: poids total de bactérie par gramme de sol * surface du modèle
	/**
	 * Total bacteria weight in the model.
	 */
	float total_initial_bacteria_weight <-  0.05*1.5#gram/(#cm*#cm)*world.shape.area;
	
	float total_CO2_produced <- 0.0;
	
	// Enzyme objectives
	MaxTotalC max_total_C;
	MaxCN max_CN;
	MaxCP max_CP;

	init {
		create MaxTotalC {
			max_total_C <- self;
		}
		create MaxCN with: (weight: 100.0) {
			max_CN <- self;
		}
		create MaxCP with: (weight: 100.0) {
			max_CP <- self;
		}
	}
}

species MicrobePopulation
{
	string bacteria_name;
	float respiration_rate <- 0.5;
	float division_enzyme_rate <- 0.0;
	
	float dividing_time <- 0.0;
	//taux labile récalcitrant
	float L_R_enzyme_rate; // labile_recalcitrante_rate
	float C_N;
	float C_P;
	
	// quantity of C N P in the colony
	// TODO: unit?
	// TODO: source?
	float C <- 60.0;
	float N <- 7.0;
	float P <- 1.0;
	
	float awake_population <- 0.5;
	float wakeup_factor <- 0.50;
	float C_actif -> {C * awake_population};
	
	// Mon estomac le cytosol 
	float cytosol_C <- 0#gram; 
	float cytosol_N <- 0#gram; 
	float cytosol_P <- 0#gram; 
	
	// ce que je percois
	float perception_C <- 0#gram; 
	float perception_N <- 0#gram;
	float perception_P <- 0#gram;
	
	
	Enzymes min_enzymes;
	Enzymes max_enzymes;
	
	// float P_C <- 8;
	// float P_N <- 12;
	// float P_P <- 12;
	
	/* 
	 R : C * (step/1)  /  (1-0.2)
	 K : C * (step/24)  / (1-0.4)
	 O : C * (step/368) / (1-0.7)s
	 */
	// Quantité de carbone voulue par la colonie
	//float C_assimilation_speed -> {C_actif * (exp(local_step / dividing_time)-1) / (division_enzyme_rate * (1 - respiration_rate))};
	float C_assimilation_speed -> {C_actif * local_step / dividing_time / (division_enzyme_rate * (1 - respiration_rate))};
	
	float CO2_produced <- 0.0;
	
	float C_wanted -> { max([0.0, C_assimilation_speed - cytosol_C]) };
	float N_wanted -> { (C+cytosol_C+C_wanted) / C_N - (cytosol_N + N)};
	float P_wanted -> { (C+cytosol_C+C_wanted) / C_P - (cytosol_P + P)};
	
	init {
		// TODO: do this in microbes
		create Enzymes with: [
			T_cellulolytic::0 #gram / #h,
			T_amino::0 #gram / #h,
			T_P::0 #gram / #h,
			T_recal::0 #gram / #h
		] {
			myself.min_enzymes <- self;
		}
		create Enzymes with: [
			T_cellulolytic::1 #gram/#h,
			T_amino::1 #gram/#h,
			T_P::1 #gram/#h,
			T_recal::1 #gram/#h
		] {
			myself.max_enzymes <- self;
		}
	}
	action respirate(float cytosol_respiration)
	{
		cytosol_C <- cytosol_C  - cytosol_respiration; 
		
		CO2_produced <- CO2_produced + cytosol_respiration;
		total_CO2_produced <- total_CO2_produced + cytosol_respiration;
	}
	
	action growth(float cytosol_division, float total_C_in_pore, float pore_carrying_capacity)
	{
		// quantity de carbone cytosol -> structurel.
		// Unconsumed C/N/P stays in the cytosol
		float new_individual <- C * (local_step/dividing_time) * (1 - total_C_in_pore / pore_carrying_capacity);
		
		//new_individual <- min([new_individual , cytosol_C]);
		new_individual <- min([new_individual , cytosol_division]);
		
		
		float assi_C_N <- cytosol_N * C_N;  //(assimilated_N=0) ? 0:(assimilated_C/assimilated_N);
		float assi_C_P <- cytosol_P * C_P; //(assimilated_P=0) ? 0:(assimilated_C/assimilated_P);
		
		new_individual <- min([new_individual, (assi_C_N) ]);
		new_individual <- min([new_individual, (assi_C_P) ]);
		
		// write bacteria_name + ": " + new_individual;
		
		float transfert_C <- new_individual;
		float transfert_N <- transfert_C / C_N;
		float transfert_P <- transfert_C / C_P;
		
		cytosol_C <- cytosol_C - transfert_C;
		cytosol_N <- cytosol_N - transfert_N;
		cytosol_P <- cytosol_P - transfert_P;
		
		
		C <- C + transfert_C;
		N <- N + transfert_N;
		P <- P + transfert_P;
	}

	action decompose(Dam dam, list<OrganicParticle> particles_to_decompose) {
		EnzymaticActivity enzymes <- produce_enzymes(dam, particles_to_decompose);
		
		float total_C_labile <- sum(particles_to_decompose collect each.C_labile);
		float total_P_labile <- sum(particles_to_decompose collect each.P_labile);
		float total_N_labile <- sum(particles_to_decompose collect each.N_labile);
		
		float total_C_recal <- sum(particles_to_decompose collect each.C_recalcitrant);
		float total_P_recal <- sum(particles_to_decompose collect each.P_recalcitrant);
		
		write "E: " + ((enzymes.X_C_cellulolytic + enzymes.X_C_amino + enzymes.X_C_P)/#gram)
			+ " / " + (total_C_labile / #gram);
					
		ask particles_to_decompose {
			// Amino activity
			float X_C_amino <- enzymes.X_C_amino * C_labile / total_C_labile; // Necessarily < to C_labile
			float X_N_amino <- enzymes.X_N_amino * N_labile / total_N_labile;
			
			// CP activity
			float X_C_P <- enzymes.X_C_P * C_labile / total_C_labile;
			float X_P_labile_to_dom <- enzymes.X_P_labile_to_dom * P_labile / total_P_labile;
			float X_P_labile_to_dim <- enzymes.X_P_labile_to_dim * P_labile / total_P_labile;
			
			// C activity
			float X_C <- enzymes.X_C_cellulolytic * C_labile / total_C_labile;
			
			// Recal activity
			float X_P_recal_to_dim <- enzymes.X_P_recal_to_dim * P_recalcitrant / total_P_recal;
			
			float X_C_recal <- enzymes.X_C_recal * C_recalcitrant / total_C_recal;
			float X_N_recal <- X_C_recal / (C_recalcitrant / N_recalcitrant);
			float X_P_recal <- X_C_recal / (C_recalcitrant / P_recalcitrant);
			
			// Commit fluxes
			C_labile <- C_labile - X_C_amino - X_C_P - X_C + X_C_recal;
			N_labile <- N_labile - X_N_amino + X_N_recal;
			P_labile <- P_labile - X_P_labile_to_dom - X_P_labile_to_dim + X_P_recal - X_P_recal_to_dim;
			
			C_recalcitrant <- C_recalcitrant - X_C_recal;
			N_recalcitrant <- N_recalcitrant - X_N_recal;
			P_recalcitrant <- P_recalcitrant - X_P_recal;
			
			ask dam {
				dim[1] <- dim[1] + X_P_recal_to_dim + X_P_labile_to_dim;
				
				dom[0] <- dom[0] + X_N_amino;
				dom[1] <- dom[1] + X_P_labile_to_dom;
				dom[2] <- dom[2] + X_C_amino + X_C_P + X_C;
			}
		}
	}
	
	EnzymaticActivity produce_enzymes(Dam dam, list<OrganicParticle> particles_to_decompose) {
		EnzymaticActivity enzymes;
		create EnzymaticActivityProblem with: [
			C_N::C_N,
			C_P::C_P,
			C_microbes::C_actif,
			dt::local_step,
			total_C_labile::sum(particles_to_decompose collect each.C_labile),
			total_N_labile::sum(particles_to_decompose collect each.N_labile),
			total_P_labile::sum(particles_to_decompose collect each.P_labile),
			total_C_recal::sum(particles_to_decompose collect each.C_recalcitrant),
			total_N_recal::sum(particles_to_decompose collect each.N_recalcitrant),
			total_P_recal::sum(particles_to_decompose collect each.P_recalcitrant),
			C_DOM::dam.dom[2],
			N_DOM::dam.dom[0],
			P_DOM::dam.dom[1],
			P_DIM::dam.dim[1],
			N_DIM::dam.dim[0],
			min_enzymes::min_enzymes,
			max_enzymes::max_enzymes
		] {
			create SimulatedAnnealing with:[
				problem::self,
				objectives::[max_total_C, max_CN, max_CP],
				max_N::1000,
				epsilon::1e-3
			] {
				do optimize;
				enzymes <- s.enzymatic_activity;
				ask s {
					ask enzymes {
						do die;
					}
					do die;
				}
				do die;
			}
			do die;
		}
		return enzymes;
	}
	
	action life(Dam dam, list<OrganicParticle> accessible_organics, float total_C_in_pore, float pore_carrying_capacity)
	{
		CO2_produced <- 0.0;
		
		if(C_wanted = 0) {
			awake_population <- 0.0;
		} else {
			awake_population <- min([(perception_C / C_wanted), 1]);
		}
		
		awake_population <- max([awake_population,wakeup_factor]);
	
		cytosol_C <- cytosol_C + perception_C;
		cytosol_N <- cytosol_N + perception_N;
		cytosol_P <- cytosol_P + perception_P;
		
		// Breathed carbon quantity from cytosol, rejected as CO2
		float cytosol_respiration <- respiration_rate * cytosol_C ;
		// Part of the left cytosol C dedicated to cell division
		float cytosol_division <- (cytosol_C-cytosol_respiration) * division_enzyme_rate;
		// Left cytosol C is used to produce enzymes
		float cytosol_enzyme <- cytosol_C - cytosol_respiration - cytosol_division;
		
		do respirate(cytosol_respiration);
		do growth(cytosol_division, total_C_in_pore, pore_carrying_capacity);
		do decompose(dam, accessible_organics);
		
		perception_C <- 0.0;
		perception_N <- 0.0;
		perception_P <- 0.0;
	}
}



experiment PopulationDecompose type: gui {
	float T_cellulolytic;
	float T_amino;
	float T_P;
	float T_recal;
	float C_microbes;
	float C_dom;
	float C_N_microbes;
	float C_P_microbes;
	
	float C_labile_init;
	float C_labile_final;
	float C_N_labile_init;
	float C_N_labile_final;
	float C_P_labile_init;
	float C_P_labile_final;
	
	float C_N_dom_init;
	float C_N_dom_final;
	float C_P_dom_init;
	float C_P_dom_final;
	
	float C_recal_init;
	float C_recal_final;
	float CP_recal <- 15.0;
	float CN_recal <- 20.0;
	
	int steps;
	
	parameter "Max cellulolytic (g/h)" category: "Constants" var: T_cellulolytic init: 1.0;
	parameter "Max amino (g/h)" category: "Constants" var: T_amino init: 1.0;
	parameter "Max P (g/h)" category: "Constants" var: T_P init: 1.0;
	parameter "Max recalcitrant (g/h)" category: "Constants" var: T_recal init: 0.0;
	parameter "Amino CN" category: "Constants" var: amino_CN init: 6.0;
	parameter "C microbes" category: "Constants" var: C_microbes init: 200#gram;
	parameter "C dom" category: "Constants" var: C_dom init: 1#gram;
	parameter "Microbe population's C/N" category: "Constants" var:C_N_microbes init:15.0;
	parameter "Microbe population's C/P" category: "Constants" var:C_P_microbes init:20.0;
	
	parameter "Initial C (labile)" category: "Labile OM" var: C_labile_init init: 10#gram;
	parameter "Final C (labile)" category: "Labile OM" var: C_labile_final init: 10#gram;
	parameter "Initial C/N (labile)" category: "Labile OM" var: C_N_labile_init init: 20.0 ;
	parameter "Final C/N (labile)" category: "Labile OM" var:C_N_labile_final init:20.0;
	parameter "Initial C/P (labile)" category: "Labile OM" var:C_P_labile_init init:20.0;
	parameter "Final C/P (labile)" category: "Labile OM" var:C_P_labile_final init:20.0;
	
	parameter "Initial C/N (DOM)" category: "DOM" var: C_N_dom_init init: 5.0;
	parameter "Final C/N (DOM)" category: "DOM" var:C_N_dom_final init:20.0;
	parameter "Initial C/P (DOM)" category: "DOM" var:C_P_dom_init init:20.0;
	parameter "Final C/P (DOM)" category: "DOM" var:C_P_dom_final init:20.0;
	
	parameter "Initial C (recalcitrant)" category: "Recalcitrant OM" var: C_recal_init init: 100#gram;
	parameter "Final C (recalcitrant)" category: "Recalcitrant OM" var: C_recal_final init: 100#gram;
	
	parameter "Count of steps" category: "Experiment" var: steps init: 1000;
	
	PopulationDecompose exp_output <- self;
	OrganicParticle organic_particle;
	Dam dam;
	MicrobePopulation microbe_population;
	
	init {
		create OrganicParticle with: (
			C_labile: init_C_labile(),
			N_labile: init_N_labile(),
			P_labile: init_P_labile(),
			C_recalcitrant: init_C_recal(),
			N_recalcitrant: init_C_recal() / CN_recal,
			P_recalcitrant: init_C_recal() / CP_recal
		) {
			myself.organic_particle <- self;
		}
		
		float N_dom <- C_dom/(C_N_dom_init + cycle * (C_N_dom_final - C_N_dom_init)/steps);
		float P_dom <- C_dom/(C_P_dom_init + cycle * (C_P_dom_final - C_P_dom_init)/steps);
		
		create Dam with: (
			dom: [N_dom, P_dom, C_dom],
			dim: [0.0, 0.0]
		) {
			myself.dam <- self;
		}
		
		PopulationDecompose exp <- self;
		create MicrobePopulation with: (
			C: C_microbes,
			C_N: C_N_microbes,
			C_P: C_P_microbes,
			awake_population: 1.0
			) {
				ask max_enzymes {
					T_cellulolytic <- exp.T_cellulolytic #gram/#h;
					T_amino <- exp.T_amino #gram/#h;
					T_P <- exp.T_P #gram/#h;
					T_recal <- exp.T_recal #gram/#h;
				}
				exp.microbe_population <- self;
		}
	}
	
	float init_C_labile {
		return C_labile_init + cycle * (C_labile_final - C_labile_init)/steps;
	}
	
	float init_N_labile {
		return init_C_labile()/(C_N_labile_init + cycle * (C_N_labile_final - C_N_labile_init)/steps);
	}
		
	float init_P_labile {
		return init_C_labile()/(C_P_labile_init + cycle * (C_P_labile_final - C_P_labile_init)/steps);
	}
	
	float init_C_recal {
		return C_recal_init + cycle * (C_recal_final - C_recal_init)/steps;
	}
	
	float init_N_dom {
		return C_dom/(C_N_dom_init + cycle * (C_N_dom_final - C_N_dom_init)/steps);
	}
	
	float init_P_dom {
		return C_dom/(C_P_dom_init + cycle * (C_P_dom_final - C_P_dom_init)/steps);
	}
	
	reflex {
		ask organic_particle {
			C_labile <- myself.init_C_labile();
			N_labile <- myself.init_N_labile();
			P_labile <- myself.init_P_labile();
			C_recalcitrant <- myself.init_C_recal();
			N_recalcitrant <- myself.init_C_recal() / myself.CN_recal;
			P_recalcitrant <- myself.init_C_recal() / myself.CP_recal;
		}
		
		ask dam {
			dom[0] <- myself.init_N_dom();
			dom[1] <- myself.init_P_dom();
			dom[2] <- myself.C_dom;
			dim[0] <- 0.0;
			dim[1] <- 0.0;
		}
		
		PopulationDecompose exp <- self;
		
		ask microbe_population {
			do decompose(exp.dam, [exp.organic_particle]);
		}
		write "Recal: " + (init_C_recal() - organic_particle.C_recalcitrant);
	}
	
	output {
		display "OM" type:2d {
			chart "OM" type:series style:line {
				data "C (labile)" value: cycle = 0 ? 0.0 : (exp_output.init_C_labile() - exp_output.organic_particle.C_labile)/#gram marker:false;
				data "N (labile)" value: cycle = 0 ? 0.0 : (exp_output.init_N_labile() - exp_output.organic_particle.N_labile)/#gram marker:false;
				data "P (labile)" value: cycle = 0 ? 0.0 : (exp_output.init_P_labile() - exp_output.organic_particle.P_labile)/#gram marker:false;
				data "C (recal)" value: cycle = 0 ? 0.0 : (exp_output.init_C_recal() - exp_output.organic_particle.C_recalcitrant)/#gram marker:false;
				data "N (recal)" value: cycle = 0 ? 0.0 : (exp_output.init_C_recal() / CN_recal - exp_output.organic_particle.N_recalcitrant)/#gram marker:false;
				data "P (recal)" value: cycle = 0 ? 0.0 : (exp_output.init_C_recal() / CP_recal - exp_output.organic_particle.P_recalcitrant)/#gram marker:false;
			}
		}
		display "C/N" type:2d {
			chart "C/N" type:series style:line {
				data "dom C/N" value: exp_output.C_dom / exp_output.init_N_dom() marker:false;
				data "available C/N" value:
				cycle = 0 ? 0.0 : 
				(exp_output.dam.available_N() > 0 ?
					exp_output.dam.available_C()/exp_output.dam.available_N() : 0.0)
				marker:false;
				data "C/N" value: exp_output.C_N_microbes marker:false;
			}
		}
		display "C/P" type:2d {
			chart "C/P" type:series style:line {
				data "dom C/P" value: exp_output.C_dom / exp_output.init_P_dom() marker:false;
				data "available C/P" value: cycle = 0 ? 0.0 :
					(exp_output.dam.available_P() > 0 ?
					exp_output.dam.available_C() / exp_output.dam.available_P() : 0.0)
				marker:false;
				data "C/P" value: exp_output.C_P_microbes marker:false;
			}
		}
	}
}