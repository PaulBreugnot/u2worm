/**
* Name: microbe_population
* Based on the internal skeleton template. 
* Author: Nicolas Marilleau, Lucas Grosjean, Paul Breugnot
* Tags: 
*/

model microbe_population

import "dam.gaml"

global {
	
	float total_CO2_produced <- 0.0;
	float enzymes_optimization_period <- 1#d;
	
	// Enzyme objectives
	MaxLabileC max_labile_C;
	MaxRecalC max_recal_C;
	MaxN max_N;
	MaxP max_P;
	MaxCN max_CN;
	MaxCP max_CP;

	init {
		create MaxLabileC with: (weight: 5.0) {
			max_labile_C <- self;
		}
		create MaxRecalC with: (weight: 1.0) {
			max_recal_C <- self;
		}
		create MaxN with: (weight: 5.0) {
			max_N <- self;
		}
		create MaxP with: (weight: 5.0) {
			max_P <- self;
		}
		create MaxCN with: (weight: 10.0) {
			max_CN <- self;
		}
		create MaxCP with: (weight: 10.0) {
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
	float C_N <- 10.0;
	float C_P <- 17.0;
	
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
	
	Enzymes enzymes;
	Enzymes min_enzymes;
	Enzymes max_enzymes;
	EnzymaticActivityProblem enzymatic_activity_problem;
	
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
		create Enzymes {
			myself.enzymes <- self;
		}
		create EnzymaticActivityProblem with: [
			dt::enzymes_optimization_period
		] {
			myself.enzymatic_activity_problem <- self;
		}
	}
	
	action set_min_max_enzymes(Enzymes _min_enzymes, Enzymes _max_enzymes) {
		self.min_enzymes <- _min_enzymes;
		self.max_enzymes <- _max_enzymes;
		ask enzymatic_activity_problem {
			min_enzymes <- _min_enzymes;
			max_enzymes <- _max_enzymes;
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
	
	action optimize_enzymes(Dam dam, list<OrganicParticle> particles_to_decompose) {
		MicrobePopulation population <- self;
		ask enzymatic_activity_problem {
			C_N <- myself.C_N;
			C_P <- myself.C_P;
			C_microbes <- myself.C_actif;
			C_labile <- sum(particles_to_decompose collect each.C_labile);
			N_labile <- sum(particles_to_decompose collect each.N_labile);
			P_labile <- sum(particles_to_decompose collect each.P_labile);
			C_recal <- sum(particles_to_decompose collect each.C_recalcitrant);
			N_recal <- sum(particles_to_decompose collect each.N_recalcitrant);
			P_recal <- sum(particles_to_decompose collect each.P_recalcitrant);
			C_DOM <- dam.dom[2];
			N_DOM <- dam.dom[0];
			P_DOM <- dam.dom[1];
			P_DIM <- dam.dim[1];
			N_DIM <- dam.dim[0];
		}
		create SimulatedAnnealing with:[
				problem::self.enzymatic_activity_problem,
				objectives::(max_enzymes.T_recal > 0 ? [max_labile_C, max_recal_C, max_N, max_P, max_CN, max_CP] : [max_labile_C, max_N, max_P, max_CN, max_CP]),
				N::1000,
				epsilon::0.0
			] {
			do optimize;
			ask population.enzymes {
				do die;
			}
			population.enzymes <- s.enzymes;
			ask s {
				ask weighted_enzymes {
					do die;
				}
				ask decomposition {
					do die;
				}
				do die;
			}
			do die;
		}
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
