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
		dividing_time <- dividing_time_A;
		carbon_use_efficiency <- carbon_use_efficiency_A;
		minimum_awake_rate <- minimum_awake_rate_A;
		
		C <- init_structural_cytosol_rate * total_C_init; //*0.2;//0.000001#gram;
		N <- C/ C_N;
		P <- C/ C_P;
		
		cytosol_C <- total_C_init - C;
		cytosol_N <- cytosol_C / C_N;
		cytosol_P <- cytosol_C / C_P;
		


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

//experiment MicrobesTestBase {
//	float init_Y_C <- 1.0;
//	float init_A_C <- 1.0;
//	float init_S_C <- 1.0;
//	float nutrient_rate <- 1.0;
//	float carrying_capacity <- 10#gram;
//	
//	bool enable_Y_Strategist <- true on_change: change_species;
//	bool enable_A_Strategist <- true on_change: change_species;
//	bool enable_S_Strategist <- true on_change: change_species;
//	
//	parameter "Enable Y Strategist" var:enable_Y_Strategist category: "Y strategist";
//	parameter "Y C init" var:init_Y_C category:"Y strategist";
//	parameter "Y dividing time" var:dividing_time_Y category:"Y strategist";
//	parameter "Y CUE" var:carbon_use_efficiency_Y category:"Y strategist";
//	
//	parameter "Enable A Strategist" var:enable_A_Strategist category: "A strategist";
//	parameter "A C init" var:init_A_C category:"A strategist";
//	parameter "A dividing time" var:dividing_time_A category:"A strategist";
//	parameter "A CUE" var:carbon_use_efficiency_A category:"A strategist";
//	
//	parameter "Enable S Strategist" var:enable_S_Strategist category: "S strategist";
//	parameter "S C init" var:init_S_C category:"S strategist";
//	parameter "S dividing time" var:dividing_time_S category:"S strategist";
//	parameter "S CUE" var:carbon_use_efficiency_S category:"S strategist";
//	
//	parameter "Rate of requested nutrients in the DAM" var:nutrient_rate;
//	parameter "Carrying capacity" var:carrying_capacity;
//	
//	
//	list<species<MicrobePopulation>> microbe_species <- init_species();
//	map<species<MicrobePopulation>, MicrobePopulation> populations;
//		
//	float C_N <- 10.0;
//	float C_P <- 17.0;
//	float total_C_dom <- 0.0;
//	
//	map<species<MicrobePopulation>, float> init_C;
//	map<species<MicrobePopulation>, float> C;
//	map<species<MicrobePopulation>, float> C_cytosol;
//	map<species<MicrobePopulation>, float> N;
//	map<species<MicrobePopulation>, float> P;
//	map<species<MicrobePopulation>, float> awake;
//	
//	action change_species {
//		microbe_species <- init_species();
//	}
//	
//	action init_species {
//		list<species<MicrobePopulation>> _microbe_species;
//		init_C <- [];
//		if enable_Y_Strategist {
//			add Y_Strategist to: _microbe_species;
//			init_C[Y_Strategist] <- init_Y_C#gram;
//		}
//		if enable_A_Strategist {
//			add A_Strategist to: _microbe_species;
//			init_C[A_Strategist] <- init_A_C#gram;
//		}
//		if enable_S_Strategist {
//			add S_Strategist to: _microbe_species;
//			init_C[S_Strategist] <- init_S_C#gram;
//		}
//		C <- [];
//		C_cytosol <- [];
//		N <- [];
//		P <- [];
//		awake <- [];
//		loop s over: _microbe_species {
//			C[s] <- init_C[s];
//			C_cytosol[s] <- 0.0;
//			N[s] <- init_C[s]/C_N;
//			P[s] <- init_C[s]/C_P;
//			awake[s] <- 1.0;
//		}
//		return _microbe_species;
//	}
//	
//	init {
//		loop s over: microbe_species {
//			create s with: [C::init_C[s], C_N::C_N, C_P::C_P] {
//				myself.populations[s] <- self;
//			}
//		}
//	}
//	
//	action update_output_data {
//		loop s over:microbe_species {
//			C[s] <- populations[s].C;
//			C_cytosol[s] <- populations[s].cytosol_C;
//			N[s] <- populations[s].N;
//			P[s] <- populations[s].P;
//			awake[s] <- populations[s].awake_population;
//		}
//	}
//	
//	output {
//		display "C" {
//			chart "Microbes populations" {
//				if enable_Y_Strategist {
//					data "Y strategists (g)" value:C[Y_Strategist]/#gram marker: false;
//				}
//				if enable_Y_Strategist {
//					data "A strategists (g)" value:C[A_Strategist]/#gram marker: false;
//				}
//				if enable_Y_Strategist {
//					data "S strategists (g)" value:C[S_Strategist]/#gram marker: false;
//				}
//			}
//		}
//		
//		display "C cytosol" {
//			chart "Microbes populations cytosol" {
//				if enable_Y_Strategist {
//					data "Y strategists (g)" value:C_cytosol[Y_Strategist]/#gram marker: false;
//				}
//				if enable_A_Strategist {
//					data "A strategists (g)" value:C_cytosol[A_Strategist]/#gram marker: false;
//				}
//				if enable_S_Strategist {
//					data "S strategists (g)" value:C_cytosol[S_Strategist]/#gram marker: false;
//				}
//				data "CO2" value: microbe_CO2_emissions marker:false;
//			}
//			
//		}
//		
//		display "Awake rate" {
//			chart "Awake rate" {
//				if enable_Y_Strategist {
//					data "Y strategist (%)" value:awake[Y_Strategist] marker: false;
//				}
//				if enable_A_Strategist {
//					data "A strategist (%)" value:awake[A_Strategist] marker: false;
//				}
//				if enable_S_Strategist {
//					data "S strategist (%)" value:awake[S_Strategist] marker: false;
//				}
//			}
//		}
//		
//		display "C conservation" {
//			chart "C conservation" {
//				data "Conservation"
//					value:(total_C_dom + sum(init_C.values))
//					/ (sum(Dam collect each.dom[2]) + sum(C_cytosol.values) + sum(C.values) + microbe_CO2_emissions)
//					marker: false;
//			}
//		}
//		display "C/N" {
//			chart "C/N" {
//				if enable_Y_Strategist {
//					data "Y strategists (g)" value:C[Y_Strategist]/N[Y_Strategist] marker: false;
//				}
//				if enable_A_Strategist {
//					data "A strategists (g)" value:C[A_Strategist]/N[A_Strategist] marker: false;
//				}
//				if enable_S_Strategist {
//					data "S strategists (g)" value:C[S_Strategist]/N[S_Strategist] marker: false;
//				}
//			}
//		}
//		display "C/P" {
//			chart "C/P" {
//				if enable_Y_Strategist {
//					data "Y strategists (g)" value:C[Y_Strategist]/P[Y_Strategist] marker: false;
//				}
//				if enable_A_Strategist {
//					data "A strategists (g)" value:C[A_Strategist]/P[A_Strategist] marker: false;
//				}
//				if enable_S_Strategist {
//					data "S strategists (g)" value:C[S_Strategist]/P[S_Strategist] marker: false;
//				}
//			}
//		}
//	}
//}
//
//experiment IndividualMicrobesGrowth parent:MicrobesTestBase {
//	map<species<MicrobePopulation>, Dam> dams;
//	
//	init {
////		do init_populations;
//		loop s over: microbe_species {
//			create Dam {
//				myself.dams[s] <- self;
//			}
//			
//			create PoreParticle with: [
//				dam: dams[s],
//				carrying_capacity: carrying_capacity
//			] {
//				// Each population is associated to its own PoreParticle and Dam
//				add myself.populations[s] to: populations;
//			}
//		}
//	}
//	
//	reflex {
//		ask populations.values {
//			do update(C, myself.carrying_capacity);
//		}
//		do feed_dams;
//		
//		ask PoreParticle {
//			do microbe_life;
//		}
//		do update_output_data;
//	}
//	
//	action feed_dams {
//		loop s over: microbe_species {
//			ask dams[s] {
//				myself.total_C_dom <- myself.total_C_dom + myself.nutrient_rate * myself.populations[s].requested_C;
//				// Each dam is fed with the nutrients requested by the population
//				dom[2] <- dom[2] + myself.nutrient_rate * myself.populations[s].requested_C;
//				dom[0] <- dom[0] + myself.nutrient_rate * myself.populations[s].requested_N;
//				dom[1] <- dom[1] + myself.nutrient_rate * myself.populations[s].requested_P;
//			}
//		}
//	}
//}
//
//
//experiment CollectiveMicrobesGrowth parent:MicrobesTestBase {
//	Dam dam;
//	init {
////		do init_populations;
//		create Dam {
//			myself.dam <- self;
//		}
//		
//		create PoreParticle with: [
//			dam: dam,
//			carrying_capacity: carrying_capacity
//		] {
//			// The three populations are added to the same PoreParticle
//			loop s over:myself.microbe_species {
//				add myself.populations[s] to: populations;
//			}
//		}
//	}
//	reflex {
//		ask populations.values {
//			do update(sum(myself.populations.values collect each.C), myself.carrying_capacity);
//		}
//		do feed_dam;
//		do update_output_data;
//		ask PoreParticle {
//			do microbe_life;
//		}
//	}
//	
//	action feed_dam {
//		ask Dam {
//			myself.total_C_dom <- myself.total_C_dom + myself.nutrient_rate * sum(myself.microbe_species collect myself.populations[each].requested_C);
//			// The dam is fed with the nutrients requested by all populations
//			dom[2] <- myself.nutrient_rate * sum(myself.microbe_species collect myself.populations[each].requested_C);
//			dom[0] <- myself.nutrient_rate * sum(myself.microbe_species collect myself.populations[each].requested_N);
//			dom[1] <- myself.nutrient_rate * sum(myself.microbe_species collect myself.populations[each].requested_P);
//		}
//	}
//}
