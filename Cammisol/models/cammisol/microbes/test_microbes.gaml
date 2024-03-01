/**
* Name: testmicrobes
* This wizard creates a new test experiment 
* Author: pbreugno
* Tags: 
*/

model test_microbes

import "microbes.gaml"

experiment MicrobesTestBase {
	float init_Y_C <- 1.0 on_change: change_species;
	float init_A_C <- 1.0 on_change: change_species;
	float init_S_C <- 1.0 on_change: change_species;
	float carrying_capacity <- 10#gram;
	
	bool enable_Y_Strategist <- true on_change: change_species;
	bool enable_A_Strategist <- true on_change: change_species;
	bool enable_S_Strategist <- true on_change: change_species;
	
	bool only_one_population <- false;
	
	parameter "Enable Y Strategist" var:enable_Y_Strategist category: "Y strategist";
	parameter "Y C init (g)" var:init_Y_C category:"Y strategist";
	parameter "Y dividing time" var:dividing_time_Y category:"Y strategist";
	parameter "Y CUE" var:carbon_use_efficiency_Y category:"Y strategist";
	
	parameter "Enable A Strategist" var:enable_A_Strategist category: "A strategist";
	parameter "A C init (g)" var:init_A_C category:"A strategist";
	parameter "A dividing time" var:dividing_time_A category:"A strategist";
	parameter "A CUE" var:carbon_use_efficiency_A category:"A strategist";
	
	parameter "Enable S Strategist" var:enable_S_Strategist category: "S strategist";
	parameter "S C init (g)" var:init_S_C category:"S strategist";
	parameter "S dividing time" var:dividing_time_S category:"S strategist";
	parameter "S CUE" var:carbon_use_efficiency_S category:"S strategist";
	
	parameter "Carrying capacity" var:carrying_capacity;
	
	
	list<species<MicrobePopulation>> microbe_species <- init_species();
	map<species<MicrobePopulation>, MicrobePopulation> populations;
		
	float C_N <- 10.0;
	float C_P <- 17.0;
	float total_C_dom <- 0.0;
	float total_N_dom <- 0.0;
	float total_P_dom <- 0.0;
	float total_C_labile <- 0.0;
	float total_N_labile <- 0.0;
	float total_P_labile <- 0.0;
	float total_C_recal <- 0.0;
	float total_N_recal <- 0.0;
	float total_P_recal <- 0.0;
	
	map<species<MicrobePopulation>, float> init_C;
	map<species<MicrobePopulation>, float> C;
	map<species<MicrobePopulation>, float> N;
	map<species<MicrobePopulation>, float> P;
	map<species<MicrobePopulation>, float> C_cytosol;
	map<species<MicrobePopulation>, float> N_cytosol;
	map<species<MicrobePopulation>, float> P_cytosol;
	map<species<MicrobePopulation>, float> awake;
	
	action change_species {
		microbe_species <- init_species();
	}
	
	action init_species {
		list<species<MicrobePopulation>> _microbe_species;
		init_C <- [];
		if enable_Y_Strategist {
			add Y_Strategist to: _microbe_species;
			init_C[Y_Strategist] <- init_Y_C#gram;
		}
		if enable_A_Strategist {
			add A_Strategist to: _microbe_species;
			init_C[A_Strategist] <- init_A_C#gram;
		}
		if enable_S_Strategist {
			add S_Strategist to: _microbe_species;
			init_C[S_Strategist] <- init_S_C#gram;
		}
		C <- [];
		C_cytosol <- [];
		N_cytosol <- [];
		P_cytosol <- [];
		N <- [];
		P <- [];
		awake <- [];
		loop s over: _microbe_species {
			C[s] <- init_C[s];
			C_cytosol[s] <- 0.0;
			N_cytosol[s] <- 0.0;
			P_cytosol[s] <- 0.0;
			N[s] <- init_C[s]/C_N;
			P[s] <- init_C[s]/C_P;
			awake[s] <- 1.0;
		}
		only_one_population <- int(enable_Y_Strategist) + int(enable_A_Strategist) + int(enable_S_Strategist) = 1;
		return _microbe_species;
	}
	
	init {
		loop s over: microbe_species {
			create s with: [C::init_C[s], C_N::C_N, C_P::C_P] {
				myself.populations[s] <- self;
			}
		}
	}
	
	action update_output_data {
		loop s over:microbe_species {
			C[s] <- populations[s].C;
			C_cytosol[s] <- populations[s].cytosol_C;
			N_cytosol[s] <- populations[s].cytosol_N;
			P_cytosol[s] <- populations[s].cytosol_P;
			N[s] <- populations[s].N;
			P[s] <- populations[s].P;
			awake[s] <- populations[s].awake_population;
		}
	}
	
	action update_dams {
		float prev_C <- sum(Dam collect each.dom[2]);
		float prev_N <- sum(Dam collect each.dom[2]);
		float prev_P <- sum(Dam collect each.dom[2]);
		do feed_dams;
		float current_C <- sum(Dam collect each.dom[2]);
		float current_N <- sum(Dam collect each.dom[0]);
		float current_P <- sum(Dam collect each.dom[1]);
		total_C_dom <- total_C_dom + (current_C - prev_C);
		total_N_dom <- total_N_dom + (current_N - prev_N);
		total_P_dom <- total_P_dom + (current_P - prev_P);
	}
	
	action feed_dams virtual:true;
	
	output {
		display "C" {
			chart "Microbes populations" {
				if enable_Y_Strategist {
					data "Y strategists (g)" value:C[Y_Strategist]/#gram marker: false;
					if only_one_population {
						data "N Y strategists (g)" value:N[Y_Strategist]/#gram marker: false;
						data "P Y strategists (g)" value:P[Y_Strategist]/#gram marker: false;
					}
				}
				if enable_A_Strategist {
					data "A strategists (g)" value:C[A_Strategist]/#gram marker: false;
					if only_one_population {
						data "N A strategists (g)" value:N[A_Strategist]/#gram marker: false;
						data "P A strategists (g)" value:P[A_Strategist]/#gram marker: false;
					}
				}
				if enable_S_Strategist {
					data "S strategists (g)" value:C[S_Strategist]/#gram marker: false;
					if only_one_population {
						data "N S strategists (g)" value:N[S_Strategist]/#gram marker: false;
						data "P S strategists (g)" value:P[S_Strategist]/#gram marker: false;
					}
				}
			}
		}
		
		display "CO2" {
			chart "CO2" {
				data "CO2" value: microbe_CO2_emissions/#gram marker:false;
			}
		}
		
		display "C cytosol" {
			chart "Microbes populations cytosol" {
				if enable_Y_Strategist {
					data "Y strategists (g)" value:C_cytosol[Y_Strategist]/#gram marker: false;
				}
				if enable_A_Strategist {
					data "A strategists (g)" value:C_cytosol[A_Strategist]/#gram marker: false;
				}
				if enable_S_Strategist {
					data "S strategists (g)" value:C_cytosol[S_Strategist]/#gram marker: false;
				}
			}
			
		}
		
		display "Awake rate" {
			chart "Awake rate" {
				if enable_Y_Strategist {
					data "Y strategist (%)" value:awake[Y_Strategist] marker: false;
				}
				if enable_A_Strategist {
					data "A strategist (%)" value:awake[A_Strategist] marker: false;
				}
				if enable_S_Strategist {
					data "S strategist (%)" value:awake[S_Strategist] marker: false;
				}
			}
		}
		
		display "Mass conservation" {
			chart "Mass conservation" {
				data "C"
					value:(total_C_dom + total_C_labile + total_C_recal + sum(init_C.values))
					/ (
						sum(Dam collect each.dom[2]) +
						sum(OrganicParticle collect each.C_labile) +
						sum(OrganicParticle collect each.C_recalcitrant) +
						sum(C_cytosol.values) + sum(C.values) + microbe_CO2_emissions
					)
					marker: false;
				data "N"
					value:(total_N_dom + total_N_labile + total_N_recal + sum(init_C.values collect (each/C_N)))
					/ (
						sum(Dam collect (each.dom[0] + each.dim[0])) +
						sum(OrganicParticle collect each.N_labile) +
						sum(OrganicParticle collect each.N_recalcitrant) +
						sum(N_cytosol.values) + sum(N.values)
					)
					marker: false;
				data "P"
					value:(total_P_dom + total_P_labile + total_P_recal + sum(init_C.values collect (each/C_P)))
					/ (
						sum(Dam collect (each.dom[1] + each.dim[1])) +
						sum(OrganicParticle collect each.P_labile) +
						sum(OrganicParticle collect each.P_recalcitrant) +
						sum(P_cytosol.values) + sum(P.values)
					)
					marker: false;
			}
		}
		display "C/N" {
			chart "C/N" {
				if enable_Y_Strategist {
					data "Y strategists (g)" value:C[Y_Strategist]/N[Y_Strategist] marker: false;
				}
				if enable_A_Strategist {
					data "A strategists (g)" value:C[A_Strategist]/N[A_Strategist] marker: false;
				}
				if enable_S_Strategist {
					data "S strategists (g)" value:C[S_Strategist]/N[S_Strategist] marker: false;
				}
			}
		}
		display "C/P" {
			chart "C/P" {
				if enable_Y_Strategist {
					data "Y strategists (g)" value:C[Y_Strategist]/P[Y_Strategist] marker: false;
				}
				if enable_A_Strategist {
					data "A strategists (g)" value:C[A_Strategist]/P[A_Strategist] marker: false;
				}
				if enable_S_Strategist {
					data "S strategists (g)" value:C[S_Strategist]/P[S_Strategist] marker: false;
				}
			}
		}
	}
}

experiment IndividualMicrobesGrowth parent:MicrobesTestBase {
	map<species<MicrobePopulation>, PoreParticle> pores;
	map<species<MicrobePopulation>, Dam> dams;
	
	init {
		loop s over: microbe_species {
			create Dam {
				myself.dams[s] <- self;
			}
			
			create PoreParticle with: [
				dam: dams[s],
				carrying_capacity: carrying_capacity
			] {
				myself.pores[s] <- self;
				// Each population is associated to its own PoreParticle and Dam
				add myself.populations[s] to: populations;
			}
		}
	}
	
	action feed_dams {
		loop s over:microbe_species {
			Dam dam <- dams[s];
			MicrobePopulation population <- populations[s];
			write "" + dam + " " + (dam = nil);
			write "s arg: " + s;
			do feed_dam(population, dam);
		}
	}
	
	
	action update_populations {
		ask populations.values {
			do update(C, myself.carrying_capacity);
		}
	}
	
	action feed_dam(MicrobePopulation population, Dam dam) virtual:true;
}

experiment IndividualMicrobesGrowth_InfiniteNutrients parent:IndividualMicrobesGrowth {
	float nutrient_rate <- 1.0;
	parameter "Rate of requested nutrients in the DAM" var:nutrient_rate;
	
	reflex {
		do update_populations;
		do update_dams;
		
		ask PoreParticle {
			do microbe_life;
		}
		do update_output_data;
	}
	
	action feed_dam(MicrobePopulation population, Dam dam) {
		write "" + dam + " " + (dam = nil);
		// Each dam is fed with the nutrients requested by the population
		ask dam {
			dom[2] <- dom[2] + myself.nutrient_rate * population.requested_C;
			dom[0] <- dom[0] + myself.nutrient_rate * population.requested_N;
			dom[1] <- dom[1] + myself.nutrient_rate * population.requested_P;			
		}
	}
}

experiment IndividualMicrobesGrowth_FixedNutrients parent:IndividualMicrobesGrowth {
	float C_dom <- 15#gram;
	map<species<MicrobePopulation>, float> C_dom_output <- [Y_Strategist::C_dom, A_Strategist::C_dom, S_Strategist::C_dom];
	
	parameter "C dom" var:C_dom;
	
	init {
		do update_populations;
		do update_dams;
	}
	
	reflex {
		do update_populations;
		ask PoreParticle {
			do microbe_life;
		}
		do update_output_data;
		loop s over:microbe_species {
			C_dom_output[s] <- dams[s].dom[2];
		}
	}
	
	action feed_dam(MicrobePopulation population, Dam dam) {
		ask dam {
			// Each dam is fed with the nutrients requested by the population
			dom[2] <- dom[2] + myself.C_dom;
			dom[0] <- dom[0] + myself.C_dom / population.requested_C_N;
			dom[1] <- dom[1] + myself.C_dom / population.requested_C_P;
		}
	}
	
	output {
		display "C dom" {
			chart "C dom" {
				if enable_Y_Strategist {
					data "C dom Y strategists (g)" value:C_dom_output[Y_Strategist]/#gram marker: false;
				}
				if enable_A_Strategist {
					data "C dom A strategists (g)" value:C_dom_output[A_Strategist]/#gram marker: false;
				}
				if enable_S_Strategist {
					data "C dom S strategists (g)" value:C_dom_output[S_Strategist]/#gram marker: false;
				}
			}
		}
	}
}

experiment CollectiveMicrobesGrowth parent:MicrobesTestBase {
	Dam dam;
	
	init {
		create Dam {
			myself.dam <- self;
		}
		
		create PoreParticle with: [
			dam: dam,
			carrying_capacity: carrying_capacity
		] {
			// The three populations are added to the same PoreParticle
			loop s over:myself.microbe_species {
				add myself.populations[s] to: populations;
			}
		}
	}
	
	action update_populations {
		ask populations.values {
			do update(sum(myself.populations collect each.C), myself.carrying_capacity);
		}
	}
	
	action feed_dams {
		do feed_dam;
	}
	
	action feed_dam virtual:true;
}

experiment CollectiveMicrobesGrowth_InfiniteNutrients parent:CollectiveMicrobesGrowth {
	float nutrient_rate <- 1.0;
	parameter "Rate of requested nutrients in the DAM" var:nutrient_rate;
	
	reflex {
		do update_populations;
		do update_dams;
		
		ask PoreParticle {
			do microbe_life;
		}
		
		do update_output_data;
	}
	
	action feed_dam {
		ask Dam {
			// The dam is fed with the nutrients requested by all populations
			dom[2] <- myself.nutrient_rate * sum(myself.microbe_species collect myself.populations[each].requested_C);
			dom[0] <- myself.nutrient_rate * sum(myself.microbe_species collect myself.populations[each].requested_N);
			dom[1] <- myself.nutrient_rate * sum(myself.microbe_species collect myself.populations[each].requested_P);
		}
	}
}

experiment CollectiveMicrobesGrowth_FixedNutrients parent:CollectiveMicrobesGrowth {
	float C_dom <- 15#gram;
	map<species<MicrobePopulation>, float> C_dom_output <- [Y_Strategist::C_dom, A_Strategist::C_dom, S_Strategist::C_dom];
	
	parameter "C dom" var:C_dom;
	
	init {
		do update_populations;
		do update_dams;
	}
	
	reflex {
		do update_populations;
		
		ask PoreParticle {
			do microbe_life;
		}
		
		do update_output_data;
	}
	
	action feed_dam {
		ask Dam {
			// The dam is fed with the nutrients requested by all populations
			dom[2] <- dom[2] + myself.C_dom;
			dom[0] <- dom[0] + myself.C_dom / (min(myself.populations collect each.requested_C_N));
			dom[1] <- dom[1] + myself.C_dom / (min(myself.populations collect each.requested_C_P));
		}
	}
}

experiment IndividualMicrobesMetabolism parent:IndividualMicrobesGrowth {
	float init_C_labile <- 1.0;
	float C_N_labile <- 20.0;
	float C_P_labile <- 34.0;
	
	float init_C_recal <- 10.0;
	float C_N_recal <- 20.0;
	float C_P_recal <- 34.0;
	
	parameter "C labile (g)" var:init_C_labile;
	parameter "C/N labile (g)" var:C_N_labile;
	parameter "C/P labile (g)" var:C_P_labile;
	parameter "C recal (g)" var:init_C_recal;
	parameter "C/N recal (g)" var:C_N_recal;
	parameter "C/P recal (g)" var:C_P_recal;
	
	map<species<MicrobePopulation>, float> C_dom_output <- [Y_Strategist::0.0, A_Strategist::0.0, S_Strategist::0.0];
	map<species<MicrobePopulation>, float> N_dom_output <- [Y_Strategist::0.0, A_Strategist::0.0, S_Strategist::0.0];
	map<species<MicrobePopulation>, float> P_dom_output <- [Y_Strategist::0.0, A_Strategist::0.0, S_Strategist::0.0];
	map<species<MicrobePopulation>, float> N_dim_output <- [Y_Strategist::0.0, A_Strategist::0.0, S_Strategist::0.0];
	map<species<MicrobePopulation>, float> P_dim_output <- [Y_Strategist::0.0, A_Strategist::0.0, S_Strategist::0.0];
	
	map<species<MicrobePopulation>, float> C_labile_output <- [Y_Strategist::0.0, A_Strategist::0.0, S_Strategist::0.0];
	map<species<MicrobePopulation>, float> N_labile_output <- [Y_Strategist::0.0, A_Strategist::0.0, S_Strategist::0.0];
	map<species<MicrobePopulation>, float> P_labile_output <- [Y_Strategist::0.0, A_Strategist::0.0, S_Strategist::0.0];
	
	map<species<MicrobePopulation>, list<float>> enzymes_output <- [
		Y_Strategist::[0.0, 0.0, 0.0, 0.0],
		A_Strategist::[0.0, 0.0, 0.0, 0.0],
		S_Strategist::[0.0, 0.0, 0.0, 0.0]
	];
	
	action feed_dam(MicrobePopulation population, Dam dam) {
		// Nothing to do
	}
	
	init {
		loop s over: microbe_species {
			create OrganicParticle with: [
				C_labile: init_C_labile#gram,
				N_labile: init_C_labile#gram/C_N_labile,
				P_labile: init_C_labile#gram/C_P_labile,
				C_recalcitrant: init_C_recal#gram,
				N_recalcitrant: init_C_recal#gram/C_N_recal,
				P_recalcitrant: init_C_recal#gram/C_P_recal
			] {
				myself.total_C_labile <- myself.total_C_labile + C_labile;
				myself.total_N_labile <- myself.total_N_labile + N_labile;
				myself.total_P_labile <- myself.total_P_labile + P_labile;
				myself.total_C_recal <- myself.total_C_recal + C_recalcitrant;
				myself.total_N_recal <- myself.total_N_recal + N_recalcitrant;
				myself.total_P_recal <- myself.total_P_recal + P_recalcitrant;
				ask myself.pores[s] {
					add myself to: self.accessible_organics;
				}
			}
		}
		ask populations.values {
			do update(C, myself.carrying_capacity);
			write "Enzyme optimization (" + species(self) +"):";
			write "  DOM: " + myself.dams[species(self)].dom[2] + ", " + myself.dams[species(self)].dom[0] + ", " + myself.dams[species(self)].dom[1];
			write "  Labile: " + sum(myself.pores[species(self)].accessible_organics collect each.C_labile)/#gram + ", "
				+ sum(myself.pores[species(self)].accessible_organics collect each.N_labile)/#gram + ", "
				+ sum(myself.pores[species(self)].accessible_organics collect each.P_labile)/#gram;
			do optimize_enzymes(myself.dams[species(self)], myself.pores[species(self)].accessible_organics);
		}
	}
	
	reflex {
		ask PoreParticle {
			do decompose;
			do microbe_life;
		}
		do update_output_data;
		loop s over:microbe_species {
			C_dom_output[s] <- dams[s].dom[2];
			N_dom_output[s] <- dams[s].dom[0];
			P_dom_output[s] <- dams[s].dom[1];
			N_dim_output[s] <- dams[s].dim[0];
			P_dim_output[s] <- dams[s].dim[1];
			C_labile_output[s] <- sum(pores[s].accessible_organics collect each.C_labile);
			N_labile_output[s] <- sum(pores[s].accessible_organics collect each.N_labile);
			P_labile_output[s] <- sum(pores[s].accessible_organics collect each.P_labile);
			enzymes_output[s] <- [
				populations[s].enzymes.T_cellulolytic,
				populations[s].enzymes.T_amino,
				populations[s].enzymes.T_P,
				populations[s].enzymes.T_recal
				];
		}
	}
		
	output {
		display "C/N/P compartments" {
			chart "C/N/P compartments" {
				if enable_Y_Strategist {
					data "C dom Y strategists (g)" value:C_dom_output[Y_Strategist]/#gram marker: false;
					if only_one_population {
						data "N dom Y strategists (g)" value:N_dom_output[Y_Strategist]/#gram marker: false;
						data "P dom Y strategists (g)" value:P_dom_output[Y_Strategist]/#gram marker: false;
						data "N dim Y strategists (g)" value:N_dim_output[Y_Strategist]/#gram marker: false;
						data "P dim Y strategists (g)" value:P_dim_output[Y_Strategist]/#gram marker: false;
					}
					data "C labile Y strategists (g)" value:C_labile_output[Y_Strategist]/#gram marker: false;
					if only_one_population {
						data "N labile Y strategists (g)" value:N_labile_output[Y_Strategist]/#gram marker: false;
						data "P labile Y strategists (g)" value:P_labile_output[Y_Strategist]/#gram marker: false;
					}
				}
				if enable_A_Strategist {
					data "C dom A strategists (g)" value:C_dom_output[A_Strategist]/#gram marker: false;
					if only_one_population {
						data "N dom A strategists (g)" value:N_dom_output[A_Strategist]/#gram marker: false;
						data "P dom A strategists (g)" value:P_dom_output[A_Strategist]/#gram marker: false;
						data "N dim Y strategists (g)" value:N_dim_output[Y_Strategist]/#gram marker: false;
						data "P dim Y strategists (g)" value:P_dim_output[Y_Strategist]/#gram marker: false;
					}
					data "C labile A strategists (g)" value:C_labile_output[A_Strategist]/#gram marker: false;
					if only_one_population {
						data "N labile A strategists (g)" value:N_labile_output[A_Strategist]/#gram marker: false;
						data "P labile A strategists (g)" value:P_labile_output[A_Strategist]/#gram marker: false;
					}
				}
				if enable_S_Strategist {
					data "C dom S strategists (g)" value:C_dom_output[S_Strategist]/#gram marker: false;
					if only_one_population {
						data "N dom S strategists (g)" value:N_dom_output[S_Strategist]/#gram marker: false;
						data "P dom S strategists (g)" value:P_dom_output[S_Strategist]/#gram marker: false;
						data "N dim Y strategists (g)" value:N_dim_output[Y_Strategist]/#gram marker: false;
						data "P dim Y strategists (g)" value:P_dim_output[Y_Strategist]/#gram marker: false;
					}
					data "C labile S strategists (g)" value:C_labile_output[S_Strategist]/#gram marker: false;
					if only_one_population {
						data "N labile S strategists (g)" value:N_labile_output[S_Strategist]/#gram marker: false;
						data "P labile S strategists (g)" value:P_labile_output[S_Strategist]/#gram marker: false;
					}
				}
			}
		}
		display "Enzymes" {
			chart "Enzymes (Y Strategists)" {
				if enable_Y_Strategist {
					data "T C (gS/gM/h)" value:enzymes_output[Y_Strategist][0]/(#gram/#gram/#h) marker: false;
					data "T N (gS/gM/h)" value:enzymes_output[Y_Strategist][1]/(#gram/#gram/#h) marker: false;
					data "T P (gS/gM/h)" value:enzymes_output[Y_Strategist][2]/(#gram/#gram/#h) marker: false;
					data "T r (gS/gM/h)" value:enzymes_output[Y_Strategist][3]/(#gram/#gram/#h) marker: false;
				}
			}
		}
	}
}
