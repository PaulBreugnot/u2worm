/**
* Name: testmicrobes
* This model allows to test and calibrate the behavior of standalone microbe
* populations. It does not add any feature to the CAMMISOL model.
* 
* Author: pbreugno 
*/

model test_microbes

import "microbes.gaml"

experiment MicrobesTestBase {
	float init_O_C <- 1.0 on_change: change_species;
	float init_F_C <- 1.0 on_change: change_species;
	float init_M_C <- 1.0 on_change: change_species;
	float carrying_capacity <- 10#gram;
	
	bool enable_O_Strategist <- true on_change: change_species;
	bool enable_F_Strategist <- true on_change: change_species;
	bool enable_M_Strategist <- true on_change: change_species;
	
	bool only_one_population <- false;
	
	parameter "Enable O Strategist" var:enable_O_Strategist category: "O strategist";
	parameter "O C init (g)" var:init_O_C category:"O strategist";
	parameter "O dividing time" var:dividing_time_O category:"O strategist";
	parameter "O CUE" var:carbon_use_efficiency_O category:"O strategist";
	
	parameter "Enable F Strategist" var:enable_F_Strategist category: "F strategist";
	parameter "F C init (g)" var:init_F_C category:"F strategist";
	parameter "F dividing time" var:dividing_time_F category:"F strategist";
	parameter "F CUE" var:carbon_use_efficiency_F category:"F strategist";
	
	parameter "Enable M Strategist" var:enable_M_Strategist category: "M strategist";
	parameter "M C init (g)" var:init_M_C category:"M strategist";
	parameter "M dividing time" var:dividing_time_M category:"M strategist";
	parameter "M CUE" var:carbon_use_efficiency_M category:"M strategist";
	
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
	
	string enzymes_display <- "Budget" among:["Budget", "Absolute"] on_change:update_enzyme_unit;
	string enzyme_unit <- "gS/gM/h";
	map<species<MicrobePopulation>, list<float>> enzymes_output <- [
		O_Strategist::[0.0, 0.0, 0.0, 0.0],
		F_Strategist::[0.0, 0.0, 0.0, 0.0],
		M_Strategist::[0.0, 0.0, 0.0, 0.0]
	];
	
	action update_enzyme_unit {
		if enzymes_display = "Budget" {
			enzyme_unit <- "gS/gM/h";
		} else if enzymes_display = "Absolute" {
			enzyme_unit <- "gS/h";
		}
	}
	
	action change_species {
		microbe_species <- init_species();
	}
	
	action init_species {
		list<species<MicrobePopulation>> _microbe_species;
		init_C <- [];
		if enable_O_Strategist {
			add O_Strategist to: _microbe_species;
			init_C[O_Strategist] <- init_O_C#gram;
		}
		if enable_F_Strategist {
			add F_Strategist to: _microbe_species;
			init_C[F_Strategist] <- init_F_C#gram;
		}
		if enable_M_Strategist {
			add M_Strategist to: _microbe_species;
			init_C[M_Strategist] <- init_M_C#gram;
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
		only_one_population <- int(enable_O_Strategist) + int(enable_F_Strategist) + int(enable_M_Strategist) = 1;
		return _microbe_species;
	}
	
	init {
		ask simulation {
			do init_enzymatic_optimisation;
			do init_enzymes;
		}
		loop s over: microbe_species {
			create s with: [C::init_C[s], C_N::C_N, C_P::C_P] {
				write N/#gram;
				write P/#gram;
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
			awake[s] <- populations[s].active_rate;
			
			if enzymes_display = "Budget" {
				enzymes_output[s] <- [
					populations[s].enzymes.T_CNP/(#gram/#gram/#h),
					populations[s].enzymes.T_N/(#gram/#gram/#h),
					populations[s].enzymes.T_P/(#gram/#gram/#h),
					populations[s].enzymes.T_recal/(#gram/#gram/#h)
					];
			} else if enzymes_display = "Absolute" {
				enzymes_output[s] <- [
					populations[s].enzymes.T_CNP * C[s] * awake[s] / (#gram/#h),
					populations[s].enzymes.T_N * C[s] * awake[s] / (#gram/#h),
					populations[s].enzymes.T_P * C[s] * awake[s] / (#gram/#h),
					populations[s].enzymes.T_recal * C[s] * awake[s] / (#gram/#h)
					];
			}
		}
	}
	
	action update_dams {
		float prev_C <- sum(Dam collect each.dom[2]);
		float prev_N <- sum(Dam collect each.dom[0]);
		float prev_P <- sum(Dam collect each.dom[1]);
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
				if enable_O_Strategist {
					data "C (O, g)" value:C[O_Strategist]/#gram marker: false;
					if only_one_population {
						data "N (O, g)" value:N[O_Strategist]/#gram marker: false;
						data "P (O, g)" value:P[O_Strategist]/#gram marker: false;
					}
				}
				if enable_F_Strategist {
					data "C (F, g)" value:C[F_Strategist]/#gram marker: false;
					if only_one_population {
						data "N (F, g)" value:N[F_Strategist]/#gram marker: false;
						data "P (F, g)" value:P[F_Strategist]/#gram marker: false;
					}
				}
				if enable_M_Strategist {
					data "C (M, g)" value:C[M_Strategist]/#gram marker: false;
					if only_one_population {
						data "N (M, g)" value:N[M_Strategist]/#gram marker: false;
						data "P (M, g)" value:P[M_Strategist]/#gram marker: false;
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
				if enable_O_Strategist {
					data "C (O, g)" value:C_cytosol[O_Strategist]/#gram marker: false;
					if only_one_population {
						data "N (O, g)" value:N_cytosol[O_Strategist]/#gram marker: false;
						data "P (O, g)" value:P_cytosol[O_Strategist]/#gram marker: false;
					}
				}
				if enable_F_Strategist {
					data "C (F, g)" value:C_cytosol[F_Strategist]/#gram marker: false;
					if only_one_population {
						data "N (F, g)" value:N_cytosol[F_Strategist]/#gram marker: false;
						data "P (F, g)" value:P_cytosol[F_Strategist]/#gram marker: false;
					}				
				}
				if enable_M_Strategist {
					data "C (M, g)" value:C_cytosol[M_Strategist]/#gram marker: false;
					if only_one_population {
						data "N (M, g)" value:N_cytosol[M_Strategist]/#gram marker: false;
						data "P (M, g)" value:P_cytosol[M_Strategist]/#gram marker: false;
					}
				}
			}
			
		}
		
		display "Awake rate" {
			chart "Awake rate" {
				if enable_O_Strategist {
					data "O (%)" value:awake[O_Strategist] marker: false;
				}
				if enable_F_Strategist {
					data "F (%)" value:awake[F_Strategist] marker: false;
				}
				if enable_M_Strategist {
					data "M (%)" value:awake[M_Strategist] marker: false;
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
				if enable_O_Strategist {
					data "O" value:C[O_Strategist]/N[O_Strategist] marker: false;
				}
				if enable_F_Strategist {
					data "F" value:C[F_Strategist]/N[F_Strategist] marker: false;
				}
				if enable_M_Strategist {
					data "M" value:C[M_Strategist]/N[M_Strategist] marker: false;
				}
			}
		}
		display "C/P" {
			chart "C/P" {
				if enable_O_Strategist {
					data "O" value:C[O_Strategist]/P[O_Strategist] marker: false;
				}
				if enable_F_Strategist {
					data "F" value:C[F_Strategist]/P[F_Strategist] marker: false;
				}
				if enable_M_Strategist {
					data "M" value:C[M_Strategist]/P[M_Strategist] marker: false;
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
			do feed_dam(population, dam);
		}
	}
	
	
	action update_populations {
		ask populations.values {
			do update;
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
	map<species<MicrobePopulation>, float> C_dom_output <- [O_Strategist::C_dom, F_Strategist::C_dom, M_Strategist::C_dom];
	
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
				if enable_O_Strategist {
					data "C dom (O, g)" value:C_dom_output[O_Strategist]/#gram marker: false;
				}
				if enable_F_Strategist {
					data "C dom (F, g)" value:C_dom_output[F_Strategist]/#gram marker: false;
				}
				if enable_M_Strategist {
					data "C dom (M, g)" value:C_dom_output[M_Strategist]/#gram marker: false;
				}
			}
		}
	}
}

experiment CollectiveMicrobesGrowth parent:MicrobesTestBase {
	Dam dam;
	PoreParticle pore;
	
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
			myself.pore <- self;
		}
	}
	
	action update_populations {
		ask populations.values {
			do update;
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
	map<species<MicrobePopulation>, float> C_dom_output <- [O_Strategist::C_dom, F_Strategist::C_dom, M_Strategist::C_dom];
	
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

experiment Sporulation parent:IndividualMicrobesGrowth {
	float phase_1_duration <- 24#h;
	float phase_2_duration <- 100#h;
	
	parameter "Duration of phase 1 (nutrients)" var:phase_1_duration;
	parameter "Duration of phase 2 (no nutrients)" var:phase_2_duration;
	parameter "Sporulation time" var:microbes_sporulation_time;
	parameter "Germination time" var:microbes_germination_time;
	
	reflex {
		do update_populations;
		do update_dams;
		
		ask PoreParticle {
			do microbe_life;
		}
		do update_output_data;
	}
	action feed_dam(MicrobePopulation population, Dam dam) {
		if local_time < phase_1_duration {
			dam.dom[2] <- population.requested_C;
			dam.dom[0] <- population.requested_N;
			dam.dom[1] <- population.requested_P;
		} else if local_time > (phase_2_duration + phase_1_duration) {
			dam.dom[2] <- population.requested_C / population.active_rate;
			dam.dom[0] <- population.requested_N / population.active_rate;
			dam.dom[1] <- population.requested_P / population.active_rate;
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
	parameter "Enzymes display" var:enzymes_display;
	
	map<species<MicrobePopulation>, float> C_dom_output <- [O_Strategist::0.0, F_Strategist::0.0, M_Strategist::0.0];
	map<species<MicrobePopulation>, float> N_dom_output <- [O_Strategist::0.0, F_Strategist::0.0, M_Strategist::0.0];
	map<species<MicrobePopulation>, float> P_dom_output <- [O_Strategist::0.0, F_Strategist::0.0, M_Strategist::0.0];
	map<species<MicrobePopulation>, float> N_dim_output <- [O_Strategist::0.0, F_Strategist::0.0, M_Strategist::0.0];
	map<species<MicrobePopulation>, float> P_dim_output <- [O_Strategist::0.0, F_Strategist::0.0, M_Strategist::0.0];
	
	map<species<MicrobePopulation>, float> C_labile_output <- [O_Strategist::0.0, F_Strategist::0.0, M_Strategist::0.0];
	map<species<MicrobePopulation>, float> N_labile_output <- [O_Strategist::0.0, F_Strategist::0.0, M_Strategist::0.0];
	map<species<MicrobePopulation>, float> P_labile_output <- [O_Strategist::0.0, F_Strategist::0.0, M_Strategist::0.0];
	
	map<species<MicrobePopulation>, float> C_recal_output <- [O_Strategist::0.0, F_Strategist::0.0, M_Strategist::0.0];
	map<species<MicrobePopulation>, float> N_recal_output <- [O_Strategist::0.0, F_Strategist::0.0, M_Strategist::0.0];
	map<species<MicrobePopulation>, float> P_recal_output <- [O_Strategist::0.0, F_Strategist::0.0, M_Strategist::0.0];
	action feed_dam(MicrobePopulation population, Dam dam) {
		// Nothing to do
		dam.dom[2] <- dam.dom[2] + population.requested_C;
		dam.dom[0] <- dam.dom[2] + population.requested_N;
		dam.dom[1] <- dam.dom[2] + population.requested_P;
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
			do update;
		}
		do update_dams;
	}
	
	reflex {
		ask populations.values {
			write "Enzyme optimization (" + species(self) +"):";
			write "  DOM: " + myself.dams[species(self)].dom[2] + ", " + myself.dams[species(self)].dom[0] + ", " + myself.dams[species(self)].dom[1];
			write "  Labile: " + sum(myself.pores[species(self)].accessible_organics collect each.C_labile)/#gram + ", "
				+ sum(myself.pores[species(self)].accessible_organics collect each.N_labile)/#gram + ", "
				+ sum(myself.pores[species(self)].accessible_organics collect each.P_labile)/#gram;
			do optimize_enzymes(myself.dams[species(self)], myself.pores[species(self)].accessible_organics);
		}
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
			C_recal_output[s] <- sum(pores[s].accessible_organics collect each.C_recalcitrant);
			N_recal_output[s] <- sum(pores[s].accessible_organics collect each.N_recalcitrant);
			P_recal_output[s] <- sum(pores[s].accessible_organics collect each.P_recalcitrant);
		}
	}
	
	output {
		display "C/N/P compartments" {
			chart "C/N/P compartments" {
				if enable_O_Strategist {
					data "C dom (O, g)" value:C_dom_output[O_Strategist]/#gram marker: false;
					if only_one_population {
						data "N dom (O, g)" value:N_dom_output[O_Strategist]/#gram marker: false;
						data "P dom (O, g)" value:P_dom_output[O_Strategist]/#gram marker: false;
						data "N dim (O, g)" value:N_dim_output[O_Strategist]/#gram marker: false;
						data "P dim (O, g)" value:P_dim_output[O_Strategist]/#gram marker: false;
					}
					data "C labile (O, g)" value:C_labile_output[O_Strategist]/#gram marker: false;
					if only_one_population {
						data "N labile (O, g)" value:N_labile_output[O_Strategist]/#gram marker: false;
						data "P labile (O, g)" value:P_labile_output[O_Strategist]/#gram marker: false;
					}
					data "C recal (O, g)" value:C_recal_output[O_Strategist]/#gram marker: false;
					if only_one_population {
						data "N recal (O, g)" value:N_recal_output[O_Strategist]/#gram marker: false;
						data "P recal (O, g)" value:P_recal_output[O_Strategist]/#gram marker: false;
					}
				}
				if enable_F_Strategist {
					data "C dom (F, g)" value:C_dom_output[F_Strategist]/#gram marker: false;
					if only_one_population {
						data "N dom (F, g)" value:N_dom_output[F_Strategist]/#gram marker: false;
						data "P dom (F, g)" value:P_dom_output[F_Strategist]/#gram marker: false;
						data "N dim (F, g)" value:N_dim_output[F_Strategist]/#gram marker: false;
						data "P dim (F, g)" value:P_dim_output[F_Strategist]/#gram marker: false;
					}
					data "C labile (F, g)" value:C_labile_output[F_Strategist]/#gram marker: false;
					if only_one_population {
						data "N labile (F, g)" value:N_labile_output[F_Strategist]/#gram marker: false;
						data "P labile (F, g)" value:P_labile_output[F_Strategist]/#gram marker: false;
					}
					data "C recal (F, g)" value:C_recal_output[F_Strategist]/#gram marker: false;
					if only_one_population {
						data "N recal (F, g)" value:N_recal_output[F_Strategist]/#gram marker: false;
						data "P recal (F, g)" value:P_recal_output[F_Strategist]/#gram marker: false;
					}
				}
				if enable_M_Strategist {
					data "C dom (M, g)" value:C_dom_output[M_Strategist]/#gram marker: false;
					if only_one_population {
						data "N dom (M, g)" value:N_dom_output[M_Strategist]/#gram marker: false;
						data "P dom (M, g)" value:P_dom_output[M_Strategist]/#gram marker: false;
						data "N dim (M, g)" value:N_dim_output[M_Strategist]/#gram marker: false;
						data "P dim (M, g)" value:P_dim_output[M_Strategist]/#gram marker: false;
					}
					data "C labile (M, g)" value:C_labile_output[M_Strategist]/#gram marker: false;
					if only_one_population {
						data "N labile (M, g)" value:N_labile_output[M_Strategist]/#gram marker: false;
						data "P labile (M, g)" value:P_labile_output[M_Strategist]/#gram marker: false;
					}
					data "C recal (M, g)" value:C_recal_output[M_Strategist]/#gram marker: false;
					if only_one_population {
						data "N recal (M, g)" value:N_recal_output[M_Strategist]/#gram marker: false;
						data "P recal (M, g)" value:P_recal_output[M_Strategist]/#gram marker: false;
					}
				}
			}
		}
		
		display "Enzymes" {
			chart "Enzymes" {
				if enable_O_Strategist {
					data "T C (O, " + enzyme_unit + ")" value:enzymes_output[O_Strategist][0] marker: false;
					data "T N (O, " + enzyme_unit + ")" value:enzymes_output[O_Strategist][1] marker: false;
					data "T P (O, " + enzyme_unit + ")" value:enzymes_output[O_Strategist][2] marker: false;
					data "T r (O, " + enzyme_unit + ")" value:enzymes_output[O_Strategist][3] marker: false;
				}
				if enable_F_Strategist {
					data "T C (F, " + enzyme_unit + ")" value:enzymes_output[F_Strategist][0] marker: false;
					data "T N (F, " + enzyme_unit + ")" value:enzymes_output[F_Strategist][1] marker: false;
					data "T P (F, " + enzyme_unit + ")" value:enzymes_output[F_Strategist][2] marker: false;
					data "T r (F, " + enzyme_unit + ")" value:enzymes_output[F_Strategist][3] marker: false;
				}
				if enable_M_Strategist {
					data "T C (M, " + enzyme_unit + ")" value:enzymes_output[M_Strategist][0] marker: false;
					data "T N (M, " + enzyme_unit + ")" value:enzymes_output[M_Strategist][1] marker: false;
					data "T P (M, " + enzyme_unit + ")" value:enzymes_output[M_Strategist][2] marker: false;
					data "T r (M, " + enzyme_unit + ")" value:enzymes_output[M_Strategist][3] marker: false;
				}
			}
		}
	}
}

experiment CollectiveMicrobesMetabolism parent:CollectiveMicrobesGrowth {
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
	parameter "Enzymes display" var:enzymes_display;
	
	float C_dom_output <- 0.0;
	float N_dom_output <- 0.0;
	float P_dom_output <- 0.0;
	float N_dim_output <- 0.0;
	float P_dim_output <- 0.0;
	
	float C_labile_output <- 0.0;
	float N_labile_output <- 0.0;
	float P_labile_output <- 0.0;
	
	float C_recal_output <- 0.0;
	float N_recal_output <- 0.0;
	float P_recal_output <- 0.0;
	action feed_dam {
		// Nothing to do
	}
	
	init {
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
			ask PoreParticle {
				add myself to: self.accessible_organics;
			}
		}
		ask populations.values {
			do update;
		}
	}
	
	reflex {
		ask populations.values {
			write "Enzyme optimization (" + species(self) +"):";
			write "  DOM: " + myself.dam.dom[2] + ", " + myself.dam.dom[0] + ", " + myself.dam.dom[1];
			write "  Labile: " + sum(myself.pore.accessible_organics collect each.C_labile)/#gram + ", "
				+ sum(myself.pore.accessible_organics collect each.N_labile)/#gram + ", "
				+ sum(myself.pore.accessible_organics collect each.P_labile)/#gram;
			do optimize_enzymes(myself.dam, myself.pore.accessible_organics);
		}
		ask PoreParticle {
			do decompose;
			do microbe_life;
		}
		do update_output_data;
		C_dom_output <- dam.dom[2];
		N_dom_output <- dam.dom[0];
		P_dom_output <- dam.dom[1];
		N_dim_output <- dam.dim[0];
		P_dim_output <- dam.dim[1];
		C_labile_output <- sum(pore.accessible_organics collect each.C_labile);
		N_labile_output <- sum(pore.accessible_organics collect each.N_labile);
		P_labile_output <- sum(pore.accessible_organics collect each.P_labile);
		C_recal_output <- sum(pore.accessible_organics collect each.C_recalcitrant);
		N_recal_output <- sum(pore.accessible_organics collect each.N_recalcitrant);
		P_recal_output <- sum(pore.accessible_organics collect each.P_recalcitrant);
	}
	
	output {
		display "C/N/P compartments" {
			chart "C/N/P compartments" {
				data "C dom (g)" value:C_dom_output/#gram marker: false;
				data "N dom (g)" value:N_dom_output/#gram marker: false;
				data "P dom (g)" value:P_dom_output/#gram marker: false;
				data "N dim (g)" value:N_dim_output/#gram marker: false;
				data "P dim (g)" value:P_dim_output/#gram marker: false;
				data "C labile (g)" value:C_labile_output/#gram marker: false;
				data "N labile (g)" value:N_labile_output/#gram marker: false;
				data "P labile (g)" value:P_labile_output/#gram marker: false;
				data "C recal (g)" value:C_recal_output/#gram marker: false;
				data "N recal (g)" value:N_recal_output/#gram marker: false;
				data "P recal (g)" value:P_recal_output/#gram marker: false;
			}
		}
		
		display "Enzymes" {
			chart "Enzymes" {
				if enable_O_Strategist {
					data "T C (O, " + enzyme_unit + ")" value:enzymes_output[O_Strategist][0] marker: false;
					data "T N (O, " + enzyme_unit + ")" value:enzymes_output[O_Strategist][1] marker: false;
					data "T P (O, " + enzyme_unit + ")" value:enzymes_output[O_Strategist][2] marker: false;
					data "T r (O, " + enzyme_unit + ")" value:enzymes_output[O_Strategist][3] marker: false;
				}
				if enable_F_Strategist {
					data "T C (F, " + enzyme_unit + ")" value:enzymes_output[F_Strategist][0] marker: false;
					data "T N (F, " + enzyme_unit + ")" value:enzymes_output[F_Strategist][1] marker: false;
					data "T P (F, " + enzyme_unit + ")" value:enzymes_output[F_Strategist][2] marker: false;
					data "T r (F, " + enzyme_unit + ")" value:enzymes_output[F_Strategist][3] marker: false;
				}
				if enable_M_Strategist {
					data "T C (M, " + enzyme_unit + ")" value:enzymes_output[M_Strategist][0] marker: false;
					data "T N (M, " + enzyme_unit + ")" value:enzymes_output[M_Strategist][1] marker: false;
					data "T P (M, " + enzyme_unit + ")" value:enzymes_output[M_Strategist][2] marker: false;
					data "T r (M, " + enzyme_unit + ")" value:enzymes_output[M_Strategist][3] marker: false;
				}
			}
		}
	}
}
