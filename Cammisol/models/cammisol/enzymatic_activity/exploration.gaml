/**
 * The purpose of this model is to explore the reactions of the enzymatic activity in different scenarios, that correspond to different environmental conditions.
 * 
 * It a standalone and experiment only model that does not add any feature to the cammisol or enzymatic activity models.
 */

model enzymatic_activity_exploration

import "../microbes/microbes.gaml"

global {
	string x_var;
	string output_file;
	map<string, species<ExplorationParameterGetter>> getter_types <- [
		"CN labile"::CN_LabileGetter,
		"CP labile"::CP_LabileGetter,
		"CN dom"::CN_DOMGetter,
		"CP dom"::CP_DOMGetter
	];
	list<species<MicrobePopulation>> populations <- [Y_Strategist, A_Strategist, S_Strategist];
	
	/* 
	 * Default values for all parameters to explore.
	 * Values are set and explored using experiment parameters.
	 */
	float C_labile <- 0.1#gram;
	float CN_labile <- 17.0;
	float CP_labile <- 31.0;
	float C_recalcitrant <- 10#gram;
	float CN_recalcitrant <- 17.0;
	float CP_recalcitrant <- 31.0;
	float C_dom <- 0#gram;
	float CN_dom <- 10.0;
	float CP_dom <- 17.0;
	
	/**
	 * A single organic particle is decomposed.
	 */
	OrganicParticle organic_particle;
	/**
	 * A single dam is considered in the environment.
	 */
	PoreParticle pore_particle;
	
	/**
	 * Output for each species (Y-A-S strategies).
	 */
	map<species<MicrobePopulation>, DataOutput> output;
	map<species<MicrobePopulation>, string> labels <- [Y_Strategist::"Y", A_Strategist::"A", S_Strategist::"S"];
	
	ExplorationParameterGetter param_getter;
	
	init {
		create getter_types[x_var] {
			param_getter <- self;
		}

		do init_enzymes;
		do init_enzymatic_optimisation;
		create OrganicParticle with: [
			C_labile: world.C_labile,
			N_labile: world.C_labile / CN_labile,
			P_labile: world.C_labile / CP_labile,
			C_recalcitrant: world.C_recalcitrant,
			N_recalcitrant: world.C_recalcitrant / CN_recalcitrant,
			P_recalcitrant: world.C_recalcitrant / CP_recalcitrant
			]{
			world.organic_particle <- self;
		}
		create Dam with: [
			dom::[world.C_dom / CN_dom, world.C_dom / CP_dom, world.C_dom]
			]{
			create PoreParticle with: [dam::self] {
				world.pore_particle <- self;
			}
		}
		loop s over:populations {
			create s with: [
				C: 1#gram,
				awake_population: 1.0
			] {
				do update;
			}
		}
		
		loop s over: populations {
			create DataOutput {
				output[s] <- self;
			}
		}
	}
	
	reflex {
		ask Y_Strategist + A_Strategist + S_Strategist {
			// Optimization
			do optimize_enzymes(world.pore_particle.dam, [world.organic_particle]);
			// Computes output
			create Decomposition {
				WeightedEnzymes weighted_enzymes;
				create WeightedEnzymes with: [
					T_C::myself.enzymes.T_C*myself.active_C,
					T_N::myself.enzymes.T_N*myself.active_C,
					T_P::myself.enzymes.T_P*myself.active_C,
					T_recal::myself.enzymes.T_recal*myself.active_C
				] {
					weighted_enzymes <- self;
				}
				
				ask decomposition_problem {
					do decomposition(weighted_enzymes, myself);
				}
				
				output[species(myself)].T_r <- myself.enzymes.T_recal;
				output[species(myself)].T_C <- myself.enzymes.T_C;
				output[species(myself)].T_N <- myself.enzymes.T_N;
				output[species(myself)].T_P <- myself.enzymes.T_P;
				
				output[species(myself)].X_C <- (X_C_eC + X_C_eN + X_C_eP)/#gram;
				output[species(myself)].X_N <- X_N_eN/#gram;
				output[species(myself)].X_P <- (X_P_eP_labile_to_dom + X_P_eP_labile_to_dim)/#gram;
				output[species(myself)].X_recal <- X_C_eR/#gram;
				
				float C_avail; float N_avail; float P_avail;
				ask decomposition_problem {
					C_avail <- C_avail_final(myself);
					N_avail <- N_avail_final(myself);
					P_avail <- P_avail_final(myself);
				}
				
				output[species(myself)].C_recal <- (world.C_recalcitrant - X_C_eR)/#gram;
				
				output[species(myself)].C_labile <- (world.C_labile + X_C_eR - (X_C_eC + X_C_eN + X_C_eP))/#gram;
				output[species(myself)].N_labile <- (world.C_labile/world.CN_labile + X_N_eR - X_N_eN)/#gram;
				output[species(myself)].P_labile <- (world.C_labile/world.CP_labile + X_P_eR_recal_to_labile - X_P_eP_labile_to_dom - X_P_eP_labile_to_dim)/#gram;
				
				output[species(myself)].C_avail <- (world.C_dom + X_C_eC + X_C_eN + X_C_eP)/#gram;
				output[species(myself)].N_avail <- (world.C_dom/world.CN_dom + X_N_eN)/#gram;
				output[species(myself)].P_avail <- (world.C_dom/world.CP_dom + X_P_eP_labile_to_dom + X_P_eP_labile_to_dim + X_P_eR_recal_to_dim)/#gram;
				
				ask weighted_enzymes {
					do die;
				}
				do die;
			}
		}
	}
	
	reflex {
		list<float> data <- [];
		map<species<MicrobePopulation>, int> data_index;
		loop s over: [Y_Strategist, A_Strategist, S_Strategist] {
			data_index[s] <- length(data);
			data <- data + [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
		}
		
		ask Y_Strategist + A_Strategist + S_Strategist {
			int current_index <- data_index[species(self)];
			data[current_index] <- self.enzymes.T_recal/(#gram/#gram/#d);
			data[current_index+1] <- self.enzymes.T_C/(#gram/#gram/#d);
			data[current_index+2] <- self.enzymes.T_N/(#gram/#gram/#d);
			data[current_index+3] <- self.enzymes.T_P/(#gram/#gram/#d);

			data[current_index+4] <- myself.output[species(self)].X_C;
			data[current_index+5] <- myself.output[species(self)].X_N;
			data[current_index+6] <- myself.output[species(self)].X_P;
			data[current_index+7] <- myself.output[species(self)].X_recal;
				
			data[current_index+8] <- myself.output[species(self)].C_recal;
			data[current_index+9] <- myself.output[species(self)].C_labile;
			data[current_index+10] <- myself.output[species(self)].N_labile;
			data[current_index+11] <- myself.output[species(self)].P_labile;
			data[current_index+12] <- myself.output[species(self)].C_avail;
			data[current_index+13] <- myself.output[species(self)].N_avail;
			data[current_index+14] <- myself.output[species(self)].P_avail;
		}
		save [self.index, param_getter.param_value()] + data header: false rewrite:false to: output_file + "_data.csv" format: "csv";
	}
}

species ExplorationParameterGetter {
	float param_value virtual:true;
}

species CN_LabileGetter parent:ExplorationParameterGetter {
	float param_value {
		return world.CN_labile;
	}
}

species CP_LabileGetter parent:ExplorationParameterGetter {
	float param_value {
		return world.CP_labile;
	}
}

species CN_DOMGetter parent:ExplorationParameterGetter {
	float param_value {
		return world.CN_dom;
	}
}

species CP_DOMGetter parent:ExplorationParameterGetter {
	float param_value {
		return world.CP_dom;
	}
}

/**
 * Data structure used to store enzymatic activities and decomposition results for each species.
 */
species DataOutput schedules:[] {
	float T_r;
	float T_C;
	float T_N;
	float T_P;
	
	float X_C;
	float X_N;
	float X_P;
	float X_recal;
	
	float C_recal;
	float C_labile;
	float N_labile;
	float P_labile;
	float C_avail;
	float N_avail;
	float P_avail;
}

/**
 * Base virtual experiment that define outputs and displays.
 */
experiment Explore {
	method exploration;
	/**
	 * Name of the experimentation used to name the CSV output file.
	 */
	string exp_name;
	/**
	 * Parameters used to plot P or N decomposition, depending on the parameter
	 * currently explored.
	 */
	string plot <- "N" among: ["N", "P"];

	/*
	 * Values of parameters used to perform config output.
	 */
	bool exp_enable_Y_strategist;
	bool exp_enable_A_strategist;
	bool exp_enable_S_strategist;
	float exp_C_labile;
	float exp_CN_labile;
	float exp_CP_labile;
	float exp_C_recalcitrant;
	float exp_CN_recalcitrant;
	float exp_CP_recalcitrant;
	float exp_C_dom;
	float exp_CN_dom;
	float exp_CP_dom;
		
	/* 
	 * Default values for model parameters. Those values are set as default
	 * read_only experiment parameters below (init facet), or superseded for
	 * exploration in child experiment parameters. In any case, experiment
	 * parameters values are used to initialise corresponding model attributes
	 * (init section below).
	 */
	parameter "Enable Y strategist" var:exp_enable_Y_strategist init:true read_only:true;
	parameter "Enable A strategist" var:exp_enable_A_strategist init:true read_only:true;
	parameter "Enable S strategist" var:exp_enable_S_strategist init:true read_only:true;
	parameter "Labile C" var:exp_C_labile init:0.1#gram read_only:true;
	parameter "Labile C/N" var:exp_CN_labile init:17.0 read_only:true;
	parameter "Labile C/P" var:exp_CP_labile init:31.0 read_only:true;
	parameter "Recalcitrant C" var:exp_C_recalcitrant init:10#gram read_only:true;
	parameter "Recalcitrant C/N" var:exp_CN_recalcitrant init:17.0 read_only:true;
	parameter "Recalcitrant C/P" var:exp_CP_recalcitrant init:31.0 read_only:true;
	parameter "C DOM" var:exp_C_dom init:0#gram read_only:true;
	parameter "C/N DOM" var:exp_CN_dom init:10.0 read_only:true;
	parameter "C/P DOM" var:exp_CP_dom init:17.0 read_only:true;
	
	action write_config_csv {
		list<string> header <- [
			"Period (days)",
			"Labile C (g)",
			"Labile C/N",
			"Labile C/P",
			"Recalcitrant C (g)",
			"Recalcitrant C/N",
			"Recalcitrant C/P",
			"C DOM (g)",
			"C/N DOM",
			"C/P DOM"
		];
		list<float> config <- [
			enzymes_optimization_period/#d,
			exp_C_labile/#gram,
			exp_CN_labile,
			exp_CP_labile,
			exp_C_recalcitrant/#gram,
			exp_CN_recalcitrant,
			exp_CP_recalcitrant,
			exp_C_dom/#gram,
			exp_CN_dom,
			exp_CP_dom
		];
		
		loop s over: [Y_Strategist, A_Strategist, S_Strategist] {
			add "Requested C/N (" + labels[s] + ")" to: header;
			add "Requested C/P (" + labels[s] + ")" to: header;
			add "Min T_C (" + labels[s] + ")" to: header;
			add "Min T_N (" + labels[s] + ")" to: header;
			add "Min T_P (" + labels[s] + ")" to: header;
			add "Min T_recal (" + labels[s] + ")" to: header;
			add "Max T_C (" + labels[s] + ")" to: header;
			add "Max T_N (" + labels[s] + ")" to: header;
			add "Max T_P (" + labels[s] + ")" to: header;
			add "Max T_recal (" + labels[s] + ")" to: header;
		}
		add microbes_CN / carbon_use_efficiency_Y to: config;
		add microbes_CP / carbon_use_efficiency_Y to: config;
		add min_T_C_Y/(#gram/#gram/#d) to: config;
		add min_T_N_Y/(#gram/#gram/#d) to: config;
		add min_T_P_Y/(#gram/#gram/#d) to: config;
		add min_T_recal_Y/(#gram/#gram/#d) to: config;
		add max_T_C_Y/(#gram/#gram/#d) to: config;
		add max_T_N_Y/(#gram/#gram/#d) to: config;
		add max_T_P_Y/(#gram/#gram/#d) to: config;
		add max_T_recal_Y/(#gram/#gram/#d) to: config;
		
		add microbes_CN / carbon_use_efficiency_A to: config;
		add microbes_CP / carbon_use_efficiency_A to: config;
		add min_T_C_A/(#gram/#gram/#d) to: config;
		add min_T_N_A/(#gram/#gram/#d) to: config;
		add min_T_P_A/(#gram/#gram/#d) to: config;
		add min_T_recal_A/(#gram/#gram/#d) to: config;
		add max_T_C_A/(#gram/#gram/#d) to: config;
		add max_T_N_A/(#gram/#gram/#d) to: config;
		add max_T_P_A/(#gram/#gram/#d) to: config;
		add max_T_recal_A/(#gram/#gram/#d) to: config;
		
		add microbes_CN / carbon_use_efficiency_S to: config;
		add microbes_CP / carbon_use_efficiency_S to: config;
		add min_T_C_S/(#gram/#gram/#d) to: config;
		add min_T_N_S/(#gram/#gram/#d) to: config;
		add min_T_P_S/(#gram/#gram/#d) to: config;
		add min_T_recal_S/(#gram/#gram/#d) to: config;
		add max_T_C_S/(#gram/#gram/#d) to: config;
		add max_T_N_S/(#gram/#gram/#d) to: config;
		add max_T_P_S/(#gram/#gram/#d) to: config;
		add max_T_recal_S/(#gram/#gram/#d) to: config;
			
		save header header: false rewrite:true to: exp_name + "_config.csv" format: "csv";
		save config header: false rewrite:false to: exp_name + "_config.csv" format: "csv";
	}
	
	/**
	 * Initialise CSV headers.
	 */
	action init_data_csv {
		list<string> header <- ["n", x_label()];
		loop strategy over: ["Y", "A", "S"] {
			header <- header + [
				"T recal (" + strategy + ")",
				"T C (" + strategy + ")",
				"T N (" + strategy + ")",
				"T P (" + strategy + ")",
				"X C (" + strategy + ")",
				"X N (" + strategy + ")",
				"X P (" + strategy + ")",
				"X C recal (" + strategy + ")",
				"C recal (" + strategy + ")",
				"C labile (" + strategy + ")",
				"N labile (" + strategy + ")",
				"P labile (" + strategy + ")",
				"C avail (" + strategy + ")",
				"N avail (" + strategy + ")",
				"P avail (" + strategy + ")"
			];
		}
		save header header:false rewrite: true to: exp_name + "_data.csv" format: "csv";
	}
	
	/**
	 * Name of the explored parameter (e.g. C/N or C/P).
	 */
	string x_label virtual: true;
	
	init {
		list<species<MicrobePopulation>> exp_populations;

		if exp_enable_Y_strategist {
			add Y_Strategist to: exp_populations;
		}
		if exp_enable_A_strategist {
			add A_Strategist to: exp_populations;
		}
		if exp_enable_S_strategist {
			add S_Strategist to: exp_populations;
		}
		populations <- exp_populations;
		C_labile <- exp_C_labile;
		CN_labile <- exp_CN_labile;
		CP_labile <- exp_CP_labile;
		C_recalcitrant <- exp_C_recalcitrant;
		CN_recalcitrant <- exp_CN_recalcitrant;
		CP_recalcitrant <- exp_CP_recalcitrant;
		C_dom <- exp_C_dom;
		CN_dom <- exp_CN_dom;
		CP_dom <- exp_CP_dom;
	}
	
	permanent {
		display "Enzymatic activity (Y)" {
			chart "Enzymatic activity (Y)" type:xy x_label: x_label() y_label: "Enzymatic activity (gC/gM/d)" {
				if exp_enable_Y_strategist {
					data "T_P" value: {mean(simulations collect each.param_getter.param_value()), mean(simulations collect (each.output[Y_Strategist].T_P / (#gram / #gram / #d)))};
					data "T_C" value: {mean(simulations collect each.param_getter.param_value()), mean(simulations collect (each.output[Y_Strategist].T_C / (#gram / #gram / #d)))};
					data "T_N" value: {mean(simulations collect each.param_getter.param_value()), mean(simulations collect (each.output[Y_Strategist].T_N / (#gram / #gram / #d)))};
					data "T_recal" value: {mean(simulations collect each.param_getter.param_value()), mean(simulations collect (each.output[Y_Strategist].T_r / (#gram / #gram / #d)))};
				}
			}
		}
			
		display "Enzymatic activity (A)" {
			chart "Enzymatic activity (A)" type:xy x_label: x_label() y_label: "Enzymatic activity (gC/gM/d)" {
				if exp_enable_A_strategist {
					data "T_P" value: {mean(simulations collect each.param_getter.param_value()), mean(simulations collect (each.output[A_Strategist].T_P / (#gram / #gram / #d)))};
					data "T_C" value: {mean(simulations collect each.param_getter.param_value()), mean(simulations collect (each.output[A_Strategist].T_C / (#gram / #gram / #d)))};
					data "T_N" value: {mean(simulations collect each.param_getter.param_value()), mean(simulations collect (each.output[A_Strategist].T_N / (#gram / #gram / #d)))};
					data "T_recal" value: {mean(simulations collect each.param_getter.param_value()), mean(simulations collect (each.output[A_Strategist].T_r / (#gram / #gram / #d)))};
				}
			}
		}
			
		display "Enzymatic activity (S)" {
			chart "Enzymatic activity (S)" type:xy x_label: x_label() y_label: "Enzymatic activity (gC/gM/d)" {
				if exp_enable_S_strategist {
					data "T_P" value: {mean(simulations collect each.param_getter.param_value()), mean(simulations collect (each.output[S_Strategist].T_P / (#gram / #gram / #d)))};
					data "T_C" value: {mean(simulations collect each.param_getter.param_value()), mean(simulations collect (each.output[S_Strategist].T_C / (#gram / #gram / #d)))};
					data "T_N" value: {mean(simulations collect each.param_getter.param_value()), mean(simulations collect (each.output[S_Strategist].T_N / (#gram / #gram / #d)))};
					data "T_recal" value: {mean(simulations collect each.param_getter.param_value()), mean(simulations collect (each.output[S_Strategist].T_r / (#gram / #gram / #d)))};
				}
			}
		}
		
		display "Decomposition (Y)" {
			chart "Decomposition (Y)" type:xy x_label: x_label() y_label: "Decomposition (g)" {
				if exp_enable_Y_strategist {
					if (plot = "N") {
						data "N labile (final)" value: {simulations mean_of each.param_getter.param_value(), simulations mean_of each.output[Y_Strategist].N_labile};
						data "N avail" value: {simulations mean_of each.param_getter.param_value(), simulations mean_of each.output[Y_Strategist].N_avail};
						data "N labile (init)" value: {simulations mean_of each.param_getter.param_value(), simulations mean_of (each.C_labile / each.CN_labile) / #gram} marker: false;
					} else if (plot = "P") {
						data "P labile (final)" value: {simulations mean_of each.param_getter.param_value(), simulations mean_of each.output[Y_Strategist].P_labile};
						data "P avail" value: {simulations mean_of each.param_getter.param_value(), simulations mean_of each.output[Y_Strategist].P_avail};
						data "P labile (init)" value: {simulations mean_of each.param_getter.param_value(), simulations mean_of (each.C_labile / each.CP_labile) / #gram} marker: false;
					}
				}
			}
		}
		display "Decomposition (A)" {
			chart "Decomposition (A)" type:xy x_label: x_label() y_label: "Decomposition (g)" {
				if exp_enable_A_strategist {
					if (plot = "N") {
						data "N labile (final)" value: {simulations mean_of each.param_getter.param_value(), simulations mean_of each.output[A_Strategist].N_labile};
						data "N avail" value: {simulations mean_of each.param_getter.param_value(), simulations mean_of each.output[A_Strategist].N_avail};
						data "N labile (init)" value: {simulations mean_of each.param_getter.param_value(), simulations mean_of (each.C_labile / each.CN_labile) / #gram} marker: false;
					} else if (plot = "P") {
						data "P labile (final)" value: {simulations mean_of each.param_getter.param_value(), simulations mean_of each.output[A_Strategist].P_labile};
						data "P avail" value: {simulations mean_of each.param_getter.param_value(), simulations mean_of each.output[A_Strategist].P_avail};
						data "P labile (init)" value: {simulations mean_of each.param_getter.param_value(), simulations mean_of (each.C_labile / each.CP_labile) / #gram} marker: false;
					}
				}
			}
		}
		display "Decomposition (S)" {
			chart "Decomposition (S)" type:xy x_label: x_label() y_label: "Decomposition (g)" {
				if exp_enable_S_strategist {
					if (plot = "N") {
						data "N labile (final)" value: {simulations mean_of each.param_getter.param_value(), simulations mean_of each.output[S_Strategist].N_labile};
						data "N avail" value: {simulations mean_of each.param_getter.param_value(), simulations mean_of each.output[S_Strategist].N_avail};
						data "N labile (init)" value: {simulations mean_of each.param_getter.param_value(), simulations mean_of (each.C_labile / each.CN_labile) / #gram} marker: false;
					} else if (plot = "P") {
						data "P labile (final)" value: {simulations mean_of each.param_getter.param_value(), simulations mean_of each.output[S_Strategist].P_labile};
						data "P avail" value: {simulations mean_of each.param_getter.param_value(), simulations mean_of each.output[S_Strategist].P_avail};
						data "P labile (init)" value: {simulations mean_of each.param_getter.param_value(), simulations mean_of (each.C_labile / each.CP_labile) / #gram} marker: false;
					}
				}
			}
		}
	}
}

/**
 * Base virtual experiment used to explore the variation of the labile C/N rate, with an empty dam.
 */
experiment ExploreCN_labile parent:Explore {
	parameter "X var" var:x_var init:"CN labile" read_only:true;
	
	float min_CN_labile <- 4.0;
	float max_CN_labile <- 70.0;
	parameter "Labile C/N" var: CN_labile min: min_CN_labile max: max_CN_labile step: 1.0;
	
	init {
		plot <- "N";
	}
	string x_label {
		return "C/N (labile)";
	}
}

// Scenario 1
/**
 * Experiment used to explore the labile C/N rate with an "high" labile C availability of 1 gram.
 */
experiment ExploreCN_High_C_labile type:batch parent:ExploreCN_labile until: cycle>0 repeat: 5 {
	parameter "Output file" var:output_file init: "CN_high_C_labile" read_only:true;
	parameter "Labile C" var:exp_C_labile init:3#gram read_only:true;
	
	init {
		exp_name <- "CN_high_C_labile";
		do write_config_csv;
		do init_data_csv;
	}
}

// Scenario 2
/**
 * Experiment used to explore the labile C/N rate with a "low" labile C availability of 0.1 gram.
 */
experiment ExploreCN_Low_C_labile type:batch parent:ExploreCN_labile until: cycle>0 repeat: 5 {
	parameter "Output file" var:output_file init: "CN_low_C_labile";
	parameter "Labile C" var:exp_C_labile init:0.1#gram read_only:true;
	
	init {
		exp_name <- "CN_low_C_labile";
		do write_config_csv;
		do init_data_csv;
	}
}

/**
 * Base virtual experiment used to explore the variation of the labile C/P rate, with an empty dam.
 */
experiment ExploreCP_labile parent:Explore {
	parameter "X var" var:x_var init:"CP labile" read_only:true;
	
	float min_CP_labile <- CP_labile;
	float max_CP_labile <- CP_labile;
	parameter "Labile C/P" var: CP_labile min: min_CP_labile max: max_CP_labile step: 1.0;
	
	init {
		plot <- "P";
	}
	string x_label {
		return "C/P (labile)";
	}
}

// Scenario 3
/**
 * Experiment used to explore the labile C/P rate with an "high" labile C availability of 1 gram.
 */
experiment ExploreCP_High_C_labile type:batch parent:ExploreCP_labile until: cycle>0 repeat: 5 {
	parameter "Output file" var:output_file init: "CP_high_C_labile";
	parameter "Labile C" var:exp_C_labile init:3#gram read_only:true;
	
	init {
		min_CP_labile <- 1.0;
		max_CP_labile <- 60.0;
		exp_name <- "CP_high_C_labile";
		do write_config_csv;
		do init_data_csv;
	}
}

// Scenario 4
/**
 * Experiment used to explore the labile C/P rate with a "low" labile C availability of 0.1 gram.
 */
experiment ExploreCP_Low_C_labile type:batch parent:ExploreCP_labile until: cycle>0 repeat: 5 {
	parameter "Output file" var:output_file init: "CP_low_C_labile";
	parameter "Labile C" var:exp_C_labile init:0.1#gram read_only:true;
	
	init {
		min_CP_labile <- 1.0;
		max_CP_labile <- 60.0;
		exp_name <- "CP_low_C_labile";
		do write_config_csv;
		do init_data_csv;
	}
}
/**
 * Base virtual experiment used to explore the variation of the dom C/N rate.
 */
experiment ExploreCN_dom parent:Explore type:batch until:cycle > 0 repeat: 5 {
	parameter "X var" var:x_var init:"CN dom" read_only:true;
	parameter "Output file" var:output_file init: "CN_dom";
	
	float min_CN_dom <- CN_dom;
	float max_CN_dom <- CN_dom;
	parameter "C/N DOM" var: CN_dom min: min_CN_dom max: max_CN_dom step: 1.0;
	parameter "Labile C" var:exp_C_labile init:1#gram read_only:true;
	parameter "C DOM" var:exp_C_dom init:0.3#gram read_only:true;
	
	string x_label {
		return "C/N (dom)";
	}
	
	init {
		min_CN_dom <- 1.0;
		max_CN_dom <- 40.0;
		exp_name <- "CN_dom";
		do write_config_csv;
		do init_data_csv;
	}
}

/**
 * Base virtual experiment used to explore the variation of the dom C/P rate.
 */
experiment ExploreCP_dom parent:Explore type:batch until:cycle > 0 repeat: 5 {
	parameter "X var" var:x_var init:"CP dom" read_only:true;
	parameter "Output file" var:output_file init: "CP_dom";
	
	float min_CP_dom <- CP_dom;
	float max_CP_dom <- CP_dom;
	parameter "C/P DOM" var: CP_dom min: min_CP_dom max: max_CP_dom step: 1.0;
	parameter "Labile C" var:exp_C_labile init:1#gram read_only:true;
	parameter "C DOM" var:exp_C_dom init:0.3#gram read_only:true;
	
	string x_label {
		return "C/P (dom)";
	}
	
	init {
		min_CP_dom <- 1.0;
		max_CP_dom <- 60.0;
		exp_name <- "CP_dom";
		do write_config_csv;
		do init_data_csv;
	}
}
