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
//		write "Init model parameters";
//		// Init model parameters
//		enzymes_optimization_period <- 1#d;
//		
//		max_T_cellulolytic_Y <- 1 #gram/ #gram / #d;
//		max_T_amino_Y <- 0.2 #gram/ #gram / #d;
//		max_T_P_Y <- 0.1 #gram/ #gram / #d;
//		max_T_recal_Y <- 0.00001 #gram/ #gram / #d;
//
//		max_T_cellulolytic_A <- 4 #gram/ #gram / #d;
//		max_T_amino_A <- 6 #gram/ #gram / #d;
//		max_T_P_A <- 1 #gram/ #gram / #d;
//		max_T_recal_A <- 2 #gram/ #gram / #d;
//		
//		max_T_cellulolytic_S <- 0.2 #gram/ #gram / #d;
//		max_T_amino_S <- 0.05 #gram/ #gram / #d;
//		max_T_P_S <- 0.01 #gram/ #gram / #d;
//		max_T_recal_S <- 0.02 #gram/ #gram / #d;
//		
//		// carbon_use_efficiency_Y <- 1.0;
//		carbon_use_efficiency_A <- 1.0;
//		carbon_use_efficiency_S <- 1.0;
//		
//		microbes_CN <- 10.0;
//		microbes_CP <- 17.0;
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

		create Y_Strategist with: [
			C: 1#gram,
			awake_population: 1.0
		] {
			do update;
		}
		
		create A_Strategist with: [
			C: 1#gram,
			awake_population: 1.0
		] {
			do update;
		}
		
		create S_Strategist with: [
			C: 1#gram,
			awake_population: 1.0
		] {
			do update;
		}
		
		loop s over: [Y_Strategist, A_Strategist, S_Strategist] {
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
					T_cellulolytic::myself.enzymes.T_cellulolytic*myself.active_C,
					T_amino::myself.enzymes.T_amino*myself.active_C,
					T_P::myself.enzymes.T_P*myself.active_C,
					T_recal::myself.enzymes.T_recal*myself.active_C
				] {
					weighted_enzymes <- self;
				}
				
				ask decomposition_problem {
					do decomposition(weighted_enzymes, myself);
				}
				
				output[species(myself)].T_r <- myself.enzymes.T_recal;
				output[species(myself)].T_C <- myself.enzymes.T_cellulolytic;
				output[species(myself)].T_N <- myself.enzymes.T_amino;
				output[species(myself)].T_P <- myself.enzymes.T_P;
				
				output[species(myself)].X_C <- (X_C_cellulolytic + X_C_amino + X_C_P)/#gram;
				output[species(myself)].X_N <- X_N_amino/#gram;
				output[species(myself)].X_P <- (X_P_labile_to_dom + X_P_labile_to_dim)/#gram;
				output[species(myself)].X_recal <- X_C_recal/#gram;
				
				float C_avail; float N_avail; float P_avail;
				ask decomposition_problem {
					C_avail <- C_avail_final(myself);
					N_avail <- N_avail_final(myself);
					P_avail <- P_avail_final(myself);
				}
				
				output[species(myself)].C_recal <- (world.C_recalcitrant - X_C_recal)/#gram;
				
				output[species(myself)].C_labile <- (world.C_labile + X_C_recal - (X_C_cellulolytic + X_C_amino + X_C_P))/#gram;
				output[species(myself)].N_labile <- (world.C_labile/world.CN_labile + X_N_recal - X_N_amino)/#gram;
				output[species(myself)].P_labile <- (world.C_labile/world.CP_labile + X_P_recal_to_labile - X_P_labile_to_dom - X_P_labile_to_dim)/#gram;
				
				output[species(myself)].C_avail <- (world.C_dom + X_C_cellulolytic + X_C_amino + X_C_P)/#gram;
				output[species(myself)].N_avail <- (world.C_dom/world.CN_dom + X_N_amino)/#gram;
				output[species(myself)].P_avail <- (world.C_dom/world.CP_dom + X_P_labile_to_dom + X_P_labile_to_dim + X_P_recal_to_dim)/#gram;
				
				ask weighted_enzymes {
					do die;
				}
				do die;
			}
		}
	}
	
	reflex {
		list<float> data <- [];
		ask Y_Strategist + A_Strategist + S_Strategist {
			data <- data + [
				self.enzymes.T_recal/(#gram/#gram/#d),
				self.enzymes.T_cellulolytic/(#gram/#gram/#d),
				self.enzymes.T_amino/(#gram/#gram/#d),
				self.enzymes.T_P/(#gram/#gram/#d)
			];

			data <- data + [myself.output[species(self)].X_C, myself.output[species(self)].X_N, myself.output[species(self)].X_P, myself.output[species(self)].X_recal];
				
			data <- data + [
				myself.output[species(self)].C_recal,
				myself.output[species(self)].C_labile, myself.output[species(self)].N_labile, myself.output[species(self)].P_labile,
				myself.output[species(self)].C_avail, myself.output[species(self)].N_avail, myself.output[species(self)].P_avail
			];
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
	 * Parameters used to plot P or N decomposition, depending on the parameter currently explored.
	 */
	string plot <- "N" among: ["N", "P"];
	
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
			C_labile/#gram,
			CN_labile,
			CP_labile,
			C_recalcitrant/#gram,
			CN_recalcitrant,
			CP_recalcitrant,
			C_dom/#gram,
			CN_dom,
			CP_dom
		];
		
		loop s over: [Y_Strategist, A_Strategist, S_Strategist] {
			add "Requested C/N (" + labels[s] + ")" to: header;
			add "Requested C/P (" + labels[s] + ")" to: header;
			add "Max T_cellulolytic (" + labels[s] + ")" to: header;
			add "Max T_amino (" + labels[s] + ")" to: header;
			add "Max T_P (" + labels[s] + ")" to: header;
			add "Max T_recal (" + labels[s] + ")" to: header;
		}
		add microbes_CN / carbon_use_efficiency_Y to: config;
		add microbes_CP / carbon_use_efficiency_Y to: config;
		add max_T_cellulolytic_Y/(#gram/#gram/#d) to: config;
		add max_T_amino_Y/(#gram/#gram/#d) to: config;
		add max_T_P_Y/(#gram/#gram/#d) to: config;
		add max_T_recal_Y/(#gram/#gram/#d) to: config;
		
		add microbes_CN / carbon_use_efficiency_A to: config;
		add microbes_CP / carbon_use_efficiency_A to: config;
		add max_T_cellulolytic_A/(#gram/#gram/#d) to: config;
		add max_T_amino_A/(#gram/#gram/#d) to: config;
		add max_T_P_A/(#gram/#gram/#d) to: config;
		add max_T_recal_A/(#gram/#gram/#d) to: config;
		
		add microbes_CN / carbon_use_efficiency_S to: config;
		add microbes_CP / carbon_use_efficiency_S to: config;
		add max_T_cellulolytic_S/(#gram/#gram/#d) to: config;
		add max_T_amino_S/(#gram/#gram/#d) to: config;
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
				"T cellulolytic (" + strategy + ")",
				"T amino (" + strategy + ")",
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
	
	permanent {
		display "Enzymatic activity (Y)" {
			chart "Enzymatic activity (Y)" type:xy x_label: x_label() y_label: "Enzymatic activity (gC/gM/d)" {
				data "T_P" value:{mean(simulations collect each.param_getter.param_value()), mean(simulations collect (each.output[Y_Strategist].T_P/(#gram/#gram/#d)))};
				data "T_cellulolytic" value:{mean(simulations collect each.param_getter.param_value()), mean(simulations collect (each.output[Y_Strategist].T_C/(#gram/#gram/#d)))};
				data "T_amino" value:{mean(simulations collect each.param_getter.param_value()), mean(simulations collect (each.output[Y_Strategist].T_N/(#gram/#gram/#d)))};
				data "T_recal" value:{mean(simulations collect each.param_getter.param_value()), mean(simulations collect (each.output[Y_Strategist].T_r/(#gram/#gram/#d)))};
			}
		}
			
		display "Enzymatic activity (A)" {
			chart "Enzymatic activity (A)" type:xy x_label: x_label() y_label: "Enzymatic activity (gC/gM/d)" {
				data "T_P" value:{mean(simulations collect each.param_getter.param_value()), mean(simulations collect (each.output[A_Strategist].T_P/(#gram/#gram/#d)))};
				data "T_cellulolytic" value:{mean(simulations collect each.param_getter.param_value()), mean(simulations collect (each.output[A_Strategist].T_C/(#gram/#gram/#d)))};
				data "T_amino" value:{mean(simulations collect each.param_getter.param_value()), mean(simulations collect (each.output[A_Strategist].T_N/(#gram/#gram/#d)))};
				data "T_recal" value:{mean(simulations collect each.param_getter.param_value()), mean(simulations collect (each.output[A_Strategist].T_r/(#gram/#gram/#d)))};
			}
		}
			
		display "Enzymatic activity (S)" {
			chart "Enzymatic activity (S)" type:xy x_label: x_label() y_label: "Enzymatic activity (gC/gM/d)" {
				data "T_P" value:{mean(simulations collect each.param_getter.param_value()), mean(simulations collect (each.output[S_Strategist].T_P/(#gram/#gram/#d)))};
				data "T_cellulolytic" value:{mean(simulations collect each.param_getter.param_value()), mean(simulations collect (each.output[S_Strategist].T_C/(#gram/#gram/#d)))};
				data "T_amino" value:{mean(simulations collect each.param_getter.param_value()), mean(simulations collect (each.output[S_Strategist].T_N/(#gram/#gram/#d)))};
				data "T_recal" value:{mean(simulations collect each.param_getter.param_value()), mean(simulations collect (each.output[S_Strategist].T_r/(#gram/#gram/#d)))};
			}
		}
		
		display "Decomposition (Y)" {
			chart "Decomposition (Y)" type:xy x_label: x_label() y_label: "Decomposition (g)" {
				if(plot = "N") {
					data "N labile (final)" value:{simulations mean_of each.param_getter.param_value(), simulations mean_of each.output[Y_Strategist].N_labile};
					data "N avail" value:{simulations mean_of each.param_getter.param_value(), simulations mean_of each.output[Y_Strategist].N_avail};
					data "N labile (init)" value:{simulations mean_of each.param_getter.param_value(), simulations mean_of (each.C_labile/each.CN_labile)/#gram} marker:false;
				} else if (plot = "P") {
					data "P labile (final)" value:{simulations mean_of each.param_getter.param_value(), simulations mean_of each.output[Y_Strategist].P_labile};
					data "P avail" value:{simulations mean_of each.param_getter.param_value(), simulations mean_of each.output[Y_Strategist].P_avail};
					data "P labile (init)" value:{simulations mean_of each.param_getter.param_value(), simulations mean_of (each.C_labile/each.CP_labile)/#gram} marker:false;	
				}
			}
		}
		display "Decomposition (A)" {
			chart "Decomposition (A)" type:xy x_label: x_label() y_label: "Decomposition (g)" {
				if(plot = "N") {
					data "N labile (final)" value:{simulations mean_of each.param_getter.param_value(), simulations mean_of each.output[A_Strategist].N_labile};
					data "N avail" value:{simulations mean_of each.param_getter.param_value(), simulations mean_of each.output[A_Strategist].N_avail};
					data "N labile (init)" value:{simulations mean_of each.param_getter.param_value(), simulations mean_of (each.C_labile/each.CN_labile)/#gram} marker:false;
				} else if (plot = "P") {
					data "P labile (final)" value:{simulations mean_of each.param_getter.param_value(), simulations mean_of each.output[A_Strategist].P_labile};
					data "P avail" value:{simulations mean_of each.param_getter.param_value(), simulations mean_of each.output[A_Strategist].P_avail};
					data "P labile (init)" value:{simulations mean_of each.param_getter.param_value(), simulations mean_of (each.C_labile/each.CP_labile)/#gram} marker:false;	
				}
			}
		}
		display "Decomposition (S)" {
			chart "Decomposition (S)" type:xy x_label: x_label() y_label: "Decomposition (g)" {
				if(plot = "N") {
					data "N labile (final)" value:{simulations mean_of each.param_getter.param_value(), simulations mean_of each.output[S_Strategist].N_labile};
					data "N avail" value:{simulations mean_of each.param_getter.param_value(), simulations mean_of each.output[S_Strategist].N_avail};
					data "N labile (init)" value:{simulations mean_of each.param_getter.param_value(), simulations mean_of (each.C_labile/each.CN_labile)/#gram} marker:false;
				} else if (plot = "P") {
					data "P labile (final)" value:{simulations mean_of each.param_getter.param_value(), simulations mean_of each.output[S_Strategist].P_labile};
					data "P avail" value:{simulations mean_of each.param_getter.param_value(), simulations mean_of each.output[S_Strategist].P_avail};
					data "P labile (init)" value:{simulations mean_of each.param_getter.param_value(), simulations mean_of (each.C_labile/each.CP_labile)/#gram} marker:false;	
				}
			}
		}
	}
}

/**
 * Base virtual experiment used to explore the variation of the labile C/N rate, with an empty dam.
 */
experiment ExploreCN_labile parent:Explore {
	parameter "X var" var:x_var init:"CN labile";
	
	float min_CN_labile <- CN_labile;
	float max_CN_labile <- CN_labile;
	parameter "CN labile" var: CN_labile min: min_CN_labile max: max_CN_labile step: 1.0;
	
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
	parameter "Output file" var:output_file init: "CN_high_C_labile";
	init {
		min_CN_labile <- 4.0;
		max_CN_labile <- 70.0;
		exp_name <- "CN_high_C_labile";
		C_labile <- 8#gram;
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
	init {
		min_CN_labile <- 4.0;
		max_CN_labile <- 70.0;
		exp_name <- "CN_low_C_labile";
		C_labile <- 0.1#gram;
		do write_config_csv;
		do init_data_csv;
	}
}

/**
 * Base virtual experiment used to explore the variation of the labile C/P rate, with an empty dam.
 */
experiment ExploreCP_labile parent:Explore {
	parameter "X var" var:x_var init:"CP labile";
	
	float min_CP_labile <- CP_labile;
	float max_CP_labile <- CP_labile;
	parameter "CP" var: CP_labile min: min_CP_labile max: max_CP_labile step: 1.0;
	
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
	
	init {
		min_CP_labile <- 1.0;
		max_CP_labile <- 60.0;
		exp_name <- "CP_high_C_labile";
		C_labile <- 8#gram;
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
	
	init {
		min_CP_labile <- 1.0;
		max_CP_labile <- 60.0;
		exp_name <- "CP_low_C_labile";
		C_labile <- 0.1#gram;
		do write_config_csv;
		do init_data_csv;
	}
}
/**
 * Base virtual experiment used to explore the variation of the dom C/N rate.
 */
experiment ExploreCN_dom parent:Explore type:batch until:cycle > 0 repeat: 5 {
	parameter "X var" var:x_var init:"CN dom";
	parameter "Output file" var:output_file init: "CN_dom";
	
	float min_CN_dom <- CN_dom;
	float max_CN_dom <- CN_dom;
	parameter "CN dom" var: CN_dom min: min_CN_dom max: max_CN_dom step: 1.0;
	
	string x_label {
		return "C/N (dom)";
	}
	
	init {
		min_CN_dom <- 1.0;
		max_CN_dom <- 40.0;
		exp_name <- "CN_dom";
		C_labile <- 1#gram;
		C_dom <- 0.3#gram;
		do write_config_csv;
		do init_data_csv;
	}
}

/**
 * Base virtual experiment used to explore the variation of the dom C/P rate.
 */
experiment ExploreCP_dom parent:Explore type:batch until:cycle > 0 repeat: 5 {
	parameter "X var" var:x_var init:"CP dom";
	parameter "Output file" var:output_file init: "CP_dom";
	
	float min_CP_dom <- CP_dom;
	float max_CP_dom <- CP_dom;
	parameter "CP dom" var: CP_dom min: min_CP_dom max: max_CP_dom step: 1.0;
	
	string x_label {
		return "C/P (dom)";
	}
	
	init {
		min_CP_dom <- 1.0;
		max_CP_dom <- 60.0;
		exp_name <- "CP_dom";
		C_labile <- 1#gram;
		C_dom <- 0.3#gram;
		do write_config_csv;
		do init_data_csv;
	}
}
