/**
* Name: enzymesactivityexploration
* Based on the internal skeleton template. 
* Author: pbreugno
* Tags: 
*/

model enzymatic_activity_exploration

import "../microbes.gaml"

global {
	/* Parameters to explore */
	float C_labile <- 0.01#gram;
	float CN_labile <- 20.0;
	float CP_labile <- 34.0;
	float C_recalcitrant <- 1#gram;
	float CN_recalcitrant <- 10.0;
	float CP_recalcitrant <- 20.0;
	float C_dom <- 0#gram;
	float CN_dom <- 10.0;
	float CP_dom <- 17.0;
	
	OrganicParticle organic_particle;
	PoreParticle pore_particle;
	/** Insert the global definitions, variables and actions here */
	
	
	map<species<MicrobePopulation>, DataOutput> output;
	
	init {
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
			C: 1#gram
		] {
			C_N <- 10.0;
			C_P <- 17.0;
			awake_population <- 1.0;
		}
		
		create A_Strategist with: [
			C: 1#gram
		] {
			C_N <- 10.0;
			C_P <- 17.0;
			awake_population <- 1.0;
		}
		
		create S_Strategist with: [
			C: 1#gram
		] {
			C_N <- 10.0;
			C_P <- 17.0;
			awake_population <- 1.0;
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
			create EnzymaticActivity {
				WeightedEnzymes weighted_enzymes;
				create WeightedEnzymes with: [
					T_cellulolytic::myself.enzymes.T_cellulolytic*myself.C_actif,
					T_amino::myself.enzymes.T_amino*myself.C_actif,
					T_P::myself.enzymes.T_P*myself.C_actif,
					T_recal::myself.enzymes.T_recal*myself.C_actif
				] {
					weighted_enzymes <- self;
				}
				
				do compute_activity(weighted_enzymes, myself.enzymatic_activity_problem);
				
				output[species(myself)].X_C <- (X_C_cellulolytic + X_C_amino + X_C_P)/#gram;
				output[species(myself)].X_N <- X_N_amino/#gram;
				output[species(myself)].X_P <- (X_P_labile_to_dom + X_P_labile_to_dim)/#gram;
				output[species(myself)].X_recal <- X_C_recal/#gram;
				
				float C_avail; float N_avail; float P_avail;
				ask myself.enzymatic_activity_problem {
					C_avail <- C_avail(myself);
					N_avail <- N_avail(myself);
					P_avail <- P_avail(myself);
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
}

species DataOutput {
	// EnzymaticActivity enzymatic_activity;
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

experiment Explore {
	method exploration;
	string exp_name;
	string plot <- "N" among: ["N", "P", "C"];
	
	action init_header {
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
		save header header:false rewrite: true to: exp_name + ".csv" format: "csv";
	}
	
	string x_label virtual: true;
	float x_param(enzymatic_activity_exploration_model sim) virtual:true;
	
	reflex {
		int i <- 0;
		ask simulations {
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
			save [i, myself.x_param(self)] + data header: false rewrite:false to: myself.exp_name + ".csv" format: "csv";
			i <- i+1;
		}
	}
	
	permanent {
		display "Enzymatic activity (Y)" {
			chart "Enzymatic activity (Y)" type:xy x_label: x_label() y_label: "Enzymatic activity (gC/gM/d)" {
				data "T_P" value:{simulations mean_of x_param(each), simulations mean_of(sum(Y_Strategist collect (each.enzymes.T_P/(#gram/#gram/#d))))};
				data "T_cellulolytic" value:{simulations mean_of x_param(each), simulations mean_of(sum(Y_Strategist collect (each.enzymes.T_cellulolytic/(#gram/#gram/#d))))};
				data "T_amino" value:{simulations mean_of x_param(each), simulations mean_of(sum(Y_Strategist collect (each.enzymes.T_amino/(#gram/#gram/#d))))};
				data "T_recal" value:{simulations mean_of x_param(each), simulations mean_of(sum(Y_Strategist collect (each.enzymes.T_recal/(#gram/#gram/#d))))};
			}
		}
			
		display "Enzymatic activity (A)" {
			chart "Enzymatic activity (A)" type:xy x_label: x_label() y_label: "Enzymatic activity (gC/gM/d)" {
				data "T_P" value:{simulations mean_of x_param(each), simulations mean_of(sum(A_Strategist collect (each.enzymes.T_P/(#gram/#gram/#d))))};
				data "T_cellulolytic" value:{simulations mean_of x_param(each), simulations mean_of(sum(A_Strategist collect (each.enzymes.T_cellulolytic/(#gram/#gram/#d))))};
				data "T_amino" value:{simulations mean_of x_param(each), simulations mean_of(sum(A_Strategist collect (each.enzymes.T_amino/(#gram/#gram/#d))))};
				data "T_recal" value:{simulations mean_of x_param(each), simulations mean_of(sum(A_Strategist collect (each.enzymes.T_recal/(#gram/#gram/#d))))};
			}
		}
			
		display "Enzymatic activity (S)" {
			chart "Enzymatic activity (S)" type:xy x_label: x_label() y_label: "Enzymatic activity (gC/gM/d)" {
				data "T_P" value:{simulations mean_of x_param(each), simulations mean_of(sum(S_Strategist collect (each.enzymes.T_P/(#gram/#gram/#d))))};
				data "T_cellulolytic" value:{simulations mean_of x_param(each), simulations mean_of(sum(S_Strategist collect (each.enzymes.T_cellulolytic/(#gram/#gram/#d))))};
				data "T_amino" value:{simulations mean_of x_param(each), simulations mean_of(sum(S_Strategist collect (each.enzymes.T_amino/(#gram/#gram/#d))))};
				data "T_recal" value:{simulations mean_of x_param(each), simulations mean_of(sum(S_Strategist collect (each.enzymes.T_recal/(#gram/#gram/#d))))};
			}
		}
		
//		display "Decomposition (Y)" {
//			chart "Decomposition (Y)" type:xy x_label: x_label() y_label: "Decomposition (g)" {
//				data "X C" value:{simulations mean_of x_param(each), simulations mean_of each.output[Y_Strategist].X_C};
//				data "X N" value:{simulations mean_of x_param(each), simulations mean_of each.output[Y_Strategist].X_N};
//				data "X P" value:{simulations mean_of x_param(each), simulations mean_of each.output[Y_Strategist].X_P};
//				data "X recal" value:{simulations mean_of x_param(each), simulations mean_of each.output[Y_Strategist].X_recal};
//				data "C" value:{simulations mean_of x_param(each), simulations mean_of each.C_labile/#gram} marker:false;
//				data "N" value:{simulations mean_of x_param(each), simulations mean_of (each.C_labile/each.CN_labile)/#gram} marker:false;
//			}
//		}
//		display "Decomposition (A)" {
//			chart "Decomposition (A)" type:xy x_label: x_label() y_label: "Decomposition (g)" {
//				data "X C" value:{simulations mean_of x_param(each), simulations mean_of each.output[A_Strategist].X_C};
//				data "X N" value:{simulations mean_of x_param(each), simulations mean_of each.output[A_Strategist].X_N};
//				data "X P" value:{simulations mean_of x_param(each), simulations mean_of each.output[A_Strategist].X_P};
//				data "X recal" value:{simulations mean_of x_param(each), simulations mean_of each.output[A_Strategist].X_recal};
//				data "C" value:{simulations mean_of x_param(each), simulations mean_of each.C_labile/#gram} marker:false;
//				data "N" value:{simulations mean_of x_param(each), simulations mean_of (each.C_labile/each.CN_labile)/#gram} marker:false;
//			}
//		}
//		display "Decomposition (S)" {
//			chart "Decomposition (S)" type:xy x_label: x_label() y_label: "Decomposition (g)" {
//				data "X C" value:{simulations mean_of x_param(each), simulations mean_of each.output[S_Strategist].X_C};
//				data "X N" value:{simulations mean_of x_param(each), simulations mean_of each.output[S_Strategist].X_N};
//				data "X P" value:{simulations mean_of x_param(each), simulations mean_of each.output[S_Strategist].X_P};
//				data "X recal" value:{simulations mean_of x_param(each), simulations mean_of each.output[S_Strategist].X_recal};
//				data "C" value:{simulations mean_of x_param(each), simulations mean_of each.C_labile/#gram} marker:false;
//				data "N" value:{simulations mean_of x_param(each), simulations mean_of (each.C_labile/each.CN_labile)/#gram} marker:false;
//			}
//		}
		display "Decomposition (Y)" {
			chart "Decomposition (Y)" type:xy x_label: x_label() y_label: "Decomposition (g)" {
				if(plot = "N") {
					data "N labile (final)" value:{simulations mean_of x_param(each), simulations mean_of each.output[Y_Strategist].N_labile};
					data "N avail" value:{simulations mean_of x_param(each), simulations mean_of each.output[Y_Strategist].N_avail};
					data "N labile (init)" value:{simulations mean_of x_param(each), simulations mean_of (each.C_labile/each.CN_labile)/#gram} marker:false;
				} else if (plot = "P") {
					data "P labile (final)" value:{simulations mean_of x_param(each), simulations mean_of each.output[Y_Strategist].P_labile};
					data "P avail" value:{simulations mean_of x_param(each), simulations mean_of each.output[Y_Strategist].P_avail};
					data "P labile (init)" value:{simulations mean_of x_param(each), simulations mean_of (each.C_labile/each.CP_labile)/#gram} marker:false;	
				}
//				data "C recal" value:{simulations mean_of x_param(each), simulations mean_of each.output[Y_Strategist].C_recal};
//				data "C labile" value:{simulations mean_of x_param(each), simulations mean_of each.output[Y_Strategist].C_labile};
//				data "C avail" value:{simulations mean_of x_param(each), simulations mean_of each.output[Y_Strategist].C_avail};
//				data "C" value:{simulations mean_of x_param(each), simulations mean_of each.C_labile/#gram} marker:false;
			}
		}
		display "Decomposition (A)" {
			chart "Decomposition (A)" type:xy x_label: x_label() y_label: "Decomposition (g)" {
				if(plot = "N") {
					data "N labile (final)" value:{simulations mean_of x_param(each), simulations mean_of each.output[A_Strategist].N_labile};
					data "N avail" value:{simulations mean_of x_param(each), simulations mean_of each.output[A_Strategist].N_avail};
					data "N labile (init)" value:{simulations mean_of x_param(each), simulations mean_of (each.C_labile/each.CN_labile)/#gram} marker:false;
				} else if (plot = "P") {
					data "P labile (final)" value:{simulations mean_of x_param(each), simulations mean_of each.output[A_Strategist].P_labile};
					data "P avail" value:{simulations mean_of x_param(each), simulations mean_of each.output[A_Strategist].P_avail};
					data "P labile (init)" value:{simulations mean_of x_param(each), simulations mean_of (each.C_labile/each.CP_labile)/#gram} marker:false;	
				}
//				data "C recal" value:{simulations mean_of x_param(each), simulations mean_of each.output[A_Strategist].C_recal};
//				data "C labile" value:{simulations mean_of x_param(each), simulations mean_of each.output[A_Strategist].C_labile};
//				data "C avail" value:{simulations mean_of x_param(each), simulations mean_of each.output[A_Strategist].C_avail};
//				data "C" value:{simulations mean_of x_param(each), simulations mean_of each.C_labile/#gram} marker:false;
			}
		}
		display "Decomposition (S)" {
			chart "Decomposition (S)" type:xy x_label: x_label() y_label: "Decomposition (g)" {
				if(plot = "N") {
					data "N labile (final)" value:{simulations mean_of x_param(each), simulations mean_of each.output[S_Strategist].N_labile};
					data "N avail" value:{simulations mean_of x_param(each), simulations mean_of each.output[S_Strategist].N_avail};
					data "N labile (init)" value:{simulations mean_of x_param(each), simulations mean_of (each.C_labile/each.CN_labile)/#gram} marker:false;
				} else if (plot = "P") {
					data "P labile (final)" value:{simulations mean_of x_param(each), simulations mean_of each.output[S_Strategist].P_labile};
					data "P avail" value:{simulations mean_of x_param(each), simulations mean_of each.output[S_Strategist].P_avail};
					data "P labile (init)" value:{simulations mean_of x_param(each), simulations mean_of (each.C_labile/each.CP_labile)/#gram} marker:false;	
				}
//				data "C recal" value:{simulations mean_of x_param(each), simulations mean_of each.output[S_Strategist].C_recal};
//				data "C labile" value:{simulations mean_of x_param(each), simulations mean_of each.output[S_Strategist].C_labile};
//				data "C avail" value:{simulations mean_of x_param(each), simulations mean_of each.output[S_Strategist].C_avail};
//				data "P avail" value:{simulations mean_of x_param(each), simulations mean_of each.output[S_Strategist].P_avail};
			}
		}
	}
}

experiment ExploreCP_labile parent:Explore type:batch until: cycle>0 repeat: 10 {
	float min_CP_labile <- CP_labile;
	float max_CP_labile <- CP_labile;
	parameter "CP" var: CP_labile min: min_CP_labile max: max_CP_labile step: 1.0;
	
	init {
		plot <- "P";
	}
	string x_label {
		return "C/P (labile)";
	}
	
	float x_param(enzymatic_activity_exploration_model sim) {
		return sim.CP_labile;
	}
}

experiment ExploreCP_High_C_labile type:batch parent:ExploreCP_labile until: cycle>0 repeat: 20 {
	init {
		min_CP_labile <- 1.0;
		max_CP_labile <- 30.0;
		exp_name <- "explore_CP_high_C_labile";
		C_labile <- 1#gram;
		do init_header;
	}
}

experiment ExploreCP_Low_C_labile type:batch parent:ExploreCP_labile until: cycle>0 repeat: 20 {
	init {
		min_CP_labile <- 1.0;
		max_CP_labile <- 60.0;
		exp_name <- "explore_CP_low_C_labile";
		C_labile <- 0.1#gram;
		do init_header;
	}
}

experiment ExploreCN_labile parent:Explore type:batch until:true {
	float min_CN_labile <- CN_labile;
	float max_CN_labile <- CN_labile;
	
	parameter "CN labile" var: CN_labile min: min_CN_labile max: max_CN_labile step: 1.0;
	
	init {
		plot <- "N";
	}
	string x_label {
		return "C/N (labile)";
	}

	float x_param(enzymatic_activity_exploration_model sim) {
		return sim.CN_labile;
	}
}

experiment ExploreCN_High_C_labile type:batch parent:ExploreCN_labile until: cycle>0 repeat: 20 {
	init {
		min_CN_labile <- 4.0;
		max_CN_labile <- 70.0;
		exp_name <- "explore_CN_high_C_labile";
		C_labile <- 1#gram;
		do init_header;
	}
}

experiment ExploreCN_Low_C_labile type:batch parent:ExploreCN_labile until: cycle>0 repeat: 20 {
	init {
		min_CN_labile <- 4.0;
		max_CN_labile <- 70.0;
		exp_name <- "explore_CN_low_C_labile";
		C_labile <- 0.1#gram;
		do init_header;
	}
}

experiment ExploreCN_dom parent:Explore type:batch until:cycle > 0 repeat: 1 {
	float min_CN_dom <- CN_dom;
	float max_CN_dom <- CN_dom;
	
	parameter "CN dom" var: CN_dom min: min_CN_dom max: max_CN_dom step: 1.0;
	
	string x_label {
		return "C/N (dom)";
	}

	float x_param(enzymatic_activity_exploration_model sim) {
		return sim.CN_dom;
	}
	
	init {
		min_CN_dom <- 1.0;
		max_CN_dom <- 60.0;
		exp_name <- "explore_CN_dom";
		C_labile <- 1#gram;
		C_dom <- 0.1#gram;
		CP_dom <- 22.0;
		do init_header;
	}
}
