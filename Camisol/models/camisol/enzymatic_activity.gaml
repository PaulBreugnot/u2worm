/**
* Name: enzymes
* Based on the internal skeleton template. 
* Author: pbreugno
* Tags: 
*/

model enzyme_activity

global {
	/** Insert the global definitions, variables and actions here */
	float amino_CN <- 6.0;
}

species Enzymes {
	// Enzymatic action rates
	float T_cellulolytic <- 0.0#gram/#gram/#h;
	float T_amino <- 0.0#gram/#gram/#h;
	float T_P <- 0.0#gram/#gram/#h;
	float T_recal <- 0.0#gram/#gram/#h;
}

species EnzymaticActivity {
	float X_C_cellulolytic <- 0.0;
	float X_C_amino <- 0.0;
	float X_N_amino <- 0.0;
	float X_C_P <- 0.0;
	float X_P_labile_to_dom <- 0.0;
	float X_P_labile_to_dim <- 0.0;
	
	float X_C_recal <- 0.0;
	float X_P_recal_to_dim <- 0.0;
	
		
	action compute_activity(
		Enzymes enzymes, EnzymaticActivityProblem problem
	) {
		ask problem {
			if(total_C_labile > 0) {
				// Cellulolytic action
				float E_C_cellulolytic <- min(
					(C_microbes * dt * enzymes.T_cellulolytic) * total_C_labile / (K_cellulolytic + total_C_labile),
					total_C_labile
				);
		
				float E_C_P <- min(
					(C_microbes * dt * enzymes.T_P) * total_C_labile / (K_P + total_C_labile),
					total_C_labile
				);
				
				float expected_C_amino <- (
					C_microbes * dt * enzymes.T_amino
				) * total_C_labile / (K_amino + total_C_labile);
				float expected_N_amino <- expected_C_amino / amino_CN;
				float limiting_amino <- expected_C_amino > 0 ? min(
					1.0, total_C_labile / expected_C_amino, total_N_labile / expected_N_amino 
				) : 0.0;
				// total_C_labile <- total_C_labile - expected_C_amino * limiting_amino;
				// total_N_labile <- total_N_labile - expected_N_amino * limiting_amino;
				float E_C_amino <- expected_C_amino * limiting_amino;
				float E_N_amino <- expected_N_amino * limiting_amino;
				
				myself.X_C_cellulolytic <- E_C_cellulolytic
					- beta_cellulolytic_amino * E_C_amino * E_C_cellulolytic / total_C_labile
					- beta_cellulolytic_P * E_C_P * E_C_cellulolytic / total_C_labile;
					
				myself.X_C_amino <- E_C_amino
					- (1-beta_cellulolytic_amino) * E_C_amino * E_C_cellulolytic / total_C_labile
					- beta_amino_P * E_C_amino * E_C_P / total_C_labile;
				myself.X_C_P <- E_C_P
					- (1-beta_cellulolytic_P) * E_C_P * E_C_cellulolytic / total_C_labile
					- (1-beta_amino_P) * E_C_amino * E_C_P / total_C_labile;
					
				myself.X_N_amino <- myself.X_C_amino / amino_CN;
				if(total_P_labile > 0) {
					float total_X_P <- myself.X_C_P / (total_C_labile / total_P_labile);
					
					myself.X_P_labile_to_dim <- alpha_P_labile_to_dim * total_X_P;
							
					// Only a fraction goes to the dom, the rest is left in the labile section
					myself.X_C_P <- alpha_CP * myself.X_C_P;
					myself.X_P_labile_to_dom <- alpha_CP * (total_X_P-myself.X_P_labile_to_dim);
				}
			}
	
			if(total_C_recal > 0) {
				float E_P_om_recal <- 0.0;
				float E_N_om_recal <- 0.0;
				float X_C_c_recal <- 0.0;
				float X_C_om_recal <- 0.0;
				
				// Recalcitrant C decomposition
				myself.X_C_recal <- min(
					C_microbes * dt * enzymes.T_recal * total_C_recal / (K_recal + total_C_recal),
					total_C_recal
				);
				if(total_P_recal > 0) {
					myself.X_P_recal_to_dim <-
						alpha_P_recal_to_dim * myself.X_C_recal / (total_C_recal / total_P_recal);
				}
			}
		}
	}
}

species SimulatedAnnealingState {
	EnzymaticActivityProblem problem;
	Enzymes enzymes;
	EnzymaticActivity enzymatic_activity;
	list<float> budget;


	init {
		enzymes <- build_enzyme();
		create EnzymaticActivity {
			myself.enzymatic_activity <- self;
			do compute_activity(myself.enzymes, myself.problem);
		}
	}
	
	Enzymes build_enzyme {
		Enzymes new_enzymes;
		ask problem {
			create Enzymes with: [
				T_cellulolytic: min_enzymes.T_cellulolytic + myself.budget[0]*(max_enzymes.T_cellulolytic - min_enzymes.T_cellulolytic),
				T_amino: min_enzymes.T_amino + myself.budget[1]*(max_enzymes.T_amino - min_enzymes.T_amino),
				T_P: min_enzymes.T_P + myself.budget[2]*(max_enzymes.T_P - min_enzymes.T_P),
				T_recal: max_enzymes.T_recal > 0 ? min_enzymes.T_recal + myself.budget[3]*(max_enzymes.T_recal - min_enzymes.T_recal) : 0.0
			] {
				new_enzymes <- self;
			}
		}
		return new_enzymes;
	}

	float C_avail {
		return problem.C_DOM + enzymatic_activity.X_C_cellulolytic + enzymatic_activity.X_C_amino + enzymatic_activity.X_C_P;
	}
	
	float N_avail {
		return problem.N_DOM + problem.N_DIM + enzymatic_activity.X_N_amino;
	}
	
	float P_avail {
		return problem.P_DOM + problem.P_DIM + enzymatic_activity.X_P_labile_to_dom + enzymatic_activity.X_P_labile_to_dim;
	}
}

species Objective {
	float weight <- 1.0;
	
	action value(SimulatedAnnealingState state) virtual: true type: float;
}

species EnzymaticActivityProblem {
	float C_N;
	float C_P;
	float C_microbes;
	float dt;
	
	float total_C_labile;
	float total_N_labile;
	float total_P_labile;
	float total_C_recal;
	float total_N_recal;
	float total_P_recal;
	
	float C_DOM;
	float P_DOM;
	float N_DOM;
	float P_DIM;
	float N_DIM;
	
	float beta_cellulolytic_amino <- 0.5;
	float beta_cellulolytic_P <- 0.5;
	float beta_amino_P <- 0.5;
	
	float alpha_CP <- 0.0;
	float alpha_P_recal_to_dim <- 1e-3;
	float alpha_P_labile_to_dim <- 0.9;
	
	Enzymes min_enzymes;
	Enzymes max_enzymes;
	
	// The smallest K is, the more powerful the enzyme is. When K=0, the enzyme instantly imposes the limiting reaction rate, whatever the substrate concentration is.
	// Concentration at which the reaction rate from [labile C] to [C dom] is half of the limiting reaction rate.
	float K_cellulolytic <- 0.0;
	// Concentration at which the reaction rate from [labile C] to [CN dom] is half of the limiting reaction rate.
	float K_amino <- 0.0;
	// Concentration at which the reaction rate from [labile C] to [P dim] is half of the limiting reaction rate.
	float K_P <- 0.0;
	// Concentration at which the reaction rate from [recal C/N/P] to [labile C/N/P] is half of the limiting reaction rate.
	float K_recal <- 0.0;
}

species SimulatedAnnealing {
	EnzymaticActivityProblem problem;
	
	int N <- 1000;
	float epsilon <- 0.0;
	
	float T_init;
	float T;
	float e;
	SimulatedAnnealingState s;
	list<Objective> objectives;
	
	list<float> init_budget <- [1.0/4, 1.0/4, 1.0/4];
	
	init {
		T_init <- float(length(objectives));
		T <- T_init;
		if(problem.max_enzymes.T_recal > 0.0) {
			add 1.0/4 to: init_budget;
		}
		
		create SimulatedAnnealingState with:[problem::problem, budget::init_budget] {
			myself.s <- self;
		}
		e <- E(s);
	}
	
	action step {
		SimulatedAnnealingState s_new <- neighbour();
		float e_new <- E(s_new);
		if e_new < e or rnd(1.0) < exp(-(e_new - e) / T) {
			ask s {
				ask enzymes {
					do die;
				}
				ask enzymatic_activity {
					do die;
				}
				do die;
			}
			s <- s_new;
			e <- e_new;
		} else {
			ask s_new {
				ask enzymes {
					do die;
				}
				ask enzymatic_activity {
					do die;
				}
				do die;
			}
		}
		T <- 0.99 * T;
	}
	
	action optimize {
		int i <- 0;
		loop while: i < N and e > epsilon {
			i <- i+1;
			do step;
		}
	}
	
	SimulatedAnnealingState neighbour {
		SimulatedAnnealingState result;
		
		list<float> _T;
		loop item over: s.budget {
			add item to: _T;
		}

		// The indexes list is shuffled to determine in which order each component will vary.
		// This the total budget (1.0) is fixed, this allows a fair access to the budget for
		// each component (i.e. the last component does not always get the last part left after
		// all other get their own part).
		list<int> indexes;
		loop i from: 0 to: length(_T)-1 {
			add i to: indexes;
		}
		indexes <- shuffle(indexes);
		float delta <- 0.1;
		float range <- 1.0;
		loop i from: 0 to: length(indexes)-2 {
			_T[indexes[i]] <- max(0, min(range, gauss(_T[indexes[i]], delta)));
			range <- range - _T[indexes[i]];
		}
		_T[indexes[length(indexes)-1]] <- range;
		
		create SimulatedAnnealingState with:[problem::problem, budget::_T] {
			result <- self;
		}
		return result;
	}
	
	float E(
		SimulatedAnnealingState state
	) {
		return sum(objectives collect (each.weight * each.value(state)))/sum(objectives collect each.weight);
	}
}

species ExactCN parent: Objective {
	action value(SimulatedAnnealingState state) type: float {
		float N_avail <- state.N_avail();
		float C_avail <- state.C_avail();
		if(N_avail > 0 and C_avail > 0) {
			return (1.0 - state.problem.C_N / (C_avail / N_avail)) ^ 2;
		}
		return #infinity;
	}
}

species MaxCN parent: Objective {
	action value(SimulatedAnnealingState state) type: float {
		float N_avail <- state.N_avail();
		float C_avail <- state.C_avail();
		if(N_avail > 0 and C_avail > 0) {
			return (1.0 - state.problem.C_N / max(C_avail / N_avail, state.problem.C_N)) ^ 2;
		}
		return #infinity;
	}
}

species ExactCP parent: Objective {
	action value(SimulatedAnnealingState state) type: float {
		float P_avail <- state.P_avail();
		float C_avail <- state.C_avail();
		if(P_avail > 0 and C_avail > 0) {
			return (1.0 - state.problem.C_P / (C_avail / P_avail)) ^ 2;
		}
		return #infinity;
	}
}

species MaxCP parent: Objective {
	action value(SimulatedAnnealingState state) type: float {
		float P_avail <- state.P_avail();
		float C_avail <- state.C_avail();
		if(P_avail > 0 and C_avail > 0) {
			return (1.0 - state.problem.C_P / max((C_avail / P_avail), state.problem.C_P)) ^ 2;	
		}
		return #infinity;
	}
}

species MaxLabileC parent: Objective {
	action value(SimulatedAnnealingState state) type: float {
		float max_C_decomposition;
		ask state.problem {
			// TODO: define this at problem level
			max_C_decomposition <- min(
				max(
					max_enzymes.T_cellulolytic,
					max_enzymes.T_amino, // TODO: limiting factor
					alpha_CP * max_enzymes.T_P
				) * state.problem.dt * state.problem.C_microbes,
				state.problem.total_C_labile
			);
		}
		return (1.0 - (state.enzymatic_activity.X_C_cellulolytic + state.enzymatic_activity.X_C_amino + state.enzymatic_activity.X_C_P)
				/ max_C_decomposition
			) ^ 2;
	}
}

species MaxP parent: Objective {
	action value(SimulatedAnnealingState state) type: float {
		float max_P_decomposition; 
		
		ask state.problem {
		max_P_decomposition <- max(
				min(
					total_C_labile > 0 and total_P_labile > 0 ?
					dt * C_microbes
						* max_enzymes.T_P / (total_C_labile / total_P_labile) : 0.0,
					total_P_labile
				),
				min(
					total_C_recal > 0 and total_P_recal > 0 ?
					alpha_P_recal_to_dim * dt * C_microbes
						* max_enzymes.T_recal / (total_C_recal / total_P_recal) : 0.0,
					total_P_recal
				)
			);
		}
		return max_P_decomposition > 0 ?
				(1.0 - (state.enzymatic_activity.X_P_labile_to_dom + state.enzymatic_activity.X_P_labile_to_dim + state.enzymatic_activity.X_P_recal_to_dim)
					/ max_P_decomposition
				) ^ 2
			: #infinity;
	}
}

species MaxN parent: Objective {
	action value(SimulatedAnnealingState state) type: float {
		float max_N_decomposition;
		ask state.problem {
			max_N_decomposition <- total_C_labile > 0 and total_N_labile > 0 ?
				dt * C_microbes * max_enzymes.T_amino / (total_C_labile / total_N_labile) : 0.0;
		}
		return max_N_decomposition > 0 ?
			(1.0 - state.enzymatic_activity.X_N_amino / max_N_decomposition) ^ 2
			: #infinity;
	}
}

species MaxRecalC parent: Objective {
	action value(SimulatedAnnealingState state) type: float {
		float max_recal_C_decomposition;
		ask state.problem {
			max_recal_C_decomposition <- min(
				max_enzymes.T_recal * dt * C_microbes,
				total_C_recal
				);
		}
		return max_recal_C_decomposition > 0 ? 
			(1.0 - state.enzymatic_activity.X_C_recal / max_recal_C_decomposition) ^ 2
			: (state.problem.max_enzymes.T_recal - state.enzymes.T_recal > 0 ?
			exp(state.enzymes.T_recal / (state.problem.max_enzymes.T_recal - state.enzymes.T_recal))
			: #infinity);
	}
}
	
species MaxTotalC parent: Objective {
	MaxLabileC max_C;
	MaxRecalC max_recal_C;
	
	init {
		create MaxLabileC {
			myself.max_C <- self;
		}
		create MaxRecalC {
			myself.max_recal_C <- self;
		}
	}
	action value(SimulatedAnnealingState state) type: float {
		float result;
		ask state.problem {
			float max_labile_C_enzyme <- max(
				max_enzymes.T_cellulolytic,
				max_enzymes.T_amino, // TODO: limiting factor
				alpha_CP * max_enzymes.T_P
			);
			if(max_enzymes.T_recal > 0 and max_labile_C_enzyme > 0) {
				float max_C_decomposition <- min(
					max_labile_C_enzyme * dt * C_microbes,
					total_C_labile
				);
				float max_C_recal_decomposition <- min(
					max_enzymes.T_recal * dt * C_microbes,
					total_C_recal
				);
				// It is possible to produced recalcitrant enzymes
				if(total_C_recal > 0 and total_C_labile > 0) {
					result <- (1.0 - (
						state.enzymatic_activity.X_C_recal +
						state.enzymatic_activity.X_C_cellulolytic +
						state.enzymatic_activity.X_C_amino +
						state.enzymatic_activity.X_C_P
					) / (max_C_decomposition + max_C_recal_decomposition)) ^ 2;
				} else if (total_C_labile > 0) {
					// Maximize labile C decomposition while prohibiting recalcitrant enzyme production (since total_C_recal = 0).
					result <- (max_enzymes.T_recal - state.enzymes.T_recal) > 0 ?
						(1.0 - (
							state.enzymatic_activity.X_C_cellulolytic +
							state.enzymatic_activity.X_C_amino +
							state.enzymatic_activity.X_C_P
						) / (max_C_decomposition)) ^ 2
							* exp(state.enzymes.T_recal / (max_enzymes.T_recal - state.enzymes.T_recal))
						: #infinity;
				} else if (total_C_recal > 0) {
					// Maximize recalcitrant C decomposition while prohibiting labile enzyme production (since total_C_labile = 0).
					float T_labile <- state.enzymes.T_cellulolytic + state.enzymes.T_amino + alpha_CP * state.enzymes.T_P;
					float T_max_labile <- max_enzymes.T_cellulolytic + max_enzymes.T_amino
						+ alpha_CP * max_enzymes.T_P;
					result <- (T_max_labile - T_labile) > 0 ?
							(1.0 - (state.enzymatic_activity.X_C_cellulolytic + state.enzymatic_activity.X_C_amino + state.enzymatic_activity.X_C_P) / (max_C_decomposition)) ^ 2
								* exp(T_labile / (T_max_labile - T_labile))
						: #infinity;
				} else {
					result <- 0.0;
				}
			} else if (max_labile_C_enzyme > 0) {
				// It is not possible to produce recalcitrant enzyme, so act like if they do not exist
				result <- myself.max_C.value(state);
			} else if(max_enzymes.T_recal > 0) {
				// It is not possible to produce labile enzyme, so act like if they do not exist
				result <- myself.max_recal_C.value(state);
			} else {
				// No C enzyme can be produced.
				result <- 0.0;
			}
		}
		return result;
	}
}


experiment TestSimulatedAnnealing type: gui {
	SimulatedAnnealing simulated_annealing;
	
	init {
		float C_N_microbes <- 10.0;
		float C_N_labile <- 10.0;
		float C_P_microbes <- 20.0;
		float C_P_labile <- 20.0;
		
		TestSimulatedAnnealing current_experiment <- self;
		
		Enzymes min_enzymes;
		create Enzymes with: [
			T_cellulolytic::0.0,
			T_amino::0.0,
			T_P::0.0,
			T_recal::0.0
		] {
			min_enzymes <- self;
		}
		
		Enzymes max_enzymes;
		create Enzymes with: [
			T_cellulolytic::1 #gram / #gram / #d,
			T_amino::1 #gram / #gram / #d,
			T_P::1 #gram / #gram / #d,
			T_recal::0.2 #gram / #gram / #d
		] {
			max_enzymes <- self;
		}
		
		Objective objective;
		create MaxTotalC {
			objective <- self;
		}
		
		create EnzymaticActivityProblem with: [
			C_N::C_N_microbes,
			C_P::C_P_microbes,
			C_microbes::10#gram,
			dt::1#d,
			total_C_labile::200#gram,
			total_N_labile::(200#gram/C_N_labile),
			total_P_labile::(200#gram/C_P_labile),
			total_C_recal::0.0,
			total_N_recal::0.0,
			total_P_recal::0.0,
			C_DOM::0.0,
			P_DOM::0.0,
			N_DOM::0.0,
			P_DIM::0.0,
			N_DIM::0.0,
			min_enzymes::min_enzymes,
			max_enzymes::max_enzymes
		] {
			create SimulatedAnnealing with:[
				problem::self,
				objectives::[objective]
			] {
				current_experiment.simulated_annealing <- self;
			}
		}
		
		ask simulated_annealing {
			do step;
		}
	}
	
	reflex {
		ask simulated_annealing {
			do step;
			float C_avail <- s.C_avail();
			float N_avail <- s.N_avail();
			float P_avail <- s.P_avail();
			write "";
			write "s: " + s;
			write "T_cellulolytic: " + s.enzymes.T_cellulolytic / (#gram / #gram / #h);
			write "T_amino: " + s.enzymes.T_amino / (#gram / #gram / #h);
			write "T_P: " + s.enzymes.T_P / (#gram / #gram / #h);
			write "T_recal: " + s.enzymes.T_recal / (#gram / #gram / #h);
			write "C/N: " + (N_avail > 0 ? C_avail / N_avail : 0.0);
			write "C/P: " + (P_avail > 0 ? C_avail / P_avail : 0.0);
			write "e: " + E(s);

		}
	}
	
	output {
		display Energy type:2d{
			chart "E" type:series {
				data "E" value: sum(SimulatedAnnealing collect each.e) marker:false;
			}
		}
	}
}

experiment ExpEnzymaticActivity type: gui {
	string C_N_objective;
	float C_N_objective_weight;
	string C_P_objective;
	float C_P_objective_weight;
	string C_labile_objective;
	float C_labile_objective_weight;
	string C_recal_objective;
	float C_recal_objective_weight;
	string C_total_objective;
	float C_total_objective_weight;
	string N_objective;
	float N_objective_weight;
	string P_objective;
	float P_objective_weight;
	list<Objective> objectives;

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
	
	int steps;
	
	parameter "C/N" category: "Objectives" var: C_N_objective init: "Max C/N" among: ["none", "Exact C/N", "Max C/N"];
	parameter "C/N weight" category: "Objectives" var: C_N_objective_weight init: 10.0;
	parameter "C/P" category: "Objectives" var: C_P_objective init: "Max C/P" among: ["none", "Exact C/P", "Max C/P"];
	parameter "C/P weight" category: "Objectives" var: C_P_objective_weight init: 10.0;
	parameter "C labile" category: "Objectives" var: C_labile_objective init: "Max labile C" among: ["none", "Max labile C"];
	parameter "C labile weight" category: "Objectives" var: C_labile_objective_weight init: 10.0;
	parameter "C recal" category: "Objectives" var: C_recal_objective init: "Max recal C" among: ["none", "Max recal C"];
	parameter "C recal weight" category: "Objectives" var: C_recal_objective_weight init: 1.0;
	parameter "C total" category: "Objectives" var: C_total_objective init: "none" among: ["none", "Max total C"];
	parameter "C total weight" category: "Objectives" var: C_total_objective_weight init: 1.0;
	parameter "N" category: "Objectives" var: N_objective init: "none" among: ["none", "Max N"];
	parameter "N weight" category: "Objectives" var: N_objective_weight init: 1.0;
	parameter "P" category: "Objectives" var: P_objective init: "none" among: ["none", "Max P"];
	parameter "P weight" category: "Objectives" var: P_objective_weight init: 1.0;
	
	parameter "Max cellulolytic (gC/gM/d)" category: "Constants" var: T_cellulolytic init: 0.1;
	parameter "Max amino (gC/gM/d)" category: "Constants" var: T_amino init: 0.1;
	parameter "Max P (gC/gM/d)" category: "Constants" var: T_P init: 0.1;
	parameter "Max recalcitrant (gC/gM/d)" category: "Constants" var: T_recal init: 0.0;
	parameter "Amino CN" category: "Constants" var: amino_CN init: 6.0;
	parameter "C microbes (g)" category: "Constants" var: C_microbes init: 1.0;
	parameter "C dom (g)" category: "Constants" var: C_dom init: 1.0;
	parameter "Microbe population's C/N" category: "Constants" var:C_N_microbes init:10.0;
	parameter "Microbe population's C/P" category: "Constants" var:C_P_microbes init:20.0;
	
	parameter "Initial C (labile, g)" category: "Labile OM" var: C_labile_init init: 1.0;
	parameter "Final C (labile, g)" category: "Labile OM" var: C_labile_final init: 1.0;
	parameter "Initial C/N (labile)" category: "Labile OM" var: C_N_labile_init init: 20.0 ;
	parameter "Final C/N (labile)" category: "Labile OM" var:C_N_labile_final init:20.0;
	parameter "Initial C/P (labile)" category: "Labile OM" var:C_P_labile_init init:20.0;
	parameter "Final C/P (labile)" category: "Labile OM" var:C_P_labile_final init:20.0;
	
	parameter "Initial C/N (DOM)" category: "DOM" var: C_N_dom_init init: 0.5 ;
	parameter "Final C/N (DOM)" category: "DOM" var:C_N_dom_final init:20.0;
	parameter "Initial C/P (DOM)" category: "DOM" var:C_P_dom_init init:20.0;
	parameter "Final C/P (DOM)" category: "DOM" var:C_P_dom_final init:20.0;
	
	parameter "Initial C (recalcitrant, g)" category: "Recalcitrant OM" var: C_recal_init init: 0.0;
	parameter "Final C (recalcitrant, g)" category: "Recalcitrant OM" var: C_recal_final init: 0.0;
	
	SimulatedAnnealing simulated_annealing;
	parameter "Count of steps" category: "Experiment" var: steps init: 100;
	
	init {
//		write "Labile C threshold: " + C_microbes * T_cellulolytic #gram/#h * 1#h / #gram;
//		write "Recalcitrant C threshold: " + C_microbes * T_recal #gram/#h * 1#h / #gram;
		Objective _C_N_objective;
		if(C_N_objective = "Exact C/N") {
			create ExactCN {
				_C_N_objective <- self;
			}
		} else if (C_N_objective = "Max C/N") {
			create MaxCN {
				_C_N_objective <- self;
			}
		}
		if _C_N_objective != nil {
			ask _C_N_objective {
				weight <- myself.C_N_objective_weight;
				add self to: myself.objectives;
			}	
		}
		
		Objective _C_P_objective;
		if(C_P_objective = "Exact C/P") {
			create ExactCP {
				_C_P_objective <- self;
			}
		} else if (C_P_objective = "Max C/P") {
			create MaxCP {
				_C_P_objective <- self;		
			}
		}
		if _C_P_objective != nil {
			ask _C_P_objective {
				weight <- myself.C_P_objective_weight;
				add self to: myself.objectives;
			}	
		}
		
		if(C_labile_objective = "Max labile C") {
			create MaxLabileC with: (weight: C_labile_objective_weight) {
				add self to: myself.objectives;
			}
		}
		
		if(C_recal_objective = "Max recalcitrant C") {
			create MaxRecalC with: (weight: C_recal_objective_weight) {
				add self to: myself.objectives;
			}
		}
		
		if(C_total_objective = "Max total C") {
			create MaxTotalC with: (weight: C_total_objective_weight) {
				add self to: myself.objectives;
			}
		}
		
		if(N_objective = "Max N") {
			create MaxN with: (weight: N_objective_weight) {
				add self to: myself.objectives;
			}
		}
				
		if(P_objective = "Max P") {
			create MaxP with: (weight: P_objective_weight) {
				add self to: myself.objectives;
			}
		}
		
	}
	
	reflex {
		ask SimulatedAnnealing {
			ask problem {
				ask min_enzymes {
					do die;
				}
				ask max_enzymes {
					do die;
				}
				do die;
			}
			ask problem {
				do die;
			}
			ask s {
				ask enzymes {
					do die;
				}
				ask enzymatic_activity {
					do die;
				}
				do die;
			}
			do die;
		}
		
		Enzymes min_enzymes;
		create Enzymes with: [
			T_cellulolytic::0.0,
			T_amino::0.0,
			T_P::0.0,
			T_recal::0.0
		] {
			min_enzymes <- self;
		}
		Enzymes max_enzymes;
		create Enzymes with: [
			T_cellulolytic::T_cellulolytic #gram/#gram/#d,
			T_amino::T_amino #gram/#gram/#d,
			T_P::T_P #gram/#gram/#d,
			T_recal::T_recal #gram/#gram/#d
		] {
			max_enzymes <- self;
		}
		ExpEnzymaticActivity exp <- self;
		
		float total_C_labile <- (exp.C_labile_init + cycle * (exp.C_labile_final - exp.C_labile_init)/exp.steps)#gram;
		float total_N_labile <- total_C_labile/(exp.C_N_labile_init + cycle * (exp.C_N_labile_final - exp.C_N_labile_init)/exp.steps);
		float total_P_labile <- total_C_labile/(exp.C_P_labile_init + cycle * (exp.C_P_labile_final - exp.C_P_labile_init)/exp.steps);
		float total_C_recal <- (exp.C_recal_init + cycle * (exp.C_recal_final - exp.C_recal_init)/exp.steps)#gram;
		create EnzymaticActivityProblem with: [
			C_N::exp.C_N_microbes,
			C_P::exp.C_P_microbes,
			C_microbes::exp.C_microbes#gram,
			dt::1#d,
			total_C_labile::total_C_labile,
			total_N_labile::total_N_labile,
			total_P_labile::total_P_labile,
			total_C_recal::total_C_recal,
			total_N_recal::0.0,
			total_P_recal::0.0,
			C_DOM::exp.C_dom#gram,
			N_DOM::(exp.C_dom/(exp.C_N_dom_init + cycle * (exp.C_N_dom_final - exp.C_N_dom_init)/exp.steps)),
			P_DOM::(exp.C_dom/(exp.C_P_dom_init + cycle * (exp.C_P_dom_final - exp.C_P_dom_init)/exp.steps)),
			P_DIM::0.0,
			N_DIM::0.0,
			min_enzymes::min_enzymes,
			max_enzymes::max_enzymes
		] {
			
			create SimulatedAnnealing with:[
				problem::self,
				objectives::exp.objectives,
				N::1000,
				epsilon::1e-3
			] {
				do optimize;
				
				float C_avail <- s.C_avail();
				float N_avail <- s.N_avail();
				float P_avail <- s.P_avail();
				write "";
				write "s: " + s;
				write "T_cellylolytic: " + s.enzymes.T_cellulolytic / (#gram / #gram / #d);
				write "T_amino: " + s.enzymes.T_amino / (#gram / #gram / #d);
				write "T_P: " + s.enzymes.T_P / (#gram / #gram / #d);
				write "C/N: " + (N_avail > 0 ? C_avail / N_avail : 0.0);
				write "C/P: " + (P_avail > 0 ? C_avail / P_avail : 0.0);
				write "e: " + E(s);
							
				ask s.enzymatic_activity {
					write "X_C_cellulolytic: " + X_C_cellulolytic;
					write "X_C_amino: " + X_C_amino;
					write "X_N_amino: " + X_N_amino;
					write "X_C_P: " + X_C_P;
					write "X_P_labile_to_dom:" + X_P_labile_to_dom;
					write "X_P_labile_to_dim: " + X_P_labile_to_dim;
					write "X_C_recal: " + X_C_recal;
					write "X_P_recal_to_dim:" + X_P_recal_to_dim;
				}
				exp.simulated_annealing <- self;
			}
		}
	}
	
	reflex when: cycle=steps {
		ask simulation {
			do pause;
		}
	}
	
	output {
		display "OM" type:2d {
			chart "OM" type:series style:line {
				data "C (labile)" value: sum(SimulatedAnnealing collect each.problem.total_C_labile)/#gram marker:false;
				data "N (labile)" value: sum(SimulatedAnnealing collect each.problem.total_N_labile)/#gram marker:false;
				data "P (labile)" value: sum(SimulatedAnnealing collect each.problem.total_P_labile)/#gram marker:false;
				data "C (recal)" value: sum(SimulatedAnnealing collect each.problem.total_C_recal)/#gram marker:false;
				data "N (recal)" value: sum(SimulatedAnnealing collect each.problem.total_N_recal)/#gram marker:false;
				data "P (recal)" value: sum(SimulatedAnnealing collect each.problem.total_P_recal)/#gram marker:false;
				data "C (avail)" value: sum(SimulatedAnnealing collect each.s.C_avail())/#gram marker:false;
				data "N (avail)" value: sum(SimulatedAnnealing collect each.s.N_avail())/#gram marker:false;
				data "P (avail)" value: sum(SimulatedAnnealing collect each.s.P_avail())/#gram marker:false;
			}
		}
		display "Enzymatic rates" type:2d {
			chart "Enzymatic rates" type:series style:line {
				data "T_cellulolytic" value: sum(SimulatedAnnealing collect (each.s.enzymes.T_cellulolytic / (#gram / #gram / #d))) marker:false;
				data "T_amino" value: sum(SimulatedAnnealing collect (each.s.enzymes.T_amino / (#gram / #gram / #d))) marker:false;
				data "T_P" value: sum(SimulatedAnnealing collect (each.s.enzymes.T_P / (#gram / #gram / #d))) marker:false;
				data "T_recal" value: sum(SimulatedAnnealing collect (each.s.enzymes.T_recal / (#gram / #gram / #d))) marker:false;
				data "T_cellulolytic (max)" value: sum(SimulatedAnnealing collect (each.problem.max_enzymes.T_cellulolytic / (#gram / #gram / #d))) marker:false;
				data "T_amino (max)" value: sum(SimulatedAnnealing collect (each.problem.max_enzymes.T_amino / (#gram / #gram / #d))) marker:false;
				data "T_P (max)" value: sum(SimulatedAnnealing collect (each.problem.max_enzymes.T_P / (#gram / #gram / #d))) marker:false;
				data "T_recal (max)" value: sum(SimulatedAnnealing collect (each.problem.max_enzymes.T_recal / (#gram / #gram / #d))) marker:false;
			}
		}
		display "C/N" type:2d {
			chart "C/N" type:series style:line {
				data "organic C/N" value: sum(SimulatedAnnealing collect (
					each.problem.total_N_labile > 0 ? each.problem.total_C_labile/each.problem.total_N_labile : 0.0
				)) marker:false;
				data "dom C/N" value: sum(SimulatedAnnealing collect (
					each.problem.N_DOM > 0 ? each.problem.C_DOM/each.problem.N_DOM : 0.0
				)) marker:false;
				data "available C/N" value: sum(SimulatedAnnealing collect (
					each.s.N_avail() > 0 ?
					each.s.C_avail()/each.s.N_avail() : 0.0
				)) marker:false;
				data "C/N" value: sum(SimulatedAnnealing collect each.problem.C_N) marker:false;
			}
		}
		display "C/P" type:2d {
			chart "C/P" type:series style:line {
				data "organic C/P" value: sum(SimulatedAnnealing collect (
					each.problem.total_P_labile > 0 ? each.problem.total_C_labile/each.problem.total_P_labile : 0.0
				)) marker:false;
				data "dom C/P" value: sum(SimulatedAnnealing collect (
					each.problem.P_DOM > 0 ? each.problem.C_DOM/each.problem.P_DOM : 0.0
				)) marker:false;
				data "available C/P" value: sum(SimulatedAnnealing collect (
					each.s.P_avail() > 0 ?
					each.s.C_avail() / each.s.P_avail() : 0.0
				)) marker:false;
				data "C/P" value: sum(SimulatedAnnealing collect each.problem.C_P) marker:false;
			}
		}
	}
}
