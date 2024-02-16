/**
* Name: enzymes
* Based on the internal skeleton template. 
* Author: pbreugno
* Tags: 
*/

model enzyme_activity

global {
	int T_CELLULOLYTIC_BUDGET <- 0;
	int T_AMINO_BUDGET <- 1;
	int T_P_BUDGET <- 2;
	int T_RECALCITRANT_BUDGET <- 3;
}

/**
 * Enzymatic action rates.
 */
species Enzymes schedules: [] {
	/**
	 * Cellulolytic action.
	 */
	float T_cellulolytic <- 0.0#gram/#gram/#d;
	/**
	 * Amino acid production action.
	 */
	float T_amino <- 0.0#gram/#gram/#d;
	/**
	 * P mineralisation action.
	 */
	float T_P <- 0.0#gram/#gram/#d;
	/**
	 * Recalcitrant cleavage action.
	 */
	float T_recal <- 0.0#gram/#gram/#d;
}

/**
 * Enzymatic action rates, multiplied by the active microbe biomass.
 * Notice the difference of units.
 */
species WeightedEnzymes schedules: [] {
	float T_cellulolytic <- 0.0#gram/#d;
	float T_amino <- 0.0#gram/#d;
	float T_P <- 0.0#gram/#d;
	float T_recal <- 0.0#gram/#d;
}

/**
 * Exchanges of matter between nutrient compartments.
 */
species Decomposition schedules: [] {
	/**
	 * Quantity of matter removed from C labile and added to C dom due to
	 * cellulolytic action.
	 */
	float X_C_cellulolytic <- 0.0;
	/**
	 * Quantity of matter removed from C labile and added to C dom due to amino
	 * action.
	 */
	float X_C_amino <- 0.0;
	/**
	 * Quantity of matter removed from C labile and added to N dom due to amino
	 * action.
	 */
	float X_N_amino <- 0.0;
	/**
	 * Quantity of matter removed from C labile and added to C dom due to P
	 * mineralisation action.
	 */
	float X_C_P <- 0.0;
	/**
	 * Quantity of matter removed from P labile and added to P dom due to P
	 * mineralisation action.
	 */
	float X_P_labile_to_dom <- 0.0;
	/**
	 * Quantity of matter removed from P labile and added to P dim due to P
	 * mineralisation action.
	 */
	float X_P_labile_to_dim <- 0.0;
	
	/**
	 * Quantity of matter removed from C recalcitrant and added to C labile.
	 */
	float X_C_recal <- 0.0;
	/**
	 * Quantity of matter removed from N recalcitrant and added to N labile.
	 */
	float X_N_recal <- 0.0;
	/**
	 * Quantity of matter removed from P recalcitrant and added to P labile.
	 */
	float X_P_recal_to_labile <- 0.0;
	/**
	 * Quantity of matter removed from P recalcitrant and added to P dim.
	 */
	float X_P_recal_to_dim <- 0.0;
}

/**
 * The decomposition problem embeds all parameters required to compute the
 * decomposition of matter due given a set of enzymatic activities.
 *
 * It is notably used to estimate the quality of an enzymatic budget allocation
 * in the EnzymaticActivityProblem, but also standalone to compute the
 * decomposition of organic particles at each time step, without optimisation.
 */
species DecompositionProblem schedules: [] {
	/**
	 * Time step duration, over which the enzymatic activity is optimized.
	 * Also used to compute the decomposition.
	 */
	float dt;
	
	/**
	 * Initial C labile that can be decomposed.
	 */
	float C_labile;
	/**
	 * Initial N labile that can be decomposed.
	 */
	float N_labile;
	/**
	 * Initial P labile that can be decomposed.
	 */
	float P_labile;
	
	/**
	 * Initial C recalcitrant that can be decomposed.
	 */
	float C_recal;
	/**
	 * Initial N recalcitrant that can be decomposed.
	 */
	float N_recal;
	/**
	 * Initial P recalcitrant that can be decomposed.
	 */
	float P_recal;
	
	/**
	 * Initial C available in the DOM.
	 */
	float C_DOM;
	/**
	 * Initial N available in the DOM.
	 */
	float N_DOM;
	/**
	 * Initial P available in the DOM.
	 */
	float P_DOM;
	/**
	 * Initial N available in the DIM.
	 */
	float N_DIM;
	/**
	 * Initial P available in the DIM.
	 */
	float P_DIM;
	
	/**
	 * Fixed CN rates of amino acids.
	 * Represent the distribution of amino acids.
	 */
	float amino_CN <- 4.0;
	
	/**
	 * Rate of C decomposed by the cellulolytic action if only the cellulolytic
	 * and amino actions request all the labile C.
	 */
	float beta_cellulolytic_amino <- 0.5;
	/**
	 * Rate of C decomposed by the cellulolytic action if only the cellulolytic
	 * and P mineralisation actions request all the labile C.
	 */
	float beta_cellulolytic_P <- 0.5;
	/**
	 * Rate of C decomposed by the amino action if only the amino and P
	 * mineralisation actions request all the labile C.
	 */
	float beta_amino_P <- 0.5;
	
	/**
	 * Rate of organic C sent from labile to DOM among the C decomposed by the P
	 * mineralisation.
	 */
	float alpha_C_e_P <- 0.1;
	/**
	 * Rate of organic P sent from labile to DIM among the P decomposed by the P
	 * mineralisation.
	 */
	float alpha_P_e_P <- 0.9;
	/**
	 * Rate of organic P sent from recalcitrant to DIM among the P decomposed by
	 * recalcitrant cleavage. Represents the the action of phytases.
	 */
	float alpha_P_e_r <- 1e-3;
	
	/**
	 * Quantity of carbon requested by the recalcitrant cleavage action.
	 */
	float d_C_recal(WeightedEnzymes enzymes) {
		return min(
					dt * enzymes.T_recal,
					C_recal
				);
	}
	
	/**
	 * Quantity of carbon requested by the cellulolytic action.
	 */
	float d_C_cellulolytic(WeightedEnzymes enzymes) {
		return min(
					dt * enzymes.T_cellulolytic,
					C_labile
				);
	}
	
	/**
	 * Quantity of carbon requested by the amino acid production action.
	 */
	float d_C_amino(WeightedEnzymes enzymes) {
		float expected_C_amino <- dt * enzymes.T_amino;
		float expected_N_amino <- expected_C_amino / amino_CN;
		float limiting_amino <- expected_C_amino > 0 ?
			min(1.0, C_labile / expected_C_amino, N_labile / expected_N_amino)
			: 0.0;
		return expected_C_amino * limiting_amino;
	}
	
	/**
	 * Quantity of carbon requested by the P mineralisation action.
	 */
	float d_C_P(WeightedEnzymes enzymes) {
		float P_attacked <- min(
					dt * enzymes.T_P,
					P_labile
				);
		return P_attacked * C_labile / P_labile;
	}
	
	/**
	 * Computes the available C that results from the provided decomposition.
	 */
	float C_avail(Decomposition decomposition) {
		return C_DOM + decomposition.X_C_cellulolytic + decomposition.X_C_amino + decomposition.X_C_P;
	}
	
	/**
	 * Computes the available N that results from the provided decomposition.
	 */
	float N_avail(Decomposition decomposition) {
		return N_DOM + N_DIM + decomposition.X_N_amino;
	}
	
	/**
	 * Computes the available P that results from the provided decomposition.
	 */
	float P_avail(Decomposition decomposition) {
		return P_DOM + P_DIM + decomposition.X_P_labile_to_dom + decomposition.X_P_labile_to_dim + decomposition.X_P_recal_to_dim;
	}
	
	/**
	 * Computes the decomposition obtained from the specified enzymes.
	 * 
	 * The result is stored in the provided Decomposition instance.
	 */
	action decomposition(
		WeightedEnzymes enzymes,
		Decomposition result
	) {
		if(C_labile > 0) {
			// Requested labile C by each action
			float d_C_cellulolytic <- d_C_cellulolytic(enzymes);
			float d_C_P <- d_C_P(enzymes);
			float d_C_amino <- d_C_amino(enzymes);
			
			result.X_C_cellulolytic <- d_C_cellulolytic
				- beta_cellulolytic_amino * d_C_amino * d_C_cellulolytic / C_labile
				- beta_cellulolytic_P * d_C_P * d_C_cellulolytic / C_labile;
				
			result.X_C_amino <- d_C_amino
				- (1-beta_cellulolytic_amino) * d_C_amino * d_C_cellulolytic / C_labile
				- beta_amino_P * d_C_amino * d_C_P / C_labile;

			float D_C_P <- d_C_P
				- (1-beta_cellulolytic_P) * d_C_P * d_C_cellulolytic / C_labile
				- (1-beta_amino_P) * d_C_amino * d_C_P / C_labile;
				
			result.X_N_amino <- result.X_C_amino / amino_CN;
			if(P_labile > 0) {
				float D_P <- D_C_P * (P_labile / C_labile);
				
				result.X_P_labile_to_dim <- alpha_P_e_P * D_P;
						
				// Only a fraction goes to the dom, the rest is left in the labile section
				result.X_C_P <- alpha_C_e_P * D_C_P;
				result.X_P_labile_to_dom <- alpha_C_e_P * (D_P-result.X_P_labile_to_dim);
			} else {
				result.X_P_labile_to_dim <- 0.0;
				result.X_C_P <- 0.0;
				result.X_P_labile_to_dom <- 0.0;	
			}
		} else {
			result.X_C_cellulolytic <- 0.0;
			result.X_C_amino <- 0.0;
			result.X_N_amino <- 0.0;
			result.X_P_labile_to_dim <- 0.0;
			result.X_C_P <- 0.0;
			result.X_P_labile_to_dom <- 0.0;
		}

		if(C_recal > 0) {
			result.X_C_recal <- d_C_recal(enzymes);
			if(N_recal > 0) {
				result.X_N_recal <- result.X_C_recal * (N_recal / C_recal);
			} else {
				result.X_N_recal <- 0.0;
			}
			if(P_recal > 0) {
				float D_P_recal <- result.X_C_recal * (P_recal / C_recal);
				result.X_P_recal_to_dim <-
					alpha_P_e_r * D_P_recal;
				result.X_P_recal_to_labile <-
					(1-alpha_P_e_r) * D_P_recal;
			} else {
				result.X_P_recal_to_dim <- 0.0;
				result.X_P_recal_to_labile <- 0.0;
			}
		} else {
			result.X_C_recal <- 0.0;
			result.X_N_recal <- 0.0;
			result.X_P_recal_to_dim <- 0.0;
			result.X_P_recal_to_labile <- 0.0;
		}
	}
}

/**
 * Represents the state of the problem optimised by the simulated annealing
 * algorithm.
 */
species SimulatedAnnealingState schedules: [] {
	/**
	 * Current problem to solve.
	 */
	EnzymaticActivityProblem problem;
	/**
	 * Current enzymatic budget allocation, used to represent
	 */
	list<float> budget;
	
	/**
	 * Current enzymatic activities, computed from the budget.
	 */
	Enzymes enzymes;
	/**
	 * Current enzymatic activities, weighted by the active microbe biomass,
	 * computed from the budget.
	 */
	WeightedEnzymes weighted_enzymes;
	/**
	 * Estimated decomposition, computed from the current problem and enzymes.
	 */
	Decomposition decomposition;

	init {
		SimulatedAnnealingState current_state <- self;
		ask problem {
			create Enzymes with: [
				T_cellulolytic: min_enzymes.T_cellulolytic + myself.budget[T_CELLULOLYTIC_BUDGET]*(max_enzymes.T_cellulolytic - min_enzymes.T_cellulolytic),
				T_amino: min_enzymes.T_amino + myself.budget[T_AMINO_BUDGET]*(max_enzymes.T_amino - min_enzymes.T_amino),
				T_P: min_enzymes.T_P + myself.budget[T_P_BUDGET]*(max_enzymes.T_P - min_enzymes.T_P),
				T_recal: min_enzymes.T_recal + myself.budget[T_RECALCITRANT_BUDGET]*(max_enzymes.T_recal - min_enzymes.T_recal)
			] {
				current_state.enzymes <- self;
			}
		}
		
		create WeightedEnzymes with: [
			T_cellulolytic::problem.C_microbes*enzymes.T_cellulolytic,
			T_amino::problem.C_microbes*enzymes.T_amino,
			T_P::problem.C_microbes*enzymes.T_P,
			T_recal::problem.C_microbes*enzymes.T_recal
		] {
			current_state.weighted_enzymes <- self;
		}
		
		create Decomposition {
			current_state.decomposition <- self;
		}
		
		ask problem.decomposition_problem {
			do decomposition(current_state.weighted_enzymes, current_state.decomposition);
		}
	}

	/**
	 * Computes the estimation of available C obtained in the current state.
	 */
	float C_avail {
		return problem.decomposition_problem.C_avail(decomposition);
	}
	/**
	 * Computes the estimation of available N obtained in the current state.
	 */
	float N_avail {
		return problem.decomposition_problem.N_avail(decomposition);
	}
	/**
	 * Computes the estimation of available P obtained in the current state.
	 */
	float P_avail {
		return problem.decomposition_problem.P_avail(decomposition);
	}
}

/**
 * Definition of a generic objective, that will be optimized by the simulated
 * annealing algorithm.
 */
species Objective schedules: [] {
	/**
	 * Relative weight of the objective used in multi-objective optimization.
	 */
	float weight <- 1.0;
	
	/**
	 * Returns the value of the objective for the specified state.
	 */
	action value(SimulatedAnnealingState state) virtual: true type: float;
}

/**
 * The enzymatic activity problem embeds all parameters required to optimize the
 * enzymatic budget allocation.
 * 
 * The environment state on which enzymatic activity is optimized is given by
 * the embedded decomposition problem, used to estimate decomposition for each
 * state in the optimization process.
 */
species EnzymaticActivityProblem schedules: [] {
	DecompositionProblem decomposition_problem;
	
	/**
	 * Microbe population active structural C.
	 */
	float C_microbes;
	/**
	 * Requested microbe CN rate.
	 */
	float C_N;
	/**
	 * Requested microbe CP rate.
	 */
	float C_P;
	/**
	 * Constitutive enzymatic rates.
	 */
	Enzymes min_enzymes;
	/**
	 * Maximum enzymatic rates.
	 */
	Enzymes max_enzymes;
	
	// Useful values used in the computation of objectives
	/**
	 * Estimation of the maximum quantity of C that can be added to C labile due
	 * to the recalcitrant cleavage action.
	 */
	float X_C_labile_max;
	/**
	 * Estimation of the maximum quantity of C that can be added to the dam due
	 * to the cellulolytic, amino acid production and P mineralisation actions.
	 */
	float X_C_dam_max;
	/**
	 * Estimation of the maximum quantity of N that can be added to the dam due
	 * to the amino acid production action. 
	 */
	float X_N_dam_max;
	/**
	 * Estimation of the maximum quantity of P that can be added to the dam due
	 * to the recalcitrant cleavage and P mineralisation actions.
	 */
	float X_P_dam_max;
	
	action update_X_C_labile_max {
		create WeightedEnzymes with: [
			T_recal::C_microbes * max_enzymes.T_recal
		] {
			myself.X_C_labile_max <- myself.decomposition_problem.d_C_recal(self);
			do die;
		}
	}
	
	action update_X_C_dam_max {
		WeightedEnzymes enzymes;
		create WeightedEnzymes with: [
			T_cellulolytic::C_microbes*max_enzymes.T_cellulolytic,
			T_amino::C_microbes*max_enzymes.T_amino,
			T_P::C_microbes*max_enzymes.T_P,
			T_recal::0.0
		] {
			// The current weighted enzymes is impossible in terms of budget, but it is a way
			// to easily compute each d_* max since:
			// - the maximum of d_C_cellulolytic is only computed from the max_enzymes.T_cellulolytic (assuming other T_* = 0)
			// - the maximum of d_C_amino is only computed from the max_enzymes.T_amino (assuming other T_* = 0)
			// - ... idem form d_C_P
			enzymes <- self;
		}
		ask decomposition_problem {
			myself.X_C_dam_max <- max(
				d_C_cellulolytic(enzymes),
				d_C_amino(enzymes),
				alpha_C_e_P * d_C_P(enzymes)
			);
		}
		ask enzymes {
			do die;
		}
	}
	
	action update_X_N_dam_max {
		WeightedEnzymes enzymes;
		create WeightedEnzymes with: [
			T_amino::C_microbes*max_enzymes.T_amino
		] {
			enzymes <- self;
		}
		ask decomposition_problem {
			myself.X_N_dam_max <- d_C_amino(enzymes) / amino_CN;
		}
		ask enzymes {
			do die;
		}
	}
		
	action update_X_P_dam_max {
		WeightedEnzymes enzymes;
		create WeightedEnzymes with: [
			T_P::C_microbes*max_enzymes.T_P,
			T_recal::C_microbes*max_enzymes.T_recal
		] {
			enzymes <- self;
		}
		ask decomposition_problem {
			myself.X_P_dam_max <- max(
				alpha_P_e_r * d_C_recal(enzymes),
				(alpha_C_e_P + alpha_P_e_P - alpha_C_e_P * alpha_P_e_P)
					* d_C_P(enzymes)
			);
		}
		ask enzymes {
			do die;
		}
	}
}

/**
 * Generic simulated annealing algorithm.
 */
species SimulatedAnnealing schedules: [] {
	/**
	 * Problem to solve.
	 */
	EnzymaticActivityProblem problem;
	
	/**
	 * Maximum count of steps.
	 */
	int N <- 1000;
	/**
	 * If the value of the objective is below this threshold, the current state is returned, even if N has not been reached yet.
	 */
	float epsilon <- 0.0;
	
	/**
	 * Initial temperature. Must be in the order of the values of the objective function.
	 */
	float T_init;
	/**
	 * Current temperature.
	 */
	float T;
	/**
	 * Current energy to minimize.
	 */
	float e;
	/**
	 * Current state.
	 */
	SimulatedAnnealingState s;
	/**
	 * List of objectives to optimize.
	 */
	list<Objective> objectives;
	
	/**
	 * Initial enzymatic budget allocation.
	 */
	list<float> init_budget;
	int budget_size;
	
	// TODO: define indexes used in neighbor here to allow T_amino = 0 and others
	
	init {
		T_init <- float(length(objectives));
		T <- T_init;
		if(problem.max_enzymes.T_recal > 0.0) {
			init_budget <- [1.0/4, 1.0/4, 1.0/4, 1.0/4];
			budget_size <- 4;
		} else {
			init_budget <- [1.0/3, 1.0/3, 1.0/3, 0.0];
			budget_size <- 3; // So the last budget stays to 0
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
				ask weighted_enzymes {
					do die;
				}
				ask decomposition {
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
				ask weighted_enzymes {
					do die;
				}
				ask decomposition {
					do die;
				}
				do die;
			}
		}
		T <- 0.99 * T;
	}
	
	action optimize {
		ask problem {
			do update_X_C_labile_max;
			do update_X_C_dam_max;
			do update_X_N_dam_max;
			do update_X_P_dam_max;
		}
		
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
		loop i from: 0 to: budget_size-1 {
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
		return sum(objectives collect (each.weight * each.value(state)^2));
	}
}

species ExactCN parent: Objective schedules: [] {
	action value(SimulatedAnnealingState state) type: float {
		float N_avail <- state.N_avail();
		float C_avail <- state.C_avail();
		if(C_avail > 0) {
			return 1.0 - state.problem.C_N * (N_avail / C_avail);
		}
		if(state.problem.decomposition_problem.C_labile > 0.0) {
			// Max out the value, since C can be decomposed but no C is decomposed.
			return 1.0;
		}
		// No C/N/P can be decomposed anyway
		return 0.0;
	}
}

species MaxCN parent: Objective schedules: [] {
	action value(SimulatedAnnealingState state) type: float {
		float N_avail <- state.N_avail();
		float C_avail <- state.C_avail();
		if(C_avail > 0) {
			// If N_avail/C_avail > required N/C, then N_avail is in excess and its not a problem so value=0 in this case.
			return 1.0 - state.problem.C_N * min(N_avail / C_avail, 1.0/state.problem.C_N);
		}
		if(state.problem.decomposition_problem.C_labile > 0.0) {
			// Max out the value, since C can be decomposed but no C is decomposed.
			return 1.0;
		}
		// No C/N/P can be decomposed anyway
		return 0.0;
	}
}

species ExactCP parent: Objective schedules: [] {
	action value(SimulatedAnnealingState state) type: float {
		float P_avail <- state.P_avail();
		float C_avail <- state.C_avail();
		if(C_avail > 0) {
			return 1.0 - state.problem.C_P * (P_avail / C_avail);
		}
		if(state.problem.decomposition_problem.C_labile > 0.0) {
			// Max out the value, since C can be decomposed but no C is decomposed.
			return 1.0;
		}
		// No C/N/P can be decomposed anyway
		return 0.0;
	}
}

species MaxCP parent: Objective schedules: [] {
	action value(SimulatedAnnealingState state) type: float {
		float P_avail <- state.P_avail();
		float C_avail <- state.C_avail();
		if(C_avail > 0) {
			// If P_avail/C_avail > required P/C, then P_avail is in excess and its not a problem so value=0 in this case.
			return 1.0 - state.problem.C_P * min((P_avail / C_avail), 1.0/state.problem.C_P);	
		}
		if(state.problem.decomposition_problem.C_labile > 0.0) {
			// Max out the value, since C can be decomposed but no C is decomposed.
			return 1.0;
		}
		// No C/N/P can be decomposed anyway
		return 0.0;
	}
}

species MaxRecalC parent: Objective schedules: [] {
	action value(SimulatedAnnealingState state) type: float {
		// Assumes max_enzymes.T_recal > 0.0. Otherwise, this objective should not be used.
		float result;
		ask state {
			if(problem.X_C_labile_max > 0.0) {
				result <- 1.0 - decomposition.X_C_recal / problem.X_C_labile_max;
			} else {
				// Forces T_recal to tend to 0 if max_recal_C_decomposition = 0.0
				result <- exp(ln(2) * budget[T_RECALCITRANT_BUDGET]) - 1.0;	
			}
		}
		return result;
	}
}

species MaxLabileC parent: Objective schedules: [] {
	action value(SimulatedAnnealingState state) type: float {
		float result;
		ask state {
			if(problem.X_C_dam_max > 0.0) {
				result <- 1.0 - (decomposition.X_C_cellulolytic + decomposition.X_C_amino + decomposition.X_C_P)
						/ problem.X_C_dam_max;
			} else {
				// Forces T_cellulolytic, T_amino and T_P budgets to 0 since no labile C/N/P can be decomposed anyway.
				result <- exp(ln(2) * (state.budget[T_CELLULOLYTIC_BUDGET] + budget[T_AMINO_BUDGET] + budget[T_P_BUDGET])/3) - 1.0;
			}
		}
		return result;
	}
}

species MaxP parent: Objective schedules: [] {
	action value(SimulatedAnnealingState state) type: float {
		float result;
		ask state {
			if(problem.X_P_dam_max > 0.0) {
				result <- 1.0 - (decomposition.X_P_recal_to_dim + decomposition.X_P_labile_to_dom + decomposition.X_P_labile_to_dim)
						/ problem.X_P_dam_max;
			} else {
				// No P can be decomposed anyway
				result <- exp(ln(2) * (budget[T_P_BUDGET] + budget[T_RECALCITRANT_BUDGET])/2) - 1.0;
			}
		}
		return result;
	}
}

species MaxN parent: Objective schedules: [] {
	action value(SimulatedAnnealingState state) type: float {
		float result;
		ask state {
			if(problem.X_N_dam_max > 0) {
				result <- 1.0 - decomposition.X_N_amino / problem.X_N_dam_max;
			} else {
				// max_N_decomposition = 0, so we force T_amino budget to 0
				result <- exp(ln(2) * state.budget[T_AMINO_BUDGET]) - 1.0;		
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
		create MaxLabileC {
			objective <- self;
		}
		create DecompositionProblem with: [
			dt::1#d,
			C_labile::200#gram,
			N_labile::(200#gram/C_N_labile),
			P_labile::(200#gram/C_P_labile),
			C_recal::0.0,
			N_recal::0.0,
			P_recal::0.0,
			C_DOM::0.0,
			P_DOM::0.0,
			N_DOM::0.0,
			P_DIM::0.0,
			N_DIM::0.0
		] {
			create EnzymaticActivityProblem with: [
				decomposition_problem::self,
				C_N::C_N_microbes,
				C_P::C_P_microbes,
				C_microbes::1#gram,
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
	string N_objective;
	float N_objective_weight;
	string P_objective;
	float P_objective_weight;
	list<Objective> objectives;

	float T_cellulolytic;
	float T_amino;
	float T_P;
	float T_recal;
	float amino_CN;
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
		float total_N_labile;
		if(exp.C_N_labile_final < #infinity and exp.C_N_labile_init < #infinity) {
			total_N_labile <- total_C_labile/(exp.C_N_labile_init + cycle * (exp.C_N_labile_final - exp.C_N_labile_init)/exp.steps);
		} else {
			total_N_labile <- 0.0;	
		}
		float total_P_labile;
		if(exp.C_P_labile_final < #infinity and exp.C_P_labile_init < #infinity) {
			total_P_labile <- total_C_labile/(exp.C_P_labile_init + cycle * (exp.C_P_labile_final - exp.C_P_labile_init)/exp.steps);
		} else {
			total_P_labile <- 0.0;
		}
		float total_C_recal <- (exp.C_recal_init + cycle * (exp.C_recal_final - exp.C_recal_init)/exp.steps)#gram;
		
		create DecompositionProblem with: [
			dt::1#d,
			amino_CN::amino_CN,
			C_labile::total_C_labile,
			N_labile::total_N_labile,
			P_labile::total_P_labile,
			C_recal::total_C_recal,
			N_recal::0.0,
			P_recal::0.0,
			C_DOM::exp.C_dom#gram,
			N_DOM::(exp.C_dom#gram/(exp.C_N_dom_init + cycle * (exp.C_N_dom_final - exp.C_N_dom_init)/exp.steps)),
			P_DOM::(exp.C_dom#gram/(exp.C_P_dom_init + cycle * (exp.C_P_dom_final - exp.C_P_dom_init)/exp.steps)),
			P_DIM::0.0,
			N_DIM::0.0
		] {
			create EnzymaticActivityProblem with: [
				decomposition_problem::self,
				C_N::exp.C_N_microbes,
				C_P::exp.C_P_microbes,
				C_microbes::exp.C_microbes#gram,
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
					write "C_DOM: " + myself.decomposition_problem.C_DOM;
					write "N_DOM: " + myself.decomposition_problem.N_DOM;
					write "e: " + E(s);
								
					ask s.decomposition {
						write "X_C_cellulolytic: " + X_C_cellulolytic / #gram;
						write "X_C_amino: " + X_C_amino / #gram;
						write "X_N_amino: " + X_N_amino / #gram;
						write "X_C_P: " + X_C_P / #gram;
						write "X_P_labile_to_dom:" + X_P_labile_to_dom / #gram;
						write "X_P_labile_to_dim: " + X_P_labile_to_dim / #gram;
						write "X_C_recal: " + X_C_recal / #gram;
						write "X_P_recal_to_dim:" + X_P_recal_to_dim / #gram;
					}
					exp.simulated_annealing <- self;
				}
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
//				data "C (labile)" value: sum(SimulatedAnnealing collect each.problem.total_C_labile)/#gram marker:false;
//				data "N (labile)" value: sum(SimulatedAnnealing collect each.problem.total_N_labile)/#gram marker:false;
//				data "P (labile)" value: sum(SimulatedAnnealing collect each.problem.total_P_labile)/#gram marker:false;
//				data "C (recal)" value: sum(SimulatedAnnealing collect each.problem.total_C_recal)/#gram marker:false;
//				data "N (recal)" value: sum(SimulatedAnnealing collect each.problem.total_N_recal)/#gram marker:false;
//				data "P (recal)" value: sum(SimulatedAnnealing collect each.problem.total_P_recal)/#gram marker:false;
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
					each.problem.decomposition_problem.N_labile > 0 ?
						each.problem.decomposition_problem.C_labile/each.problem.decomposition_problem.N_labile : 0.0
				)) marker:false;
				data "dom C/N" value: sum(SimulatedAnnealing collect (
					each.problem.decomposition_problem.N_DOM > 0 ?
					each.problem.decomposition_problem.C_DOM/each.problem.decomposition_problem.N_DOM : 0.0
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
					each.problem.decomposition_problem.P_labile > 0 ?
					each.problem.decomposition_problem.C_labile/each.problem.decomposition_problem.P_labile : 0.0
				)) marker:false;
				data "dom C/P" value: sum(SimulatedAnnealing collect (
					each.problem.decomposition_problem.P_DOM > 0 ?
					each.problem.decomposition_problem.C_DOM/each.problem.decomposition_problem.P_DOM : 0.0
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
