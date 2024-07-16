/**
* Name: enzymes
* Author: Paul Breugnot
* 
* Implementation of the buzy-pop model.
*/

model enzyme_activity

global {
	/**
	 * Index of the budget allocated to CNP extraction.
	 */
	int T_CNP_BUDGET <- 0;
	/**
	 * Index of the budget allocated to N extraction.
	 */
	int T_N_BUDGET <- 1;
	/**
	 * Index of the budget allocated to P extraction.
	 */
	int T_P_BUDGET <- 2;
	/**
	 * Index of the budget allocated to recalcitrant cleavage.
	 */
	int T_RECALCITRANT_BUDGET <- 3;
}

/**
 * Enzymatic action rates.
 */
species Enzymes schedules: [] {
	/**
	 * CNP extraction.
	 * 
	 * Decomposition of sugars, carbohydrates, proteins and fatty acids.reduces
	 * the size of carbon based polymers (depolymerisation) using endocleaving
	 * and exocleaving enzymes. Produces C/N/P compounds with C:N and C:P ratios
	 * corresponding to the C:N and C:P ratios of labile OM.
	 */
	float T_CNP <- 0.0#gram/#gram/#d;

	/**
	 * N extraction.
	 * 
	 * Decomposition of polypeptides and amino sugars. Extract assimilable amino
	 * acids from polypeptides or N-acetyl glucosamine - which is the building
	 * block of chitin molecules - from amino sugars. Represents the grouped
	 * actions of aminopeptidases and N-acetyl glucosamidase. The C:N ratio
	 * of products is equal to extracted_CN.
	 */
	float T_N <- 0.0#gram/#gram/#d;

	/**
	 * P extraction.
	 * 
	 * Mineralisation of P. Extracts P into the DIM. Represents the
	 * decomposition of phosphate groups by the combined action of
	 * phosphomonoesterases and phosphodiesterases.
	 */
	float T_P <- 0.0#gram/#gram/#d;

	/**
	 * Recalcitrant cleavage.
	 * 
	 * Breaks coarse recalcitrant C based molecules into smaller and more labile
	 * compounds, but not yet assimilable. This action brings together the
	 * activities of numerous enzymes such as laccases, peroxidases, cellulases,
	 * tannases, lipases, endoproteases, endochitinases, and endonucleases.
	 * Doing so, trapped N and P is released in the labile OM in addition to C
	 * according to the C:N and C:P ratios of the recalcitrant OM.
	 * 
	 * In order to model the special action of phytases, a small fraction
	 * alpha_P_er of the released P is directly sent to the DIM, modelling the
	 * first P extracted from recalcitrant phytic acid molecules.
	 */
	float T_recal <- 0.0#gram/#gram/#d;
}

/**
 * Enzymatic action rates, multiplied by the active microbe biomass.
 * Notice the difference of units compared to Enzymes.
 */
species WeightedEnzymes schedules: [] {
	float T_CNP <- 0.0#gram/#d;
	float T_N <- 0.0#gram/#d;
	float T_P <- 0.0#gram/#d;
	float T_recal <- 0.0#gram/#d;
}

/**
 * Exchanges of matter between nutrient compartments.
 */
species Decomposition schedules: [] {
	/**
	 * Quantity of C removed from C labile and added to C dom due to
	 * CNP extraction.
	 */
	float X_C_eCNP <- 0.0;

	/**
	 * Quantity of N removed from N labile and added to N dom due to
	 * CNP extraction.
	 */
	float X_N_eCNP <- 0.0;

	/**
	 * Quantity of P removed from P labile and added to P dom due to
	 * CNP extraction.
	 */
	float X_P_eCNP <- 0.0;
	
	/**
	 * Quantity of C removed from C labile and added to C dom due to N
	 * extraction.
	 */
	float X_C_eN <- 0.0;
	
	/**
	 * Quantity of N removed from N labile and added to N dom due to N
	 * extraction.
	 */
	float X_N_eN <- 0.0;
	
	/**
	 * Quantity of P removed from P labile and added to P dom due to P
	 * extraction.
	 */
	float X_P_eP <- 0.0;
	
	/**
	 * Quantity of C removed from C recalcitrant and added to C labile due to
	 * recalcitrant cleavage.
	 */
	float X_C_er <- 0.0;

	/**
	 * Quantity of N removed from N recalcitrant and added to N labile due to
	 * recalcitrant cleavage.
	 */
	float X_N_er <- 0.0;

	/**
	 * Quantity of P removed from P recalcitrant and added to P labile due to
	 * recalcitrant cleavage.
	 */
	float X_P_er_recal_to_labile <- 0.0;

	/**
	 * Quantity of P removed from P recalcitrant and added to P dim due to
	 * recalcitrant cleavage. (specific action of phytases)
	 */
	float X_P_er_recal_to_dim <- 0.0;
}

/**
 * The decomposition problem embeds all parameters required to compute the
 * decomposition of matter given a set of enzymatic activities.
 *
 * It is notably used to estimate the quality of an enzymatic budget allocation
 * in the EnzymaticActivityProblem, but also standalone to compute the
 * decomposition of organic particles at each time step.
 */
species DecompositionProblem schedules: [] {
	/**
	 * Time step duration, over which the enzymatic activity is optimized.
	 * Also used to compute the decomposition.
	 */
	float dt;
	
	/**
	 * Initial C recalcitrant that can be decomposed.
	 */
	float C_recal_init;
	/**
	 * Initial N recalcitrant that can be decomposed.
	 */
	float N_recal_init;
	/**
	 * Initial P recalcitrant that can be decomposed.
	 */
	float P_recal_init;
	
	/**
	 * Initial C labile that can be decomposed.
	 */
	float C_labile_init;
	/**
	 * Initial N labile that can be decomposed.
	 */
	float N_labile_init;
	/**
	 * Initial P labile that can be decomposed.
	 */
	float P_labile_init;
	
	/**
	 * Initial C available in the DOM.
	 */
	float C_DOM_init;
	/**
	 * Initial N available in the DOM.
	 */
	float N_DOM_init;
	/**
	 * Initial P available in the DOM.
	 */
	float P_DOM_init;
	/**
	 * Initial N available in the DIM.
	 */
	float N_DIM_init;
	/**
	 * Initial P available in the DIM.
	 */
	float P_DIM_init;
	
	/**
	 * Fixed C:N rate of matter produced by the N extraction action.
	 * 
	 * Represents the distribution of amino acids and N-acetyl glucosamine in
	 * the substrate.
	 */
	float extracted_CN <- 5.0;

	/**
	 * Count of C in the polymer section bound to each P and attacked by P
	 * extraction. This quantity of C is requested, but not included in products
	 * as only the phosphate is extracted to the DIM.
	 */
	float phosphatase_CP <- 1.0;

	/**
	 * Rate of organic P sent from recalcitrant to DIM among the P decomposed by
	 * recalcitrant cleavage. Represents the the action of phytases.
	 */
	float alpha_P_er <- 1e-3;

	/**
	 * Rate of C decomposed by C extraction if only C and N extractions request
	 * all the labile C.
	 */
	float beta_CNP_N <- 0.5;
	/**
	 * Rate of C decomposed by C extraction if only C and P extractions request
	 * all the labile C.
	 */
	float beta_CNP_P <- 0.5;
	/**
	 * Rate of C decomposed by N extraction if only N and P extractions request
	 * all the labile C.
	 */
	float beta_N_P <- 0.5;

	/**
	 * Quantity of substrate that can be requested by the recalcitrant cleavage
	 * action.
	 * 
	 * It is expressed as a quantity of recalcitrant C.
	 */
	float recal_substrate {
		return C_recal_init;
	}
	
	/**
	 * Quantity of substrate that can be requested by C extraction.
	 * 
	 * It is expressed as a quantity of labile C.
	 */
	float CNP_extraction_substrate {
		return C_labile_init;
	}
	
	/**
	 * Quantity of substrate that can be requested by N extraction.
	 * 
	 * It is expressed as a quantity of labile C included in polypeptides and
	 * chitin (amino acids and N-acetyl glucosamine based molecules).
	 */
	float N_extraction_substrate {
		return min(C_labile_init, N_labile_init * extracted_CN);
	}
	
	/**
	 * Quantity of substrate that can be requested by P extraction.
	 * 
	 * It is expressed as a quantity of labile C bound to phosphatase.
	 */
	float P_extraction_substrate {
		return min(C_labile_init, P_labile_init * phosphatase_CP);
	}
	
	/**
	 * Quantity of carbon requested by the recalcitrant cleavage action.
	 */
	float d_C_recal(WeightedEnzymes enzymes) {
		return min(
					dt * enzymes.T_recal,
					recal_substrate()
				);
	}
	
	/**
	 * Quantity of carbon requested by CNP extraction.
	 */
	float d_CNP_extraction(WeightedEnzymes enzymes) {
		return min(
					dt * enzymes.T_CNP,
					CNP_extraction_substrate()
				);
	}
	
	/**
	 * Quantity of carbon requested by N extraction.
	 */
	float d_N_extraction(WeightedEnzymes enzymes) {
		return min(
			dt * enzymes.T_N,
			N_extraction_substrate()
		);
	}
	
	/**
	 * Quantity of carbon requested by P extraction.
	 */
	float d_P_extraction(WeightedEnzymes enzymes) {
		return min(
					dt * enzymes.T_P,
					P_extraction_substrate()
				);
	}
	
	/**
	 * Computes the quantity of recalcitrant C resulting from the provided
	 * decomposition.
	 */
	float C_recal_final(Decomposition decomposition) {
		return C_recal_init - decomposition.X_C_er;
	}
	
	/**
	 * Computes the quantity of recalcitrant N resulting from the provided
	 * decomposition.
	 */
	float N_recal_final(Decomposition decomposition) {
		return N_recal_init - decomposition.X_N_er;
	}
	
	/**
	 * Computes the quantity of recalcitrant P resulting from the provided
	 * decomposition.
	 */
	float P_recal_final(Decomposition decomposition) {
		return P_recal_init - decomposition.X_P_er_recal_to_labile - decomposition.X_P_er_recal_to_dim;
	}
	
	/**
	 * Represents the quantity of labile C resulting from the provided
	 * decomposition, without products of the decomposition of recalcitrant OM.
	 */
	float C_labile_decomp(Decomposition decomposition) {
		return C_labile_init - decomposition.X_C_eCNP - decomposition.X_C_eN;
	}
	
	/**
	 * Represents the quantity of labile N resulting from the provided
	 * decomposition, without products of the decomposition of recalcitrant OM.
	 */
	float N_labile_decomp(Decomposition decomposition) {
		return N_labile_init - decomposition.X_N_eCNP - decomposition.X_N_eN;
	}
	
	/**
	 * Represents the quantity of labile P resulting from the provided
	 * decomposition, without products of the decomposition of recalcitrant OM.
	 */
	float P_labile_decomp(Decomposition decomposition) {
		return P_labile_init - decomposition.X_P_eCNP - decomposition.X_P_eP;
	}
	
	/**
	 * Computes the quantity of labile C resulting from the provided
	 * decomposition, including products of the decomposition of recalcitrant
	 * OM.
	 */
	float C_labile_final(Decomposition decomposition) {
		return C_labile_decomp(decomposition) + decomposition.X_C_er;
	}
	
	/**
	 * Computes the quantity of labile N resulting from the provided
	 * decomposition, including products of the decomposition of recalcitrant
	 * OM.
	 */
	float N_labile_final(Decomposition decomposition) {
		return N_labile_decomp(decomposition) + decomposition.X_N_er;
	}
	
	/**
	 * Computes the quantity of labile P resulting from the provided
	 * decomposition, including products of the decomposition of recalcitrant
	 * OM.
	 */
	float P_labile_final(Decomposition decomposition) {
		return P_labile_decomp(decomposition) + decomposition.X_P_er_recal_to_labile;
	}
	
	/**
	 * Computes the quantity of C in the DOM resulting from the provided
	 * decomposition.
	 */
	float C_DOM_final(Decomposition decomposition) {
		return C_DOM_init + decomposition.X_C_eCNP + decomposition.X_C_eN;
	}
	
	/**
	 * Computes the quantity of N in the DOM resulting from the provided
	 * decomposition.
	 */
	float N_DOM_final(Decomposition decomposition) {
		return N_DOM_init + decomposition.X_N_eN + decomposition.X_N_eCNP;
	}
	
	/**
	 * Computes the quantity of P in the DOM resulting from the provided
	 * decomposition.
	 */
	float P_DOM_final(Decomposition decomposition) {
		return P_DOM_init + decomposition.X_P_eCNP;
	}
	
	/**
	 * Computes the quantity of N in the DIM resulting from the provided
	 * decomposition.
	 */
	float N_DIM_final(Decomposition decomposition) {
		return N_DIM_init;
	}
	
	/**
	 * Computes the quantity of P in the DIM resulting from the provided
	 * decomposition.
	 */
	float P_DIM_final(Decomposition decomposition) {
		return P_DIM_init + decomposition.X_P_eP + decomposition.X_P_er_recal_to_dim;
	}
	
	/**
	 * Computes the quantity of available C resulting in the DAM from the
	 * provided decomposition.
	 */
	float C_avail_final(Decomposition decomposition) {
		return C_DOM_final(decomposition);
	}
	
	/**
	 * Computes the quantity available N resulting in the DAM from the provided
	 * decomposition.
	 */
	float N_avail_final(Decomposition decomposition) {
		return N_DOM_final(decomposition) + N_DIM_final(decomposition);
	}
		
	/**
	 * Computes the quantity of available P resulting in the DAM from the
	 * provided decomposition.
	 */
	float P_avail_final(Decomposition decomposition) {
		return P_DOM_final(decomposition) + P_DIM_final(decomposition);
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
		if(C_labile_init > 0) {
			// Requested labile C by each action
			float d_CNP <- d_CNP_extraction(enzymes);
			float d_P <- d_P_extraction(enzymes);
			float d_N <- d_N_extraction(enzymes);

			result.X_C_eCNP <- d_CNP
				- beta_CNP_N * d_N * d_CNP / C_labile_init
				- beta_CNP_P * d_P * d_CNP / C_labile_init;

			result.X_C_eN <- d_N
				- (1-beta_CNP_N) * d_N * d_CNP / C_labile_init
				- beta_N_P * d_N * d_P / C_labile_init;

			float D_C_P <- d_P
				- (1-beta_CNP_P) * d_P * d_CNP / C_labile_init
				- (1-beta_N_P) * d_N * d_P / C_labile_init;
	
			result.X_N_eN <- result.X_C_eN / extracted_CN;
			
			result.X_P_eP <- D_C_P / phosphatase_CP;
			
			float left_C <- C_labile_init - result.X_C_eN - D_C_P;
			if left_C > 0 {
				result.X_N_eCNP <- result.X_C_eCNP * (N_labile_init-result.X_N_eN) / left_C;
				result.X_P_eCNP <- result.X_C_eCNP * (P_labile_init-result.X_P_eP) / left_C;	
			}
		} else {
			// C/N/P extraction
			result.X_C_eCNP <- 0.0;
			result.X_N_eCNP <- 0.0;
			result.X_P_eCNP <- 0.0;
			// N extraction
			result.X_C_eN <- 0.0;
			result.X_N_eN <- 0.0;
			// P extraction
			result.X_P_eP <- 0.0;
		}

		if(C_recal_init > 0) {
			result.X_C_er <- d_C_recal(enzymes);
			result.X_N_er <- result.X_C_er * (N_recal_init / C_recal_init);
			
			float D_P_recal <- result.X_C_er * (P_recal_init / C_recal_init);
			result.X_P_er_recal_to_dim <-
				alpha_P_er * D_P_recal;
			result.X_P_er_recal_to_labile <-
				(1-alpha_P_er) * D_P_recal;
		} else {
			result.X_C_er <- 0.0;
			result.X_N_er <- 0.0;
			result.X_P_er_recal_to_dim <- 0.0;
			result.X_P_er_recal_to_labile <- 0.0;
		}
	}
}

/**
 * Represents the state of the problem optimized by the simulated annealing
 * algorithm.
 */
species SimulatedAnnealingState schedules: [] {
	/**
	 * Current problem to solve.
	 */
	EnzymaticActivityProblem problem;
	/**
	 * Current enzymatic budget allocation, used to represent the part of the
	 * budget allocated to each action.
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
				T_CNP: min_enzymes.T_CNP + myself.budget[T_CNP_BUDGET]*(max_enzymes.T_CNP - min_enzymes.T_CNP),
				T_N: min_enzymes.T_N + myself.budget[T_N_BUDGET]*(max_enzymes.T_N - min_enzymes.T_N),
				T_P: min_enzymes.T_P + myself.budget[T_P_BUDGET]*(max_enzymes.T_P - min_enzymes.T_P),
				T_recal: min_enzymes.T_recal + myself.budget[T_RECALCITRANT_BUDGET]*(max_enzymes.T_recal - min_enzymes.T_recal)
			] {
				current_state.enzymes <- self;
			}
		}
		
		create WeightedEnzymes with: [
			T_CNP::problem.C_microbes*enzymes.T_CNP,
			T_N::problem.C_microbes*enzymes.T_N,
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
	/**
	 * Decomposition problem used to estimate the quality of each budget.
	 */
	DecompositionProblem decomposition_problem;
	
	/**
	 * Active structural C of the microbe population.
	 */
	float C_microbes;
	/**
	 * CN ratio requested by the microbe population.
	 */
	float C_N;
	/**
	 * CP ratio requested by the microbe population.
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
	 * Estimation of the maximum quantity of labile C that can be produced by
	 * the recalcitrant cleavage action.
	 * 
	 * The value is computed by the update_X_C_labile_max() method.
	 */
	float X_C_labile_max;
	
	/**
	 * Estimation of the maximum quantity of available C that can be produced by
	 * CNP and N extraction actions.
	 * 
	 * The value is computed by the update_X_C_dam_max() method.
	 */
	float X_C_dam_max;
	
	/**
	 * Estimation of the maximum quantity of available N that can be produced by
	 * CNP and N extraction.
	 * 
	 * The value is computed by the update_X_N_dam_max() method.
	 */
	float X_N_dam_max;
	
	/**
	 * Estimation of the maximum quantity of available P that can be produced by
	 * the recalcitrant cleavage, CNP and P extraction actions.
	 * 
	 * The value is computed by the update_X_P_dam_max() method.
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
			T_CNP::C_microbes*max_enzymes.T_CNP,
			T_N::C_microbes*max_enzymes.T_N
		] {
			// The current weighted enzymes is impossible in terms of budget, but it is a way
			// to easily compute each d_* max since:
			// - the maximum of d_CNP_extraction is only computed from
			// max_enzymes.T_CNP (assuming other T_* = 0)
			// - the maximum of d_N_extraction is only computed from
			// max_enzymes.T_N (assuming other T_* = 0)
			// - ... idem for d_C_P
			enzymes <- self;
		}
		ask decomposition_problem {
			myself.X_C_dam_max <- max(
				d_CNP_extraction(enzymes),
				d_N_extraction(enzymes)
			);
		}
		ask enzymes {
			do die;
		}
	}
	
	action update_X_N_dam_max {
		WeightedEnzymes enzymes;
		create WeightedEnzymes with: [
			T_CNP::C_microbes*max_enzymes.T_CNP,
			T_N::C_microbes*max_enzymes.T_N
		] {
			enzymes <- self;
		}
		ask decomposition_problem {
			myself.X_N_dam_max <- C_labile_init > 0 ? max(
				d_CNP_extraction(enzymes) * N_labile_init / C_labile_init,
				d_N_extraction(enzymes) / extracted_CN
				) : 0.0;
		}
		ask enzymes {
			do die;
		}
	}
		
	action update_X_P_dam_max {
		WeightedEnzymes enzymes;
		create WeightedEnzymes with: [
			T_CNP::C_microbes*max_enzymes.T_CNP,
			T_P::C_microbes*max_enzymes.T_P,
			T_recal::C_microbes*max_enzymes.T_recal
		] {
			enzymes <- self;
		}
		ask decomposition_problem {
			myself.X_P_dam_max <- max(
				C_recal_init > 0.0 ? alpha_P_er * d_C_recal(enzymes) * P_recal_init / C_recal_init : 0.0,
				C_labile_init > 0.0 ? d_P_extraction(enzymes) / phosphatase_CP : 0.0,
				C_labile_init > 0.0 ? d_CNP_extraction(enzymes) * P_labile_init / C_labile_init : 0.0
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
	 * Maximum count of steps.
	 */
	int N <- 1000;
	/**
	 * If the value of the objective is below this threshold, the current state
	 * is returned, even if N has not been reached yet.
	 */
	float epsilon <- 0.0;
	
	/**
	 * Initial temperature. Must be in the order of the values of the objective
	 * function.
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

	/**
	 * Indexes of actions used in the optimization (when the enzymatic activity
	 * of an action is fixed, it is ignored in the optimization process).
	 */
	list<int> budget_indexes;
	
	action init_state(EnzymaticActivityProblem problem) {
		T_init <- float(length(objectives));
		T <- T_init;
		budget_indexes <- [];
		if(problem.max_enzymes.T_CNP - problem.min_enzymes.T_CNP > 0.0) {
			add T_CNP_BUDGET to: budget_indexes;
		}
		if(problem.max_enzymes.T_N - problem.min_enzymes.T_N > 0.0) {
			add T_N_BUDGET to: budget_indexes;
		}
		if(problem.max_enzymes.T_P - problem.min_enzymes.T_P > 0.0) {
			add T_P_BUDGET to: budget_indexes;
		}
		if(problem.max_enzymes.T_recal - problem.min_enzymes.T_recal > 0.0) {
			add T_RECALCITRANT_BUDGET to: budget_indexes;
		}
		init_budget <- [0.0, 0.0, 0.0, 0.0];
		loop i over: budget_indexes {
			// The budget is initially distributed equally between all possible
			// actions
			init_budget[i] <- 1.0 / length(budget_indexes);
		}

		create SimulatedAnnealingState with: [problem::problem, budget::init_budget] {
			myself.s <- self;
		}
		e <- E(s);
	}
	
	action step {
		SimulatedAnnealingState s_new <- neighbour();
		float e_new <- E(s_new);
		if e_new < e or rnd(1.0) < exp(-(e_new - e) / T) {
			ask s {
				// Clears current state
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
			// Update current state
			s <- s_new;
			e <- e_new;
		} else {
			ask s_new {
				// Discards new state
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
	
	action optimize(EnzymaticActivityProblem problem) {
		do init_state(problem);
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
	
	/**
	 * Generates a random neighbour of the current state.
	 */
	SimulatedAnnealingState neighbour {
		SimulatedAnnealingState result;
		
		list<float> _T;
		loop item over: s.budget {
			add item to: _T;
		}

		// The indexes list is shuffled to determine in which order each
		// component will vary. Since the total budget (1.0) is fixed, this
		// allows a fair access to the budget for each action (i.e. the last
		// component does not always get the last part left after all other get
		// their own part).
		list<int> indexes <- shuffle(budget_indexes);
		float delta <- 0.1;
		float range <- 1.0;
		loop i from: 0 to: length(indexes)-2 {
			_T[indexes[i]] <- max(0, min(range, gauss(_T[indexes[i]], delta)));
			range <- range - _T[indexes[i]];
		}
		_T[indexes[length(indexes)-1]] <- range;
		
		create SimulatedAnnealingState with:[problem::s.problem, budget::_T] {
			result <- self;
		}
		return result;
	}
	
	/**
	 * Energy of the system to minimize.
	 */
	float E(
		SimulatedAnnealingState state
	) {
		// Multi-objective function to minimize
		return sum(objectives collect (each.weight * each.value(state)^2));
	}
}

/**
 * Objective that is minimized when the available CN rate is exactly equal to
 * the requested CN rate.
 */
species ExactCN parent: Objective schedules: [] {
	action value(SimulatedAnnealingState state) type: float {
		float N_avail <- state.problem.decomposition_problem.N_avail_final(state.decomposition);
		float C_avail <- state.problem.decomposition_problem.C_avail_final(state.decomposition);
		if(C_avail > 0) {
			return 1.0 - state.problem.C_N * (N_avail / C_avail);
		}
		if(state.problem.decomposition_problem.C_labile_init > 0.0) {
			// Max out the value, since C can be decomposed but no C is decomposed.
			return 1.0;
		}
		// No C/N/P can be decomposed anyway
		return 0.0;
	}
}

/**
 * Objective that is minimized when the available CN rate is at most equal to
 * the requested CN rate.
 */
species CapCN parent: Objective schedules: [] {
	action value(SimulatedAnnealingState state) type: float {
		float N_avail <- state.problem.decomposition_problem.N_avail_final(state.decomposition);
		float C_avail <- state.problem.decomposition_problem.C_avail_final(state.decomposition);
		if(C_avail > 0) {
			// If N_avail/C_avail > required N/C, then N_avail is in excess and its not a problem so value=0 in this case.
			return 1.0 - state.problem.C_N * min(N_avail / C_avail, 1.0/state.problem.C_N);
		}
		if(state.problem.decomposition_problem.C_labile_init > 0.0) {
			// Max out the value, since C can be decomposed but no C is decomposed.
			return 1.0;
		}
		// No C/N/P can be decomposed anyway
		return 0.0;
	}
}

/**
 * Objective that is minimized when the available CP rate is exactly equal to
 * the requested CP rate.
 */
species ExactCP parent: Objective schedules: [] {
	action value(SimulatedAnnealingState state) type: float {
		float P_avail <- state.problem.decomposition_problem.P_avail_final(state.decomposition);
		float C_avail <- state.problem.decomposition_problem.C_avail_final(state.decomposition);
		if(C_avail > 0) {
			return 1.0 - state.problem.C_P * (P_avail / C_avail);
		}
		if(state.problem.decomposition_problem.C_labile_init > 0.0) {
			// Max out the value, since C can be decomposed but no C is decomposed.
			return 1.0;
		}
		// No C/N/P can be decomposed anyway
		return 0.0;
	}
}

/**
 * Objective that is minimized when the available CP rate is at most equal to
 * the requested CP rate.
 */
species CapCP parent: Objective schedules: [] {
	action value(SimulatedAnnealingState state) type: float {
		float P_avail <- state.problem.decomposition_problem.P_avail_final(state.decomposition);
		float C_avail <- state.problem.decomposition_problem.C_avail_final(state.decomposition);
		if(C_avail > 0) {
			// If P_avail/C_avail > required P/C, then P_avail is in excess and its not a problem so value=0 in this case.
			return 1.0 - state.problem.C_P * min((P_avail / C_avail), 1.0/state.problem.C_P);	
		}
		if(state.problem.decomposition_problem.C_labile_init > 0.0) {
			// Max out the value, since C can be decomposed but no C is decomposed.
			return 1.0;
		}
		// No C/N/P can be decomposed anyway
		return 0.0;
	}
}

/**
 * Objective minimized when the quantity of produced labile C is maximized.
 * 
 * Labile C is notably produced from recalcitrant cleavage action.
 */
species MaxLabileC parent: Objective schedules: [] {
	action value(SimulatedAnnealingState state) type: float {
		// Assumes max_enzymes.T_recal > 0.0. Otherwise, this objective should not be used.
		float result;
		ask state {
			if(problem.X_C_labile_max > 0.0) {
				result <- 1.0 - decomposition.X_C_er / problem.X_C_labile_max;
			} else {
				// Forces T_recal to tend to 0 if max_recal_C_decomposition = 0.0
				result <- exp(ln(2) * budget[T_RECALCITRANT_BUDGET]) - 1.0;	
			}
		}
		return result;
	}
}

/**
 * Objective minimized when the quantity of produced available C is maximized.
 * 
 * Available C is produced from the C, N and P extraction actions.
 */
species MaxC parent: Objective schedules: [] {
	action value(SimulatedAnnealingState state) type: float {
		float result;
		ask state {
			if(problem.X_C_dam_max > 0.0) {
				result <- 1.0 - (decomposition.X_C_eN + decomposition.X_C_eCNP)
						/ problem.X_C_dam_max;
			} else {
				// Forces T_C, T_N and T_P budgets to 0 since no labile C/N/P can be decomposed anyway.
				result <- exp(ln(2) * (state.budget[T_CNP_BUDGET] + budget[T_N_BUDGET] + budget[T_P_BUDGET])/3) - 1.0;
			}
		}
		return result;
	}
}

/**
 * Objective minimized when the quantity of produced available N is maximized.
 * 
 * Available N is produced from N extraction.
 */
species MaxN parent: Objective schedules: [] {
	action value(SimulatedAnnealingState state) type: float {
		float result;
		ask state {
			if(problem.X_N_dam_max > 0) {
				result <- 1.0 - (decomposition.X_N_eN + decomposition.X_N_eCNP) / problem.X_N_dam_max;
			} else {
				// max_N_decomposition = 0, so we force T_N budget to 0
				result <- exp(ln(2) * state.budget[T_N_BUDGET]) - 1.0;		
			}
		}
		return result;
	}
}

/**
 * Objective minimized when the quantity of produced available P is maximized.
 * 
 * Available P is produced from the P extraction.
 */
species MaxP parent: Objective schedules: [] {
	action value(SimulatedAnnealingState state) type: float {
		float result;
		ask state {
			if(problem.X_P_dam_max > 0.0) {
				result <- 1.0 - (decomposition.X_P_er_recal_to_dim + decomposition.X_P_eP + decomposition.X_P_eCNP)
						/ problem.X_P_dam_max;
			} else {
				// No P can be decomposed anyway
				result <- exp(ln(2) * (budget[T_P_BUDGET] + budget[T_RECALCITRANT_BUDGET])/2) - 1.0;
			}
		}
		return result;
	}
}

/**
 * An experiment with many user defined parameters that allows to test the
 * behaviour of the enzymatic activity model in specific configurations.
 */
experiment EnzymaticActivityWorkbench type: gui {
	string C_N_objective;
	float C_N_objective_weight;
	string C_P_objective;
	float C_P_objective_weight;
	string C_labile_objective;
	float C_labile_objective_weight;
	string C_objective;
	float C_objective_weight;
	string N_objective;
	float N_objective_weight;
	string P_objective;
	float P_objective_weight;
	list<Objective> objectives;

	float min_T_CNP;
	float min_T_N;
	float min_T_P;
	float min_T_recal;
	float max_T_CNP;
	float max_T_N;
	float max_T_P;
	float max_T_recal;
	
	float C_microbes;
	float C_N_microbes;
	float C_P_microbes;
	float CUE;
	
	float dt;
	float extracted_CN;
	float phosphatase_CP;
	float alpha_P_er;
	float beta_CNP_N;
	float beta_CNP_P;
	float beta_N_P;
	
	float C_labile_init;
	float C_labile_final;
	float C_N_labile_init;
	float C_N_labile_final;
	float C_P_labile_init;
	float C_P_labile_final;
	
	float C_dom_init;
	float C_dom_final;
	float C_N_dom_init;
	float C_N_dom_final;
	float C_P_dom_init;
	float C_P_dom_final;
	
	float C_recal_init;
	float C_recal_final;
	float C_N_recal_init;
	float C_N_recal_final;
	float C_P_recal_init;
	float C_P_recal_final;
	
	int N;
	float epsilon;
	
	int steps;
	bool show_max_rates;
	
	parameter "C/N" category: "Objectives" var: C_N_objective init: "Cap C/N" among: ["none", "Exact C/N", "Cap C/N"];
	parameter "C/N weight" category: "Objectives" var: C_N_objective_weight init: 10.0;
	parameter "C/P" category: "Objectives" var: C_P_objective init: "Cap C/P" among: ["none", "Exact C/P", "Cap C/P"];
	parameter "C/P weight" category: "Objectives" var: C_P_objective_weight init: 10.0;
	parameter "C labile" category: "Objectives" var: C_labile_objective init: "Max labile C" among: ["none", "Max labile C"];
	parameter "C labile weight" category: "Objectives" var: C_labile_objective_weight init: 1.0;
	parameter "C" category: "Objectives" var: C_objective init: "Max C" among: ["none", "Max C"];
	parameter "C weight" category: "Objectives" var: C_objective_weight init: 5.0;
	parameter "N" category: "Objectives" var: N_objective init: "none" among: ["none", "Max N"];
	parameter "N weight" category: "Objectives" var: N_objective_weight init: 5.0;
	parameter "P" category: "Objectives" var: P_objective init: "none" among: ["none", "Max P"];
	parameter "P weight" category: "Objectives" var: P_objective_weight init: 5.0;
	
	parameter "Min CNP extraction (gC/gM/d)" category: "Microbe population" var: min_T_CNP init: 0.0;
	parameter "Min N extraction (gC/gM/d)" category: "Microbe population" var: min_T_N init: 0.0;
	parameter "Min P extraction (gC/gM/d)" category: "Microbe population" var: min_T_P init: 0.0;
	parameter "Min recalcitrant (gC/gM/d)" category: "Microbe population" var: min_T_recal init: 0.0;
	parameter "Max CNP extraction (gC/gM/d)" category: "Microbe population" var: max_T_CNP init: 3.0;
	parameter "Max N extraction (gC/gM/d)" category: "Microbe population" var: max_T_N init: 1.0;
	parameter "Max P extraction (gC/gM/d)" category: "Microbe population" var: max_T_P init: 0.2;
	parameter "Max recalcitrant (gC/gM/d)" category: "Microbe population" var: max_T_recal init: 1.0;

	parameter "C microbes (g)" category: "Microbe population" var: C_microbes init: 1.0;
	parameter "Microbe population's C/N" category: "Microbe population" var:C_N_microbes init:10.0;
	parameter "Microbe population's C/P" category: "Microbe population" var:C_P_microbes init:17.0;
	parameter "Microbe population's CUE" category: "Microbe population" var:CUE init:0.3;
	
	parameter "dt" category: "Constants" var:dt init: 1#d;
	parameter "Extracted CN" category: "Constants" var: extracted_CN init: 5.0;
	parameter "Phosphatase CP" category: "Constants" var: phosphatase_CP init: 1.0;
	parameter "Alpha P er" category: "Constants" var: alpha_P_er init: 0.001;
	parameter "Beta CNP N" category: "Constants" var: beta_CNP_N init: 0.5;
	parameter "Beta CNP P" category: "Constants" var: beta_CNP_P init: 0.5;
	parameter "Beta N P" category: "Constants" var: beta_N_P init: 0.5;
		
	parameter "Initial C dom (g)" category: "DOM" var: C_dom_init init: 0.0;
	parameter "Final C dom (g)" category: "DOM" var: C_dom_final init: 0.0;
	parameter "Initial C/N dom" category: "DOM" var: C_N_dom_init init: 10.0;
	parameter "Final C/N dom" category: "DOM" var:C_N_dom_final init:10.0;
	parameter "Initial C/P dom" category: "DOM" var:C_P_dom_init init:17.0;
	parameter "Final C/P dom" category: "DOM" var:C_P_dom_final init:17.0;
	
	parameter "Initial C (labile, g)" category: "Labile OM" var: C_labile_init init: 1.0;
	parameter "Final C (labile, g)" category: "Labile OM" var: C_labile_final init: 1.0;
	parameter "Initial C/N (labile)" category: "Labile OM" var: C_N_labile_init init: 20.0;
	parameter "Final C/N (labile)" category: "Labile OM" var:C_N_labile_final init:20.0;
	parameter "Initial C/P (labile)" category: "Labile OM" var:C_P_labile_init init:20.0;
	parameter "Final C/P (labile)" category: "Labile OM" var:C_P_labile_final init:20.0;
	
	parameter "Initial C (recalcitrant, g)" category: "Recalcitrant OM" var: C_recal_init init: 1.0;
	parameter "Final C (recalcitrant, g)" category: "Recalcitrant OM" var: C_recal_final init: 1.0;
	parameter "Initial C/N (recalcitrant)" category: "Recalcitrant OM" var: C_N_recal_init init: 17.0;
	parameter "Final C/N (recalcitrant)" category: "Recalcitrant OM" var:C_N_recal_final init:17.0;
	parameter "Initial C/P (recalcitrant)" category: "Recalcitrant OM" var:C_P_recal_init init:31.0;
	parameter "Final C/P (recalcitrant)" category: "Recalcitrant OM" var:C_P_recal_final init:31.0;
	
	parameter "Maximum N" category: "Simulated annealing" var:N init: 1000;
	parameter "Epsilon" category: "Simulated annealing" var:epsilon init: 1e-3;
	
	parameter "Count of steps" category: "Experiment" var: steps init: 100;
	parameter "Show max enzymatic rates" category: "Experiment" var:show_max_rates init: false;
	
	SimulatedAnnealing simulated_annealing;
	DecompositionProblem decomposition_problem;
	EnzymaticActivityProblem problem;
	
	init {
		Objective _C_N_objective;
		if(C_N_objective = "Exact C/N") {
			create ExactCN {
				_C_N_objective <- self;
			}
		} else if (C_N_objective = "Cap C/N") {
			create CapCN {
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
		} else if (C_P_objective = "Cap C/P") {
			create CapCP {
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
		
		if(C_objective = "Max C") {
			create MaxC with: (weight: C_objective_weight) {
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
			ask s.problem {
				ask decomposition_problem {
					do die;
				}
				ask min_enzymes {
					do die;
				}
				ask max_enzymes {
					do die;
				}
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
			T_CNP::min_T_CNP #gram/#gram/#d,
			T_N::min_T_N #gram/#gram/#d,
			T_P::min_T_P #gram/#gram/#d,
			T_recal::min_T_recal #gram/#gram/#d
		] {
			min_enzymes <- self;
		}
		Enzymes max_enzymes;
		create Enzymes with: [
			T_CNP::max_T_CNP #gram/#gram/#d,
			T_N::max_T_N #gram/#gram/#d,
			T_P::max_T_P #gram/#gram/#d,
			T_recal::max_T_recal #gram/#gram/#d
		] {
			max_enzymes <- self;
		}
		EnzymaticActivityWorkbench exp <- self;
		
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
		float total_N_recal;
		if(exp.C_N_recal_final < #infinity and exp.C_N_recal_init < #infinity) {
			total_N_recal <- total_C_recal/(exp.C_N_recal_init + cycle * (exp.C_N_recal_final - exp.C_N_recal_init)/exp.steps);
		} else {
			total_N_recal <- 0.0;	
		}
		float total_P_recal;
		if(exp.C_P_recal_final < #infinity and exp.C_P_recal_init < #infinity) {
			total_P_recal <- total_C_recal/(exp.C_P_recal_init + cycle * (exp.C_P_recal_final - exp.C_P_recal_init)/exp.steps);
		} else {
			total_P_recal <- 0.0;
		}
		
		float C_dom <- (exp.C_dom_init + cycle * (exp.C_dom_final - exp.C_dom_init)/exp.steps)#gram;
		float N_dom;
		if(exp.C_N_dom_final < #infinity and exp.C_N_dom_init < #infinity) {
			N_dom <- C_dom/(exp.C_N_dom_init + cycle * (exp.C_N_dom_final - exp.C_N_dom_init)/exp.steps);
		} else {
			N_dom <- 0.0;	
		}
		float P_dom;
		if(exp.C_P_dom_final < #infinity and exp.C_P_dom_init < #infinity) {
			P_dom <- C_dom/(exp.C_P_dom_init + cycle * (exp.C_P_dom_final - exp.C_P_dom_init)/exp.steps);
		} else {
			P_dom <- 0.0;
		}
		
		create DecompositionProblem with: [
			dt::dt,
			extracted_CN::extracted_CN,
			phosphatase_CP::phosphatase_CP,
			alpha_P_er::alpha_P_er,
			beta_CNP_N::beta_CNP_N,
			beta_CNP_P::beta_CNP_P,
			beta_N_P::beta_N_P,
			C_labile_init::total_C_labile,
			N_labile_init::total_N_labile,
			P_labile_init::total_P_labile,
			C_recal_init::total_C_recal,
			N_recal_init::total_N_recal,
			P_recal_init::total_P_recal,
			C_DOM_init::C_dom,
			N_DOM_init::N_dom,
			P_DOM_init::P_dom,
			P_DIM_init::0.0,
			N_DIM_init::0.0
		] {
			myself.decomposition_problem <- self;
		}
		create EnzymaticActivityProblem with: [
			decomposition_problem::decomposition_problem,
			C_N::exp.C_N_microbes/exp.CUE,
			C_P::exp.C_P_microbes/exp.CUE,
			C_microbes::exp.C_microbes#gram,
			min_enzymes::min_enzymes,
			max_enzymes::max_enzymes
		] {
			myself.problem <- self;
		}
		
		create SimulatedAnnealing with:[
			objectives::exp.objectives,
			N::N,
			epsilon::epsilon
		] {
			do optimize(myself.problem);
			
			float C_avail <- s.problem.decomposition_problem.C_avail_final(s.decomposition);
			float N_avail <- s.problem.decomposition_problem.N_avail_final(s.decomposition);
			float P_avail <- s.problem.decomposition_problem.P_avail_final(s.decomposition);
			write "";
			write "s: " + s;
			write "T_CNP: " + s.enzymes.T_CNP / (#gram / #gram / #d);
			write "T_N: " + s.enzymes.T_N / (#gram / #gram / #d);
			write "T_P: " + s.enzymes.T_P / (#gram / #gram / #d);
			write "C/N: " + (N_avail > 0 ? C_avail / N_avail : 0.0);
			write "C/P: " + (P_avail > 0 ? C_avail / P_avail : 0.0);
			write "C_DOM: " + myself.decomposition_problem.C_DOM_init;
			write "N_DOM: " + myself.decomposition_problem.N_DOM_init;
			write "e: " + E(s);
						
			ask s.decomposition {
				write "X_C_eCNP: " + X_C_eCNP / #gram;
				write "X_N_eCNP: " + X_N_eCNP / #gram;
				write "X_P_eCNP: " + X_P_eCNP / #gram;
				write "X_C_eN: " + X_C_eN / #gram;
				write "X_N_eN: " + X_N_eN / #gram;
				write "X_P_eP: " + X_P_eP / #gram;
				write "X_C_er: " + X_C_er / #gram;
				write "X_P_er_recal_to_dim:" + X_P_er_recal_to_dim / #gram;
			}
			exp.simulated_annealing <- self;
		}
	}
	
	reflex when: cycle=steps {
		ask simulation {
			do pause;
		}
	}
	
	output {
		display "e" type:2d {
			chart "e" type:series style:line {
				data "e" value:sum(SimulatedAnnealing collect each.E(each.s)) marker:false;
			}
		}
		display "Labile" type:2d {
			chart "Labile" type:series style:line {
				data "C labile (decomp)" value: sum(SimulatedAnnealing collect each.s.problem.decomposition_problem.C_labile_decomp(each.s.decomposition))/#gram marker:false;
				data "N labile (decomp)" value: sum(SimulatedAnnealing collect each.s.problem.decomposition_problem.N_labile_decomp(each.s.decomposition))/#gram marker:false;
				data "P labile (decomp)" value: sum(SimulatedAnnealing collect each.s.problem.decomposition_problem.P_labile_decomp(each.s.decomposition))/#gram marker:false;
				data "C labile (final)" value: sum(SimulatedAnnealing collect each.s.problem.decomposition_problem.C_labile_final(each.s.decomposition))/#gram marker:false;
				data "N labile (final)" value: sum(SimulatedAnnealing collect each.s.problem.decomposition_problem.N_labile_final(each.s.decomposition))/#gram marker:false;
				data "P labile (final)" value: sum(SimulatedAnnealing collect each.s.problem.decomposition_problem.P_labile_final(each.s.decomposition))/#gram marker:false;
			}
		}
		display "DAM" type:2d {
			chart "DAM" type:series style:line {
				data "C (avail)" value: sum(SimulatedAnnealing collect each.s.problem.decomposition_problem.C_avail_final(each.s.decomposition))/#gram marker:false;
				data "N (avail)" value: sum(SimulatedAnnealing collect each.s.problem.decomposition_problem.N_avail_final(each.s.decomposition))/#gram marker:false;
				data "P (avail)" value: sum(SimulatedAnnealing collect each.s.problem.decomposition_problem.P_avail_final(each.s.decomposition))/#gram marker:false;
			}
		}
		display "Enzymatic rates" type:2d {
			chart "Enzymatic rates" type:series style:line {
				data "T_CNP" value: sum(SimulatedAnnealing collect (each.s.enzymes.T_CNP / (#gram / #gram / #d))) marker:false;
				data "T_N" value: sum(SimulatedAnnealing collect (each.s.enzymes.T_N / (#gram / #gram / #d))) marker:false;
				data "T_P" value: sum(SimulatedAnnealing collect (each.s.enzymes.T_P / (#gram / #gram / #d))) marker:false;
				data "T_recal" value: sum(SimulatedAnnealing collect (each.s.enzymes.T_recal / (#gram / #gram / #d))) marker:false;
				if show_max_rates {
					data "T_CNP (max)" value: problem != nil ? problem.max_enzymes.T_CNP / (#gram / #gram / #d) : 0.0 marker:false;
					data "T_N (max)" value: problem != nil ? problem.max_enzymes.T_N / (#gram / #gram / #d) : 0.0 marker:false;
					data "T_P (max)" value: problem != nil ? problem.max_enzymes.T_P / (#gram / #gram / #d) : 0.0 marker:false;
					data "T_recal (max)" value: problem != nil ? problem.max_enzymes.T_recal / (#gram / #gram / #d) : 0.0 marker:false;			
				}
			}
		}
		display "C/N" type:2d {
			chart "C/N" type:series style:line {
				data "available C/N" value: sum(SimulatedAnnealing collect (
					each.s.problem.decomposition_problem.N_avail_final(each.s.decomposition) > 0 ?
					each.s.problem.decomposition_problem.C_avail_final(each.s.decomposition)
						/ each.s.problem.decomposition_problem.N_avail_final(each.s.decomposition) : 0.0
				)) marker:false;
				data "C/N" value: problem != nil ? problem.C_N : 0.0 marker:false;
			}
		}
		display "C/P" type:2d {
			chart "C/P" type:series style:line {
				data "available C/P" value: sum(SimulatedAnnealing collect (
					each.s.problem.decomposition_problem.P_avail_final(each.s.decomposition) > 0 ?
					each.s.problem.decomposition_problem.C_avail_final(each.s.decomposition)
						/ each.s.problem.decomposition_problem.P_avail_final(each.s.decomposition) : 0.0
				)) marker:false;
				data "C/P" value: problem != nil ? problem.C_P : 0.0 marker:false;
			}
		}
	}
}
