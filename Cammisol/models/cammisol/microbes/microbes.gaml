/**
* Name: microbes
* Author: Paul Breugnot
* 
* Instantiation of Y-A-S Strategists as microbe populations.
*/

model microbes

import "../environment/pore_particle.gaml"

/**
 * Parameter specialization for Y-A-S strategies. See the documentation of the
 * microbe_population module for the meaning of each parameter.
 */
global {
	list<species<MicrobePopulation>> bacteria_types <- [Y_Strategist, A_Strategist, S_Strategist];

	float dividing_time_Y <- 1#h;
	float dividing_time_A <- 24#h;
	float dividing_time_S <- 368#h;
	
	float carbon_use_efficiency_Y <- 0.7;
	float carbon_use_efficiency_A <- 0.3;
	float carbon_use_efficiency_S <- 0.5;
	
	float minimum_active_rate_Y <- 0.002;
	float minimum_active_rate_A <- 0.018;
	float minimum_active_rate_S <- 1.0;
	
	float min_T_C_Y <- 0.0;
	float min_T_N_Y <- 0.0;
	float min_T_P_Y <- 0.0;
	float min_T_recal_Y <- 0.0;
	float max_T_C_Y <- 1 #gram / #gram / #d;
	float max_T_N_Y <- 0.2 #gram / #gram / #d;
	float max_T_P_Y <- 0.04 #gram / #gram / #d;
	float max_T_recal_Y <- 1e-6 #gram / #gram / #d;
	
	float min_T_C_A <- 0.0;
	float min_T_N_A <- 0.0;
	float min_T_P_A <- 0.0;
	float min_T_recal_A <- 0.0;
	float max_T_C_A <- 5 #gram/ #gram / #d;
	float max_T_N_A <- 1 #gram / #gram / #d;
	float max_T_P_A <- 0.2 #gram / #gram / #d;
	float max_T_recal_A <- 1 #gram / #gram / #d;
	
	float min_T_C_S <- 0.04 #gram/ #gram / #d;
	float min_T_N_S <- 0.05 #gram/ #gram / #d;
	float min_T_P_S <- 0.005 #gram/ #gram / #d;
	float min_T_recal_S <- 0.0;
	float max_T_C_S <- 0.14 #gram/ #gram / #d;
	float max_T_N_S <- 0.07 #gram / #gram / #d;
	float max_T_P_S <- 0.009 #gram / #gram / #d;
	float max_T_recal_S <- 1e-6 #gram / #gram / #d;
	
	Enzymes min_enzymes_Y;
	Enzymes max_enzymes_Y;
	Enzymes min_enzymes_A;
	Enzymes max_enzymes_A;
	Enzymes min_enzymes_S;
	Enzymes max_enzymes_S;
	
	/**
	 * This method must be called from the main global method to allow the usage of Y-A-S_Strategists.
	 */
	action init_enzymes {
		create Enzymes with: [
			name::"Min enzymes (Y)",
			T_C::min_T_C_Y,
			T_N::min_T_N_Y,
			T_P::min_T_P_Y,
			T_recal::min_T_recal_Y
		] {
			myself.min_enzymes_Y <- self;
		}
		
		create Enzymes with: [
			name::"Max enzymes (Y)",
			T_C::max_T_C_Y,
			T_N::max_T_N_Y,
			T_P::max_T_P_Y,
			T_recal::max_T_recal_Y
		] {
			myself.max_enzymes_Y <- self;
		}
		
		create Enzymes with: [
			name::"Min enzymes (A)",
			T_C::min_T_C_A,
			T_N::min_T_N_A,
			T_P::min_T_P_A,
			T_recal::min_T_recal_A
		] {
			myself.min_enzymes_A <- self;
		}
		
		create Enzymes with: [
			name::"Max enzymes (A)",
			T_C::max_T_C_A,
			T_N::max_T_N_A,
			T_P::max_T_P_A,
			T_recal::max_T_recal_A
		] {
			myself.max_enzymes_A <- self;
		}	
		
		create Enzymes with: [
			name::"Min enzymes (S)",
			T_C::min_T_C_S,
			T_N::min_T_N_S,
			T_P::min_T_P_S,
			T_recal::min_T_recal_S
		] {
			myself.min_enzymes_S <- self;
		}
		
		create Enzymes with: [
			name::"Max enzymes (S)",
			T_C::max_T_C_S,
			T_N::max_T_N_S,
			T_P::max_T_P_S,
			T_recal::max_T_recal_S
		] {
			myself.max_enzymes_S <- self;
		}	
	}
}
/**
 * Y strategists (Yield) can decompose labile OM, but not recalcitrant OM, as
 * they invest most of their metabolism in growth rather than in complex
 * enzymatic systems. While enough nutrients are available, they grow rapidly
 * and exponentially. But once not enough nutrients are available, they rapidly
 * go dormant (sporulation).
 */
species Y_Strategist parent:MicrobePopulation schedules:[] {
	init {
		dividing_time <- dividing_time_Y;
		carbon_use_efficiency <- carbon_use_efficiency_Y;
		minimum_active_rate <- minimum_active_rate_Y;
		
		cytosol_mineralization_rate <- 1.0;
		
		min_enzymes <- min_enzymes_Y;
		max_enzymes <- max_enzymes_Y;
	}
}

/**
 * A stategists (resource Acquisition) grow slowly than Y strategists, to invest
 * more of their metabolism in enzymatic systems. In consequence, if not enough
 * labile OM is available, they have the ability to also decompose recalcitrant
 * OM, but only if it is necessary. They also go dormant once not enough
 * nutrients are available.
 */
species A_Strategist parent:MicrobePopulation schedules:[] {
	init {
		dividing_time <- dividing_time_A;
		carbon_use_efficiency <- carbon_use_efficiency_A;
		minimum_active_rate <- minimum_active_rate_A;
		
		cytosol_mineralization_rate <- 0.0;
		
		min_enzymes <- min_enzymes_A;
		max_enzymes <- max_enzymes_A;
	}
}

/**
 * S strategists (Stress tolerant) have a very slow metabolism to survive in
 * poor environments. In order to limit the energy requirement of their
 * metabolism, they do not have the ability to produce complex enzymes as it
 * would require complex genomes that are costly to maintain: they can decompose
 * simple labile compounds, but not recalcitrant OM. For the same reason, they
 * do not have much adaptation capacity to the environment: they survive in any
 * circumstance and never go dormant, but they have a limited capacity to adapt
 * to the environment. In consequence, they will never grow rapidly, even if
 * lots of nutrients are available.
 */
species S_Strategist parent:MicrobePopulation schedules:[] {
	init {
		dividing_time <- dividing_time_S;
		carbon_use_efficiency <- carbon_use_efficiency_S;
		minimum_active_rate <- minimum_active_rate_S;

		cytosol_mineralization_rate <- 0.0;
		
		min_enzymes <- min_enzymes_S;
		max_enzymes <- max_enzymes_S;
	}
}
