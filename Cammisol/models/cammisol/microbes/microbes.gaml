/**
* Name: microbes
* Author: Paul Breugnot
* 
* Instantiation of O-F-M Strategists (Opportunists, Foragers, Minimalists) as
* microbe populations.
*/

model microbes

import "../environment/pore_particle.gaml"

/**
 * Parameter specialization for O-F-M (Opportunists, Foragers, Minimalists)
 * strategies, adapted from the Y-A-S Strategies [1]. See the documentation of
 * each class for a detailed description of each strategist.
 *
 * See the documentation of the microbe_population module for the meaning of
 * each parameter.
 * 
 * [1] Malik, A.A., Martiny, J.B.H., Brodie, E.L., Martiny, A.C., Treseder,
 * K.K., Allison, S.D., 2020. Defining trait-based microbial strategies with
 * consequences for soil carbon cycling under climate change.
 * https://doi.org/10.1038/s41396-019-0510-0
 */
global {
	list<species<MicrobePopulation>> bacteria_types <- [O_Strategist, F_Strategist, M_Strategist];

	float dividing_time_O <- 1#h;
	float dividing_time_F <- 24#h;
	float dividing_time_M <- 368#h;
	
	float carbon_use_efficiency_O <- 0.7;
	float carbon_use_efficiency_F <- 0.3;
	float carbon_use_efficiency_M <- 0.5;
	
	float minimum_active_rate_O <- 0.002;
	float minimum_active_rate_F <- 0.018;
	float minimum_active_rate_M <- 1.0;
	
	float min_T_CNP_O <- 0.0;
	float min_T_N_O <- 0.0;
	float min_T_P_O <- 0.0;
	float min_T_recal_O <- 0.0;
	float max_T_CNP_O <- 0.6 #gram / #gram / #d;
	float max_T_N_O <- 0.2 #gram / #gram / #d;
	float max_T_P_O <- 0.04 #gram / #gram / #d;
	float max_T_recal_O <- 1e-6 #gram / #gram / #d;
	
	float min_T_CNP_F <- 0.0;
	float min_T_N_F <- 0.0;
	float min_T_P_F <- 0.0;
	float min_T_recal_F <- 0.0;
	float max_T_CNP_F <- 3 #gram / #gram / #d;
	float max_T_N_F <- 1 #gram / #gram / #d;
	float max_T_P_F <- 0.2 #gram / #gram / #d;
	float max_T_recal_F <- 1 #gram / #gram / #d;
	
	float min_T_CNP_M <- 0.2 #gram/ #gram / #d;
	float min_T_N_M <- 0.008 #gram/ #gram / #d;
	float min_T_P_M <- 0.0016 #gram/ #gram / #d;
	float min_T_recal_M <- 0.0;
	float max_T_CNP_M <- 0.25 #gram / #gram / #d;
	float max_T_N_M <- 0.028 #gram / #gram / #d;
	float max_T_P_M <- 0.0056 #gram / #gram / #d;
	float max_T_recal_M <- 1e-6 #gram / #gram / #d;
	
	Enzymes min_enzymes_O;
	Enzymes max_enzymes_O;
	Enzymes min_enzymes_F;
	Enzymes max_enzymes_F;
	Enzymes min_enzymes_M;
	Enzymes max_enzymes_M;
	
	/**
	 * This method must be called from the main global method to allow the usage of Y-A-S_Strategists.
	 */
	action init_enzymes {
		create Enzymes with: [
			name::"Min enzymes (O)",
			T_CNP::min_T_CNP_O,
			T_N::min_T_N_O,
			T_P::min_T_P_O,
			T_recal::min_T_recal_O
		] {
			myself.min_enzymes_O <- self;
		}
		
		create Enzymes with: [
			name::"Max enzymes (O)",
			T_CNP::max_T_CNP_O,
			T_N::max_T_N_O,
			T_P::max_T_P_O,
			T_recal::max_T_recal_O
		] {
			myself.max_enzymes_O <- self;
		}
		
		create Enzymes with: [
			name::"Min enzymes (F)",
			T_CNP::min_T_CNP_F,
			T_N::min_T_N_F,
			T_P::min_T_P_F,
			T_recal::min_T_recal_F
		] {
			myself.min_enzymes_F <- self;
		}
		
		create Enzymes with: [
			name::"Max enzymes (F)",
			T_CNP::max_T_CNP_F,
			T_N::max_T_N_F,
			T_P::max_T_P_F,
			T_recal::max_T_recal_F
		] {
			myself.max_enzymes_F <- self;
		}	
		
		create Enzymes with: [
			name::"Min enzymes (M)",
			T_CNP::min_T_CNP_M,
			T_N::min_T_N_M,
			T_P::min_T_P_M,
			T_recal::min_T_recal_M
		] {
			myself.min_enzymes_M <- self;
		}
		
		create Enzymes with: [
			name::"Max enzymes (M)",
			T_CNP::max_T_CNP_M,
			T_N::max_T_N_M,
			T_P::max_T_P_M,
			T_recal::max_T_recal_M
		] {
			myself.max_enzymes_M <- self;
		}	
	}
}

/**
 * Opportunists (O strategists) correspond to the _Yield_ strategy of the
 * Y-A-S classification. They can decompose labile OM, but not recalcitrant OM,
 * as they invest most of their metabolism in growth (assimilation and fixation)
 * rather than in complex enzymatic systems. While enough nutrients are
 * available, they grow rapidly and exponentially. But once not enough nutrients
 * are available, they rapidly go dormant (sporulation). Their focus on growth
 * is such that they cannot decompose as much nutrients as they can assimilate:
 * as opportunists, they benefit from the decomposition of nutrients by other
 * populations and the amendment of ready to assimilate available nutrients.
 */
species O_Strategist parent:MicrobePopulation schedules:[] {
	init {
		dividing_time <- dividing_time_O;
		carbon_use_efficiency <- carbon_use_efficiency_O;
		minimum_active_rate <- minimum_active_rate_O;
		
		cytosol_mineralization_rate <- 1.0;
		
		min_enzymes <- min_enzymes_O;
		max_enzymes <- max_enzymes_O;
	}
}

/**
 * Foragers (F strategists) correspond to the _resource Acquisition_ strategy of
 * the Y-A-S classification. They grow slower than O strategists, to invest more
 * of their metabolism in enzymatic systems. In consequence, if not enough
 * labile OM is available, they have the ability to forage recalcitrant OM, but
 * only if it is necessary to maintain growth. They also go dormant once not
 * enough nutrients are available.
 */
species F_Strategist parent:MicrobePopulation schedules:[] {
	init {
		dividing_time <- dividing_time_F;
		carbon_use_efficiency <- carbon_use_efficiency_F;
		minimum_active_rate <- minimum_active_rate_F;
		
		cytosol_mineralization_rate <- 0.0;
		
		min_enzymes <- min_enzymes_F;
		max_enzymes <- max_enzymes_F;
	}
}

/**
 * Minimalists (M strategists) only partly correspond to the _Stress tolerant_
 * strategy of the Y-A-S decomposition. Indeed, we think that stress tolerance
 * includes many different strategies also shared by Y and A classes. We finally
 * prefer to dedicate the Minimalist strategy to strong oligotrophy, with an
 * emphasis on a low cost metabolism. Such slow metabolism allow _Minimalists_
 * to survive in poor environments, reducing competition for the available
 * substrate with an high uptake efficiency optimised for environments with a
 * low nutrient concentration [2]. In order to limit the energy requirement of
 * their metabolism, they do not have the ability to produce complex enzymes as
 * it would require complex genomes that are costly to maintain: they can
 * decompose simple labile compounds, but not recalcitrant OM. For the same
 * reason, they do not have much adaptation capacity to the environment: they
 * survive in any circumstance and never go dormant, but they have a limited
 * capacity to adapt to the environment. In consequence, they will never grow
 * rapidly, even if lots of nutrients are available.
 * 
 * [2] Coche, A., Babey, T., Rapaport, A., Gonod, L.V., Garnier, P., Nunan, N.,
 * de Dreuzy, J.-R., 2022. Competition within low-density bacterial populations
 * as an unexpected factor regulating carbon decomposition in bulk soil.
 * https://doi.org/10.1016/j.soilbio.2021.108423
 */
species M_Strategist parent:MicrobePopulation schedules:[] {
	init {
		dividing_time <- dividing_time_M;
		carbon_use_efficiency <- carbon_use_efficiency_M;
		minimum_active_rate <- minimum_active_rate_M;

		cytosol_mineralization_rate <- 0.0;
		
		min_enzymes <- min_enzymes_M;
		max_enzymes <- max_enzymes_M;
	}
}
