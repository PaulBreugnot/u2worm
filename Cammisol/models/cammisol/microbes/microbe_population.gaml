/**
* Name: microbe_population
* Based on the internal skeleton template. 
* Author: Nicolas Marilleau, Lucas Grosjean, Paul Breugnot
* Tags: 
*/

model microbe_population

import "../environment/dam.gaml"

global {
	float microbe_CO2_emissions <- 0.0#gram;
	float enzymes_optimization_period <- 1#d;
	
	float microbes_CN <- 10.0;
	float microbes_CP <- 17.0;
	
	// Enzyme objectives
	MaxLabileC max_labile_C;
	MaxC max_C;
	CapCN max_CN;
	CapCP max_CP;
	
	DecompositionProblem decomposition_problem;
	EnzymaticActivityProblem enzymatic_activity_problem;
	SimulatedAnnealing simulated_annealing;
	
	action init_enzymatic_optimisation {
		create MaxLabileC with: (weight: 1.0) {
			max_labile_C <- self;
		}
		create MaxC with: (weight: 5.0) {
			max_C <- self;
		}
		create CapCN with: (weight: 10.0) {
			max_CN <- self;
		}
		create CapCP with: (weight: 10.0) {
			max_CP <- self;
		}
		
		create DecompositionProblem with: [
			dt::enzymes_optimization_period
		] {
			decomposition_problem <- self;
		}
		create EnzymaticActivityProblem with: [
			decomposition_problem: decomposition_problem
		] {
			enzymatic_activity_problem <- self;
		}
		create SimulatedAnnealing with:[
			objectives::[max_labile_C, max_C, max_CN, max_CP],
			N::1000,
			epsilon::0.0
		] {
			simulated_annealing <- self;	
		}
	}
}

/**
 * Microbe populations represent a spatially confined aggregate of microbes with
 * its own functional characteristics. Microbes in a population do not
 * necessarily belong to a common species, but rather share a common functional
 * behavior, represented by parameters defined at the scale of each population. 
 */
species MicrobePopulation schedules:[]
{
	float carbon_use_efficiency;
	/**
	 * Relative growth rate coefficients used in classical growth models.
	 * https://en.wikipedia.org/wiki/Population_growth#Logistic_equation
	 */
	float dividing_time;
	/**
	 * Structural C:N rate.
	 */
	float C_N <- microbes_CN;
	/**
	 * Structural C:P rate.
	 */
	float C_P <- microbes_CP;
	
	/**
	 * Structural C.
	 * 
	 * The structural C represents C fixed by the population.
	 */
	float C <- 0#gram;
	/**
	 * Structural N.
	 */
	float N <- 0#gram;
	/**
	 * Structural P.
	 */
	float P <- 0#gram;
	
	/**
	 * Rate of the population that is currently awoken and active.
	 */
	float active_rate <- 1.0;
	/**
	 * Minimum rate of the population that is always awoken and active.
	 */
	float minimum_active_rate;
	/**
	 * Characteristic sporulation time (time to go dormant or inactive).
	 */
	float sporulation_time <- 10#h;
	/**
	 * Characteristic germination time (time to go awoken or active).
	 */
	float germination_time <- 1#h;
	
	/**
	 * Cytosol C.
	 * 
	 * The cytosol is a nutrient storage that is not yet assimilated by the
	 * population. It can either be assimilated later or rejected in the DAM.
	 */
	float cytosol_C <- 0#gram; 
	/**
	 * Cytosol N.
	 */
	float cytosol_N; 
	/**
	 * Cytosol P.
	 */
	float cytosol_P;
	
	/**
	 * Defines the capacity of the population to store nutrients in the cytosol.
	 * At each time, after assimilation, growth of respiration, a rate
	 * cytosol_mineralization_rate of each C/N/P quantity is left in the cytosol
	 * while the rest is rejected as CO2 and in the DIM.
	 */
	float cytosol_mineralization_rate <- 0.0;
	
	/**
	 * Current enzymatic activities.
	 */
	Enzymes enzymes;
	/**
	 * Constitutive enzymatic activities.
	 */
	Enzymes min_enzymes;
	/**
	 * Maximum enzymatic activities.
	 */
	Enzymes max_enzymes;
	
	/**
	 * Currently active C.
	 */
	float active_C;
	/**
	 * C requested at the current time step to maintain an ideal exponential
	 * growth.
	 */
	float requested_C;
	/**
	 * N requested at the current time step to maintain an ideal exponential
	 * growth.
	 */
	float requested_N;
	/**
	 * P requested at the current time step to maintain an ideal exponential
	 * growth.
	 */
	float requested_P;
	/**
	 * C:N rate requested at the current time step to maintain an ideal
	 * exponential growth.
	 */
	float requested_C_N;
	/**
	 * C:P rate requested at the current time step to maintain an ideal
	 * exponential growth.
	 */
	float requested_C_P;
	
	init {
		create Enzymes {
			myself.enzymes <- self;
		}
		
		N <- C / C_N;
		P <- C / C_P;

		cytosol_N <- cytosol_C / C_N;
		cytosol_P <- cytosol_C / C_P;
	}
	
	/**
	 * Updates active C and requested C/N/P for the current time step.
	 */
	action update {
		active_C <- C * active_rate;
		requested_C <- active_C * local_step / (carbon_use_efficiency * dividing_time);
		requested_N <- carbon_use_efficiency * requested_C / C_N;
		requested_P <- carbon_use_efficiency * requested_C / C_P;
		requested_C_N <- C_N / carbon_use_efficiency;
		requested_C_P <- C_P / carbon_use_efficiency;
	}
	
	/**
	 * Assimilates nutrients from the specified dam to the cytosol.
	 */
	action assimilate(Dam dam, float perceived_rate) {
		// C/N/P is retrieved in priority from the cytosol. Is enough nutrients
		// is available in the cytosol, nothing is consumed from the perceived
		// C/N/P.
		float assimilated_C <- min(dam.dom[2] * perceived_rate, max(0.0, requested_C-cytosol_C));
		float assimilated_N <- min(dam.dom[0] * perceived_rate, max(0.0, requested_N-cytosol_N));
		float assimilated_P <- min(dam.dom[1] * perceived_rate, max(0.0, requested_P-cytosol_P));

		cytosol_C <- cytosol_C + assimilated_C;
		cytosol_N <- cytosol_N + assimilated_N;
		cytosol_P <- cytosol_P + assimilated_P;

		dam.dom[0] <- dam.dom[0] - assimilated_N;
		dam.dom[1] <- dam.dom[1] - assimilated_P;
		dam.dom[2] <- dam.dom[2] - assimilated_C;

		assimilated_N <- min(dam.dim[0] * perceived_rate, max(0.0, requested_N-cytosol_N));
		assimilated_P <- min(dam.dim[1] * perceived_rate, max(0.0, requested_P-cytosol_P));

		cytosol_N <- cytosol_N + assimilated_N;
		cytosol_P <- cytosol_P + assimilated_P;

		dam.dim[0] <- dam.dim[0] - assimilated_N;
		dam.dim[1] <- dam.dim[1] - assimilated_P;
	}
	
	/**
	 * Updates active rate according to perceived nutrients.
	 */
	action sleep(Dam dam, float perceived_rate) {
		float perceived_C_respiration <- (1-carbon_use_efficiency) * (perceived_rate * dam.dom[2]+cytosol_C);
		float perceived_C_growth <- min(
			carbon_use_efficiency * (perceived_rate * dam.dom[2]+cytosol_C),
			(perceived_rate * (dam.dom[0]+dam.dim[0])+cytosol_N) * C_N,
			(perceived_rate * (dam.dom[1]+dam.dim[1])+cytosol_P) * C_P
		);		
		float max_awake_population <- min(1.0, (perceived_C_respiration + perceived_C_growth) / (C * local_step / (carbon_use_efficiency * dividing_time)));
		
		float d_awake_population <- 1/(active_rate > max_awake_population ? sporulation_time : germination_time)
			* active_rate * (1 - active_rate / max(minimum_active_rate, max_awake_population));

		// The minimum awake rate prevents awake_population to be set to 0.0, what would prevent it to awake (see d_awake_population equation above)
		active_rate <- min(1.0, max(minimum_active_rate, active_rate + d_awake_population * local_step));
	}
	
	/**
	 * Reject CO2 from cytosol C, according to the CUE of the population.
	 */
	action respirate
	{
		float C_respiration <- (1-carbon_use_efficiency) * cytosol_C;
		cytosol_C <- cytosol_C - C_respiration; 
		
		microbe_CO2_emissions <- microbe_CO2_emissions + C_respiration;
	}
	
	/**
	 * Consumes nutrients from the cytosol to grow, fixing C/N/P in structural
	 * C/N/P.
	 * 
	 * Growth is limited by the carrying capacity.
	 */
	action grow(Dam dam, float total_C_in_pore, float carrying_capacity)
	{		
		float C_usable_for_growth <- min(
			cytosol_C,
			cytosol_N*C_N,
			cytosol_P*C_P
		);
		
		float C_used_for_growth <- max(0.0, C_usable_for_growth * (1.0 - total_C_in_pore/carrying_capacity));
		
		// C/N/P assimilated for growth but finally not used due to carrying capacity
		float C_not_used_for_growth <- C_usable_for_growth - C_used_for_growth;
		float N_not_used_for_growth <- C_not_used_for_growth / C_N;
		float P_not_used_for_growth <- C_not_used_for_growth / C_P;
		
		cytosol_C <- cytosol_C - C_not_used_for_growth;
		cytosol_N <- cytosol_N - N_not_used_for_growth;
		cytosol_P <- cytosol_P - P_not_used_for_growth;
				
		ask dam {
			dom[2] <- dom[2] + C_not_used_for_growth;
			dom[0] <- dom[0] + N_not_used_for_growth;
			dom[1] <- dom[1] + P_not_used_for_growth;
		}
		
		float N_growth <- C_used_for_growth / C_N;
		float P_growth <- C_used_for_growth / C_P;
		
		cytosol_C <- cytosol_C - C_used_for_growth;
		cytosol_N <- cytosol_N - N_growth;
		cytosol_P <- cytosol_P - P_growth;
		
		C <- C + C_used_for_growth;
		N <- N + N_growth;
		P <- P + P_growth;
	}
	
	/**
	 * Reject C/N/P left in the cytosol as CO2 and mineralized N/P, according to
	 * the cytosol_mineralization_rate. If cytosol_mineralization_rate is 1.0,
	 * everything is left in the cytosol.
	 */
	action reject(Dam dam) {	
		float left_C_in_cytosol <- cytosol_mineralization_rate * cytosol_C;
		float left_N_in_cytosol <- cytosol_mineralization_rate * cytosol_N;
		float left_P_in_cytosol <- cytosol_mineralization_rate * cytosol_P;
		
		ask dam {
			dim[0] <- dim[0] + myself.cytosol_N - left_N_in_cytosol;
			dim[1] <- dim[1] + myself.cytosol_P - left_P_in_cytosol;
		}
		microbe_CO2_emissions <- microbe_CO2_emissions + cytosol_C - left_C_in_cytosol;
				
		cytosol_C <- left_C_in_cytosol;
		cytosol_N <- left_N_in_cytosol;
		cytosol_P <- left_P_in_cytosol;
	}
	
	/**
	 * Complete microbe metabolism.
	 */
	action life(
		Dam dam, float perceived_rate, float total_C_in_pore, float pore_carrying_capacity
	) {	
		do sleep(dam, perceived_rate);
		do assimilate(dam, perceived_rate);
		do respirate;
		do grow(dam, total_C_in_pore, pore_carrying_capacity);
		do reject(dam);
	}
	
	/**
	 * Optimizes the current enzymatic activity according to the composition of
	 * the environment, using the buzy-pop model implemented in the
	 * enzymatic-activity module.
	 */
	action optimize_enzymes(Dam dam, list<OrganicParticle> particles_to_decompose) {
		MicrobePopulation population <- self;
		ask decomposition_problem {
			// Updates initial substrate and available matter
			// Organic particles around are considered as a single and aggregated organic compartment
			C_labile_init <- sum(particles_to_decompose collect each.C_labile);
			N_labile_init <- sum(particles_to_decompose collect each.N_labile);
			P_labile_init <- sum(particles_to_decompose collect each.P_labile);
			C_recal_init <- sum(particles_to_decompose collect each.C_recalcitrant);
			N_recal_init <- sum(particles_to_decompose collect each.N_recalcitrant);
			P_recal_init <- sum(particles_to_decompose collect each.P_recalcitrant);
			C_DOM_init <- dam.dom[2];
			N_DOM_init <- dam.dom[0];
			P_DOM_init <- dam.dom[1];
			P_DIM_init <- dam.dim[1];
			N_DIM_init <- dam.dim[0];
		}
		ask enzymatic_activity_problem {
			// Update required C/N and C/P rate, and microbe biomass
			C_N <- myself.requested_C_N;
			C_P <- myself.requested_C_P;
			C_microbes <- myself.active_C;
			min_enzymes <- myself.min_enzymes;
			max_enzymes <- myself.max_enzymes;
		}
		// TODO: Optimize this.
		ask simulated_annealing {
			do optimize(enzymatic_activity_problem);
			ask population.enzymes {
				do die;
			}
			population.enzymes <- s.enzymes;
			ask s {
				ask weighted_enzymes {
					do die;
				}
				ask decomposition {
					do die;
				}
				do die;
			}
		}
	}
}
