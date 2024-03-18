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
	
	// Enzyme objectives
	MaxLabileC max_labile_C;
	MaxRecalC max_recal_C;
	MaxN max_N;
	MaxP max_P;
	MaxCN max_CN;
	MaxCP max_CP;
	
	action init_microbes_objectives {
		create MaxLabileC with: (weight: 5.0) {
			max_labile_C <- self;
		}
		create MaxRecalC with: (weight: 1.0) {
			max_recal_C <- self;
		}
		create MaxN with: (weight: 5.0) {
			max_N <- self;
		}
		create MaxP with: (weight: 5.0) {
			max_P <- self;
		}
		create MaxCN with: (weight: 10.0) {
			max_CN <- self;
		}
		create MaxCP with: (weight: 10.0) {
			max_CP <- self;
		}
	}
}


species MicrobePopulation schedules:[]
{
	float carbon_use_efficiency;
	
	float dividing_time;
	float C_N <- 10.0;
	float C_P <- 17.0;
	
	// TODO: unit?
	// TODO: source?
	float C <- 0#gram;
	float N;
	float P;
	
	float awake_population <- 1.0;
	float minimum_awake_rate;
	float sporulation_time <- 10#h;
	float germination_time <- 1#h;
	
	// Mon estomac le cytosol 
	float cytosol_C <- 0#gram; 
	float cytosol_N; 
	float cytosol_P;
	
	float cytosol_mineralization_rate <- 0.0;
	
	Enzymes enzymes;
	Enzymes min_enzymes;
	Enzymes max_enzymes;
	DecompositionProblem decomposition_problem;
	EnzymaticActivityProblem enzymatic_activity_problem;
	
	float active_C;
	float requested_C;
	float requested_N;
	float requested_P;
	float requested_C_N;
	float requested_C_P;
	
	init {
		create Enzymes {
			myself.enzymes <- self;
		}
		create DecompositionProblem with: [
			dt::enzymes_optimization_period
		] {
			myself.decomposition_problem <- self;
		}
		create EnzymaticActivityProblem with: [
			decomposition_problem: decomposition_problem
		] {
			myself.enzymatic_activity_problem <- self;
		}
		
		N <- C / C_N;
		P <- C / C_P;

		cytosol_N <- cytosol_C / C_N;
		cytosol_P <- cytosol_C / C_P;
	}
	
	action set_min_max_enzymes(Enzymes _min_enzymes, Enzymes _max_enzymes) {
		self.min_enzymes <- _min_enzymes;
		self.max_enzymes <- _max_enzymes;
		ask enzymatic_activity_problem {
			min_enzymes <- _min_enzymes;
			max_enzymes <- _max_enzymes;
		}
	}
	
	action respirate(float C_respiration)
	{
		cytosol_C <- cytosol_C - C_respiration; 
		
		microbe_CO2_emissions <- microbe_CO2_emissions + C_respiration;
	}
	
	action growth(float C_growth)
	{
		float N_growth <- C_growth / C_N;
		float P_growth <- C_growth / C_P;
		
		cytosol_C <- cytosol_C - C_growth;
		cytosol_N <- cytosol_N - N_growth;
		cytosol_P <- cytosol_P - P_growth;
		
		C <- C + C_growth;
		N <- N + N_growth;
		P <- P + P_growth;
	}
	
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
		}
		create SimulatedAnnealing with:[
				problem::self.enzymatic_activity_problem,
				objectives::(max_enzymes.T_recal > 0 ? [max_labile_C, max_recal_C, max_N, max_P, max_CN, max_CP] : [max_labile_C, max_N, max_P, max_CN, max_CP]),
				N::1000,
				epsilon::0.0
			] {
			do optimize;
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
			do die;
		}
	}
	
	action update {
		active_C <- C * awake_population;
		requested_C <- active_C * local_step / (carbon_use_efficiency * dividing_time);
		requested_N <- carbon_use_efficiency * requested_C / C_N;
		requested_P <- carbon_use_efficiency * requested_C / C_P;
		requested_C_N <- C_N / carbon_use_efficiency;
		requested_C_P <- C_P / carbon_use_efficiency;
	}
	
	action dormancy(Dam dam, float perceived_rate) {
		float perceived_C_respiration <- (1-carbon_use_efficiency) * perceived_rate * (dam.dom[2]+cytosol_C);
		float perceived_C_growth <- min(
			carbon_use_efficiency * perceived_rate * (dam.dom[2]+cytosol_C),
			perceived_rate * (dam.dom[0]+dam.dim[0]+cytosol_N) * C_N,
			perceived_rate * (dam.dom[1]+dam.dim[1]+cytosol_P) * C_P
		);		
		float max_awake_population <- min(1.0, (perceived_C_respiration + perceived_C_growth) / (C * local_step / (carbon_use_efficiency * dividing_time)));
		float d_awake_population <- 1/(awake_population > max_awake_population ? sporulation_time : germination_time)
			* awake_population * (1 - awake_population / max(minimum_awake_rate, max_awake_population));
		awake_population <- min(1.0, max(0.0, awake_population + d_awake_population * local_step));
	}
	
	action assimilate(Dam dam, float perceived_rate) {
			// If total_X_consumed = total_X_requested, X_consum = X_requested for all microbe population
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
	
	action life(
		Dam dam, list<OrganicParticle> accessible_organics, float total_C_in_pore, float pore_carrying_capacity
	)
	{	
		float C_respiration <- (1-carbon_use_efficiency) * cytosol_C;
		float C_usable_for_growth <- min(
			carbon_use_efficiency * cytosol_C,
			cytosol_N*C_N,
			cytosol_P*C_P
		);
		
		float C_used_for_growth <- max(0.0, C_usable_for_growth * (1.0 - total_C_in_pore/pore_carrying_capacity));
		
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
		
		do respirate(C_respiration);
		do growth(C_used_for_growth);

		
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
}
