/**
* Name: PoreParticle
* Based on the internal empty template. 
* Author: pbreugno
* Tags: 
*/


model pore_particle

import "microbe_population.gaml"


species PoreParticle {
	int grid_x;
	int grid_y;
	float carrying_capacity;
	
	/**
	 * OrganicParticle embedded within the PoreParticle. This is where nematodes spit out some C/N/P.
	 */
	OrganicParticle organic_particle;
	
	Dam dam;
	
	list<MicrobePopulation> populations;
	
	list<PoreParticle> pore_neighbors;
	list<OrganicParticle> accessible_organics;


	reflex decompose {
		float total_C_labile <- sum(accessible_organics collect each.C_labile);
		float total_C_recal <- sum(accessible_organics collect each.C_recalcitrant);
		
		WeightedEnzymes total_enzymes;
		create WeightedEnzymes with: [
			T_cellulolytic: sum(populations collect (each.active_C() * each.enzymes.T_cellulolytic)),
			T_amino: sum(populations collect (each.active_C() * each.enzymes.T_amino)),
			T_P: sum(populations collect (each.active_C() * each.enzymes.T_P)),
			T_recal: sum(populations collect (each.active_C() * each.enzymes.T_recal))
		] {
			total_enzymes <- self;
		}
		
		DecompositionProblem decomposition_problem;
		create DecompositionProblem with: [
			dt::local_step
			] {
			decomposition_problem <- self;
		}
		Decomposition decomposition;
		create Decomposition {
			decomposition <- self;
		}
		ask accessible_organics {
			OrganicParticle organic <- self;
			PoreParticle pore <- myself;
			
			WeightedEnzymes local_enzymes;
			create WeightedEnzymes with: [
				T_cellulolytic: total_enzymes.T_cellulolytic * C_labile / total_C_labile,
				T_amino: total_enzymes.T_amino * C_labile / total_C_labile,
				T_P: total_enzymes.T_P * C_labile / total_C_labile,
				T_recal: total_enzymes.T_recal * C_recalcitrant / total_C_recal
			] {
				local_enzymes <- self;
			}
			write "D eC: " + local_enzymes.T_cellulolytic * local_step / #gram;
			write "D eN: " + local_enzymes.T_amino * local_step / #gram;
			write "D eP: " + local_enzymes.T_P * local_step / #gram;
			write "D er: " + local_enzymes.T_recal * local_step / #gram;
			
			ask decomposition_problem {
				C_recal_init <- organic.C_recalcitrant;
				N_recal_init <- organic.N_recalcitrant;
				P_recal_init <- organic.P_recalcitrant;
				C_labile_init <- organic.C_labile;
				N_labile_init <- organic.N_labile;
				P_labile_init <- organic.P_labile;
				C_DOM_init <- pore.dam.dom[2];
				N_DOM_init <- pore.dam.dom[0];
				P_DOM_init <- pore.dam.dom[1];
				N_DIM_init <- pore.dam.dim[0];
				P_DIM_init <- pore.dam.dim[1];
				
				do decomposition(local_enzymes, decomposition);
				
				organic.C_recalcitrant <- C_recal_final(decomposition);
				organic.N_recalcitrant <- N_recal_final(decomposition);
				organic.P_recalcitrant <- P_recal_final(decomposition);
				
				organic.C_labile <- C_labile_final(decomposition);
				organic.N_labile <- N_labile_final(decomposition);
				organic.P_labile <- P_labile_final(decomposition);
				
				pore.dam.dom[2] <- C_DOM_final(decomposition);
				pore.dam.dom[0] <- N_DOM_final(decomposition);
				pore.dam.dom[1] <- P_DOM_final(decomposition);
				
				pore.dam.dim[0] <- N_DIM_final(decomposition);
				pore.dam.dim[1] <- P_DIM_final(decomposition);
			}
			
			ask local_enzymes {
				do die;
			}
		}
		ask total_enzymes {
			do die;
		}
		ask decomposition_problem {
			do die;
		}
		ask decomposition {
			do die;
		}
	}

	reflex {
		do microbe_life;
	}
	
	action microbe_life {
		float logistic_factor <- max(0.0, 1 - sum(populations collect each.C) / carrying_capacity);
		
		float total_C_wanted <- sum(populations collect(each.requested_C())) * logistic_factor;
		float total_N_wanted <- sum(populations collect(each.requested_N())) * logistic_factor;
		float total_P_wanted <- sum(populations collect(each.requested_P())) * logistic_factor;
		
		float total_C_consumed <- min([dam.dom[2], total_C_wanted]);
		float total_N_consumed <- min([dam.dom[0], total_N_wanted]);
		float total_P_consumed <- min([dam.dom[1], total_P_wanted]);
		
		float total_bacteria_C <- sum(populations collect (each.C + each.cytosol_C));
		
		ask shuffle(populations) 
		{
			// TODO: optimise calls to resquested_X()
			
			// Proportion of the total_*_consumed that will be consumed by the current microbe population
			float C_rate  <- total_C_wanted > 0 ? (requested_C() / total_C_wanted) : 0;
			float N_rate  <- total_N_wanted > 0 ? (requested_N() / total_N_wanted) : 0;
			float P_rate  <- total_P_wanted > 0 ? (requested_P() / total_P_wanted) : 0;
			
			// If total_X_consumed = total_X_wanted, X_consum = X_wanted for all microbe population
			float assimilated_C <- total_C_consumed * C_rate;
			float assimilated_N <- total_N_consumed * N_rate;
			float assimilated_P <- total_P_consumed * P_rate;
			
//			write "";
//			write "C: " + (requested_C() > 0 ? 100*assimilated_C / requested_C() : -1);
//			write "N: " + (requested_N() > 0 ? 100*assimilated_N / requested_N() : -1);
//			write "P: " + (requested_P() > 0 ? 100*assimilated_P / requested_P() : -1);

			myself.dam.dom[0] <- max([0.0, myself.dam.dom[0] - assimilated_N]);
			myself.dam.dom[1] <- max([0.0, myself.dam.dom[1] - assimilated_P]);
			myself.dam.dom[2] <- max([0.0, myself.dam.dom[2] - assimilated_C]);

			do life(
				myself.dam, myself.accessible_organics,
				assimilated_C, assimilated_N, assimilated_P,
				total_bacteria_C, myself.carrying_capacity
			);
		}	
	}
}

