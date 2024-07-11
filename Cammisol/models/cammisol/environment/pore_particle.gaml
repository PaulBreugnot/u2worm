/**
* Name: PoreParticle
* Based on the internal empty template. 
* Author: pbreugno
* Tags: 
*/


model pore_particle

import "../microbes/microbe_population.gaml"

global {
	DecompositionProblem pore_decomposition_problem;
	Decomposition pore_decomposition;
	WeightedEnzymes total_enzymes_in_pore;
	WeightedEnzymes local_enzymes_in_pore;
	
	init {
		create DecompositionProblem with: [
			dt::local_step
			] {
			pore_decomposition_problem <- self;
		}
		create Decomposition {
			pore_decomposition <- self;
		}
		create WeightedEnzymes {
			total_enzymes_in_pore <- self;
		}
		create WeightedEnzymes {
			local_enzymes_in_pore <- self;
		}
	}
}

species PoreParticle schedules:[] {
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


	action decompose {
		float total_C_labile <- sum(accessible_organics collect each.C_labile);
		float total_C_recal <- sum(accessible_organics collect each.C_recalcitrant);
		
		ask total_enzymes_in_pore {
			T_CNP <- sum(myself.populations collect (each.active_C * each.enzymes.T_CNP));
			T_N <- sum(myself.populations collect (each.active_C * each.enzymes.T_N));
			T_P <- sum(myself.populations collect (each.active_C * each.enzymes.T_P));
			T_recal <- sum(myself.populations collect (each.active_C * each.enzymes.T_recal));
		}
		
		ask accessible_organics {
			OrganicParticle organic <- self;
			PoreParticle pore <- myself;
			
			ask local_enzymes_in_pore {
				T_CNP <- total_C_labile > 0.0 ? total_enzymes_in_pore.T_CNP * myself.C_labile / total_C_labile : 0.0;
				T_N <- total_C_labile > 0.0 ? total_enzymes_in_pore.T_N * myself.C_labile / total_C_labile : 0.0;
				T_P <- total_C_labile > 0.0 ? total_enzymes_in_pore.T_P * myself.C_labile / total_C_labile : 0.0;
				T_recal <- total_C_recal > 0.0 ? total_enzymes_in_pore.T_recal * myself.C_recalcitrant / total_C_recal : 0.0;
			}
			
			ask pore_decomposition_problem {
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
				
				do decomposition(local_enzymes_in_pore, pore_decomposition);
				
				organic.C_recalcitrant <- C_recal_final(pore_decomposition);
				organic.N_recalcitrant <- N_recal_final(pore_decomposition);
				organic.P_recalcitrant <- P_recal_final(pore_decomposition);
				
				organic.C_labile <- C_labile_final(pore_decomposition);
				organic.N_labile <- N_labile_final(pore_decomposition);
				organic.P_labile <- P_labile_final(pore_decomposition);
				
				pore.dam.dom[2] <- C_DOM_final(pore_decomposition);
				pore.dam.dom[0] <- N_DOM_final(pore_decomposition);
				pore.dam.dom[1] <- P_DOM_final(pore_decomposition);
				
				pore.dam.dim[0] <- N_DIM_final(pore_decomposition);
				pore.dam.dim[1] <- P_DIM_final(pore_decomposition);
			}
		}
	}
	
	action microbe_life {
		ask populations {
			do update;
		}
		
		float total_microbes_C <- sum(populations collect each.C);
		
		ask shuffle(populations) {
			float perceived_rate <- C / total_microbes_C;
			do life(
				myself.dam, perceived_rate, sum(myself.populations collect each.C), myself.carrying_capacity
			);
		}	
	}
}

