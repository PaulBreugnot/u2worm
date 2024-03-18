/**
* Name: PoreParticle
* Based on the internal empty template. 
* Author: pbreugno
* Tags: 
*/


model pore_particle

import "../microbes/microbe_population.gaml"


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
		
		WeightedEnzymes total_enzymes;
		create WeightedEnzymes with: [
			T_cellulolytic: sum(populations collect (each.active_C * each.enzymes.T_cellulolytic)),
			T_amino: sum(populations collect (each.active_C * each.enzymes.T_amino)),
			T_P: sum(populations collect (each.active_C * each.enzymes.T_P)),
			T_recal: sum(populations collect (each.active_C * each.enzymes.T_recal))
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
				T_cellulolytic: total_C_labile > 0.0 ? total_enzymes.T_cellulolytic * C_labile / total_C_labile : 0.0,
				T_amino: total_C_labile > 0.0 ? total_enzymes.T_amino * C_labile / total_C_labile : 0.0,
				T_P: total_C_labile > 0.0 ? total_enzymes.T_P * C_labile / total_C_labile : 0.0,
				T_recal: total_C_recal > 0.0 ? total_enzymes.T_recal * C_recalcitrant / total_C_recal : 0.0
			] {
				local_enzymes <- self;
			}
			
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
	
	action microbe_life {
		ask populations {
			do update;
		}
		
		// TODO: Allocate C/N/P depending on the size of the population, not on the requested nutrients
//		float total_requested_C <- sum(populations collect max(0.0, each.requested_C-each.cytosol_C));
//		float total_requested_N_from_dom <- sum(populations collect max(0.0, each.requested_N-each.cytosol_N));
//		float total_requested_P_from_dom <- sum(populations collect max(0.0, each.requested_P-each.cytosol_P));
		
//		float total_C_consumed <- min(dam.dom[2], total_requested_C);
//		float total_N_consumed_from_dom <- min(dam.dom[0], total_requested_N_from_dom);
//		float total_P_consumed_from_dom <- min(dam.dom[1], total_requested_P_from_dom);
		
		float total_microbes_C <- sum(populations collect each.C);
		
		ask shuffle(populations) {
			float perceived_rate <- C / total_microbes_C;
			do dormancy(myself.dam, perceived_rate);
			do assimilate(myself.dam, perceived_rate);
			do life(
				myself.dam, myself.accessible_organics, sum(myself.populations collect each.C), myself.carrying_capacity
			);
		}	
	}
}

