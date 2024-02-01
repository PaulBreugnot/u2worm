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


	action decompose {
		float _total_C_labile <- sum(accessible_organics collect each.C_labile);
		float _total_P_labile <- sum(accessible_organics collect each.P_labile);
		float _total_N_labile <- sum(accessible_organics collect each.N_labile);
		
		float _total_C_recal <- sum(accessible_organics collect each.C_recalcitrant);
		float _total_N_recal <- sum(accessible_organics collect each.N_recalcitrant);
		float _total_P_recal <- sum(accessible_organics collect each.P_recalcitrant);
		
		WeightedEnzymes enzymes;
		create WeightedEnzymes with: [
			T_cellulolytic: sum(populations collect (each.C_actif * each.enzymes.T_cellulolytic)),
			T_amino: sum(populations collect (each.C_actif * each.enzymes.T_amino)),
			T_P: sum(populations collect (each.C_actif * each.enzymes.T_P)),
			T_recal: sum(populations collect (each.C_actif * each.enzymes.T_recal))
		] {
			enzymes <- self;
		}
		
		EnzymaticActivity enzymatic_activity;
		create EnzymaticActivityProblem with:[
			dt::local_step,
			total_C_labile::_total_C_labile,
			total_P_labile::_total_P_labile,
			total_N_labile::_total_N_labile,
			total_C_recal::_total_C_recal,
			total_N_recal::_total_N_recal,
			total_P_recal::_total_P_recal,
			C_DOM::dam.dom[2],
			N_DOM::dam.dom[0],
			P_DOM::dam.dom[1],
			P_DIM::dam.dim[1],
			N_DIM::dam.dim[0]
		] {
			create EnzymaticActivity {
				do compute_activity(enzymes, myself);
				enzymatic_activity <- self;
				ask enzymes {
					do die;
				}
			}
			do die;
		}
					
		ask accessible_organics {
			// Amino activity
			float X_C_amino <- 0.0; // Necessarily < to C_labile
			float X_N_amino <- 0.0;
			
			// CP activity
			float X_C_P <- 0.0;
			float X_P_labile_to_dom <- 0.0;
			float X_P_labile_to_dim <- 0.0;
			
			// C activity
			float X_C <- 0.0;
			
			// Recal activity
			
			float X_C_recal <- 0.0;
			float X_N_recal <- 0.0;
			float X_P_recal_to_dim <- 0.0;
			float X_P_recal_to_labile <- 0.0;
			
			if(_total_C_labile > 0.0) {
				// Amino activity
				X_C_amino <- enzymatic_activity.X_C_amino * C_labile / _total_C_labile; // Necessarily < to C_labile
				if(_total_N_labile > 0.0) {
					X_N_amino <-enzymatic_activity.X_N_amino * N_labile / _total_N_labile;
				}
				
				// CP activity
				X_C_P <- enzymatic_activity.X_C_P * C_labile / _total_C_labile;
				if(_total_P_labile > 0.0) {
					X_P_labile_to_dom <- enzymatic_activity.X_P_labile_to_dom * P_labile / _total_P_labile;
					X_P_labile_to_dim <- enzymatic_activity.X_P_labile_to_dim * P_labile / _total_P_labile;	
				}
				
				// C activity
				X_C <- enzymatic_activity.X_C_cellulolytic * C_labile / _total_C_labile;
			}
			
			// Recal activity
			if(_total_P_recal > 0.0) {
				// Phytase action
				X_P_recal_to_dim <- enzymatic_activity.X_P_recal_to_dim * P_recalcitrant / _total_P_recal;
			}
			
			// TODO: check the consistency of the recal decomposition. What if all C in an organic particle and all N in an other?
			if(_total_C_recal > 0.0) {
				X_C_recal <- enzymatic_activity.X_C_recal * C_recalcitrant / _total_C_recal;
				// Note: X_N_recal=0 and X_P_recal=0 anyway if there is no recalcitrant C.
				if(_total_N_recal > 0.0) {
					X_N_recal <- enzymatic_activity.X_N_recal * N_recalcitrant / _total_N_recal;
				}
				if(_total_P_recal > 0.0) {
					X_P_recal_to_dim <- enzymatic_activity.X_P_recal_to_dim * P_recalcitrant / _total_P_recal;
					X_P_recal_to_labile <- enzymatic_activity.X_P_recal_to_labile * P_recalcitrant / _total_P_recal;			
				}
			}
			
			// Commit fluxes
			C_labile <- C_labile + X_C_recal - X_C_amino - X_C_P - X_C;
			N_labile <- N_labile + X_N_recal - X_N_amino;
			P_labile <- P_labile + X_P_recal_to_labile - X_P_labile_to_dom - X_P_labile_to_dim - X_P_recal_to_dim;
			
			C_recalcitrant <- C_recalcitrant - X_C_recal;
			N_recalcitrant <- N_recalcitrant - X_N_recal;
			P_recalcitrant <- P_recalcitrant - X_P_recal_to_labile - X_P_recal_to_dim;
			
			ask myself.dam {
				dim[1] <- dim[1] + X_P_recal_to_dim + X_P_labile_to_dim;
				
				dom[0] <- dom[0] + X_N_amino;
				dom[1] <- dom[1] + X_P_labile_to_dom;
				dom[2] <- dom[2] + X_C_amino + X_C_P + X_C;
			}
		}
		ask enzymatic_activity {
			do die;
		}
	}

	reflex microbe_life {
		float total_C_wanted <- sum(populations collect(each.C_wanted));
		float total_N_wanted <- sum(populations collect(each.N_wanted));
		float total_P_wanted <- sum(populations collect(each.P_wanted));
		
		float total_C_consumed <- min([dam.dom[2], total_C_wanted]);
		float total_N_consumed <- min([dam.dom[0], total_N_wanted]);
		float total_P_consumed <- min([dam.dom[1], total_P_wanted]);
		
		float total_bacteria_C <- sum(populations collect (each.C + each.cytosol_C));
		
		ask shuffle(populations) 
		{
			// Proportion of the total_*_consumed that will be consumed by the current microbe population
			float C_rate  <- total_C_wanted > 0 ? (C_wanted / total_C_wanted) : 0;
			float N_rate  <- total_N_wanted > 0 ? (N_wanted / total_N_wanted) : 0;
			float P_rate  <- total_P_wanted > 0 ? (P_wanted / total_P_wanted) : 0;
			
			// If total_X_consumed = total_X_wanted, X_consum = X_wanted for all microbe population
			float C_consumed <- total_C_consumed * C_rate;
			float N_consumed <- total_N_consumed * N_rate;
			float P_consumed <- total_P_consumed * P_rate;
			
			// Makes the C/N/P available for the microbe population
			perception_C <- C_consumed;
			perception_N <- N_consumed;
			perception_P <- P_consumed;

			myself.dam.dom[0] <- max([0.0, myself.dam.dom[0] - N_consumed]);
			myself.dam.dom[1] <- max([0.0, myself.dam.dom[1] - P_consumed]);
			myself.dam.dom[2] <- max([0.0, myself.dam.dom[2] - C_consumed]);
			
			do life(myself.dam, myself.accessible_organics, total_bacteria_C, myself.carrying_capacity);
		}	
	}
}

