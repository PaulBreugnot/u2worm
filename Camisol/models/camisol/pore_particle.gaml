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

	reflex decompe_dam {
		ask dam {
			do decompe(myself.accessible_organics);
		}
	}
	
	reflex microbe_eat {
		float total_C_wanted <- sum(populations collect(each.C_wanted));
		float total_N_wanted <- sum(populations collect(each.N_wanted));
		float total_P_wanted <- sum(populations collect(each.P_wanted));
		
		float total_C_consumed <- min([dam.available_C(), total_C_wanted]);
		float total_N_consumed <- min([dam.available_N(), total_N_wanted]);
		float total_P_consumed <- min([dam.available_P(), total_P_wanted]);
		
		float total_bacteria_C <- sum(populations collect each.C);
		
		ask shuffle(populations) 
		{
			// Proportion of the total_*_consumed that will be consumed by the current microbe population
			float C_rate  <- total_C_wanted > 0 ? (C_wanted / total_C_wanted) : 0;
			float N_rate  <- total_N_wanted > 0 ? (N_wanted / total_N_wanted) : 0;
			float P_rate  <- total_P_wanted > 0 ? (P_wanted / total_P_wanted) : 0;
			
			// If total_*_consumed = total_*_wanted, *_consum = *_wanted for all microbe population
			float C_consumed <- total_C_consumed * C_rate;
			float N_consumed <- total_N_consumed * N_rate;
			float P_consumed <- total_P_consumed * P_rate;
			
			// Makes the C/N/P available for the microbe population
			perception_C <- C_consumed;
			perception_N <- N_consumed;
			perception_P <- P_consumed;
			
			// TODO: check if this is really necessary
			C_consumed <- min([myself.dam.dom[2],C_consumed]);
			P_consumed <- min([myself.dam.dom[1],P_consumed]);
			N_consumed <- min([myself.dam.dom[0],N_consumed]);

			myself.dam.dom[0] <- myself.dam.dom[0] - N_consumed;
			myself.dam.dom[1] <- myself.dam.dom[1] - P_consumed;
			myself.dam.dom[2] <- myself.dam.dom[2] - C_consumed;
			
			do life(total_bacteria_C, myself.carrying_capacity);
		}	
	}
}
