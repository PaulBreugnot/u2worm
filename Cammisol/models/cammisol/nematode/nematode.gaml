/**
* Name: nematode
* Based on the internal empty template. 
* Author: pbreugno
* Tags: 
*/


model nematode

import "../microbes/microbes.gaml"

global {
	float nematode_CO2_emissions <- 0.0#gram;
}

species Nematode schedules:[]
{
	
	float C_N <- 6.0; // TODO: source?
	float C_P <- 47.0; // TODO: source?
	float carbon_use_efficiency <- 0.16; // TODO: source?
	float rejected_microbe_necromass_rate <- 0.5; // TODO: à vérifier
	
	// quantity of C N P
	float C <- (470#gram * 10^-6); // TODO: source?
	float N <- (470#gram * 10^-6)/C_N;
	float P <- (470#gram * 10^-6)/C_P;
	
	float stomack_C;
	float stomack_N;
	float stomack_P;
	
	PoreParticle current_pore;
	
	bool awake <- true;
	float predation_rate <- 1e5 * 1e-12 #gram/#day; // About 100k bacterias (1e-12#gram) each day
	float threshold_DOC <- 0.1; // If less than 10% of maximum predation rate, sleep
//	float predation_rate <- 100*0.5e-9#gram/#mn;
//	float predation_rate <- 0.00001#gram/#day; // 0.5E-9#gram/#mn;
	
	//un nematode mange entre 10K et 100k bacteries par jour -> moins de 10k il dort (ref article de colman)
	
	action life
	{
		float requested_C <- predation_rate * local_step;
		float perceived_C <- perceive();
		awake <- (perceived_C / requested_C) >= threshold_DOC;
		
		if(awake)
		{
			do move;
			do predation(perceived_C, requested_C);
			do respirate;	
			do anabolize;
		}
	}
	
	float perceive
	{
		return sum(current_pore.populations collect(each.C + each.cytosol_C));
	}
	
	action move
	{
		list<PoreParticle> neighbors_pores <- current_pore.pore_neighbors;
		
		float max_perceive <- sum(current_pore.populations collect (each.C + each.cytosol_C));
		PoreParticle most_attractive <- current_pore;
		ask neighbors_pores
		{
			float perceived_C <- sum(self.populations collect(each.C + each.cytosol_C));
			if(perceived_C > max_perceive) {
				max_perceive <- perceived_C;
				most_attractive <- self;
			}
		}
		current_pore <- most_attractive;
		location <- any_location_in(current_pore);
	}
	
	action respirate
	{
		float C_respiration <- (1-carbon_use_efficiency) * stomack_C;
		stomack_C <- stomack_C  - C_respiration;
		nematode_CO2_emissions <- nematode_CO2_emissions + C_respiration;
	}
	
	action anabolize
	{
		float C_assimilated <- min(stomack_C, stomack_N * C_N, stomack_P * C_P);
		stomack_C <- stomack_C  - C_assimilated; 
		
		float N_assimilated <- C_assimilated/C_N;
		float P_assimilated <- C_assimilated/C_P;
		
		// Assimilated C/N/P is actually rejected to labile as nematode necromass, to ensure nematode biomass stability
		do reject_to_labile(C_assimilated, N_assimilated, P_assimilated);
		
		stomack_N <- stomack_N - N_assimilated;
		stomack_P <- stomack_P - P_assimilated;
		
		// Excess N/P is rejected in DIM
		current_pore.dam.dim[0] <- current_pore.dam.dim[0] + stomack_N;
		current_pore.dam.dim[1] <- current_pore.dam.dim[1] + stomack_P;
		
		stomack_N <- 0.0;
		stomack_P <- 0.0;		
	}
	
	action reject_to_labile(float C_rejected, float N_rejected, float P_rejected)
	{
		ask current_pore.organic_particle {
			C_labile <- C_labile + C_rejected;
			N_labile <- N_labile + N_rejected;
			P_labile <- P_labile + P_rejected;
		}
	}
	
 	action predation(float perceived_C, float requested_C)
 	{
		float total_C_predated <- min([perceived_C, requested_C]);
		// Total C predated is computed from global perceptions, but N and P predated must be computed from the local state of each population.
		float total_N_predated <- 0.0;
		float total_P_predated <- 0.0;
		
		ask current_pore.populations
		{	
			// perceived_C = sum(C + cytocol_C) on all the populations. So, if total_C_predated = perceived_C,
			// all the C is consumed from each population. Else if total_C_predated = requested_C, an equal proportion
			// of total_C_predated/perceived_C is consumed from each population so that the total C consumed is equal to
			// requested_C.
			float total_C_predated_from_this_population <- (total_C_predated / perceived_C) * (C + cytosol_C);
			
			if (total_C_predated_from_this_population > 0.0) {
				float C_predated_from_structural <- C/(C+cytosol_C) * total_C_predated_from_this_population;
				float C_predated_from_cytosol <- total_C_predated_from_this_population - C_predated_from_structural;
				
				float N_predated_from_structural <- C_predated_from_structural * N/C;
				float N_predated_from_cytosol <- C_predated_from_cytosol * cytosol_N / cytosol_C;
				
				float P_predated_from_structural <- C_predated_from_structural * P/C;
				float P_predated_from_cytosol <- C_predated_from_cytosol * cytosol_P / cytosol_C;
				
				C <- C - C_predated_from_structural;
				N <- N - N_predated_from_structural;
				P <- P - P_predated_from_structural;
				
				cytosol_C <- cytosol_C - C_predated_from_cytosol;
				cytosol_N <- cytosol_N - N_predated_from_cytosol;
				cytosol_P <- cytosol_P - P_predated_from_cytosol;
				
				total_N_predated <- total_N_predated + N_predated_from_structural + N_predated_from_cytosol;
				total_P_predated <- total_P_predated + P_predated_from_structural + P_predated_from_cytosol;
			}
		}
		
		float microbe_necromass_C <- total_C_predated*rejected_microbe_necromass_rate;
		float microbe_necromass_N <- total_N_predated*rejected_microbe_necromass_rate;
		float microbe_necromass_P <- total_P_predated*rejected_microbe_necromass_rate;
		
		do reject_to_labile(microbe_necromass_C, microbe_necromass_N, microbe_necromass_P);
		
		stomack_C <- stomack_C + total_C_predated - microbe_necromass_C;
		stomack_N <- stomack_N + total_N_predated - microbe_necromass_N;
		stomack_P <- stomack_P + total_P_predated - microbe_necromass_P;
 	}
 	
 	aspect red_dot {
 		if(awake) {
 			draw circle(0.00005) at: location color: #red;
 		} else {
 			draw circle(0.00005) at: location color: #brown;
 		}
 	}
}
