/**
* Name: microbe_population
* Based on the internal skeleton template. 
* Author: Nicolas Marilleau, Lucas Grosjean, Paul Breugnot
* Tags: 
*/

model microbe_population

import "dam.gaml"

global {
	// 5E8 -> 5E9 bacterie / gramme de sol
	/*
	 * A modifier................................................ 
	 */
	 // TODO: poids total de bactérie par gramme de sol * surface du modèle
	/**
	 * Total bacteria weight in the model.
	 */
	float total_initial_bacteria_weight <-  0.05*1.5#gram/(#cm*#cm)*world.shape.area;
	
	float total_CO2_produced <- 0.0;
}

species MicrobePopulation
{
	string bacteria_name;
	float taux_respiration <- 0.5;
	float dividing_time <- 0.0;
	//taux labile récalcitrant
	float L_R_rate; // labile_recalcitrante_rate
	float C_N;
	float C_P;
	
	// quantity of C N P in the colony
	// TODO: unit?
	// TODO: source?
	float C <- 60.0;
	float N <- 7.0;
	float P <- 1.0;
	
	float awake_population <- 0.5;
	float wakeup_factor <- 0.50;
	float C_actif -> {C * awake_population};
	
	// Mon estomac le cytosol 
	float cytosol_C <- 0#gram; 
	float cytosol_N <- 0#gram; 
	float cytosol_P <- 0#gram; 
	
	// ce que je percois
	float perception_C <- 0#gram; 
	float perception_N <- 0#gram;
	float perception_P <- 0#gram;

	/* 
	 R : C * (step/1)  /  (1-0.2)
	 K : C * (step/24)  / (1-0.4)
	 O : C * (step/368) / (1-0.7)s
	 */
	// Quantité de carbone voulue par la colonie
	float C_assimilation_speed -> {C_actif * (step / dividing_time) * (1 - taux_respiration)};
	
	
	float CO2_produced <- 0.0;
	
	float C_wanted -> { C_assimilation_speed };
	float N_wanted -> { (C+cytosol_C+C_wanted) / C_N - (cytosol_N + N)};
	float P_wanted -> { (C+cytosol_C+C_wanted) / C_P - (cytosol_P + P)};
	
	float division_enzyme_rate <- 0.0;
	float cytosol_enzyme_facteur <- 0.0;
	
	Dam dam;
	
	action respirate(float cytosol_respiration) // pas de perte c
	{
		cytosol_C <- cytosol_C  - cytosol_respiration; 
		
		CO2_produced <- CO2_produced + cytosol_respiration;
		total_CO2_produced <- total_CO2_produced + cytosol_respiration;
	}
	
	action growth(float cytosol_division, float total_C_in_pore, float pore_carrying_capacity) // pas de perte c
	{
		// quantity de carbone cytosol -> structurel
		float new_individual <- C * (step / dividing_time) * (1 - total_C_in_pore / pore_carrying_capacity);
		
		//new_individual <- min([new_individual , cytosol_C]);
		new_individual <- min([new_individual , cytosol_division]);
		
		
		float assi_C_N <- cytosol_N * C_N;  //(assimilated_N=0) ? 0:(assimilated_C/assimilated_N);
		float assi_C_P <- cytosol_P * C_P; //(assimilated_P=0) ? 0:(assimilated_C/assimilated_P);
		
		new_individual <- min([new_individual, (assi_C_N) ]);
		new_individual <- min([new_individual, (assi_C_P) ]);
		
		// write bacteria_name + ": " + new_individual;
		
		float transfert_C <- new_individual;
		float transfert_N <- transfert_C / C_N;
		float transfert_P <- transfert_C / C_P;
		
		cytosol_C <- cytosol_C - transfert_C;
		cytosol_N <- cytosol_N - transfert_N;
		cytosol_P <- cytosol_P - transfert_P;
		
		
		C <- C + transfert_C;
		N <- N + transfert_N;
		P <- P + transfert_P;
	}
	
	action decompose(float cytosol_enzyme) // perte c
	{
		cytosol_C <- cytosol_C - cytosol_enzyme;
		
		float wanted_P <- cytosol_enzyme / C_P;
		float wanted_N <- cytosol_enzyme / C_N;
		
		ask dam
		{
			do inject_enzim(myself.L_R_rate*cytosol_enzyme, (1-myself.L_R_rate)*cytosol_enzyme, wanted_P, wanted_N);
		}
	}
	
	action life(float total_C_in_pore, float pore_carrying_capacity)
	{
		CO2_produced <- 0.0;
		
		if(C_wanted = 0) {
			awake_population <- 0.0;
		} else {
			awake_population <- min([(perception_C / C_wanted), 1]);
		}
		
		awake_population <- max([awake_population,wakeup_factor]);
		
	
		cytosol_C <- cytosol_C + perception_C * awake_population;
		cytosol_N <- cytosol_N + perception_N * awake_population;
		cytosol_P <- cytosol_P + perception_P * awake_population;
		
		// Breathed carbon quantity from cytosol, rejected as CO2
		float cytosol_respiration <- (taux_respiration) * cytosol_C ;
		// Part of the left cytosol C dedicated to cell division
		float cytosol_division <- (cytosol_C-cytosol_respiration) * division_enzyme_rate;
		// Left cytosol C is used to produce enzymes
		float cytosol_enzyme <- cytosol_C - cytosol_respiration - cytosol_division;
		
		cytosol_enzyme_facteur <- cytosol_enzyme;
		
		
		do respirate(cytosol_respiration);
		do growth(cytosol_division, total_C_in_pore, pore_carrying_capacity);
		do decompose(cytosol_enzyme);
		
	
		ask dam
		{
			// TODO: explanations
			dom[2]<- dom[2] + myself.perception_C * (1 - myself.awake_population);
			dom[1]<- dom[1] + myself.perception_P * (1 - myself.awake_population);
			dom[0]<- dom[0] + myself.perception_N * (1 - myself.awake_population);
		}
		
		perception_C <- 0.0;
		perception_N <- 0.0;
		perception_P <- 0.0;
	}
}
