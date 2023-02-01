/**
* Name: nematode
* Based on the internal empty template. 
* Author: pbreugno
* Tags: 
*/


model nematode

import "microbes.gaml"

/* Insert your model definition here */

species Nematode
{
	float threshold_DOC <- 0.5; // TODO: à vérifier
	
	float C_N <- 6.0; // TODO: source?
	float C_P <- 47.0; // TODO: source?
	float CO2_efficiency <- 0.84; // TODO: source?
	float t1 <- 0.5; // TODO: à vérifier
	
	// quantity of C N P
	float C <- (470#gram * 10^-6); // TODO: source?
	float N <- (470#gram * 10^-6)/C_N;
	float P <- (470#gram * 10^-6)/C_P;
	
	float stomack_C;
	float stomack_N;
	float stomack_P;
	
	float CO2_production <- 0.0 ;
	PoreParticle current_pore;
	
	bool awake <- true;
//	float predation_rate <- 0.5e-9#gram/#mn;
	float predation_rate <- 0.00001#gram/#day; // 0.5E-9#gram/#mn;
	
	//un nematode mange entre 10K et 100k bacteries par jour -> moins de 10k il dort (ref article de colman)
	float wanted_C <- predation_rate * step;
	
	//perception of the nematode
	float available_C <- 0.0;
	float available_N <- 0.0; 
	float available_P <- 0.0; 
	
	reflex life
	{
		do perceive;
		awake <- (available_C / wanted_C) >= threshold_DOC;
		
		if(awake)
		{
			do move;
			do predating;
			do respirate;	
			do anabolize;
			
		}
	}
	
	action perceive
	{
		available_C <- sum(current_pore.populations collect(each.C + each.cytosol_C));
		available_N <- sum(current_pore.populations collect(each.N + each.cytosol_N));
		available_P <- sum(current_pore.populations collect(each.P + each.cytosol_P));
	}
	
	action move
	{
		list<PoreParticle> neighbors_pores <- current_pore.pore_neighbors;
		
		float max_perceive <- sum(neighbors_pores[0].populations collect(each.C + each.cytosol_C));
		PoreParticle most_attractive <- neighbors_pores[0];
		ask neighbors_pores
		{
			float tmp <- sum(self.populations collect(each.C + each.cytosol_C));
			if(tmp >= max_perceive)
			{
				max_perceive <- tmp;
				most_attractive <- self;
			}
		}
		current_pore <- most_attractive;
		location <- current_pore.location;
	}
	
	action respirate
	{
		float carbone_respiration <- CO2_efficiency * stomack_C;
		stomack_C <- stomack_C  - carbone_respiration; 
		
		CO2_production <- CO2_production + carbone_respiration;
		total_CO2_produced <- total_CO2_produced + carbone_respiration;
	}
	
	action anabolize
	{
		float carbone_assimilation <- (1 - CO2_efficiency) * stomack_C;
		stomack_C <- stomack_C  - carbone_assimilation; 
		
		float qty_N_necromasse <- carbone_assimilation/C_N;
		float qty_P_necromasse <- carbone_assimilation/C_P;
		
		// necromasse nematode
		do reject(carbone_assimilation, qty_N_necromasse, qty_P_necromasse);
		
		stomack_N <- stomack_N - qty_N_necromasse;
		stomack_P <- stomack_P - qty_P_necromasse;
		
		// le reste de l'estomac va dans la dim
		current_pore.dam.dim[0] <- current_pore.dam.dim[0] + stomack_N;
		current_pore.dam.dim[1] <- current_pore.dam.dim[1] + stomack_P;
		
		stomack_N <- 0.0;
		stomack_P <- 0.0;		
	}
	
	action reject(float carbone, float azote, float phosphore)
	{
		ask current_pore.organic_particle {
			C_labile <- C_labile + carbone;
			N <- N + azote;
			P <- P + phosphore;
		}
	}
	
 	action predating
 	{
		float C_to_catch <- min([available_C, wanted_C]);
		
		float total_C <- 0.0;
		float total_N <- 0.0;
		float total_P <- 0.0;
		
		ask shuffle(current_pore.populations)
		{	
			// available_C = sum(C + cytocol_C) on all the population. So, if C_to_catch = available_C,
			// all the C is consumed from each pore. Else if C_to_catch = wanted_C, an equal proportion
			// of wanted_C/available_C is consumed from each pore such as the total C consumed is equal to
			// wanted_C.
			float total_C_catched_from_pore <- C_to_catch / myself.available_C * (C + cytosol_C);
			
			float C_catched_from_colony;
			//élément dans la partie structure
			if(C+cytosol_C = 0){
				C_catched_from_colony <- 0.0; // No C is catched from this pore
			}else{
				// A proportion of C/total_microbe_C is catched from the colony's C deposit
				C_catched_from_colony <- C/(C+cytosol_C) * total_C_catched_from_pore;
			}
			float N_catched_from_colony <- min([N, C_catched_from_colony / C_N]);
			float P_catched_from_colony <- min([P, C_catched_from_colony / C_P]);
			
			//élément dans la partie cytosol
			float C_catched_from_cytosol <- total_C_catched_from_pore - C_catched_from_colony ;
			float N_catched_from_cytosol <- C_catched_from_cytosol / C_N;
			float P_catched_from_cytosol <- C_catched_from_cytosol / C_P;

			C <- max([0.0, C - C_catched_from_colony]);
			N <- max([0.0, N - N_catched_from_colony]);
			P <- max([0.0, P - P_catched_from_colony]);
			
			cytosol_C <- max([0.0, cytosol_C  - C_catched_from_cytosol]);
			cytosol_N <- max([0.0, cytosol_N  - N_catched_from_cytosol]);
			cytosol_P <- max([0.0, cytosol_P  - P_catched_from_cytosol]);
			
			total_C <- total_C + C_catched_from_colony + C_catched_from_cytosol;
			total_N <- total_N + N_catched_from_colony + N_catched_from_cytosol;
			total_P <- total_P + P_catched_from_colony + P_catched_from_cytosol;
		}
		
		do reject(total_C*t1,total_N*t1,total_P*t1);
		
		stomack_C <- stomack_C + ((1-t1) * total_C);
		stomack_N <- stomack_N + ((1-t1) * total_N);
		stomack_P <- stomack_P + ((1-t1) * total_P);
 	}
}
