/**
* Name: dam
* Based on the internal empty template. 
* Author: pbreugno
* Tags: 
*/


model dam

import "organic_particle.gaml"

global {
	/**
	 * How many recalcitrant C is obtained from organic particles from a unit of C enzyme.
	 */
	float enzyme_cr_factor <- 8.0;
	/**
	 * How many labile C is obtained from organic particles from a unit of C enzyme.
	 */
	float enzyme_cl_factor <- 8.0;
	/**
	 * How many N/P is obtained from organic particles from a unit of N/P enzyme.
	 */
	 float enzyme_n_p_factor <- 12.0;
}

/**
 * Dissolved Available Matter.
 * 
 * Represents the solution of accessible matter in a pore.
 */
species Dam 
{	
	//solution organique du sol
	//					N	P	C
	list<float> dom <- [0.0,0.0,0.0];
	
	//solution inorganique du sol
	//					N	 P
	list<float> dim <- [0.0,0.0];
	
//	float C_N -> {dom[2] / ((dom at 0) + (dim at 0))};
//	float C_P -> {(dom at 2) / ((dom at 1) + (dim at 1))};
	
	//---------------------- DEM ------------------------ //
	
	//les Enzymes
	//enzyme carbone labile
	float enzyme_Cl <- 0.0;
	//enzyme carbone récalcitrant
	float enzyme_Cr <- 0.0;
	// enzyme Phosphore
	float enzyme_P <- 0.0;
	// enzyme Azone
	float enzyme_N <- 0.0;
	
	
	//les Fixateurs
	
	// enzyme Acide Carboxilique (CA+)
	// float acide_C;
	// H+ pour le ph
	// float ph_h;

	
	float available_N {
		return dom[0] + dim[0];
	}
	
	float available_P {
		return dom[1] + dim[1];
	}
		
	float available_C {
		return dom[2];
	}
	
	action inject_enzim(float e_Cl,float e_Cr,float e_P,float e_N)
	{
		enzyme_Cl <- enzyme_Cl + e_Cl;
		enzyme_Cr <- enzyme_Cr + e_Cr;
		
		enzyme_P <- enzyme_P + e_P;
		enzyme_N <- enzyme_N + e_N;
	}
		
	action decompe(list<OrganicParticle> particles_to_decompose)
	{
		float add_C <- 0.0;
		float add_P <- 0.0;
		float add_N <- 0.0;
		int organic_particles_count <- length(particles_to_decompose);
		
		ask shuffle(particles_to_decompose) 
		{
			/* Decompose C_recalcitrant */
			float qty_enzyme_Cr_basic <- (myself.enzyme_Cr/organic_particles_count);
			float carbone_expected_cr <- (qty_enzyme_Cr_basic * enzyme_cr_factor);
			
			float dec_Cr <- min([C_recalcitrant, carbone_expected_cr]);
			C_recalcitrant  <- C_recalcitrant - dec_Cr;
			
			/* Decompose C_labile */
			float qty_enzyme_Cl_basic <- (myself.enzyme_Cl/organic_particles_count);
			float carbone_expected_cl <- (qty_enzyme_Cl_basic * enzyme_cl_factor);
			
			float dec_Cl <- min([C_labile, carbone_expected_cl]);
			C_labile  <- C_labile - dec_Cl;
			
			// Total C decomposed
			add_C <- add_C + dec_Cl + dec_Cr + qty_enzyme_Cl_basic + qty_enzyme_Cr_basic;
			
			/* Decompose P */
			float qty_enzyme_P <- ((myself.enzyme_P / organic_particles_count) * enzyme_n_p_factor);
			float dec_P <- min([P,qty_enzyme_P]);
			P <- P - dec_P;
			add_P <- add_P + dec_P;
			
			/* Decompose N */
			float qty_enzyme_N <- ((myself.enzyme_N / organic_particles_count) * enzyme_n_p_factor);
			float dec_N <- min([N,qty_enzyme_N]);
			N <- N - dec_N;
			add_N <- add_N + dec_N;
		}
		
		
		dom <- [ dom[0] + add_N, dom[1] + add_P, dom[2] + add_C];
		
		enzyme_Cl <- 0.0;
		enzyme_Cr <- 0.0;
		enzyme_P <- 0.0;
		enzyme_N <- 0.0;
	}
	
}
