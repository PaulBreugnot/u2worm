/**
* Name: dam
* Based on the internal empty template. 
* Author: pbreugno
* Tags: 
*/


model dam

import "organic_particle.gaml"
import "../enzymatic_activity/enzymatic_activity.gaml"

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
	
	init {
	}
	
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
	
	string state {
		return name + " [dom:("
			+ dom[2] + ", " + dom[0] + ", " + dom[1]
			+ "), dim:(" + dim[0] + ", " + dim[1] + ")]";
	}
}
