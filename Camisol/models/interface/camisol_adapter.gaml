/**
* Name: camisoladapter
* Based on the internal empty template. 
* Author: pbreugno
* Tags: 
*/


model camisoladapter

import "../camisol.gaml"

/* Insert your model definition here */
experiment Simple type:gui {
	
//	action fertilize(string fertilizer) {
//		ask world {
//			do apport_MO(fertilizer);
//		}
//		// do world.apport_MO(fertilizer);
//	}

	float get_production {
		return kilo_of_production;
	}
	
	output {
	}
}