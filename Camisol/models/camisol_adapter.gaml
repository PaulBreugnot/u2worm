/***
* Name: Camisoleadapter
* Author: Lucas GROSJEAN, Nicolas MARILLEAU
***/

model Camisoleadapter

import "camisol.gaml"

global
{
//	geometry shape <- square(1#cm);
}

experiment Simple type: gui
{	
	float kilo_of_production_adapt <- kilo_of_production;
	
	action extract_N_P(float N, float P)
	{
		float sumN <- sum(Dam collect(each.dim[0]));
		float sumP <- sum(Dam collect(each.dim[1]));
		write "NNX "+sumN+" "+N;
		
		ask Dam{
			dim[0] <- dim[0] - N*(dim[0]/sumN);
			dim[1] <- dim[1] - P*(dim[1]/sumP);
		}
		float sumN2 <- sum(Dam collect(each.dim[0]));
		float sumP2 <- sum(Dam collect(each.dim[1]));
		write "NN "+sumN+" "+sumN2;
		write "PP "+sumP+" "+sumP2;
	}
	
	float get_DIM_N
	{
		return sum(Dam collect(each.dim[0]));
	}
	
	float get_DOM_N
	{
		return sum(Dam collect(each.dom[0]));
	}
	
	float get_DOM_C
	{
		return sum(Dam collect(each.dom[2]));
	}
	
	float get_C_L
	{
		return sum(particle collect(each.C_labile));
	}
	
	float get_C_R
	{
		return sum(particle collect(each.C_recalcitrant));
	}
	
	float get_N
	{
		return sum(particle collect(each.N));
	}
	
	float get_P
	{
		return sum(particle collect(each.P));
	}
	
	float get_area
	{
		return world.shape.area;
	}
	
	float get_DOM_P
	{
		return sum(Dam collect(each.dom[1]));
	}
	
	float get_DIM_P
	{
		return sum(Dam collect(each.dim[1]));
	}
	
	float get_production_total
	{
		return kilo_of_production;
	}
	
	float get_production_without_charette
	{
		return kilo_of_production_adapt;
	}
	
	output
	{
		
	}
}

/* Insert your model definition here */

