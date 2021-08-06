/***
* Name: Camisoleadapter
* Author: nicolas
* Description: 
* Tags: Tag1, Tag2, TagN
***/

model Camisoleadapter

import "camisole_V7.gaml"
global
{
	geometry shape <- square(1#cm);
}

experiment Simple type: gui
{
	/*bool ini_new_soil(map<string, float> sc ){
		//soil_caracteristics <- sc;
		return true;
	}*/
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
	action init_new_year(float cl, float cr, float No, float Po ){
		
		int nb_om <- length(organic_particles);
		float OM <- /*000001*/ 1.5#gram/(#cm*#cm)*world.shape.area;
		ask organic_particles{
			C_labile <- C_labile +cl* OM*self.shape.area ; ///nb_om;
			C_recalcitrant <-C_recalcitrant + cr* OM*self.shape.area ; ///nb_om;
			N <- N+ No*OM*self.shape.area;
			P <- P+ Po*OM*self.shape.area;
		//	write "fdf"+C_labile+ " "+C_recalcitrant+" "+N+" "+P;
		
		}
		
		
		return true;
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
	output
	{
		
	}
}

/* Insert your model definition here */

