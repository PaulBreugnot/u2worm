/**
* Name: mockcamisoladapter
* Based on the internal empty template. 
* Author: pbreugno
* Tags: 
*/


model mockcamisoladapter

global {
	init {
		step <- 1#month;
	}
	
	float production(float N_seed, float P_seed, float N_plant, float P_plant, float harvest_index, float N_from_soil, float plot_surface) {
		int i <- rnd(3);
		list<float> intervals <- [0.0, 20.0, 100.0, 500.0, 600.0];
		return rnd(intervals[i], intervals[i+1]);
	}

	action fertilize(
		float solubles,
		float hemicellulose,
		float cellulose,
		float lignine,
		float C_N,
		float C_P,
		float sample_dose
	) {
		// Do nothing
	}
}

experiment Simple type:gui {

	output {
	}
}