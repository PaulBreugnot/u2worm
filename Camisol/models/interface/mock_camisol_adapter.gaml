/**
* Name: mockcamisoladapter
* Based on the internal empty template. 
* Author: pbreugno
* Tags: 
*/


model mockcamisoladapter

global {
	float local_step <- 10#m;
	
	float production(float N_seed, float P_seed, float N_plant, float P_plant, float harvest_index, float N_from_soil, float plot_surface) {
		list<float> intervals <- [0.0, 20.0, 100.0, 500.0];
		float mock_quantity;
		int i <- rnd(2);
		mock_quantity <- mock_quantity + rnd(1, 3) * rnd(intervals[i], intervals[i+1]);
		return mock_quantity;
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
