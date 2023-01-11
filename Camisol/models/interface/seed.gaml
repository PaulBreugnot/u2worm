/**
* Name: crop
* Based on the internal skeleton template. 
* Author: pbreugno
* Tags: 
*/

model seed

import "landscape.gaml"

/**
 * Defines features used to handle seeds.
 */
global {
	/**
	 * Images used to represent each crop.
	 */
	list<image_file> seed_images;
	
	// Note: the regular init function is not used, since the controller model
	// needs the seeds but "import seed.gaml" in controller.gaml causes the
	// seed init to be called AFTER the controller init.
	action init_seeds {
		write "Building available seeds...";
		// Initializes all the available seeds
		loop i from:1 to:19 {
			add image_file(image_path + definition + "/crops/crop_" + i + ".png") to:seed_images;
			create Seed number: 1 with:(type: i);
		}
	}
}

/**
 * A basic seed.
 */
species Seed {
	/**
	 * Type of the seed.
	 */
	int type <- 1 min: 1 max:19;

	// TODO: parameters of each seed.
	float C_N_seed <- 0.0;
	float C_P_seed <- 0.0;
	float C_N_plant <- 0.0;
	float C_P_plant <- 0.0;
	float harvest_index <- 0.0;
}

/**
 * Adapter used to visualize seeds using the corresponding icon.
 */
species SeedView {
	/**
	 * Seed to visualize.
	 */
	Seed seed;
	/**
	 * Seed icon size.
	 */
	float icon_size;

	aspect debug {
		draw circle(icon_size/2);
	}
	aspect default {
		draw seed_images[seed.type-1] size:icon_size;
	}
}
