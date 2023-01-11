/**
* Name: fertilizer
* Based on the internal empty template. 
* Author: pbreugno
* Tags: 
*/


model fertilizer

import "seed.gaml"

/**
 * Defines features used to handle fertilizers.
 */
 
global {
	/**
	 * Images used to represent each fertilizer.
	 */
	list<image_file> fertilizer_images;
	
	// Note: the regular init function is not used, since the controller model
	// needs the fertilizers but "import fertilizer.gaml" in controller.gaml causes the
	// fertilizer init to be called AFTER the controller init.
	action init_fertilizers {
		write "Building available fertilizers...";
		// Initializes all the available fertilizers
		
		// Organic fertilizers
		loop i from:1 to:7 {
			add image_file(image_path + definition + "/fertilizers/fertilizer_" + i + ".png") to:fertilizer_images;
			create OrganicFertilizer number: 1 with:(type: i) returns: new_fertilizer;
		}
		// Chemical fertilizers
		loop i from:8 to:10 {
			add image_file(image_path + definition + "/fertilizers/fertilizer_" + i + ".png") to:fertilizer_images;
			create ChemicalFertilizer number: 1 with:(type: i) returns: new_fertilizer;
		}
	}
}

/**
 * A basic fertilizer.
 */
species Fertilizer {
	/**
	 * Type of fertilizer.
	 */
	int type <- 1 min: 1 max:10;
	// TODO: fertilizer parameters
	
	float solubles <- 0.0;
	float hemicellulose <- 0.0;
	float cellulose <- 0.0;
	float lignine <- 0.0;
	float C_N <- 0.0;
	float C_P <- 0.0;
}

species OrganicFertilizer parent: Fertilizer {
	
}

species ChemicalFertilizer parent: Fertilizer {
	
}

/**
 * Adaptor used to visualize a fertilizer icon.
 */
species FertilizerView {
	/**
	 * Fertilizer to visualize.
	 */
	Fertilizer fertilizer;
	/**
	 * Fertilizer icon size.
	 */
	float icon_size;

	aspect debug {
		draw circle(icon_size/2);
	}
	aspect default {
		draw fertilizer_images[fertilizer.type-1] size:icon_size;
	}
}
