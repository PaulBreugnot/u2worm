/**
* Name: crop
* Based on the internal skeleton template. 
* Author: pbreugno
* Tags: 
*/

model seed

import "landscape.gaml"
import "../csv_loader.gaml"

/**
 * Defines features used to handle seeds.
 */
global {
	list<string> crops <- [
		"Arachide", "Ble", "Choux", "Coton", "Haricot", "Mais", "Manioc", "Mil",
		"Patate douce", "Pomme de terre", "Riz", "Soja", "Sorgho", "Tomate"
	];
	
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
		loop i from:0 to:length(crops)-1 {
			add image_file(image_path + definition + "/crops/" + crops[i] + ".png") to:seed_images;
			create Seed number: 1 with:(type: i+1, name: crops[i]);
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
	float N_seed <- 0.0;
	float P_seed <- 0.0;
	float N_plant <- 0.0;
	float P_plant <- 0.0;
	float harvest_index <- 0.0;
	float N_from_soil <- 0.0;
	
	init {
		map<string, float> crop_data <- crops_data[name];
		if(length(crop_data) = 0) {
			map<string, float> seed_crop_data <- crops_data[name + "_g"];
			map<string, float> plant_crop_data <- crops_data[name + "_p"];
			if(length(seed_crop_data) = 0 or length(plant_crop_data) = 0) {
				error "Crop " + name + " not found in the provided CSV data.";
			} else {
				// Distinct data for plant and seeds, as with corn or rice
				N_seed <- seed_crop_data[CSV_N]#kg/#ton;
				N_plant <- plant_crop_data[CSV_N]#kg/#ton;
				P_seed <- seed_crop_data[CSV_P]#kg/#ton;
				P_plant <- plant_crop_data[CSV_P]#kg/#ton;
				harvest_index <- seed_crop_data[CSV_HARVEST_INDEX];
				N_from_soil <- seed_crop_data[CSV_N_FROM_SOIL];
			}
		} else {
				N_seed <- crop_data[CSV_N]#kg/#ton;
				N_plant <- N_seed#kg/#ton;
				P_seed <- crop_data[CSV_P]#kg/#ton;
				P_plant <- P_seed#kg/#ton;
				harvest_index <- crop_data[CSV_HARVEST_INDEX];
				N_from_soil <- crop_data[CSV_N_FROM_SOIL];
		}
	}
	
	string to_string {
		return "Seed(" + name + ")";
	}
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
