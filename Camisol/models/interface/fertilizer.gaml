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
	
	list<string> organic_fertilizers <- [
		"Corne de zebu", "Cendre d'eucalyptus", "Mada Compost", "Vermicompost moyen",
		"Fumier traditionnel Itasy", "Fumier ameliore", "Fientes de volailles"
	];
	list<string> chemical_fertilizers <- [
		"Hyperfos", "Guanomad"
	];
	
	// Note: the regular init function is not used, since the controller model
	// needs the fertilizers but "import fertilizer.gaml" in controller.gaml causes the
	// fertilizer init to be called AFTER the controller init.
	action init_fertilizers {
		write "Building available fertilizers...";
		// Initializes all the available fertilizers
		// Organic fertilizers
		loop i from:0 to:length(organic_fertilizers)-1 {
			add image_file(image_path + definition + "/fertilizers/" + organic_fertilizers[i] + ".png") to:fertilizer_images;
			create OrganicFertilizer number: 1 with:(type: i+1, name: organic_fertilizers[i]) returns: new_fertilizer;
		}
		// Chemical fertilizers
		loop i from:0 to:length(chemical_fertilizers)-1 {
			add image_file(image_path + definition + "/fertilizers/" + chemical_fertilizers[i] + ".png") to:fertilizer_images;
			create ChemicalFertilizer number: 1 with:(type: i+length(organic_fertilizers)+1, name: chemical_fertilizers[i]) returns: new_fertilizer;
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
	int type;
	
	init {
		map<string, float> fertilizer_data <- fertilizers_data[name];
		if(length(fertilizer_data) = 0) {
			error "Fertilizer " + name + " not found in the provided CSV data.";
		} else {
			C_rate <- fertilizer_data[CSV_C_RATE];
			solubles <- fertilizer_data[CSV_QUANTITE_SOLUBLE];
			hemicellulose <- fertilizer_data[CSV_HEMICELLULOSE];
			cellulose <- fertilizer_data[CSV_CELLULOSE];
			lignine <- fertilizer_data[CSV_LIGNINE];
			C_N <- fertilizer_data[CSV_C_N];
			C_P <- fertilizer_data[CSV_C_P];
			sample_dose <- fertilizer_data[CSV_DOSE_ESSAI]#kg/(10000#m2);
		}
	}
	
	float C_rate <- 0.0;
	float solubles <- 0.0;
	float hemicellulose <- 0.0;
	float cellulose <- 0.0;
	float lignine <- 0.0;
	float C_N <- 0.0;
	float C_P <- 0.0;
	float sample_dose <- 0.0;
	
	string to_string {
		return "Fertilizer(" + name + ")";
	}
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
