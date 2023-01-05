/**
* Name: harvest
* Based on the internal skeleton template. 
* Author: pbreugno
* Tags: 
*/

model harvest

import "landscape.gaml"
import "seed.gaml"
import "fertilizer.gaml"

/**
 * Defines features used to represents harvests, i.e. quantities of crops
 * collected at each time step on each plot.
 * 
 * On the final interface, the harvest can only be visualized for a single epoch.
 */
global {
	/**
	 * Current time, initialize at 0. One unit of time represents a complete cycle growth.
	 */
	int current_time <- 0;

	/**
	 * Width/height pixel ratio for quantity images.
	 */
	list<float> harvest_images_ratio <- [
		231/137, // Basket
		174/209, // Bag
		364/156  // Barrow
	];
	/**
	 * Scales each quantity image.
	 */
	list<point> harvest_sizes <- [
		{cell_size, cell_size/harvest_images_ratio[0]}, // Basket
		{0.8*cell_size*harvest_images_ratio[1], 0.8*cell_size}, // Bag
		{0.8*cell_size*harvest_images_ratio[2], 0.8*cell_size}  // Barrow
	];
	/**
	 * The locations at which each harvest can be represent for each plot.
	 */
	list<point> harvest_locations <- [
		{0.5*cell_size, 2.8*cell_size}, // Plot 0
		{0.5*cell_size, 3.8*cell_size}, // Plot 1
		{0.5*cell_size, 5.2*cell_size}, // Plot 2
		{0.5*cell_size, 6.5*cell_size}  // Plot 3
	];
	/**
	 * Images used to represent a crop quantity. A quantity can only have the
	 * following possible values:
	 * - 0: a basket of crop
	 * - 1: a bag of crop
	 * - 2: a barrow of crop
	 */
	list<image_file> harvest_images <- [
		image_file("../../images/harvest/basket.png"),
		image_file("../../images/harvest/bag.png"),
		image_file("../../images/harvest/barrow.png")
	];
	
	/**
	 * Currently selected epoch to visualize.
	 */
	EpochView selected_epoch <- nil;
	
	action init_calendar {
		loop i from: 0 to: 5 {
			create Epoch number: 1 with: (time: i) returns: _epochs;
		}
	}
	
	init {
	}
}

/**
 * A basic epoch.
 */
species Epoch {
	/**
	 * Date/time of the epoch.
	 */
	int time;
}

/**
 * Adaptor used to visualize an epoch as an icon.
 */
species EpochView {
	/**
	 * Epoch to visualize.
	 */
	Epoch epoch;
	/**
	 * Epoch icon size.
	 */
	float icon_size <- cell_size;
	
	aspect default {
		// Default colors
		rgb epoch_color <- #white;
		rgb epoch_border <- #black;
		
		if current_time < epoch.time {
			// Epoch not yet simulated
			epoch_color <- epoch_color - rgb(0, 0, 0, 0.5);
			epoch_border <- #grey;
		} else if current_time = epoch.time {
			// Epoch currently simulated
			epoch_color <- #green;
		} else if self = selected_epoch {
			// Passed epoch
			epoch_border <- #blue;
		}
		draw square(icon_size) color: epoch_color border: epoch_border;
		draw string(epoch.time+1) at: location+{0, 0, 0.5} color: epoch_border anchor: #center font: font("Helvetica", 25);
	}
}

/**
 * Harvests are used to represent the simulation history, storing all the
 * parameters used for each plot at each epoch and the quantity of crops harvested.
 */
species Harvest {
	/**
	 * Time at which crops were harvested.
	 */
	int time;
	/**
	 * Plot identifier on which crops were harvested.
	 */
	int plot;
	/**
	 * Crop type harvested.
	 */
	int seed;
	/**
	 * Fertilizer type.
	 */
	int fertilizer;
	/**
	 * Quantity of crops harvested (between 1 and 3)
	 */
	int quantity;
}

/**
 * Adaptor used to visualize the harvest of all plots for a given epoch.
 */
species HarvestView {
	/**
	 * Harvest to visualize.
	 */
	Harvest harvest;
	
	init {
		location <- harvest_locations[harvest.plot];
	}
	
	aspect default {
		draw harvest_images[harvest.quantity]
			size: harvest_sizes[harvest.quantity]
			at: location+{cell_size+harvest_sizes[harvest.quantity].x/2, 0};
		draw fertilizer_images[harvest.fertilizer]
			size: 0.5*cell_size
			at: location + {0.1*cell_size, 0};
		draw seed_images[harvest.seed]
			size: 0.5*cell_size
			at: location + {0.7*cell_size, 0};
	}
}

experiment debug_harvest type: gui {

	output {
		display camisol type:opengl axes: false {
			event mouse_down {
				ask Harvest {
					do die;
				}
				ask HarvestView {
					do die;
				}

				loop i from:0 to:3 {
					create Harvest number:1 with:(plot: i, seed:i, fertilizer: i, quantity: (current_time+i+2) mod 3) returns: harvests;
					create HarvestView number:1 with:(harvest: harvests[0]);
				}
				current_time <- current_time+1;
			}
			image "landscape" position: {0, 0, -0.01} size: {1, 1} file: landscape;
			image "plots" position: {0, 0, -0.01} size: {1, 1} file: plots;
			species HarvestView aspect: default;
			species EpochView aspect: default;
			// grid HelpGrid border:#black;
		}
	}
}
