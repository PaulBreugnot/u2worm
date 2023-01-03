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

global {
	int current_time <- 0;

	/** Insert the global definitions, variables and actions here */
	list<float> harvest_images_ratio <- [
		231/137,
		174/209,
		364/156
	];
	list<point> harvest_sizes <- [
		{cell_size, cell_size/harvest_images_ratio[0]},
		{0.8*cell_size*harvest_images_ratio[1], 0.8*cell_size},
		{0.8*cell_size*harvest_images_ratio[2], 0.8*cell_size}
	];
	list<point> harvest_locations <- [
		{0.5*cell_size, 2.8*cell_size},
		{0.5*cell_size, 3.8*cell_size},
		{0.5*cell_size, 5.2*cell_size},
		{0.5*cell_size, 6.5*cell_size}
	];
	list<image_file> harvest_images <- [
		image_file("../../images/harvest/basket.png"),
		image_file("../../images/harvest/bag.png"),
		image_file("../../images/harvest/barrow.png")
	];
	
	EpochView current_epoch;
	EpochView selected_epoch <- nil;
}

species Calendar {
	init {
		create Epoch number: 1 with: (time: 0) returns: _init_epoch;
		create EpochView number: 1 with: (epoch: _init_epoch[0], location: {1.5*cell_size, 8.5*cell_size}) returns: calendar_views;
		current_epoch <- calendar_views[0];
		loop i from: 1 to: 5 {
			create Epoch number: 1 with: (time: i) returns: _epochs;
			create EpochView number: 1 with: (epoch: _epochs[0], location: {(1.5 + 2*i)*cell_size, 8.5*cell_size});
		}
	}
}

species Epoch {
	int time;
}

species EpochView {
	float size <- cell_size;
	Epoch epoch;
	
	aspect default {
		rgb epoch_color <- #white;
		rgb epoch_border <- #black;
		if current_time < epoch.time {
			epoch_color <- epoch_color - rgb(0, 0, 0, 0.5);
			epoch_border <- #grey;
		} else if current_time = epoch.time {
			epoch_color <- #green;
		} else if self = selected_epoch {
			epoch_border <- #blue;
		}
		draw square(size) color: epoch_color border: epoch_border;
		draw string(epoch.time+1) at: location+{0, 0, 0.5} color: epoch_border anchor: #center font: font("Helvetica", 25);
	}
}

species Harvest {
	int epoch;
	int plot;
	int seed;
	int fertilizer;
	int quantity;
}

species HarvestView {
	Harvest harvest;
	bool show <- false;
	
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
				ask Calendar {
					do die;
				}
				if(true) {
					create Calendar number: 1;
				}
				loop i from:0 to:3 {
					create Harvest number:1 with:(plot: i, seed:i, fertilizer: i, quantity: (current_time+i+2) mod 3) returns: harvests;
					create HarvestView number:1 with:(harvest: harvests[0]);
				}
				current_time <- current_time+1;
			}
			image "landscape" position: {0, 0, -0.001} size: {1, 1} file: landscape;
			image "plots" position: {0, 0, -0.001} size: {1, 1} file: plots;
			species HarvestView aspect: default;
			species EpochView aspect: default;
			// grid HelpGrid border:#black;
		}
	}
}
