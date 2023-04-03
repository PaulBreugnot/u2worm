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
	float fertilizer_icon_sep <- 0.8*cell_size;
	
	/**
	 * Current time, initialize at 0. One unit of time represents a complete cycle growth.
	 */
	int current_time <- 0;

	/**
	 * Width/height pixel ratio for quantity images.
	 */
	list<float> harvest_images_ratio <- [
		221.871/131.128, // Basket
		166.684/201.108, // Bag
		349.310/149.796  // Barrow
	];
	/**
	 * Scales each quantity image.
	 */
	list<point> harvest_sizes <- [
		{0.7*cell_size, 0.7*cell_size/harvest_images_ratio[0]}, // Basket
		{0.9*cell_size*harvest_images_ratio[1], 0.9*cell_size}, // Bag
		{0.8*cell_size*harvest_images_ratio[2], 0.8*cell_size}  // Barrow
	];
	/**
	 * The locations at which each harvest can be represent for each plot.
	 */
	list<point> harvest_locations <- [
		{114#m, 96#m}, // Plot 1
		{114#m, 77#m}, // Plot 2
		{114#m, 57#m}, // Plot 3
		{114#m, 43#m}  // Plot 4
	];
	/**
	 * Images used to represent a crop quantity. A quantity can only have the
	 * following possible values:
	 * - 0: a basket of crop
	 * - 1: a bag of crop
	 * - 2: a barrow of crop
	 */
	list<image_file> harvest_images <- [
		image_file(image_path + definition + "/harvest/basket.png"),
		image_file(image_path + definition + "/harvest/bag.png"),
		image_file(image_path + definition + "/harvest/barrow.png")
	];
	
	// TODO: Calibrate this
	// TODO: define as volumes? + volume by mass for each crop?
	list<float> harvest_thresholds <- [
		20.0#kg, 50.0#kg, 200.0#kg
	];
	
	/**
	 * Currently selected epoch to visualize.
	 */
	EpochView selected_epoch <- nil;
	
	action init_calendar {
		loop i from: 0 to: 5 {
			create Epoch number: 1 with: (time: i);
		}
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
	
	image epoch_image;
	/**
	 * Epoch to visualize.
	 */
	Epoch epoch;
	/**
	 * Epoch icon size.
	 */
	float icon_size <- cell_size;
	
	init {
		epoch_image <- image(image_path + definition + "/epochs/epoch_" + (epoch.time+1) + ".png");
	}
	aspect default {
		if(selected_epoch = self) {
			draw square(icon_size) color: rgb(125, 162, 87);
		}
		draw epoch_image size: icon_size;
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
	 * Surface of the plot on which crops were harvested.
	 * This is used to compute the quantity_index from harvest_thresholds.
	 */
	float plot_surface;
	/**
	 * Crop type harvested.
	 */
	Seed seed;
	/**
	 * Fertilizer type.
	 */
	list<Fertilizer> fertilizers;
	/**
	 * Quantity of crops harvested (between 1 and 3)
	 */
	float quantity;
	
	int barrows <- 0;
	int bags <- 0;
	int baskets <- 0;
	
	init {
		if(quantity > 0.0) {
			barrows <- round(quantity/harvest_thresholds[2]);
			if (barrows = 0) {
				bags <- round(quantity/harvest_thresholds[1]);
				if (bags = 0) {
					baskets <- max([1, round(quantity/harvest_thresholds[0])]);
				}
			}
		 	write "Quantity " + quantity + "kg = " + barrows + " barrows, " + bags + " bags, " + baskets + " baskets.";
		}
	}
}

/**
 * Adaptor used to visualize the harvest of all plots for a given epoch.
 */
species HarvestView {
	/**
	 * Harvest to visualize.
	 */
	Harvest harvest;
	
	point seed_icon_location;
	point harvest_icon_location;
	point fertilizer_icon_location;
	
	init {
	}
	
	aspect default {
		float total_length <- 0.0;
		list<image_file> images;
		list<point> sizes;
		list<point> locations;
		point current_location <- harvest_icon_location;
		if(harvest.barrows > 0) {
			total_length <- total_length + harvest_sizes[2].x * (0 + 0.9*(harvest.barrows-1));
			loop i from: 0 to: harvest.barrows-1 {
				add harvest_images[2] to: images;
				add harvest_sizes[2] to: sizes;
				add current_location to: locations;
				current_location <- current_location + point(0.9*harvest_sizes[2].x, 0.0, 0.0);
			}
		}
		if(harvest.bags > 0) {
			total_length <- total_length + harvest_sizes[1].x * (0 + 0.9*(harvest.bags-1));
			loop i from: 0 to: harvest.bags-1 {
				add harvest_images[1] to: images;
				add harvest_sizes[1] to: sizes;
				add current_location to: locations;
				current_location <- current_location + point(0.9*harvest_sizes[1].x, 0.0, 0.0);
			}
		}
		if(harvest.baskets > 0) {
			total_length <- total_length + harvest_sizes[0].x * (0 + 0.9*(harvest.baskets-1));
			loop i from: 0 to: harvest.baskets-1 {
				add harvest_images[0] to: images;
				add harvest_sizes[0] to: sizes;
				add current_location to: locations;
				current_location <- current_location + point(0.9*harvest_sizes[0].x, 0.0, 0.0);
			}
		}
		loop i from: 0 to: length(locations)-1 {
			locations[i] <- locations[i] - point(total_length/2, 0.0, 0.0);
		}
//		write "Images: " + images;
//		write "Locations: " + locations;
//		write "sizes: " + sizes;
		loop i from: 0 to: length(images)-1 {
			draw images[i] at: locations[i] size: sizes[i];
		}
		loop i from: 0 to: length(harvest.fertilizers)-1 {
			draw fertilizer_images[harvest.fertilizers[i].type-1]
				size: 0.9*cell_size
				at: fertilizer_icon_location + {i*fertilizer_icon_sep, 0};
		}
		if(harvest.seed != nil) {
			draw seed_images[harvest.seed.type-1]
				size: 0.9*cell_size
				at: seed_icon_location;
		}
	}
}

experiment debug_harvest type: gui {

	action _init_ {
		create simulation;
		ask simulation {
			do init_seeds();
			do init_fertilizers();
		}
	}
	output {
		display camisol type:opengl axes: false {
			event #mouse_down {
				ask Harvest {
					do die;
				}
				ask HarvestView {
					do die;
				}

				loop i from:0 to:3 {
					create Harvest number:1 with:(plot: i, seed:Seed[i], fertilizers: [Fertilizer[i]], quantity: (current_time+i+2) mod 3) returns: harvests;
					create HarvestView number:1 with:(harvest: harvests[0]);
				}
				current_time <- current_time+1;
			}
			image landscape position: {0, 0, -0.01} size: {1, 1};
			image plots position: {0, 0, -0.01} size: {1, 1};
			species HarvestView aspect: default;
			species EpochView aspect: default;
			// grid HelpGrid border:#black;
		}
	}
}
