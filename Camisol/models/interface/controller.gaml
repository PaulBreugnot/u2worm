/**
* Name: userinput
* Based on the internal skeleton template. 
* Author: pbreugno
* Tags: 
*/

model buttons

import "button.gaml"
import "plot.gaml"

/**
 * The controller model contains all the heavy application mechanic, in order to:
 * - link all species with their view adaptor, and eventually associate a button to them
 * - manage Seeds, Soils, Fertilizers selection and application to Plots
 * - launch Plots micro simulation and manage the simulation Calendar
 * - ...
 */
global {
	// Helper values that can be used to place UI elements
	float num_cell_width <- env_width / cell_size;
	float num_cell_height <- env_height / cell_size;

	/**
	 * Buttons that should be managed by the button model.
	 */
	list<species<Button>> handled_buttons <- [SoilButton, SeedButton, SeedButtonMenu, FertilizerButton, FertilizerButtonMenu, EpochButton, RunButton];

	/**
	 * Mouse move buttons event handle adaptor.
	 */
	action mouse_move_buttons {
		do mouse_move_buttons_list(handled_buttons);
	}
	
	init {
		do init_seeds;
		do init_fertilizers;
		do init_soils;
		do init_calendar;
		
		write "Building buttons...";
		int i <- 0;
		
		create SoilButtonMenu;
		
		create SeedButtonMenu with: (
			location:{(num_cell_width-4.5)*cell_size, 1*cell_size},
			button_size: 0.8*cell_size
		);
		
		create FertilizerButtonMenu with: (
			location:{(num_cell_width-4.5)*cell_size, 2*cell_size},
			button_size: 0.8*cell_size
		);
		
		ask SeedButtonMenu {
			// Shows the seed button at initialization
			do click;
		}
		
		// Init epoch views/buttons
		i <- 0;
		ask Epoch {
			create EpochView number: 1 with: (
				epoch: self,
				location: {(1.5 + 2*i)*cell_size, 8.5*cell_size}
			) returns: new_epoch_view;
			create EpochButton number:1 with: (
				epoch_view: new_epoch_view[0],
				location:new_epoch_view[0].location,
				button_size: cell_size
			);
			i <- i + 1;
		}
		create RunButton number:1 with: (
			button_size: 0.8*cell_size,
			location: {(1.5 + 2*i)*cell_size, 8.5*cell_size}
		);
	}
}

/**
 * Button used for soil selection.
 */
species SoilButton parent: Button {
	SoilView soil_view;
	
	action mouse_enter {
		soil_view.icon_size <- 0.9*cell_size;
	}
	action mouse_leave {
		soil_view.icon_size <- 0.8*cell_size;
	}
	agent click {
		if selected_item = nil {
			// Creates a new SoilView bound to the the same Soil as the button
			create SoilView number:1 returns:new_soil_icon with:(
				soil:self.soil_view.soil,
				icon_size:0.8*cell_size,
				location:#user_location
			);
			selected_soil<-new_soil_icon[0];
			write "Selected soil: " + selected_soil.soil.color;
		}
		return selected_soil;
	}
	
	action post_click {
		ask selected_soil {
			do die;
		}
		selected_soil <- nil;
	}
}

/**
 * Seed selection button.
 */
species SeedButton parent: Button {
	SeedView seed_view;
	
	action mouse_enter {
		seed_view.icon_size <- 0.9*cell_size;
	}
	action mouse_leave {
		seed_view.icon_size <- 0.8*cell_size;
	}
	agent click {
		if selected_item = nil {
			// Creates a new SeedView bound to the the same Seed as the button
			create SeedView number:1 returns:new_seed_view with:(seed:self.seed_view.seed,icon_size:0.8*cell_size,location:#user_location);
			selected_seed<-new_seed_view[0];
			write "Selected seed: " + selected_seed.seed.type;
		}
		return selected_seed;
	}
	action post_click {
		ask selected_seed {
			do die;
		}
		selected_seed <- nil;
	}
}

/**
 * Fertilizer selection button
 */
species FertilizerButton parent: Button {
	FertilizerView fertilizer_view;
	
	action mouse_enter {
		fertilizer_view.icon_size <- 0.9*cell_size;
	}
	action mouse_leave {
		fertilizer_view.icon_size <- 0.8*cell_size;
	}
	agent click {
		if selected_item = nil {
			// Creates a new FertilizerView bound to the the same Fertilizer as the button
			create FertilizerView number:1 returns:new_fertilizer_view with:(
				fertilizer:self.fertilizer_view.fertilizer,
				icon_size:0.8*cell_size,
				location:#user_location
			);
			selected_fertilizer<-new_fertilizer_view[0];
			write "Selected fertilizer: " + selected_fertilizer.fertilizer.type;
		}
		return selected_fertilizer;
	}
	action post_click {
		ask selected_fertilizer {
			do die;
		}
		selected_fertilizer <- nil;
	}
}

/**
 * Epoch selection button. Only passed epochs can be selected.
 */
species EpochButton parent: Button {
	EpochView epoch_view;
	
	init {
		shape<-square(button_size);
	}
	action mouse_enter {
		if epoch_view.epoch.time < current_time {
			epoch_view.icon_size <- 1.1*button_size;
		}
	}
	action mouse_leave {
		if epoch_view.epoch.time < current_time {
			epoch_view.icon_size <- button_size;
		}
	}
	agent click {
		if epoch_view.epoch.time < current_time {
			selected_epoch <- epoch_view;
		}
		return nil;
	}
	action post_click {
		// Nothing to do
	}
}

species RunButton parent: Button {
	image_file button_image <- image_file("../../images/play-green.png");
	float icon_size <- button_size;
	
	init {
		shape<-square(button_size);
	}
	action mouse_enter {
		icon_size <- 1.1*button_size;
	}
	action mouse_leave {
		icon_size <- button_size;
	}
	agent click {
		if current_time < 6 {
			ask Plot parallel: true{
				do grow;
			}
			current_time <- current_time+1;
		}
		return nil;
	}
	action post_click {
		
	}
	
	aspect default {
		draw button_image size: icon_size;
	}
}

species ButtonBackground {
	geometry menu_background;
	float menu_background_opacity <- 0.5;

	action compute_background(list<point> coordinates) {
		loop i from: 0 to: length(coordinates)-1 {
			coordinates[i] <- coordinates[i]*cell_size;
		}

		list<geometry> menu_shape;
		loop i from: 0 to: length(coordinates)-1 {
			point prev_point <- coordinates[(i+length(coordinates)-1) mod length(coordinates)];
			point next_point <- coordinates[(i+1) mod length(coordinates)];
			point point1 <- {coordinates[i].x, coordinates[i].y};
			point point2 <- {coordinates[i].x, coordinates[i].y};
			point center <- {coordinates[i].x, coordinates[i].y};
			float begin_angle;
			float end_angle;
			float radius <- 0.1*cell_size;
			if (prev_point.x = coordinates[i].x) {
				// By construction, next_point.y = coordinates[i].y
				if(prev_point.y < coordinates[i].y) {
					// .
					// 1
					point1 <- point1 - {0, radius};
					if(next_point.x < coordinates[i].x) {
						// 2.
						//  1
						point2 <- point2 - {radius, 0};
						center <- center + {-radius, -radius};
						begin_angle <- 0.0;
						end_angle <- 90.0;
					} else {
						// .2
						// 1
						point2 <- point2 + {radius, 0};
						center <- center + {radius, -radius};
						begin_angle <- 180.0;
						end_angle <- 90.0;
					}
				} else {
					//1
					//.
					point1 <- point1 + {0, radius};
					if(next_point.x < coordinates[i].x) {
						//  1
						// 2.
						point2 <- point2 - {radius, 0};
						center <- center + {-radius, radius};
						begin_angle <- 0.0;
						end_angle <- -90.0;
					} else {
						// 1
						// .2
						point2 <- point2 + {radius, 0};
						center <- center + {radius, radius};
						begin_angle <- 180.0;
						end_angle <- 270.0;
					}
				}
			} else {
				// By construction, prev_point.y = coordinates[i].y
				if(prev_point.x < coordinates[i].x) {
					// 1.
					point1 <- point1 - {radius, 0};
					if(next_point.y < coordinates[i].y) {
						// 1.
						//  2
						point2 <- point2 - {0, radius};
						center <- center + {-radius, -radius};
						begin_angle <- 90.0;
						end_angle <- 360.0;
					} else {
						//  2
						// 1.
						point2 <- point2 + {0, radius};
						center <- center + {-radius, radius};
						begin_angle <- 270.0;
						end_angle <- 360.0;
					}
				} else {
					// .1
					point1 <- point1 + {radius, 0};
					if(next_point.y < coordinates[i].y) {
						// .1
						// 2
						point2 <- point2 - {0, radius};
						center <- center + {radius, -radius};
						begin_angle <- 90.0;
						end_angle <- 180.0;
					} else {
						// 2
						// .1
						point2 <- point2 + {0, radius};
						center <- center + {radius, radius};
						begin_angle <- 270.0;
						end_angle <- 180.0;
					}
				}
			}
			 // add point1 to: menu_shape;
			 // add center to: menu_shape;
			 // add point2 to: menu_shape;
			 write("P: " + coordinates[i]/cell_size + ": " + begin_angle + "->" + end_angle);
			 loop j from: 0 to: 10 {
			 	float angle <- begin_angle + j * (end_angle-begin_angle)/10;
			 	write("    " + angle);
			 	add center + {radius*cos(angle), radius*sin(angle)} to: menu_shape;
			 }
		}
		menu_background <- polygon(menu_shape);
	}
}
/**
 * ButtonMenu interface.
 * 
 * The ButtonMenu is a special button that can create and display other buttons
 * when clicked.
 */
species ButtonMenu parent: Button {
	/**
	 * Size of the menu icon.
	 */
	float icon_size <- 0.8*cell_size;
	/**
	 * True iff the ButtonMenu is currently activated, i.e. it is displaying associated buttons.
	 */
	bool activated <- false;
	/**
	 * Builds the buttons coordinates organized as a grid on the right of the environment.
	 */
	list<point> button_coordinates <- build_button_coordinates();

	ButtonBackground background;

	init {
		create ButtonBackground returns: menu_background;
		background <- menu_background[0];
	}

	list<point> build_button_coordinates {
		int i_button <- 0;
		
		point init_coordinates <- {(num_cell_width-2.5)*cell_size, 1*cell_size};
		
		/*
		 * Defines the vector to add to the current button coordinates to plot buttons in the following order:
		 *  1  0  2
		 *  4  3  5
		 * ... 6 ...
		 */
		list<point> button_moves <- [
			// Go to left...
			{-cell_size, 0#m},
			// then go forward...
			{2*cell_size, 0#m},
			// then go below center
			{-cell_size, cell_size}
		];
		list<point> coordinates <- [init_coordinates];
		loop i from: 1 to: 18 {	
			add coordinates[i-1] + button_moves[i_button mod 3] to: coordinates;
			i_button <- i_button + 1;
		}
		return coordinates;
	}
	
	/**
	 * Creates and shows buttons associated to this menu.
	 */
	action show_buttons virtual: true;
	/**
	 * Hides and kill displayed buttons.
	 */
	action hide_buttons virtual: true;
	
	action mouse_enter {
		if(!activated) {
			icon_size <- 0.9*cell_size;
		}
	}
	action mouse_leave {
		icon_size <- 0.8*cell_size;
	}
	agent click {
		if !activated {
			do show_buttons;
		}
		// Nothing to return as no selection is performed
		return nil;
	}
	action post_click {
		// Nothing to do
	}
}

/**
 * A button menu used to display available Seed types, initialized in the seed model.
 */
species SeedButtonMenu parent: ButtonMenu {
	/**
	 * Image used to represent the menu button.
	 */
	image_file image <- image_file(image_path + definition + "/crops/crop.png");

	
	init {
		ask background {
			do compute_background([
				{12, 0.5}, {16, 0.5}, {16, 6.5}, {15, 6.5}, {15, 7.5},
				{14, 7.5}, {14, 6.5}, {13, 6.5}, {13, 1.5}, {12, 1.5}
			]);
		}
	}
	
	action show_buttons {
		// Activates the current menu
		activated <- true;
		// Hide the fertilizer menu
		ask FertilizerButtonMenu {
			do hide_buttons;
		}
		// Builds the seed menu
		int i <- 0;
		// For all available Seed...
		ask Seed {
			// ... creates a view
			create SeedView number: 1 with: (
				location: myself.button_coordinates[i],
				seed: self,
				icon_size: 0.8*cell_size
			) returns: new_seed_view;
			// ... associates a button to the view.
			create SeedButton number: 1 with: (
				location: myself.button_coordinates[i],
				seed_view: new_seed_view[0],
				button_size: 0.8*cell_size
			);
			i <- i+1;
		}
	}
	
	action hide_buttons {
		// Deactivates the seed menu
		activated <- false;
		// Deletes the button and views of the menu.
		// Notice that the available Seed are NOT deleted.
		ask SeedButton {
			ask seed_view {
				do die;
			}
			do die;
		}
	}
	
	aspect default {
		if(activated) {
			draw background.menu_background color: rgb(#white, background.menu_background_opacity);
		}
		draw image size:icon_size;
	}
	aspect debug {
		draw circle(0.5*icon_size) color:#red;
	}
}

/**
 * A button menu used to display available Fertilizer types, initialized in the fertilizer model.
 */
species FertilizerButtonMenu parent: ButtonMenu {

	image_file image <- image_file(image_path + definition + "/fertilizers/fertilizer.png");
	
	init {
		ask background {
			do compute_background([
				{12, 1.5}, {16, 1.5}, {16, 5.5}, {13, 5.5},
				{13, 2.5}, {12, 2.5}
			]);
		}
	}
	
	action show_buttons {
		// Activates the fertilizer menu
		activated <- true;
		// Hides the seed menu
		ask SeedButtonMenu {
			do hide_buttons;
		}
		int i <- 1;
		// First button block
		ask OrganicFertilizer {
			// Initializes a view/button for all organic fertilizers
			create FertilizerView number: 1 with: (
				location: myself.button_coordinates[3+i-1],
				fertilizer: self,
				icon_size: 0.8*cell_size
			) returns: new_fertilizer_view;
			create FertilizerButton number: 1 with: (
				location: new_fertilizer_view[0].location, // Same location as the associated view
				fertilizer_view: new_fertilizer_view[0],
				button_size: 0.8*cell_size
			) returns: new_button;
			i<-i+1;
		}
		// Second button block
		ask ChemicalFertilizer {
			// Initializes a view/button for all chemical fertilizers
			create FertilizerView number: 1 with: (
				location: myself.button_coordinates[5+i-1],
				fertilizer: self,
				icon_size: 0.8*cell_size
			) returns: new_fertilizer_view;
			create FertilizerButton number: 1 with: (
				location: new_fertilizer_view[0].location, // Same location as the associated view
				fertilizer_view: new_fertilizer_view[0],
				button_size: 0.8*cell_size
			) returns: new_button;
			i<-i+1;
		}
	}
	
	action hide_buttons {
		// Deactivates the seed menu
		activated <- false;
		// Deletes the button and views of the menu.
		// Notice that the available Fertilizers are NOT deleted.
		ask FertilizerButton {
			ask fertilizer_view {
				do die;
			}
			do die;
		}
	}
	
	aspect default {
		if(activated) {
			draw background.menu_background color: rgb(#white, background.menu_background_opacity);
		}
		draw image size:icon_size;
	}
	aspect debug {
		draw circle(0.5*icon_size) color:#blue;
	}
}

species SoilButtonMenu {
	ButtonBackground background;

	init {
		create ButtonBackground returns: button_background;
		background <- button_background[0];
		ask background {
			do compute_background([
				{6, 1.5}, {9, 1.5}, {9, 2.5}, {6, 2.5}
			]);
		}
		// Creates a soil view/button for each available Soil
		int i <- 0;
		ask Soil {
			create SoilView number:1 with: (
				soil: self,
				location:{(i+6.5)*cell_size, 2*cell_size}
			) returns: new_soil_views;
			create SoilButton number: 1 with: (
				soil_view: new_soil_views[0],
				location: new_soil_views[0].location,
				button_size: 0.8*cell_size
			);
			i <- i+1;
		}
	}
	
	aspect default {
		draw background.menu_background color: rgb(#white, background.menu_background_opacity);
	}
}

experiment debug_controller type:gui {	
	output {
		display camisol type:opengl axes: false background: #black {
			event mouse_move action:mouse_move_buttons;
			event mouse_down action:mouse_down_buttons;
			event mouse_up action:mouse_up_buttons;
			
			species SoilButtonMenu aspect: default;
			species SoilView aspect:default;
			species SeedButtonMenu aspect:default;
			species SeedView aspect:default;
			species FertilizerButtonMenu aspect:default;
			species FertilizerView aspect:default;
			species EpochView aspect: default;
			species RunButton aspect: default;
//			grid HelpGrid border:#black;
		}
	}
}
