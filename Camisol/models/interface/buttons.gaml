/**
* Name: userinput
* Based on the internal skeleton template. 
* Author: pbreugno
* Tags: 
*/

model buttons

import "plot.gaml"

global {
	// Helper values that can be used to place UI elements
	float num_cell_width <- env_width / cell_size;
	float num_cell_height <- env_height / cell_size;
	
	Button current_button_focus;
	bool disable_up <- false;
	
	action mouse_move_buttons {
		list<Button> buttons_under_mouse <- (SoilButton+SeedButton+SeedButtonMenu+FertilizerButton+FertilizerButtonMenu) overlapping #user_location;
		if length(buttons_under_mouse) = 0 or buttons_under_mouse[0] != current_button_focus {
			if (current_button_focus != nil) {
				ask current_button_focus {
					do mouse_leave;
				}
				current_button_focus <- nil;
			}
		}
		if (selected_seed = nil and selected_fertilizer = nil and selected_soil = nil) {
			if length(buttons_under_mouse) > 0 and buttons_under_mouse[0] != current_button_focus {
				current_button_focus <- buttons_under_mouse[0];
				ask buttons_under_mouse[0] {
					do mouse_enter;
				}
			}
		} else if (selected_seed != nil) {
			selected_seed.location <- #user_location;
		} else if (selected_fertilizer != nil) {
			selected_fertilizer.location <- #user_location;
		} else if (selected_soil != nil) {
			selected_soil.location <- #user_location;
		}
	}
	
	action mouse_down_buttons {
		if current_button_focus != nil and selected_seed = nil and selected_fertilizer = nil and selected_soil = nil {
			ask current_button_focus {
				do click;
				disable_up <- true;
			}
		}
	}
	
	action mouse_up_buttons {
		// Plots should handle selected items in mouse down actions, so they can be killed safely at mouse up
		if (!disable_up) {
			if (selected_seed != nil) {
				ask selected_seed {
					do die;
				}

				selected_seed <- nil;
			} else if (selected_fertilizer != nil) {
				ask selected_fertilizer {
					do die;
				}

				selected_fertilizer <- nil;
			} else if (selected_soil != nil) {
				ask selected_soil {
					do die;
				}

				selected_soil <- nil;
			}
		}
		disable_up <- false;
	}
	
	init {
		// list<point> soil_button_coordinates <- [{1, 3*cell_size}, ]
		loop i from: 0 to: 2 {
			create Soil number:1 with: (color: soil_colors[i]) returns: new_soils;
			create SoilView number:1 with: (soil: new_soils[0], location:{(i+1)*cell_size, 3.5*cell_size}) returns: new_soil_views;
			create SoilButton number: 1 with: (soil_view: new_soil_views[0], location:{(i+1)*cell_size, 3.5*cell_size});
		}
		create SeedButtonMenu number:1 with: (location:{(num_cell_width-4.5)*cell_size, 1*cell_size});
		create FertilizerButtonMenu number:1 with: (location:{(num_cell_width-4.5)*cell_size, 2*cell_size});
		
		ask SeedButtonMenu {
			do click;
		}
	}
}


species Button {
	/**
	 * Shape of the button, used to catch mouse events.
	 */
	geometry shape <- circle(0.5*0.8*cell_size); // 0.5 factor for radius
	action mouse_enter virtual: true;
	action mouse_leave virtual: true;
	action click virtual: true;
}

species SoilButton parent: Button {
	SoilView soil_view;
	
	action mouse_enter {
		soil_view.size <- 0.9*cell_size;
	}
	action mouse_leave {
		soil_view.size <- 0.8*cell_size;
	}
	action click {
		if selected_seed = nil and selected_fertilizer = nil and selected_soil = nil {
			// Creates a new SoilView bound to the the same Soil as the button
			create SoilView number:1 returns:new_soil_icon with:(soil:self.soil_view.soil,size:0.8*cell_size,location:#user_location);
			selected_soil<-new_soil_icon[0];
			write "Selected soil: " + selected_soil.soil.color;
		}
	}
}

species SeedButton parent: Button {
	SeedView seed_view;
	
	action mouse_enter {
		seed_view.size <- 0.9*cell_size;
	}
	action mouse_leave {
		seed_view.size <- 0.8*cell_size;
	}
	action click {
		if selected_seed = nil and selected_fertilizer = nil {
			// Creates a new SeedView bound to the the same Seed as the button
			create SeedView number:1 returns:new_seed_view with:(seed:self.seed_view.seed,size:0.8*cell_size,location:#user_location);
			selected_seed<-new_seed_view[0];
			write "Selected seed: " + selected_seed.seed.type;
		}
	}
}

species FertilizerButton parent: Button {
	FertilizerView fertilizer_view;
	
	action mouse_enter {
		fertilizer_view.size <- 0.9*cell_size;
	}
	action mouse_leave {
		fertilizer_view.size <- 0.8*cell_size;
	}
	action click {
		if selected_seed = nil and selected_fertilizer = nil {
			// Creates a new FertilizerView bound to the the same Fertilizer as the button
			create FertilizerView number:1 returns:new_fertilizer_view with:(fertilizer:self.fertilizer_view.fertilizer,size:0.8*cell_size,location:#user_location);
			selected_fertilizer<-new_fertilizer_view[0];
			write "Selected fertilizer: " + selected_fertilizer.fertilizer.type;
		}
	}
}

species ButtonMenu parent: Button {
	float size <- 0.8*cell_size;
	bool activated <- false;
	list<point> button_coordinates <- build_button_coordinates();
	
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
	
	action show_buttons virtual: true;
	action hide_buttons virtual: true;
	
	action mouse_enter {
		size <- 0.9*cell_size;
	}
	action mouse_leave {
		size <- 0.8*cell_size;
	}
	action click {
		if !activated {
			do show_buttons;
		}
	}
}

species SeedButtonMenu parent: ButtonMenu {
	image_file image <- image_file("../../images/crops/crop.png");
	list<Seed> seeds;
	list<SeedView> seed_views;
	list<Button> buttons;

	init {
		loop i from: 1 to: 19 {
			create Seed number: 1 with:(type: i) returns: new_seed;
			add new_seed[0] to: seeds;
		}
	}
	
	action show_buttons {
		activated <- true;
		ask FertilizerButtonMenu {
			do hide_buttons;
		}
		loop i from: 1 to: 19 {
			create SeedView number: 1 with: (location: button_coordinates[i-1], seed: seeds[i-1], size: 0.8*cell_size) returns: new_seed_view;
			create SeedButton number: 1 with: (location: button_coordinates[i-1], seed_view: new_seed_view[0]) returns: new_button;
			add new_seed_view[0] to: seed_views;
			add new_button[0] to: buttons;
		}
	}
	
	action hide_buttons {
		activated <- false;
		ask seed_views {
			do die;
		}
		ask buttons {
			do die;
		}
		seed_views <- [];
		buttons <- [];
	}
	
	aspect default {
		draw image size:size;
	}
	aspect debug {
		draw circle(0.5*size) color:#red;
	}
}

species FertilizerButtonMenu parent: ButtonMenu {
	image_file image <- image_file("../../images/fertilizers/fertilizer.png");
	list<Fertilizer> fertilizers;
	list<FertilizerView> fertilizer_views;
	list<Button> buttons;
	
	init {
		loop i from: 1 to: 10 {
			create Fertilizer number: 1 with:(type: i) returns: new_fertilizer;
			add new_fertilizer[0] to: fertilizers;
		}
	}
	
	action show_buttons {
		activated <- true;
		ask SeedButtonMenu {
			do hide_buttons;
		}
		// First button block
		loop i from: 1 to: 7 {
			create FertilizerView number: 1 with: (location: button_coordinates[3+i-1], fertilizer: fertilizers[i-1], size: 0.8*cell_size) returns: new_fertilizer_view;
			create FertilizerButton number: 1 with: (location: button_coordinates[3+i-1], fertilizer_view: new_fertilizer_view[0]) returns: new_button;
			add new_fertilizer_view[0] to: fertilizer_views;
			add new_button[0] to: buttons;
		}
		// Second button block
		loop i from: 8 to: 10 {
			create FertilizerView number: 1 with: (location: button_coordinates[5+i-1], fertilizer: fertilizers[i-1], size: 0.8*cell_size) returns: new_fertilizer_view;
			create FertilizerButton number: 1 with: (location: button_coordinates[5+i-1], fertilizer_view: new_fertilizer_view[0]) returns: new_button;
			add new_fertilizer_view[0] to: fertilizer_views;
			add new_button[0] to: buttons;
		}
	}
	
	action hide_buttons {
		activated <- false;
		ask fertilizer_views {
			do die;
		}
		ask buttons {
			do die;
		}
		fertilizer_views <- [];
		buttons <- [];
	}
	
	aspect default {
		draw image size:size;
	}
	aspect debug {
		draw circle(0.5*size) color:#blue;
	}
}


experiment debug_buttons type:gui {	
	output {
		display camisol type:opengl axes: false {
			event mouse_move action:mouse_move_buttons;
			event mouse_down action:mouse_down_buttons;
			event mouse_up action:mouse_up_buttons;
			
			species SoilView aspect:default;
			species SeedButtonMenu aspect:default;
			species SeedView aspect:default;
			species FertilizerButtonMenu aspect:default;
			species FertilizerView aspect:default;
		}
	}
}
