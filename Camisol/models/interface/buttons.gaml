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
	
	action mouse_move_buttons {
		list<Button> buttons_under_mouse <- (SeedButton+SeedButtonMenu) overlapping #user_location;
		if length(buttons_under_mouse) = 0 or buttons_under_mouse[0] != current_button_focus {
			if (current_button_focus != nil) {
				ask current_button_focus {
					do mouse_leave;
				}
				current_button_focus <- nil;
			}
		}
		if (selected_seed = nil) {
			if length(buttons_under_mouse) > 0 and buttons_under_mouse[0] != current_button_focus {
				current_button_focus <- buttons_under_mouse[0];
				ask buttons_under_mouse[0] {
					do mouse_enter;
				}
			}
		} else {
			selected_seed.location <- #user_location;
		}
	}
	
	action mouse_down_buttons {
		if (selected_seed != nil) {
			ask selected_seed {
				do die;
			}
			selected_seed <- nil;
		}
		if current_button_focus != nil {
			ask current_button_focus {
				do click;
			}
		}
	}
	
	init {
		create SeedButtonMenu number:1 with: (location:{(num_cell_width-4.5)*cell_size, 1*cell_size});
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

species SeedButton parent: Button {
	SeedView seed_view;
	
	action mouse_enter {
		seed_view.size <- 0.9*cell_size;
	}
	action mouse_leave {
		seed_view.size <- 0.8*cell_size;
	}
	action click {
		if selected_seed = nil {
			// Creates a new SeedView bound to the the same Seed as the button
			create SeedView number:1 returns:new_seed_view with:(seed:self.seed_view.seed,size:0.8*cell_size,location:#user_location);
			selected_seed<-new_seed_view[0];
			write "Selected seed: " + selected_seed.seed.type;
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
	
	list<Button> buttons;
	
	action show_buttons virtual: true;
	
	action hide_buttons {
		activated <- false;
		ask buttons {
			do die;
		}
		buttons <- [];
	}
	
	action mouse_enter {
		size <- 0.9*cell_size;
	}
	action mouse_leave {
		size <- 0.8*cell_size;
	}
	action click {
		if !activated {
			activated <- true;
			do show_buttons;
		}
	}
}

species SeedButtonMenu parent: ButtonMenu {
	list<Seed> seeds;

	init {
		loop i from: 1 to: 19 {
			create Seed number: 1 with:(type: i) returns: new_seed;
			add new_seed[0] to: seeds;
		}
	}
	
	action show_buttons {
		loop i from: 1 to: 19 {
			create SeedView number: 1 with: (location: button_coordinates[i-1], seed: seeds[i-1], size: 0.8*cell_size) returns: new_seed_view;
			create SeedButton number: 1 with: (location: button_coordinates[i-1], seed_view: new_seed_view[0]) returns: new_button;
			add new_button[0] to: buttons;
		}
	}
	
	aspect debug {
		draw circle(0.5*size) color:#red;
	}
}


experiment debug_buttons type:gui {	
	output {
		display camisol type:opengl axes: false {
			event mouse_move action:mouse_move_buttons;
			event mouse_down action:mouse_down_buttons;
			
			species SeedButtonMenu aspect:debug;
			species SeedView aspect:default;
		}
	}
}
