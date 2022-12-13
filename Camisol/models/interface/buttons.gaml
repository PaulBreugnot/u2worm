/**
* Name: userinput
* Based on the internal skeleton template. 
* Author: pbreugno
* Tags: 
*/

model buttons

import "plot.gaml"

global {	
	
	action next_button_coordinates {

	}
	list<point> button_coordinates {
		int i_button <- 0;
		float num_cell_width <- env_width / cell_size;
		float num_cell_height <- env_height / cell_size;
		
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
	
	init {
		list<point> coordinates <- button_coordinates();
		loop i from: 1 to: 19 {
			create Seed number: 1 with:(type: i) returns: new_seed;
			create SeedView number: 1 with: (location: coordinates[i-1], seed: new_seed[0], size: 0.8*cell_size) returns: new_seed_view;
			create Button number: 1 with: (location: coordinates[i-1], seed_view: new_seed_view[0]);
			do next_button_coordinates;
		}
	}
	
	Button current_button_focus;
	
	action mouse_move_seed_buttons {
		list<Button> buttons_under_mouse <- Button overlapping #user_location;
		if length(buttons_under_mouse) = 0 or buttons_under_mouse[0] != current_button_focus {
			if (current_button_focus != nil) {
				current_button_focus.seed_view.size <- 0.8*cell_size;
				current_button_focus <- nil;
			}
		}
		if (selected_seed = nil) {
			if length(buttons_under_mouse) > 0 and buttons_under_mouse[0] != current_button_focus {
				current_button_focus <- buttons_under_mouse[0];
				current_button_focus.seed_view.size <- 0.9*cell_size;
			}
		} else {
			selected_seed.location <- #user_location;
		}
	}
	
	action mouse_down_seed_buttons {
		if selected_seed = nil {
			ask Button overlapping #user_location {
				// Creates a new SeedView bound to the the same Seed as the button
				create SeedView number:1 returns:new_seed_view with:(seed:self.seed_view.seed,size:0.8*cell_size,location:#user_location);
				selected_seed<-new_seed_view[0];
				write "Selected seed: " + selected_seed.seed.type;
			}
		} else {
			ask selected_seed {
				do die;
			}
			selected_seed <- nil;
		}
	}
}


species Button {
	SeedView seed_view;
	/**
	 * Shape of the button, used to catch mouse events.
	 */
	geometry shape <- circle(0.5*0.8*cell_size); // 0.5 factor for radius
}

experiment debug_buttons type:gui {	
	output {
		display camisol type:opengl axes: false {
			event mouse_move action:mouse_move_seed_buttons;
			event mouse_down action:mouse_down_seed_buttons;
			
			species SeedView aspect:default;
		}
	}
}
