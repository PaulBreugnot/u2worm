/**
* Name: userinput
* Based on the internal skeleton template. 
* Author: pbreugno
* Tags: 
*/

model userinput

import "landscape.gaml"

global {
	
	list<string> buttons;
	int i_button <- 0;
	/**
	 * Current button coordinates.
	 */
	point button_coordinates <- init_crop_button_coordinates();
	/**
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
	
	point init_crop_button_coordinates {
		float num_cell_width <- env_width / cell_size;
		float num_cell_height <- env_height / cell_size;
		
		return {(num_cell_width-2.5)*cell_size, 0.5*cell_size};
	}
	
	action next_button_coordinates {
		button_coordinates <- button_coordinates + button_moves[i_button mod 3];
		i_button <- i_button + 1;
	}
	
	init {
		loop i from: 1 to: 19 {
			buttons <- buttons + ("crop_" + i + ".png");
		}
		loop button over: buttons {
			create Button number:1 with: (location: button_coordinates, image: image_file("../../images/crops/" + button));
			do next_button_coordinates;
		}
	}
	
	Button current_focus;
	
	action mouse_move_crop_buttons {
		list<Button> buttons_under_mouse <- Button overlapping #user_location;
		if length(buttons_under_mouse) = 0 or buttons_under_mouse[0] != current_focus {
			if (current_focus != nil) {
				current_focus.size <- 0.8*cell_size;
				current_focus <- nil;
			}
		}
		if length(buttons_under_mouse) > 0 and buttons_under_mouse[0] != current_focus {
			current_focus <- buttons_under_mouse[0];
			current_focus.size <- 0.9*cell_size;
		}
	}
}


species Button {
	image_file image;
	float size <- 0.8*cell_size;
	/**
	 * Shape of the button, used to catch mouse events.
	 */
	geometry shape <- circle(0.5*size); // 0.5 factor for radius
	
	aspect debug {
		draw shape color:#white border:#black;
	}
	aspect button_image {
		draw image size:size;
	}
}

experiment userinput type:gui {	
	output {
		display camisol type:opengl axes: false {
			event mouse_move action:mouse_move_crop_buttons;
			
			species Button aspect:button_image;
		}
	}
}
