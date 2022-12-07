/**
* Name: application
* Based on the internal skeleton template. 
* Author: pbreugno
* Tags: 
*/

model application

import "user_input.gaml"

global {
	/** Insert the global definitions, variables and actions here */
}

experiment application type:gui {	
	output {
		display camisol type:opengl axes: false fullscreen: true {
			event mouse_move action:mouse_move_crop_buttons;
			event mouse_down action:mouse_down_crop_buttons;
			event mouse_up action:mouse_up_crop_buttons;
			
			image "sun" position: {0, 0, 0} size: {1, 1} file: sun;
			image "landscape" position: {0, 0, 0} size: {1, 1} file: landscape;
			image "plots" position: {0, 0, 0} size: {1, 1} file: plots;
			// grid HelpGrid border:#black;
			
			species Button aspect:button_image;
			species Seed aspect:crop_image;
		}
	}
}