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
			
			image "landscape" position: {0, 0, 0} size: {1, 1} file: landscape;
			// grid HelpGrid border:#black;
			species Button aspect:button_image;
		}
	}
}