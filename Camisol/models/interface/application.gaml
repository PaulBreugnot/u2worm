/**
* Name: application
* Based on the internal skeleton template. 
* Author: pbreugno
* Tags: 
*/

model application

import "buttons.gaml"
import "plot.gaml"

global {
	/** Insert the global definitions, variables and actions here */
}

experiment application type:gui {	
	output {
		display camisol type:opengl axes: false fullscreen: false {
			// Button actions
			event mouse_move action:mouse_move_buttons;
			event mouse_down action:mouse_down_buttons;
			// Plot actions
			event mouse_move action:mouse_move_plots;
			event mouse_down action:mouse_down_plots;
			
			image "sun" position: {0, 0, 0} size: {1, 1} file: sun;
			image "landscape" position: {0, 0, 0} size: {1, 1} file: landscape;
			image "plots" position: {0, 0, 0} size: {1, 1} file: plots;
			// grid HelpGrid border:#black;
			
			species PlotView aspect:default;
			species SeedButtonMenu aspect:debug;
			species SeedView aspect:default;
		}
	}
}