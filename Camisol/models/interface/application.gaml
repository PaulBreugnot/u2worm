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
		display camisol type:opengl axes: false fullscreen: false show_fps:true {
			// Button actions
			event mouse_move action:mouse_move_buttons;
			event mouse_down action:mouse_down_buttons;
			event mouse_up action:mouse_up_buttons;
			// Plot actions
			event mouse_move action:mouse_move_plots;
			event mouse_down action:mouse_down_plots;
			
			image "sun" position: {0, 0, -0.001} size: {1, 1} file: sun;
			image "landscape" position: {0, 0, -0.001} size: {1, 1} file: landscape;
			// image "plots" position: {0, 0, -0.001} size: {1, 1} file: plots;
			// grid HelpGrid border:#black;
			
			species PlotView aspect:default;
			image "plants" position: {0, 0, 0} size: {1, 1} file: plants;
			species SoilView aspect:default;
			species SeedButtonMenu aspect:default;
			species SeedView aspect:default;
			species FertilizerButtonMenu aspect:default;
			species FertilizerView aspect:default;
			species EpochView aspect:default;
		}
	}
}