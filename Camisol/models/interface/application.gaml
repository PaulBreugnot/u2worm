/**
* Name: application
* Based on the internal skeleton template. 
* Author: pbreugno
* Tags: 
*/

model application

import "controller.gaml"
import "plot.gaml"

/**
 * Main application. All the features of the interface models are imported in
 * this model and activated in the application experiment.
 */
global {
	
}

experiment application type:gui {
	action _init_ {
		gama.pref_display_slice_number <- 12; /* 128 too slow ! */
		gama.pref_display_show_rotation <- false;
		gama.pref_display_show_errors <- false;
		gama.pref_errors_display <- false;
		gama.pref_errors_stop <- false;
		gama.pref_errors_in_editor <- false;
		gama.pref_display_numkeyscam <- false;
		create simulation;
	}
	
	output {
		// layout #none consoles: false tabs:false toolbars:false controls:false editors: false navigator: false tray: false background: #black;

		
		display camisol type:opengl axes: false fullscreen: true show_fps:true background: #black {
			camera #from_above dynamic:false locked:true distance:min([env_width/2/tan(45/2), env_height/2/tan(45/2)]);
//			graphics "World "{
//				ask world {draw shape color:rgb(10, 0, 0, 0.5) border: #black;}
//			}
			// Button actions
			event mouse_move action:mouse_move_buttons;
			event mouse_down action:mouse_down_buttons;
			event mouse_up action:mouse_up_buttons;
			// Plot actions
			event mouse_move action:mouse_move_plots;
			event mouse_down action:mouse_down_plots;
			
			image "sun" position: {0, 0, -0.001} size: {1, 0.33} file: sun refresh: false;
			image "landscape" position: {0, 0, -0.001} size: {1, 1} file: landscape refresh: false;
			// image "plots" position: {0, 0, -0.001} size: {1, 1} file: plots;
//			grid HelpGrid border:#black;
			
			species PlotView aspect:default;
			// image "plants" position: {0, 0, 0} size: {1, 1} file: plants refresh: false;
			species ButtonBox aspect: default;
			species SoilView aspect:default;
			species SeedButtonMenu aspect:default;
			species SeedView aspect:default;
			species FertilizerButtonMenu aspect:default;
			species FertilizerView aspect:default;
			species EpochView aspect:default;
			species RunButton aspect: default;
		}
	}
}