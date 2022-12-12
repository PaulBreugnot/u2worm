/**
* Name: plot
* Based on the internal empty template. 
* Author: pbreugno
* Tags: 
*/


model plot

import "landscape.gaml"

/* Insert your model definition here */
global {
	init {
		create Plot number:4 returns: new_plots;
		/*
		 * Polygon shapes gather manually from clicks in the debug_plots experiment.
		 */
		new_plots[0].shape <- polygon([
			{49.62765121459961,103.64986419677734,0.0},
			{55.96125030517578,100.33258056640625,0.0},
			{66.51646423339844,90.68199920654297,0.0},
			{164.2290496826172,90.07862091064453,0.0},
			{182.02249145507812,101.23721313476562,0.0},
			{155.78465270996094,103.0467758178711,0.0},
			{49.326107025146484,103.64986419677734,0.0}
			]);
		new_plots[1].shape <- polygon([
			{71.12135314941406,85.66154479980469,0.0},
			{79.50666809082031,72.84355926513672,0.0},
			{113.16856384277344,70.32786560058594,0.0},
			{158.45059204101562,76.4373779296875,0.0},
			{163.72146606445312,86.61988830566406,0.0},
			{155.2181854248047,87.3386001586914,0.0},
			{117.36358642578125,84.703125,0.0},
			{91.84754943847656,86.02091979980469,0.0}
		]);
		new_plots[2].shape <- polygon([
			{79.18619537353516,64.98503875732422,0.0},
			{83.85208129882812,53.604713439941406,0.0},
			{115.26211547851562,51.10100173950195,0.0},
			{149.17575073242188,53.604713439941406,0.0},
			{155.4347686767578,67.37484741210938,0.0},
			{126.5286636352539,64.52976989746094,0.0}
		]);
		new_plots[3].shape <- polygon([
			{84.64569854736328,47.21467208862305,0.0},
			{91.24068450927734,42.133304595947266,0.0},
			{113.29595947265625,40.94395065307617,0.0},
			{129.9456024169922,42.45760726928711,0.0},
			{138.81112670898438,42.99821853637695,0.0},
			{145.94662475585938,48.51198959350586,0.0},
			{123.56702423095703,47.106571197509766,0.0},
			{105.836181640625,47.106571197509766,0.0}		
		]);
		/*
		 * Plant growth images
		 */
		 loop i from: 0 to: 3 {
		 	new_plots[i].images[0] <- image_file("../../images/plots/plot_" + (i+1) + "_plant_1.png");
		 }
	}
	
	Plot current_plot_focus;
	bool seed_is_selected <- false;
	
	action mouse_move_plots {
		list<Plot> plots_under_mouse <- Plot overlapping #user_location;
		if length(plots_under_mouse) = 0 or plots_under_mouse[0] != current_plot_focus {
			if (current_plot_focus != nil) {
				current_plot_focus.current_image <- nil;
				current_plot_focus <- nil;
			}
		}
		if (seed_is_selected) {
			if length(plots_under_mouse) > 0 and plots_under_mouse[0] != current_plot_focus and plots_under_mouse[0].growth_state = 0 {
				current_plot_focus <- plots_under_mouse[0];
				current_plot_focus.current_image <- current_plot_focus.images[0];
			}
		}
	}
}

species Plot {
	list<image_file> images <- [nil, nil, nil, nil];
	image_file current_image;
	int growth_state <- 0 min:0 max:3;
	
	aspect debug {
		draw shape border:#black;
	}
	aspect plot_state {
		if current_image != nil {
			draw current_image at: {env_width/2, env_height/2, 0};
		}
	}
}

experiment debug_plots type:gui {
	output {
		display camisol type:opengl axes: true {
			event mouse_up {
				write #user_location;
				seed_is_selected <- !seed_is_selected;
			}
			event mouse_move action:mouse_move_plots;
			image "plots" position: {0, 0, 0} size: {1, 1} file: plots;
			// species Plot aspect:debug;
			
			species Plot aspect:plot_state;
			// grid HelpGrid border:#black;
		}
	}
}