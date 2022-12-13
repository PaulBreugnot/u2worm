/**
* Name: plot
* Based on the internal empty template. 
* Author: pbreugno
* Tags: 
*/


model plot

import "landscape.gaml"
import "seed.gaml"

/* Insert your model definition here */
global {
	init {
		create Plot number:4 returns: new_plots;
		create PlotView number: 4 returns: new_plot_views;
		/*
		 * Polygon shapes gather manually from clicks in the debug_plots experiment.
		 */
		new_plot_views[0].shape <- polygon([
			{49.62765121459961,103.64986419677734,0.0},
			{55.96125030517578,100.33258056640625,0.0},
			{66.51646423339844,90.68199920654297,0.0},
			{164.2290496826172,90.07862091064453,0.0},
			{182.02249145507812,101.23721313476562,0.0},
			{155.78465270996094,103.0467758178711,0.0},
			{49.326107025146484,103.64986419677734,0.0}
			]);
		new_plot_views[0].icon_location <- {43.637484782513866,93.46230227181701,0.0};
		
		new_plot_views[1].shape <- polygon([
			{71.12135314941406,85.66154479980469,0.0},
			{79.50666809082031,72.84355926513672,0.0},
			{113.16856384277344,70.32786560058594,0.0},
			{158.45059204101562,76.4373779296875,0.0},
			{163.72146606445312,86.61988830566406,0.0},
			{155.2181854248047,87.3386001586914,0.0},
			{117.36358642578125,84.703125,0.0},
			{91.84754943847656,86.02091979980469,0.0}
		]);
		new_plot_views[1].icon_location <- {62.510271224840764,78.69243296044709,0.0};
		
		new_plot_views[2].shape <- polygon([
			{79.18619537353516,64.98503875732422,0.0},
			{83.85208129882812,53.604713439941406,0.0},
			{115.26211547851562,51.10100173950195,0.0},
			{149.17575073242188,53.604713439941406,0.0},
			{155.4347686767578,67.37484741210938,0.0},
			{126.5286636352539,64.52976989746094,0.0}
		]);
		new_plot_views[2].icon_location <- {69.0747186830414,58.17894706231469,0.0};
		
		new_plot_views[3].shape <- polygon([
			{84.64569854736328,47.21467208862305,0.0},
			{91.24068450927734,42.133304595947266,0.0},
			{113.29595947265625,40.94395065307617,0.0},
			{129.9456024169922,42.45760726928711,0.0},
			{138.81112670898438,42.99821853637695,0.0},
			{145.94662475585938,48.51198959350586,0.0},
			{123.56702423095703,47.106571197509766,0.0},
			{105.836181640625,47.106571197509766,0.0}		
		]);
		new_plot_views[3].icon_location <- {77.68979988923124,41.3576191854663,0.0};
		/*
		 * Plant growth images
		 */
		 loop i from: 0 to: 3 {
		 	new_plots[i].number <- i+1;
		 	new_plot_views[i].plot <- new_plots[i];
		 	new_plot_views[i].images[0] <- image_file("../../images/plots/plot_" + (i+1) + "_plant_1.png");
		 }
	}
	
	PlotView current_plot_focus;
	SeedView selected_seed;
	
	action mouse_move_plots {
		list<PlotView> plots_under_mouse <- PlotView overlapping #user_location;
		if length(plots_under_mouse) = 0 or plots_under_mouse[0] != current_plot_focus {
			if (current_plot_focus != nil) {
				// If the plot has not been planted
				if(current_plot_focus.plot.seed = nil) {
					// Growth state set back to 0
					current_plot_focus.plot.growth_state <- 0;
				}
				current_plot_focus <- nil;
			}
		}
		if (selected_seed != nil) {
			if length(plots_under_mouse) > 0 and plots_under_mouse[0] != current_plot_focus {
				current_plot_focus <- plots_under_mouse[0];
				// Simulates a growth state of 1 for not planted plots, for visual purpose only
				if (current_plot_focus.plot.growth_state = 0) {
					current_plot_focus.plot.growth_state <- 1;		
				}
			}
		}
	}
	action mouse_down_plots {
		if selected_seed != nil and current_plot_focus != nil {
			// Only if the plot is not yet planted (note: the mouse move event sets the growth state to 1 for those plots)
			if current_plot_focus.plot.growth_state = 1 {
				current_plot_focus.plot.seed <- selected_seed.seed;
				write "Seed " + selected_seed.seed.type + " planted to " + current_plot_focus.plot.number;
			}
			
		}
	}
}

species Plot {
	int number;
	Seed seed;
	int growth_state <- 0 min:0 max:3;
	

}

species PlotView {
	Plot plot;
	list<image_file> images <- [nil, nil, nil, nil];
	image_file current_image;
	point icon_location;
	float icon_size <- 0.5*cell_size;
	
	aspect debug {
		draw shape border:#black;
	}
	aspect default {
		if plot.growth_state > 0 {
			draw images[plot.growth_state-1] at: {env_width/2, env_height/2, 0};
			if plot.seed != nil {
				draw seed_images[plot.seed.type-1] at: icon_location size:icon_size;
			}
		}
	}
}

experiment debug_plots type:gui {
	output {
		display camisol type:opengl axes: true {
			// For debug purpose
			event mouse_up {
				write #user_location;
				if(selected_seed = nil) {
					create Seed number:1 with: (type: rnd(19)+1)returns: new_seed;
					create SeedView number:1 with:(seed: new_seed[0]) returns:new_seed_view;
					selected_seed <- new_seed_view[0];
				} else {
					ask selected_seed {
						do die;
					}
					selected_seed <- nil;
				}
			}
			event mouse_move action:mouse_move_plots;
			event mouse_down action:mouse_down_plots;
			image "plots" position: {0, 0, 0} size: {1, 1} file: plots;
			// species Plot aspect:debug;
			
			species PlotView aspect:default;
			// grid HelpGrid border:#black;
		}
	}
}