/**
* Name: plot
* Based on the internal empty template. 
* Author: pbreugno
* Tags: 
*/


model plot

import "soil.gaml"
import "harvest.gaml"

/**
 * Defines the Plot features that allow a parcel to be planted and fertilized.
 */
 
global {
	init {
		create Plot number:4 returns: new_plots;
		create Soil number:4 returns: new_soils;
		create PlotView number: 4 returns: new_plot_views;
		
		/*
		 * Polygon shapes gather manually from clicks in the debug_plots experiment.
		 * Shapes are only used to catch mouse clicks on plots.
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
		new_plot_views[0].seed_icon_location <- {43.637484782513866,93.46230227181701,0.0};
		new_plot_views[0].fertilizer_icon_location <- {184.36069638917638,91.8210529376853,0.0};
		
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
		new_plot_views[1].seed_icon_location <- {62.510271224840764,78.69243296044709,0.0};
		new_plot_views[1].fertilizer_icon_location <- {170.00144871789786,79.92347306323204,0.0};
		
		new_plot_views[2].shape <- polygon([
			{79.18619537353516,64.98503875732422,0.0},
			{83.85208129882812,53.604713439941406,0.0},
			{115.26211547851562,51.10100173950195,0.0},
			{149.17575073242188,53.604713439941406,0.0},
			{155.4347686767578,67.37484741210938,0.0},
			{126.5286636352539,64.52976989746094,0.0}
		]);
		new_plot_views[2].seed_icon_location <- {69.0747186830414,58.17894706231469,0.0};
		new_plot_views[2].fertilizer_icon_location <- {159.74429336008703,58.99936552500823,0.0};
		
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
		new_plot_views[3].seed_icon_location <- {77.68979988923124,41.3576191854663,0.0};
		new_plot_views[3].fertilizer_icon_location <- {150.7185905138058,44.2294962136383,0.0};
		
		/*
		 * Plant growth images
		 */
		 loop i from: 0 to: 3 {
		 	new_plots[i].number <- i+1;
		 	new_plots[i].soil <- new_soils[i];
		 	
		 	new_plot_views[i].plot <- new_plots[i];
		 	new_plot_views[i].growth_images[0] <- image_file("../../images/plots/plot_" + (i+1) + "_plant_1.png");
		 }
	}
	
	/**
	 * Reference to the plot under the cursor. Used to preview seeds and soils.
	 */
	PlotView current_plot_focus;
	/**
	 * Holds a reference to the original soil so that we can get back to it after a soil preview.
	 */
	Soil current_plot_focus_old_soil;
	bool soil_changed;
	
	SoilView selected_soil;
	SeedView selected_seed;
	FertilizerView selected_fertilizer;
	
	action mouse_move_plots {
		list<PlotView> plots_under_mouse <- PlotView overlapping #user_location;
		// No plot under focus or same plot
		if length(plots_under_mouse) = 0 or plots_under_mouse[0] != current_plot_focus {
			if (current_plot_focus != nil) {
				// If the plot has not been planted
				if(current_plot_focus.plot.seed = nil) {
					// Growth state set back to 0
					current_plot_focus.plot.growth_state <- 0;
				}
				// If the soil has not changed after the preview
				if(current_plot_focus_old_soil != nil and !soil_changed) {
					current_plot_focus.plot.soil <- current_plot_focus_old_soil;
					current_plot_focus_old_soil <- nil;
				}
				current_plot_focus <- nil;
			}
		}
		// New plot under focus case
		if length(plots_under_mouse) > 0 and plots_under_mouse[0] != current_plot_focus {
			current_plot_focus <- plots_under_mouse[0];
			// Simulates a growth state of 1 for not planted plots, for visual purpose only
			if (selected_seed != nil and current_plot_focus.plot.growth_state = 0) {
				current_plot_focus.plot.growth_state <- 1;		
			} else if (selected_soil != nil) {
				// Saves current soil
				current_plot_focus_old_soil <- current_plot_focus.plot.soil;
				// Preview the selected soil
				current_plot_focus.plot.soil <- selected_soil.soil;
				// False until the soil is actually changed by a click
				soil_changed <- false;
			}
		}

		
	}
	action mouse_down_plots {
		if current_plot_focus != nil {
			if selected_seed != nil {
				// Only if the plot is not yet planted (note: the mouse move event sets the growth state to 1 for those plots)
				if current_plot_focus.plot.growth_state = 1 {
					current_plot_focus.plot.seed <- selected_seed.seed;
					write "Seed " + selected_seed.seed.type + " planted to " + current_plot_focus.plot.number;
				}
			}
			if selected_fertilizer != nil {
				// Only if the plot has not grown yet
				if current_plot_focus.plot.growth_state <= 1 {
					current_plot_focus.plot.fertilizer <- selected_fertilizer.fertilizer;
					write "Plot " + current_plot_focus.plot.number + " fertilized with " + selected_fertilizer.fertilizer.type;
				}
			}
			if selected_soil != nil {
				// TODO: copy soil parameters (other than color, already handled by the move/preview action) from the selected soil to the plot's soil
				current_plot_focus.plot.soil <- selected_soil.soil;
				write "Soil of plot " + current_plot_focus.plot.number + " set to " + selected_soil.soil.color;
				soil_changed <- true;
			}
		}
	}
}

/**
 * A parcel that can be planted and fertilized.
 */
species Plot {
	/**
	 * Plot identifier.
	 */
	int number;
	/**
	 * Define the soil color and base parameters.
	 */
	Soil soil;
	/**
	 * Seed currently planted (nil if none is planted).
	 */
	Seed seed;
	/**
	 * Fertilizer currently applied (nil if none is applied).
	 */
	Fertilizer fertilizer;
	/**
	 * The growth state:
	 * - 0: no seed
	 * - 1: seed planted
	 * - 2: start to grow
	 * - 3: ready to harvest
	 */
	int growth_state <- 0 min:0 max:3;
	/**
	 * List of harvest performed on this plot so that harvest[epoch.time]
	 * corresponds to the harvest performed at the given epoch.
	 */
	list<Harvest> harvests;
	
	int micro_model_num_steps <- 1000;
		
	action grow {
		float N<-soil.N;
		float P<-soil.P;
		
		if(fertilizer != nil) {
			N <- N + fertilizer.N;
			P <- P + fertilizer.P;
		}
		loop i from: 0 to: 1000 {
			
		}
	}
}

/**
 * Adaptor used to visualise a Plot.
 */
species PlotView {
	/**
	 * Visualised plot.
	 */
	Plot plot;
	/**
	 * Images of the plot for each growth_state level.
	 */
	list<image_file> growth_images <- [nil, nil, nil, nil];
	/**
	 * Location of the seed icon for the current plot.
	 */
	point seed_icon_location;
	/**
	 * Location of the fertilizer icon for the current plot.
	 */
	point fertilizer_icon_location;
	/**
	 * Size of seed and fertilizer icons.
	 */
	float icon_size <- 0.5*cell_size;
	
	/**
	 * Can be used to visualize the plot shape.
	 */
	aspect debug {
		draw shape border:#black;
	}
	aspect default {
		if plot.fertilizer != nil {
			draw fertilizer_images[plot.fertilizer.type-1] at: fertilizer_icon_location size:icon_size;
		}
		if plot.seed != nil {
			draw seed_images[plot.seed.type-1] at: seed_icon_location size:icon_size;
		}
		if plot.soil != nil {
			draw soil_images[plot.soil.color][plot.number-1] at: {env_width/2, env_height/2, 0};
		}
		if plot.growth_state > 0 {
			draw growth_images[plot.growth_state-1] at: {env_width/2, env_height/2, 0};
		}
	}
}

experiment debug_plots type:gui {
	output {
		display camisol type:opengl axes: true {
			// For debug purpose: click outside of the plots to select a random seed
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
			// species Plot aspect:debug;
			
			species PlotView aspect:default;
			// grid HelpGrid border:#black;
		}
	}
}