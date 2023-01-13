/**
* Name: plot
* Based on the internal empty template. 
* Author: pbreugno
* Tags: 
*/


model plot

import "harvest.gaml"
import "soil.gaml"

import "camisol_adapter.gaml" as Camisol

/**
 * Defines the Plot features that allow a parcel to be planted and fertilized.
 */
 
global {
	init {
		create Plot number:4 returns: new_plots;
		create Soil number:4 returns: new_soils;
		create PlotView number: 4 returns: new_plot_views;
		
		// Factors to convert coordinates within the SVG file in world coordinates in meters
		float SVG_to_world_x <- env_width/582.850;
		float SVG_to_world_y <- env_height/330.320;
		
		/*
		 * Polygon shapes gather manually from clicks in the debug_plots experiment.
		 * Shapes are only used to catch mouse clicks on plots.
		 */
		new_plot_views[0].shape <- polygon([
			{49.62765121459961#m, 103.64986419677734#m,0.0},
			{55.96125030517578#m, 100.33258056640625#m,0.0},
			{66.51646423339844#m, 90.68199920654297#m,0.0},
			{164.2290496826172#m, 90.07862091064453#m,0.0},
			{182.02249145507812#m, 101.23721313476562#m,0.0},
			{155.78465270996094#m, 103.0467758178711#m,0.0},
			{49.326107025146484#m, 103.64986419677734#m,0.0}
			]);
		// Computes the center coordinate of the plot in world coordinates (meters)
		new_plot_views[0].location <- {
			SVG_to_world_x * (105.002 + (412.479 - 105.002)/2),
			SVG_to_world_y * (194.409 + (228.650 - 194.409)/2)
		};
		// Image size in meter
		new_plot_views[0].plot_image_size <- {
			SVG_to_world_x * (412.479 - 105.002),
			SVG_to_world_y * (228.650 - 194.409)
		};
		new_plot_views[0].seed_icon_location <- {43.637484782513866#m, 93.46230227181701#m,0.0};
		new_plot_views[0].fertilizer_icon_location <- {184.36069638917638#m, 91.8210529376853#m,0.0};
		
		new_plot_views[1].shape <- polygon([
			{71.12135314941406#m, 85.66154479980469#m, 0.0},
			{79.50666809082031#m, 72.84355926513672#m, 0.0},
			{113.16856384277344#m, 70.32786560058594#m, 0.0},
			{158.45059204101562#m, 76.4373779296875#m, 0.0},
			{163.72146606445312#m, 86.61988830566406#m, 0.0},
			{155.2181854248047#m, 87.3386001586914#m, 0.0},
			{117.36358642578125#m, 84.703125#m, 0.0},
			{91.84754943847656#m, 86.02091979980469#m, 0.0}
		]);
		// Computes the center coordinate of the plot in world coordinates (meters)
		new_plot_views[1].location <- {
			SVG_to_world_x * (155.802 + (370.119 - 155.802)/2),
			SVG_to_world_y * (151.179 + (192.585 - 151.179)/2)
		};
		// Image size in meter
		new_plot_views[1].plot_image_size <- {
			SVG_to_world_x * (370.119 - 155.802),
			SVG_to_world_y * (192.585 - 151.179)
		};
		new_plot_views[1].seed_icon_location <- {62.510271224840764#m, 78.69243296044709#m,0.0};
		new_plot_views[1].fertilizer_icon_location <- {170.00144871789786#m, 79.92347306323204#m,0.0};
		
		new_plot_views[2].shape <- polygon([
			{79.18619537353516#m, 64.98503875732422#m, 0.0},
			{83.85208129882812#m, 53.604713439941406#m, 0.0},
			{115.26211547851562#m, 51.10100173950195#m, 0.0},
			{149.17575073242188#m, 53.604713439941406#m, 0.0},
			{155.4347686767578#m, 67.37484741210938#m, 0.0},
			{126.5286636352539#m, 64.52976989746094#m, 0.0}
		]);
		// Computes the center coordinate of the plot in world coordinates (meters)
		new_plot_views[2].location <- {
			SVG_to_world_x * (174.287 + (350.746 - 174.287)/2),
			SVG_to_world_y * (107.657 + (150.276 - 107.657)/2)
		};
		// Image size in meter
		new_plot_views[2].plot_image_size <- {
			SVG_to_world_x * (350.746 - 174.287),
			SVG_to_world_y * (150.276 - 107.657)
		};
		new_plot_views[2].seed_icon_location <- {69.0747186830414#m, 58.17894706231469#m, 0.0};
		new_plot_views[2].fertilizer_icon_location <- {159.74429336008703#m, 58.99936552500823#m, 0.0};
		
		new_plot_views[3].shape <- polygon([
			{84.64569854736328#m, 47.21467208862305#m, 0.0},
			{91.24068450927734#m, 42.133304595947266#m, 0.0},
			{113.29595947265625#m, 40.94395065307617#m, 0.0},
			{129.9456024169922#m, 42.45760726928711#m, 0.0},
			{138.81112670898438#m, 42.99821853637695#m, 0.0},
			{145.94662475585938#m, 48.51198959350586#m, 0.0},
			{123.56702423095703#m, 47.106571197509766#m, 0.0},
			{105.836181640625#m, 47.106571197509766#m, 0.0}		
		]);
		// Computes the center coordinate of the plot in world coordinates (meters)
		new_plot_views[3].location <- {
			SVG_to_world_x * (187.359 + (329.928 - 187.359)/2),
			SVG_to_world_y * (87.847 + (107.802 - 87.847)/2)
		};
		// Image size in meter
		new_plot_views[3].plot_image_size <- {
			SVG_to_world_x * (329.928 - 187.359),
			SVG_to_world_y * (107.802 - 87.847)
		};
		new_plot_views[3].seed_icon_location <- {77.68979988923124#m, 41.3576191854663#m, 0.0};
		new_plot_views[3].fertilizer_icon_location <- {150.7185905138058#m, 44.2294962136383#m, 0.0};
		
		create Camisol.Simple number: 4;
		/*
		 * Plant growth images
		 */
		 loop i from: 0 to: 3 {
		 	new_plots[i].number <- i+1;
		 	new_plots[i].soil <- new_soils[i];	 	
		 	new_plot_views[i].plot <- new_plots[i];
		 	new_plot_views[i].growth_images[0] <- image_file(image_path + definition + "/plots/plot_" + (i+1) + "_plant_1.png");
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

	reflex {
		ask Camisol.Simple[0] parallel:true {
			ask simulation {
				do _step_;
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
	int number <- 1 min:1 max: 4;
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


//	action grow {
//		write "Starting camisol on plot " + number;
//		ask Camisol.Simple[number-1] {
//			ask simulation {
//				loop i from: 0 to: 50 {
//					write i;
//					do _step_;
//					ask world {
//						cycle <- i;
//					}
//				}
//				write "Done: " + kilo_of_production;
//			}
//		}
//	}
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
	point plot_image_size;
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
			// draw soil_images[plot.soil.color][plot.number-1] at: {env_width/2, env_height/2, 0};
			draw soil_images[plot.soil.color][plot.number-1] at: location size:plot_image_size;
		}
//		if plot.growth_state > 0 {
//			draw growth_images[plot.growth_state-1] at: {env_width/2, env_height/2, 0};
//		}
	}
}

experiment debug_plots type:gui {
	
	action _init_ {
		create simulation;
		ask simulation {
			do init_soils();
		}
	}
	
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
			
			species PlotView aspect:default;
			// grid HelpGrid border:#black;
		}
	}
}
