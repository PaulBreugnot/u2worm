/**
* Name: plot
* Based on the internal empty template. 
* Author: pbreugno
* Tags: 
*/


model plot

import "button.gaml"

//import "mock_camisol_adapter.gaml" as Camisol
import "camisol_adapter.gaml" as Camisol

/**
 * Defines the Plot features that allow a parcel to be planted and fertilized.
 */
 
global {
	float epoch_duration <- 6#month;
//	float epoch_duration <- 100#h; // For interface test purpose only

	date output_time <- #now;
	string output_file_suffix <- "" + output_time.year + "" +output_time.month + "" + output_time.day + "" + output_time.hour + "" + output_time.minute;
	float fertilizer_icon_sep <- 0.8*cell_size;
	
	action init_plots {
		// Factors to convert coordinates within the SVG file in world coordinates in meters
		float SVG_to_world_x <- env_width/582.850;
		float SVG_to_world_y <- env_height/330.320;

		create Plot with: (number: 1, soil: default_soil) {
			create PlotView with: (
				plot: self,
				// Coordinates gathered manually from debug_plot
				shape: polygon([
					{63.568172454833984#m,89.06570434570312#m,0.0},
					{87.54286193847656#m,89.22566986083984#m,0.0},
					{145.56124877929688#m,89.5452880859375#m,0.0},
					{162.02391052246094#m,88.90589141845703#m,0.0},
					{180.08412170410156#m,100.41363525390625#m,0.0},
					{153.07284545898438#m,103.29051971435547#m,0.0},
					{46.46635437011719#m,103.76995086669922#m,0.0},
					{52.54005432128906#m,99.93405151367188#m,0.0},
					{57.494747161865234#m,92.90177917480469#m,0.0}
				]),
				// Computes the center coordinate of the plot in world coordinates (meters)
				plot_image_location: {
					SVG_to_world_x * (105.002 + (412.479 - 105.002)/2),
					SVG_to_world_y * (194.409 + (228.650 - 194.409)/2)
				},
				// Image size in meter
				plot_image_size: {
					SVG_to_world_x * (412.479 - 105.002),
					SVG_to_world_y * (228.650 - 194.409)
				},
				selected_plot_image_size: {
					SVG_to_world_x * 328.423,
					SVG_to_world_y * 55.188
				},
				growth_images_locations: [
					{SVG_to_world_x * (124.060 + 276.150/2), SVG_to_world_y * (198.230 + 27/2)}, // Plant 1 (seeds)
					{SVG_to_world_x * (120.790 + 282.451/2), SVG_to_world_y * (196.318 + 28.612/2)}, // Plant 2
					{SVG_to_world_x * (118.849 + 289.851/2), SVG_to_world_y * (193.130 + 30.881/2)} // Plant 3
				],
				growth_images_sizes: [
					{SVG_to_world_x * 276.150, SVG_to_world_y * 27}, // Plant 1 (seeds)
					{SVG_to_world_x * 282.451, SVG_to_world_y * 28.612}, // Plant 2
					{SVG_to_world_x * 289.851, SVG_to_world_y * 30.881} // Plant 2
				],
				seed_icon_location: {43.637484782513866#m, 93.46230227181701#m,0.0},
				fertilizer_icon_location: {184.36069638917638#m, 91.8210529376853#m,0.0},
				growth_images: [],
				selected_image: image_file(image_path + definition + "/plots/plot_1_selected.png")
			) {
				loop i from: 1 to: 3 {
					add image_file(image_path + definition + "/plots/plot_" + myself.number + "_plant_" + i + ".png") to: growth_images;
				}
			}
		}
		
		create Plot with: (number: 2, soil: default_soil) {
			create PlotView with: (
				plot: self,
				// Coordinates gathered manually from debug_plot
				shape: polygon([
					{68.88167572021484#m,85.81029510498047#m,0.0},
					{77.08274841308594#m,72.2039566040039#m,0.0},
					{100.75379943847656#m,69.40814208984375#m,0.0},
					{118.83322143554688#m,69.2217788696289#m,0.0},
					{142.31826782226562#m,71.83123016357422#m,0.0},
					{155.9242401123047#m,75.5588607788086#m,0.0},
					{161.52137192122848#m,87.30346384544302#m,0.0},
					{157.04278564453125#m,86.369384765625#m,0.0},
					{152.9420623779297#m,87.30138397216797#m,0.0},
					{120.13812255859375#m,84.8782958984375#m,0.0},
					{89.94327545166016#m,86.369384765625#m,0.0},
					{68.6953125#m,85.99665832519531#m,0.0}
				]),
				// Computes the center coordinate of the plot in world coordinates (meters)
				plot_image_location: {
					SVG_to_world_x * (155.802 + (370.119 - 155.802)/2),
					SVG_to_world_y * (151.179 + (192.585 - 151.179)/2)
				},
				// Image size in meter
				plot_image_size: {
					SVG_to_world_x * (370.119 - 155.802),
					SVG_to_world_y * (192.585 - 151.179)
				},
				selected_plot_image_size: {
					SVG_to_world_x * 239.530,
					SVG_to_world_y * 65.989
				},
				growth_images_locations: [
					{SVG_to_world_x * (161.920 + 203/2), SVG_to_world_y * (155.140 + 36.4/2)},
					{SVG_to_world_x * (160.987 + 204.835/2), SVG_to_world_y * (153.649 + 38.061/2)},
					{SVG_to_world_x * (160.190 + 206.870/2), SVG_to_world_y * (152.117 + 38.920/2)}
				],
				growth_images_sizes: [
					{SVG_to_world_x * 203, SVG_to_world_y * 36.4},
					{SVG_to_world_x * 204.835, SVG_to_world_y * 38.061},
					{SVG_to_world_x * 206.870, SVG_to_world_y * 38.920}
				],
				seed_icon_location: {62.510271224840764#m, 78.69243296044709#m,0.0},
				fertilizer_icon_location: {170.00144871789786#m, 76.92347306323204#m,0.0},
				growth_images: [],
				selected_image: image_file(image_path + definition + "/plots/plot_2_selected.png")
			) {
				loop i from: 1 to: 3 {
					add image_file(image_path + definition + "/plots/plot_" + myself.number + "_plant_" + i + ".png") to: growth_images;
				}
			}
		}
		
		create Plot with: (number: 3, soil: default_soil) {
			create PlotView with: (
				plot: self,
				// Coordinates gathered manually from debug_plot
				shape: polygon([
					{81.577392578125#m,51.642005920410156#m,0.0},
					{113.87451171875#m,49.19090270996094#m,0.0},
					{131.4650115966797#m,51.06534194946289#m,0.0},
					{146.60440063476562#m,52.07450485229492#m,0.0},
					{153.24830455490866#m,68.22027926859504#m,0.0},
					{134.64045835249556#m,66.92201832277529#m,0.0},
					{122.3814468383789#m,64.9069595336914#m,0.0},
					{76.67540740966797#m,65.62792205810547#m,0.0}
				]),
				// Computes the center coordinate of the plot in world coordinates (meters)
				plot_image_location: {
					SVG_to_world_x * (174.287 + (350.746 - 174.287)/2),
					SVG_to_world_y * (107.657 + (150.276 - 107.657)/2)
				},
				// Image size in meter
				plot_image_size: {
					SVG_to_world_x * (350.746 - 174.287),
					SVG_to_world_y * (150.276 - 107.657)
				},
				selected_plot_image_size: {
					SVG_to_world_x * 193.803,
					SVG_to_world_y * 59.963
				},
				growth_images_locations: [
					{SVG_to_world_x * (180.340 + 161.170/2), SVG_to_world_y * (113.680 + 31.36/2)},
					{SVG_to_world_x * (178.630 + 165.582/2), SVG_to_world_y * (112.933 + 31.967/2)},
					{SVG_to_world_x * (177.880 + 168.110/2), SVG_to_world_y * (108.953 + 36.047/2)}
				],
				growth_images_sizes: [
					{SVG_to_world_x * 161.170, SVG_to_world_y * 31.36},
					{SVG_to_world_x * 165.582, SVG_to_world_y * 31.967},
					{SVG_to_world_x * 168.110, SVG_to_world_y * 36.047}
				],
				seed_icon_location: {69.0747186830414#m, 58.17894706231469#m, 0.0},
				fertilizer_icon_location: {159.74429336008703#m, 58.99936552500823#m, 0.0},
				growth_images: [],
				selected_image: image_file(image_path + definition + "/plots/plot_3_selected.png")
			) {
				loop i from: 1 to: 3 {
					add image_file(image_path + definition + "/plots/plot_" + myself.number + "_plant_" + i + ".png") to: growth_images;
				}
			}
		}
		
		create Plot with: (number: 4, soil: default_soil) {
			create PlotView with: (
				plot: self,
				// Coordinates gathered manually from debug_plot
				shape: polygon([
					{82.53433990478516#m,47.36542510986328#m,0.0},
					{88.92381286621094#m,41.295223236083984#m,0.0},
					{101.54296875#m,40.815921783447266#m,0.0},
					{113.04410552978516#m,40.17705535888672#m,0.0},
					{124.70464324951172#m,40.97563552856445#m,0.0},
					{128.3787384033203#m,42.25352478027344#m,0.0},
					{137.3238067626953#m,42.413394927978516#m,0.0},
					{144.0327911376953#m,48.80302429199219#m,0.0},
					{129.17832634098838#m,48.62767017124233#m,0.0},
					{124.0657958984375#m,47.20555114746094#m,0.0},
					{92.09748645444498#m,47.98819593804603#m,0.0},
					{82.37462615966797#m,47.36542510986328#m,0.0}
				]),
				// Computes the center coordinate of the plot in world coordinates (meters)
				plot_image_location: {
					SVG_to_world_x * (187.359 + (329.928 - 187.359)/2),
					SVG_to_world_y * (87.847 + (107.802 - 87.847)/2)
				},
				// Image size in meter
				plot_image_size: {
					SVG_to_world_x * (329.928 - 187.359),
					SVG_to_world_y * (107.802 - 87.847)
				},
				selected_plot_image_size: {
					SVG_to_world_x * 158.467,
					SVG_to_world_y * 35.853
				},
				growth_images_locations: [
					{SVG_to_world_x * (196.990 + 121.990/2), SVG_to_world_y * (90.430 + 15.22/2)},
					{SVG_to_world_x * (193.560 + 130.789/2), SVG_to_world_y * (89.407 + 16.233/2)},
					{SVG_to_world_x * (192.590 + 134.940/2), SVG_to_world_y * (86.738 + 18.682/2)}
				],
				growth_images_sizes: [
					{SVG_to_world_x * 121.990, SVG_to_world_y * 15.220},
					{SVG_to_world_x * 130.789, SVG_to_world_y * 16.233},
					{SVG_to_world_x * 134.940, SVG_to_world_y * 18.682}
				],
				seed_icon_location: {77.68979988923124#m, 41.3576191854663#m, 0.0},
				fertilizer_icon_location: {150.7185905138058#m, 44.2294962136383#m, 0.0},
				growth_images: [],
				selected_image: image_file(image_path + definition + "/plots/plot_4_selected.png")
			) {
				loop i from: 1 to: 3 {
					add image_file(image_path + definition + "/plots/plot_" + myself.number + "_plant_" + i + ".png") to: growth_images;
				}
			}
		}
		
		create Camisol.Simple number: 4;
	}
	
	/**
	 * Reference to the plot under the cursor. Used to preview seeds and soils.
	 */
	PlotView current_plot_focus;
	bool move_on_current_plot_focus <- false;
	bool drag_on_current_plot_focus <- false;
	
	SoilView selected_soil;
	SeedView selected_seed;
	FertilizerView selected_fertilizer;
	
	action enter_plot(PlotView plot_view) {
		current_plot_focus <- plot_view;
		if (selected_seed != nil or selected_soil != nil or selected_fertilizer != nil) {
			current_plot_focus.selected <- true;
		}
	}
	action leave_plot {
		if (current_plot_focus != nil) {
			current_plot_focus.selected <- false;
			current_plot_focus <- nil;
			move_on_current_plot_focus <- false;
			drag_on_current_plot_focus <- false;
		}
	}
	
	bool mouse_hover_plots {
		bool new_plot_under_focus <- false;
		list<PlotView> plots_under_mouse <- PlotView overlapping #user_location;
		// No plot under focus or same plot
		if length(plots_under_mouse) = 0 or plots_under_mouse[0] != current_plot_focus {
			do leave_plot;
		}
		// New plot under focus case
		if length(plots_under_mouse) > 0 and plots_under_mouse[0] != current_plot_focus {
			do enter_plot(plots_under_mouse[0]);
			new_plot_under_focus <- true;
		}
		return new_plot_under_focus;
	}
	action mouse_move_plots {
		bool new_plot_under_focus <- mouse_hover_plots();
		if(!move_on_current_plot_focus) {
			move_on_current_plot_focus <- new_plot_under_focus;
		}
	}
	action mouse_drag_plots {
		bool new_plot_under_focus <- mouse_hover_plots();
		if(!drag_on_current_plot_focus) {
			drag_on_current_plot_focus <- new_plot_under_focus;
		}
	}
	
	action put_selected_item_on_plot {
		current_plot_focus.selected <- false;
		if selected_seed != nil {
			// Only if the plot is not grown yet
			if current_plot_focus.plot.growth_state <= 1 {
				ask current_plot_focus {
					do add_crop(selected_seed.seed);
				}
				current_plot_focus.plot.seed <- selected_seed.seed;
				write "Seed " + selected_seed.seed.type + " planted to " + current_plot_focus.plot.number;
			}
			ask selected_seed {
				do die;
			}
			selected_seed <- nil;
		}
		if selected_fertilizer != nil {
			// Only if the plot has not grown yet
			if current_plot_focus.plot.growth_state <= 1 {
				if length(current_plot_focus.plot.fertilizers) < 6 {
					ask current_plot_focus {
						do add_fertilizer(selected_fertilizer.fertilizer);
					}
				}
			}
			ask selected_fertilizer {
				do die;
			}
			selected_fertilizer <- nil;
		}
		if selected_soil != nil {
			current_plot_focus.plot.soil <- selected_soil.soil;
			write "Soil of plot " + current_plot_focus.plot.number + " set to " + selected_soil.soil.color;
			ask selected_soil {
				do die;
			}
			selected_soil <- nil;
		}
	}
	
	action mouse_down_plots {
		if current_plot_focus != nil and move_on_current_plot_focus {
			do put_selected_item_on_plot;
		}
	}
	
	action mouse_up_plots {
		if current_plot_focus != nil and drag_on_current_plot_focus {
			do put_selected_item_on_plot;
		}
	}
}

species EndThreadCallback {
	action call virtual: true;
}


species CsvOutput {
	int time;
	int plot_surface;
	string soil;
	string seed;
	string fertilizer_1;
	string fertilizer_2;
	string fertilizer_3;
	string fertilizer_4;
	string fertilizer_5;
	string fertilizer_6;
	float quantity;
}

/**
 * A parcel that can be planted and fertilized.
 */
species Plot skills: [thread] {
	EndThreadCallback end_thread_callback;
	
	/**
	 * Plot identifier.
	 */
	int number <- 1 min:1 max: 4;
	float surface <- 100#m2;
	
	/**
	 * Define the soil color and base parameters.
	 */
	Soil soil;
	/**
	 * Seed currently planted (nil if none is planted).
	 */
	Seed seed;
	/**
	 * Handfuls of fertilizers currently applied (empty if none is applied).
	 */
	list<Fertilizer> fertilizers;

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

	bool camisol_running <- false;
	int progress <- 0;
	
	string output_file <- "../output/plot_" + number + "_" + output_file_suffix + ".csv";

	action plant(Seed seed_to_plant) {
		seed<-seed_to_plant;
		growth_state<-1;
	}
	action remove_plant {
		seed <- nil;
		growth_state<-0;
	}
	
	action thread_action {
		camisol_running <- true;
		Plot current_plot <- self;
		float harvest;
		ask Camisol.Simple[number-1] {
			// Fertilize
			ask simulation {
				if current_time = 0 {
					do init_soil(
						C: current_plot.soil.C,
						N: current_plot.soil.N,
						P: current_plot.soil.P
					);
				}
				loop fertilizer over: current_plot.fertilizers {
					write "[Camisol] Fertilizes plot " + current_plot.number + " with " + fertilizer.name + ".";
					do fertilize
						C_rate: fertilizer.C_rate
						N_dim: fertilizer.N_dim
						P_dim: fertilizer.P_dim
						solubles: fertilizer.solubles
						hemicellulose: fertilizer.hemicellulose
						cellulose: fertilizer.cellulose
						lignine: fertilizer.lignine
						C_N: fertilizer.C_N
						C_P: fertilizer.C_P
						sample_dose: fertilizer.sample_dose;
				}
				
				// Simulate
				int max_step <- int(epoch_duration/local_step);
				write "[Camisol] Starting camisol simulation on plot " + current_plot.number + " for " + max_step + " steps.";
				loop i from: 1 to: max_step {
					do _step_;
					ask current_plot {
						progress <- int(100 * i / max_step);
						if(seed != nil) {
							growth_state <- 1 + int((2*i)/max_step);
						}
					}
				}
				// Harvest crops
				Seed crop <- current_plot.seed;
				if (crop != nil) {
					write "[Camisol] Producing " + current_plot.seed.name + " on plot " + current_plot.number;
					harvest <- production(
						N_seed: crop.N_seed,
						P_seed: crop.P_seed,
						N_plant: crop.N_plant,
						P_plant: crop.P_plant,
						harvest_index: crop.harvest_index,
						N_from_soil: crop.N_from_soil,
						plot_surface: current_plot.surface);
					write "[Camisol] " + harvest + "kg of " + current_plot.seed.name + " produced on plot " + current_plot.number;
				}
			}
		}
		write "[Camisol] Ending camisol simulation on plot " + number;
		
		create Harvest with: (
			time: current_time,
			plot: self.number,
			plot_surface: self.surface,
			seed: self.seed,
			fertilizers: self.fertilizers collect each, // List copy
			quantity: harvest
		) {
			add self to: myself.harvests;
		}
		
		create CsvOutput with: (
			time: current_time,
			plot_surface: self.surface,
			soil: self.soil.color,
			seed: self.seed = nil? "" : self.seed.name,
			fertilizer_1: length(self.fertilizers) > 0 ? self.fertilizers[0].name : "",
			fertilizer_2: length(self.fertilizers) > 1 ? self.fertilizers[1].name : "",
			fertilizer_3: length(self.fertilizers) > 2 ? self.fertilizers[2].name : "",
			fertilizer_4: length(self.fertilizers) > 3 ? self.fertilizers[3].name : "",
			fertilizer_5: length(self.fertilizers) > 4 ? self.fertilizers[4].name : "",
			fertilizer_6: length(self.fertilizers) > 5 ? self.fertilizers[5].name : "",
			quantity: harvest
		) {
			save [time, plot_surface, soil, seed, fertilizer_1, fertilizer_2, fertilizer_3, fertilizer_4, fertilizer_5, fertilizer_6, quantity] to: myself.output_file format:csv header:true rewrite: false;
		}
		
		camisol_running <- false;
		progress <- 0;
		if(end_thread_callback != nil) {
			// Only called from one plot, that waits until all other plots have finished
			ask end_thread_callback {
				do call;
			}			
		}
	}
}


species FertilizerHandfulButton parent: Button {
	PlotView plot_view;
	int index;
	bool hidden <- false;
	
	action mouse_enter {
		
	}
	action mouse_leave {
		
	}
	agent click {
		ask plot_view {
			do remove_fertilizer(myself);
		}
		return nil;
	}
	action post_click {
	}
	
	action hide {
		hidden <- true;
	}
	action show {
		hidden <- false;
	}
	aspect default {
		if(!hidden) {
			draw fertilizer_images[plot_view.plot.fertilizers[index].type-1] size:button_size;
		}
	}
}

species PlotCropButton parent: Button {
	PlotView plot_view;
	bool hidden <- false;
	
	action mouse_enter {
		
	}
	action mouse_leave {
		
	}
	agent click {
		ask plot_view {
			do remove_crop(myself);
		}
		return nil;
	}
	action post_click {
	}
	
	action hide {
		hidden <- true;
	}
	action show {
		hidden <- false;
	}
	aspect default {
		if(!hidden) {
			draw seed_images[plot_view.plot.seed.type-1] size:button_size;
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
	float button_box_width <- 0.0;
	ButtonBox fertilizers_button_box;
	ButtonBox crop_button_box;
	list<FertilizerHandfulButton> fertilizer_handful_buttons;
	PlotCropButton plot_crop_button;
	bool selected <- false;

	/**
	 * Images of the plot for each growth_state level.
	 */
	list<image_file> growth_images <- [nil, nil, nil, nil];
	image_file selected_image;
	
	/**
	 * Location of the seed icon for the current plot.
	 */
	point seed_icon_location;
	/**
	 * Location of the fertilizer icon for the current plot.
	 */
	point fertilizer_icon_location;
	point plot_image_location;
	list<point> growth_images_locations;
	
	point plot_image_size;
	point selected_plot_image_size;
	list<point> growth_images_sizes;
	
	
	/**
	 * Size of seed and fertilizer icons.
	 */
	float icon_size <- 0.9*cell_size;
	
	init {
		create ButtonBox with: (
			button_types: [FertilizerHandfulButton], visible: false, envelope: nil, hidden_envelope: nil, layer: -1) {
			myself.fertilizers_button_box <- self;
		}
		create ButtonBox with: (
			button_types: [PlotCropButton], visible: false, envelope: nil, hidden_envelope: nil, layer: -1) {
			myself.crop_button_box <- self;
		}
	}
	
	action add_crop(Seed crop) {
		crop_button_box.hidden_envelope <- rectangle(
			seed_icon_location + {-icon_size/2, icon_size/2},
			seed_icon_location + {icon_size/2, -icon_size/2}
		);
		create PlotCropButton with: (
			location: seed_icon_location,
			button_size: icon_size,
			plot_view: self
		) {
			myself.plot_crop_button <- self;
		}
		ask plot {
			do plant(crop);
		}
	}
	
	action remove_crop(PlotCropButton crop_button) {
		write "Crop " + plot.seed.type + " removed from plot " + plot.number;
		plot_crop_button <- nil;
		ask plot {
			do remove_plant;
		}
		ask crop_button {
			do die;
		}
	}
	
	action add_fertilizer(Fertilizer fertilizer) {
		if length(plot.fertilizers) = 0 {
			button_box_width <- button_box_width+icon_size;
		} else {
			button_box_width <- button_box_width+fertilizer_icon_sep;
		}
		fertilizers_button_box.hidden_envelope <- rectangle(
			fertilizer_icon_location + {-icon_size/2, icon_size/2},
			fertilizer_icon_location + {-icon_size/2+button_box_width, -icon_size/2}
		);
		create FertilizerHandfulButton with:(
			location:fertilizer_icon_location+{length(plot.fertilizers)*fertilizer_icon_sep, 0},
			button_size: icon_size,
			plot_view: self, index: length(plot.fertilizers)
		) {
			add self to:myself.fertilizer_handful_buttons;
		}
		add fertilizer to: plot.fertilizers;
		write "Plot " + current_plot_focus.plot.number + " fertilized with " + selected_fertilizer.fertilizer.type;	
	}
	
	action remove_fertilizer(FertilizerHandfulButton handful) {
		write "Fertilizer " + plot.fertilizers[handful.index].type + " removed from plot " + plot.number;
		remove from: plot.fertilizers index: handful.index;
		remove from: fertilizer_handful_buttons index: handful.index;
		int i <- handful.index;
		loop while: i < length(fertilizer_handful_buttons) {
			ask fertilizer_handful_buttons[i] {
				index <- index-1;
				location <- location - {fertilizer_icon_sep, 0};
			}
			i <- i+1;
		}
		ask handful {
			do die;
		}
		ask world {
			do mouse_move_buttons;
		}
		if length(plot.fertilizers) = 0 {
			button_box_width <- 0.0;
		} else {
			button_box_width <- button_box_width-fertilizer_icon_sep;
		}
		fertilizers_button_box.hidden_envelope <- rectangle(
			fertilizer_icon_location + {-icon_size/2, icon_size/2},
			fertilizer_icon_location + {-icon_size/2+button_box_width, -icon_size/2}
		);
	}
	
	/**
	 * Can be used to visualize the plot shape.
	 */
	aspect debug {
		if (selected_epoch.epoch.time = current_time) {
			draw shape border:#black;
			draw fertilizers_button_box.hidden_envelope color: fertilizers_button_box.background_color border:#black;
		}
	}
	aspect default {
		if (selected_epoch.epoch.time = current_time) {
			if selected {
				if selected_soil != nil {
					// Preview soil instead of the selection image
					draw soil_images[selected_soil.soil.color][plot.number-1] at: plot_image_location size:plot_image_size;	
				} else {
					draw selected_image at: plot_image_location size:selected_plot_image_size;
				}
			} else {
				// draw soil_images[plot.soil.color][plot.number-1] at: {env_width/2, env_height/2, 0};
				draw soil_images[plot.soil.color][plot.number-1] at: plot_image_location size:plot_image_size;
			}
			if plot.growth_state > 0 {
				draw growth_images[plot.growth_state-1] at: growth_images_locations[plot.growth_state-1] size: growth_images_sizes[plot.growth_state-1];
			}
		} else {
			draw soil_images[plot.soil.color][plot.number-1] at: plot_image_location size:plot_image_size;
		}
		if plot.camisol_running {
			draw "" + plot.progress + "%" at: (plot_image_location + point(0.0, 0.0, 0.05)) anchor: #center;
		}
	}
}

experiment debug_plots type:gui {
	
	action _init_ {
		create simulation;
		ask simulation {
			do init_soils();
			do init_plots();
			do init_calendar();
			create EpochView with: (epoch: ((Epoch where (each.time=current_time))[0])) {
				selected_epoch <- self;
			}
		}
	}
	
	reflex grow {
		ask Plot {
			growth_state <- (growth_state + 1) mod 4;
		}
	}
	
	output {
		display camisol type:opengl axes: true {
			// For debug purpose: click outside of the plots to select a random seed
			event #mouse_up {
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
			event #mouse_move { ask simulation {do mouse_move_plots;}}
			event #mouse_drag { ask simulation {do mouse_drag_plots;}}
			event #mouse_down { ask simulation {do mouse_down_plots;}}
			event #mouse_up { ask simulation {do mouse_up_plots;}}
			
			// image "plots" position: {0, 0, 0} size: {1, 1} file: plots;
			species PlotView aspect:default;
			// grid HelpGrid border:#black;
		}
	}
}
