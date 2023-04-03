/**
* Name: userinput
* Based on the internal skeleton template. 
* Author: pbreugno
* Tags: 
*/

model buttons

import "plot.gaml"

/**
 * The controller model contains all the heavy application mechanic, in order to:
 * - link all species with their view adaptor, and eventually associate a button to them
 * - manage Seeds, Soils, Fertilizers selection and application to Plots
 * - launch Plots micro simulation and manage the simulation Calendar
 * - ...
 */
global {
	float interface_minimum_cycle_duration <- 0.1#s;
	string NORMAL <- "NORMAL";
	string SIMULATION <- "SIMULATION";
	string HARVEST <- "HARVEST";
	string mode <- NORMAL among: [NORMAL, SIMULATION, HARVEST];
	
	// Helper values that can be used to place UI elements
	float num_cell_width <- env_width / cell_size;
	float num_cell_height <- env_height / cell_size;

	/**
	 * Buttons that should be managed by the button model.
	 */
	list<species<Button>> handled_buttons <- [SoilButton, SeedButton, SeedButtonMenu, FertilizerButton, FertilizerButtonMenu, EpochButton, RunButton];

	init {
		do load_csv_data;
		do init_seeds;
		do init_fertilizers;
		do init_soils;
		do init_plots;
		do init_calendar;
		
		write "Building buttons...";
		int i <- 0;
		
		create SoilButtonMenu;
		
		create SeedButtonMenu with: (
			location:{(num_cell_width-4.5)*cell_size, 1*cell_size},
			button_size: 0.8*cell_size
		);
		
		create FertilizerButtonMenu with: (
			location:{(num_cell_width-4.5)*cell_size, 2*cell_size},
			button_size: 0.8*cell_size
		);
		
		ask SeedButtonMenu {
			// Shows the seed button at initialization
			do show_buttons;
		}
		ask FertilizerButtonMenu {
			// Hide fertilizer button background
			do deactivate;
		}
		
		create EpochsMenu;
	}
	
	action mouse_down_harvest {
		if mode = HARVEST {
			if(current_plot_focus != nil) {
				ask current_plot_focus {
					if(self.plot.growth_state > 0) {
						create HarvestView with: (
							harvest: self.plot.harvests[current_time],
							seed_icon_location: self.seed_icon_location,
							harvest_icon_location: harvest_locations[self.plot.number-1],
							fertilizer_icon_location: self.fertilizer_icon_location
						);
						self.plot.growth_state <- 0;
						self.plot.seed <- nil;
					}
					self.plot.fertilizers <- [];
					self.fertilizers_button_box.hidden_envelope <- nil;
					ask fertilizer_handful_buttons {
						do die;
					}
					ask plot_crop_button {
						do die;
					}
				}
				int global_growth <- 0;
				ask Plot {
					global_growth <- global_growth + growth_state;
				}
				if(global_growth = 0) {
					// All plots are harvested
					ask RunButton {
						do go_to_next_epoch;
					}
					ask EpochButton {
						// Re-enables only EpochButtons.
						// Other buttons will be re-enabled when the current epoch is selected.
						do enable;
					}
					mode <- NORMAL;
				}
			}
		}
	}
}

/**
 * Button used for soil selection.
 */
species SoilButton parent: Button {
	SoilView soil_view;
	
	action mouse_enter {
		soil_view.icon_size <- 0.9*cell_size;
	}
	action mouse_leave {
		soil_view.icon_size <- 0.8*cell_size;
	}
	agent click {
		if selected_item = nil {
			// Creates a new SoilView bound to the the same Soil as the button
			create SoilView number:1 returns:new_soil_icon with:(
				soil:self.soil_view.soil,
				icon_size:0.8*cell_size,
				location:#user_location
			);
			selected_soil<-new_soil_icon[0];
			write "Selected soil: " + selected_soil.soil.to_string();
		}
		return selected_soil;
	}
	
	action post_click {
		if current_plot_focus = nil {
			ask selected_soil {
				do die;
			}
			selected_soil <- nil;		
		}
	}
}

/**
 * Seed selection button.
 */
species SeedButton parent: Button {
	SeedView seed_view;
	
	action mouse_enter {
		seed_view.icon_size <- 0.9*cell_size;
	}
	action mouse_leave {
		seed_view.icon_size <- 0.8*cell_size;
	}
	agent click {
		if selected_item = nil {
			// Creates a new SeedView bound to the the same Seed as the button
			create SeedView number:1 returns:new_seed_view with:(seed:self.seed_view.seed,icon_size:0.8*cell_size,location:#user_location);
			selected_seed<-new_seed_view[0];
			write "Selected seed: " + selected_seed.seed.to_string();
		}
		return selected_seed;
	}
	action post_click {
		if current_plot_focus = nil {
			ask selected_seed {
				do die;
			}
			selected_seed <- nil;
		}
	}
}

/**
 * Fertilizer selection button
 */
species FertilizerButton parent: Button {
	FertilizerView fertilizer_view;
	
	action mouse_enter {
		fertilizer_view.icon_size <- 0.9*cell_size;
	}
	action mouse_leave {
		fertilizer_view.icon_size <- 0.8*cell_size;
	}
	agent click {
		if selected_item = nil {
			// Creates a new FertilizerView bound to the the same Fertilizer as the button
			create FertilizerView number:1 returns:new_fertilizer_view with:(
				fertilizer:self.fertilizer_view.fertilizer,
				icon_size:0.8*cell_size,
				location:#user_location
			);
			selected_fertilizer<-new_fertilizer_view[0];
			write "Selected fertilizer: " + selected_fertilizer.fertilizer.to_string();
		}
		return selected_fertilizer;
	}
	action post_click {
		if current_plot_focus = nil {
			ask selected_fertilizer {
				do die;
			}
			selected_fertilizer <- nil;
		}
	}
}

/**
 * Epoch selection button. Only passed epochs can be selected.
 */
species EpochButton parent: Button {
	EpochView epoch_view;
	
	action mouse_enter {
		if(selected_epoch != epoch_view) {
			epoch_view.icon_size <- 1.1*button_size;
		}
	}
	action mouse_leave {
		if(selected_epoch != epoch_view) {
			epoch_view.icon_size <- button_size;
		}
	}
	agent click {
		if(selected_epoch != epoch_view) {
			epoch_view.icon_size <- button_size;
			selected_epoch <- epoch_view;

			if(selected_epoch.epoch.time != current_time) {

				do display_harvest;
			} else {
				ask HarvestView {
					// Clears all current harvest view, if any
					do die;
				}
				ask FertilizerHandfulButton + PlotCropButton + SeedButton + FertilizerButton {
					do enable;
				}
				ask FertilizerHandfulButton {
					do show;
				}
				ask PlotCropButton {
					do show;
				}
				if(current_time = 0) {
					ask SoilButton {
						do enable;
					}
				}
				if(current_time < 6) {
					ask RunButton {
						do enable;
					}
				}
			}
		}
		return nil;
	}
	
	action display_harvest {
		ask FertilizerHandfulButton + PlotCropButton + SeedButton + FertilizerButton + SoilButton + RunButton {
			do disable;
		}
		ask FertilizerHandfulButton {
			// Hides the current fertilizers
			do hide;
		}
		ask PlotCropButton {
			// Hides the current crop
			do hide;
		}
		ask HarvestView {
			// Clears all current harvest view, if any
			do die;
		}
		ask PlotView {
			create HarvestView with: (
				harvest: self.plot.harvests[selected_epoch.epoch.time],
				seed_icon_location: self.seed_icon_location,
				harvest_icon_location: harvest_locations[self.plot.number-1],
				fertilizer_icon_location: self.fertilizer_icon_location
			);
		}
	}
	action post_click {
		// Nothing to do
	}
}

species PlotEndThreadCallback parent: EndThreadCallback {
	RunButton run_button;
	
	action call {
		// True if at least one plot is still running Camisol
		bool camisol_running_on_plot <- true;
		loop while: camisol_running_on_plot {
			camisol_running_on_plot <- false;
			ask Plot {
				camisol_running_on_plot <- camisol_running_on_plot or self.camisol_running;
			}
		}
		ask run_button {
			do end_simulation;
		}
	}
}
species RunButton parent: Button {
	PlotEndThreadCallback plot_end_thread_callback;
	
	image_file button_image <- image_file(image_path + definition + "/epochs/play.png");
	image_file running_image <- image_file(image_path + definition + "/epochs/running.png");
	float icon_size <- button_size;
	
	bool running <- false;
	int init_cycle;
	
	init {
		create PlotEndThreadCallback with:(run_button: self) {
			myself.plot_end_thread_callback <- self;
		}
		shape<-square(button_size);
		ask Plot[0] {
			end_thread_callback <- myself.plot_end_thread_callback;
		}
	}
	action mouse_enter {
		icon_size <- 1.1*button_size;
	}
	action mouse_leave {
		icon_size <- button_size;
	}
	
	action launch_simulation {
		running <- true;
		init_cycle <- cycle; // Application cycle (not camisol cycle)
		// Runs the current epoch
		mode <- SIMULATION;
		
		ask FertilizerHandfulButton + SeedButton + FertilizerButton + SoilButton + EpochButton + RunButton {
			// Disables all buttons while the Camisol simulations are running
			do disable;
		}

		ask Plot {
			do run_thread;
		}
	}

	action end_simulation {
		int global_growth <- 0;
		ask PlotView {
			global_growth <- global_growth + plot.growth_state;
			if(plot.seed = nil) {
				// Nothing to harvest on this plot, but still needs an harvest view
				create HarvestView with: (
					harvest: self.plot.harvests[current_time],
					seed_icon_location: self.seed_icon_location,
					harvest_icon_location: harvest_locations[self.plot.number-1],
					fertilizer_icon_location: self.fertilizer_icon_location
				);
				self.plot.fertilizers <- [];
				self.fertilizers_button_box.hidden_envelope <- nil;
				ask fertilizer_handful_buttons {
					do die;
				}
			}
		}
		if(global_growth = 0) {
			// No harvest to do, all plots were empty
			ask RunButton {
				do go_to_next_epoch;
			}
			ask EpochButton {
				// Re-enables only EpochButtons.
				// Other buttons will be re-enabled when the current epoch is selected.
				do enable;
			}
			mode <- NORMAL;
		} else {
			mode <- HARVEST;
		}
		running <- false;
		write "[Camisol] Simulation duration: " + int((cycle-init_cycle)*interface_minimum_cycle_duration/#mn) + " minutes.";
	}
	
	action go_to_next_epoch {
		if current_time = 0 {
			// Soil can only be set at the first epoch
			ask SoilButton {
				do die;
			}
			ask SoilView {
				do die;
			}
			ask SoilButtonMenu {
				ask button_box {
					visible <- false;
				}
			}
		}
		current_time <- current_time+1;
		if current_time < 6 {
			self.location <- self.location + {2*cell_size, 0};
			create Epoch with: (time: current_time) {
				create EpochView with: (
					epoch: self,
					location: {(1.5 + 2*time)*cell_size, 8.5*cell_size}
				) {
					create EpochButton with: (
						epoch_view: self,
						location:self.location,
						button_size: cell_size
					);
				}
			}
		} else {
			do disable;
		}
	}
	
	agent click {
		icon_size <- button_size;
		do launch_simulation;
		return nil;
	}
	action post_click {
		
	}
	
	aspect default {
		if(running) {
			draw running_image size: cell_size rotate: (-(360 * (init_cycle-cycle)*interface_minimum_cycle_duration/2#s) mod 360);
		} else {
			draw button_image size: icon_size;
		}
	}
}


/**
 * ButtonMenu interface.
 * 
 * The ButtonMenu is a special button that can create and display other buttons
 * when clicked.
 */
species ButtonMenu parent: Button {
	/**
	 * Size of the menu icon.
	 */
	float icon_size <- 0.8*cell_size;
	/**
	 * True iff the ButtonMenu is currently activated, i.e. it is displaying associated buttons.
	 */
	bool active <- false;
	/**
	 * Builds the buttons coordinates organized as a grid on the right of the environment.
	 */
	list<point> button_coordinates <- build_button_coordinates();

	ButtonBox button_box;
	
	action activate {
		active <- true;
		ask button_box {
			visible <- true;
		}
	}
	action deactivate {
		active <- false;
		ask button_box {
			visible <- false;
		}
	}
	list<point> build_button_coordinates {
		int i_button <- 0;
		
		point init_coordinates <- {(num_cell_width-2.5)*cell_size, 1*cell_size};
		
		/*
		 * Defines the vector to add to the current button coordinates to plot buttons in the following order:
		 *  1  0  2
		 *  4  3  5
		 * ... 6 ...
		 */
		list<point> button_moves <- [
			// Go to left...
			{-cell_size, 0#m},
			// then go forward...
			{2*cell_size, 0#m},
			// then go below center
			{-cell_size, cell_size}
		];
		list<point> coordinates <- [init_coordinates];
		loop i from: 1 to: 18 {	
			add coordinates[i-1] + button_moves[i_button mod 3] to: coordinates;
			i_button <- i_button + 1;
		}
		return coordinates;
	}
	
	/**
	 * Creates and shows buttons associated to this menu.
	 */
	action show_buttons virtual: true;
	/**
	 * Hides and kill displayed buttons.
	 */
	action hide_buttons virtual: true;
	
	action mouse_enter {
		if(!active) {
			icon_size <- 0.9*cell_size;
		}
	}
	action mouse_leave {
		icon_size <- 0.8*cell_size;
	}
	agent click {
		if active {
			do hide_buttons;
		} else {
			do show_buttons;
		}
		// Nothing to return as no selection is performed
		return nil;
	}
	action post_click {
		// Nothing to do
	}
}

/**
 * A button menu used to display available Seed types, initialized in the seed model.
 */
species SeedButtonMenu parent: ButtonMenu {
	/**
	 * Image used to represent the menu button.
	 */
	image seed_image <- image(image_path + definition + "/crops/crop.png");

	
	init {
		create ButtonBox with: (
			button_types: [SeedButton, SeedButtonMenu]
		) {
			// Full crops layout
//			do compute_background([
//				{12, 0.5}, {16, 0.5}, {16, 6.5}, {15, 6.5}, {15, 7.5},
//				{14, 7.5}, {14, 6.5}, {13, 6.5}, {13, 1.5}, {12, 1.5}
//			]);

			// Real crops layout
			do compute_background([
				{12, 0.5}, {16, 0.5}, {16, 5.5}, {13, 5.5}, {13, 1.5}, {12, 1.5}
			]);
			// Envelope of the SeedButtonMenu only
			hidden_envelope <- envelope(myself);
			myself.button_box <- self;
		}
	}
	
	action show_buttons {
		// Activates the current menu
		do activate;
		// Hide the fertilizer menu
		ask FertilizerButtonMenu {
			if active {
				do hide_buttons;
			}
		}
		// Builds the seed menu
		int i <- 0;
		// For all available Seed...
		ask Seed {
			// ... creates a view
			create SeedView number: 1 with: (
				location: myself.button_coordinates[i],
				seed: self,
				icon_size: 0.8*cell_size
			) returns: new_seed_view;
			// ... associates a button to the view.
			create SeedButton number: 1 with: (
				location: myself.button_coordinates[i],
				seed_view: new_seed_view[0],
				button_size: 0.8*cell_size,
				enabled: mode=NORMAL and (selected_epoch = nil or selected_epoch.epoch.time = current_time)
			);
			i <- i+1;
		}
	}
	
	action hide_buttons {
		// Deactivates the seed menu
		do deactivate;
		// Deletes the button and views of the menu.
		// Notice that the available Seed are NOT deleted.
		ask SeedButton {
			ask seed_view {
				do die;
			}
			do die;
		}
	}
	
	aspect default {
		draw seed_image size:icon_size;
	}
	aspect debug {
		draw circle(0.5*icon_size) color:#red;
	}
}

/**
 * A button menu used to display available Fertilizer types, initialized in the fertilizer model.
 */
species FertilizerButtonMenu parent: ButtonMenu {

	image fertilizer_image <- image(image_path + definition + "/fertilizers/fertilizer.png");
	
	init {
		create ButtonBox with: (button_types: [FertilizerButton, FertilizerButtonMenu]) {
			do compute_background([
				{12, 1.5}, {16, 1.5}, {16, 5.5}, {13, 5.5},
				{13, 2.5}, {12, 2.5}
			]);
			// Envelope of the FertilizerButtonMenu only
			hidden_envelope <- envelope(myself);
			myself.button_box <- self;
		}
	}
	
	action show_buttons {
		// Activates the fertilizer menu
		do activate;
		// Hides the seed menu
		ask SeedButtonMenu {
			if active {
				do hide_buttons;	
			}
		}
		int i <- 1;
		// First button block
		ask OrganicFertilizer {
			// Initializes a view/button for all organic fertilizers
			create FertilizerView number: 1 with: (
				location: myself.button_coordinates[3+i-1],
				fertilizer: self,
				icon_size: 0.8*cell_size
			) returns: new_fertilizer_view;
			create FertilizerButton number: 1 with: (
				location: new_fertilizer_view[0].location, // Same location as the associated view
				fertilizer_view: new_fertilizer_view[0],
				button_size: 0.8*cell_size,
				enabled: mode=NORMAL and (selected_epoch = nil or selected_epoch.epoch.time = current_time)
			) returns: new_button;
			i<-i+1;
		}
		i<-8;
		// Second button block
		ask ChemicalFertilizer {
			// Initializes a view/button for all chemical fertilizers
			create FertilizerView number: 1 with: (
				location: myself.button_coordinates[5+i-1],
				fertilizer: self,
				icon_size: 0.8*cell_size
			) returns: new_fertilizer_view;
			create FertilizerButton number: 1 with: (
				location: new_fertilizer_view[0].location, // Same location as the associated view
				fertilizer_view: new_fertilizer_view[0],
				button_size: 0.8*cell_size,
				enabled: selected_epoch = nil or selected_epoch.epoch.time = current_time
			) returns: new_button;
			i<-i+1;
		}
	}
	
	action hide_buttons {
		// Deactivates the seed menu
		do deactivate;
		// Deletes the button and views of the menu.
		// Notice that the available Fertilizers are NOT deleted.
		ask FertilizerButton {
			ask fertilizer_view {
				do die;
			}
			do die;
		}
	}
	
	aspect default {
		draw fertilizer_image size:icon_size;
	}
	aspect debug {
		draw circle(0.5*icon_size) color:#blue;
	}
}

species SoilButtonMenu {
	ButtonBox button_box;
	init {
		create ButtonBox with:(button_types:[SoilButton]) {
			do compute_background([
				{6, 1.5}, {9, 1.5}, {9, 2.5}, {6, 2.5}
			]);
			myself.button_box <- self;
		}
		// Creates a soil view/button for each available Soil
		int i <- 0;
		ask Soil {
			create SoilView with: (
				soil: self,
				location:{(i+6.5)*cell_size, 2*cell_size}
			) returns: new_soil_views;
			create SoilButton with: (
				soil_view: new_soil_views[0],
				location: new_soil_views[0].location,
				button_size: 0.8*cell_size
			);
			i <- i+1;
		}
	}
}

species EpochsMenu {
	init {
		create ButtonBox with:(button_types:[EpochButton, RunButton]) {
			do compute_background([
				{1, 8}, {14, 8}, {14, 9}, {1, 9}
			]);
		}
		create Epoch with: (time: current_time) {
			create EpochView with: (
				epoch: self,
				location: {1.5*cell_size, 8.5*cell_size}
			) {
				create EpochButton with: (
					epoch_view: self,
					location:self.location,
					button_size: cell_size
				);
				selected_epoch <- self;
			}
		}
		create RunButton with: (
			button_size: 0.8*cell_size,
			location: {3.5*cell_size, 8.5*cell_size}
		);
	}
}

experiment debug_controller type:gui {	
	output {
		display camisol type:opengl axes: false background: #black {
			event #mouse_move action:mouse_move_buttons;
			event #mouse_down action:mouse_down_buttons;
			event #mouse_up action:mouse_up_buttons;
			
			species ButtonBox aspect: default;
			species SoilView aspect:default;
			species SeedButtonMenu aspect:default;
			species SeedView aspect:default;
			species FertilizerButtonMenu aspect:default;
			species FertilizerView aspect:default;
			species EpochView aspect: default;
			species RunButton aspect: default;
//			grid HelpGrid border:#black;
		}
	}
}
