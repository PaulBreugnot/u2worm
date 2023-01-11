/**
* Name: soil
* Based on the internal empty template. 
* Author: pbreugno
* Tags: 
*/


model soil

import "fertilizer.gaml"

/**
 * Defines features used to manage the different soil types.
 */
global {
	/**
	 * Available soil colors.
	 */
	list<string> soil_colors <- ["red", "brown", "black"];
	/**
	 * Icons corresponding to each color.
	 */
	map<string, image_file> soil_icons;
	/**
	 * Soil image of each plot for a given color, so that soil_images[<color>][<index>]
	 * represents the image of the plot identified by <index> with the specified <color>.
	 */
	map<string, list<image_file>> soil_images;
	
	// Note: the regular init function is not used, since the controller model
	// needs the soils but "import soil.gaml" in controller.gaml causes the
	// soil init to be called AFTER the controller init.
	action init_soils {
		// Loads all plot images
		loop soil_color over: soil_colors {
			soil_icons[soil_color] <- image_file(image_path + definition + "/plots/plot_" + soil_color + ".png");
			soil_images[soil_color] <- [];
			loop plot from: 1 to: 4 {
				add image_file(image_path + definition + "/plots/plot_" + plot + "_" + soil_color + ".png") to:soil_images[soil_color];
			}
			create Soil number:1 with: (color: soil_color);
		}
	}
}

species Soil {
	// TODO: What parameters should be associated to the soil?
	
	/**
	 * The color of the soil. Brown by default.
	 */
	string color <- "brown" among: soil_colors;
}

/**
 * Adaptor used to visualize a soil icon. The visualization of the soil on a plot is performed in
 * the Plot visualization (PlotView).
 */
species SoilView {
	/**
	 * Visualized soil.
	 */
	Soil soil;
	/**
	 * Size of the soil icon.
	 */
	float icon_size <- 0.8*cell_size;
	
	aspect default {
		draw soil_icons[soil.color] size:icon_size;
	}
}