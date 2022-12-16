/**
* Name: soil
* Based on the internal empty template. 
* Author: pbreugno
* Tags: 
*/


model soil

import "landscape.gaml"

global {
	list<string> soil_colors <- ["red", "brown", "black"];
	map<string, image_file> soil_icons;
	map<string, list<image_file>> soil_images;
	
	init {
		// Loads all plot images
		loop soil_color over: soil_colors {
			soil_icons[soil_color] <- image_file("../../images/plots/plot_" + soil_color + ".png");
			soil_images[soil_color] <- [];
			loop plot from: 1 to: 4 {
				add image_file("../../images/plots/plot_" + plot + "_" + soil_color + ".png") to:soil_images[soil_color];
			}
		}
	}
}

species Soil {
	// TODO: What parameters should be associated to the soil?
	string color <- "brown" among: soil_colors;
}

species SoilView {
	Soil soil;
	float size <- 0.8*cell_size;
	
	aspect default {
		draw soil_icons[soil.color] size:size;
	}
}