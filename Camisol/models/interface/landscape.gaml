/**
* Name: landscape
* Based on the internal skeleton template. 
* Author: pbreugno
* Tags: 
*/

model landscape

/**
 * Defines the base environment size and backgrounds.
 */
global {
	string definition <- "720p";
	string image_path <- "../../images/";
	file landscape <- file(image_path + definition + "/landscape.png");
	image_file sun <- image_file(image_path + "common/sun.png");
	image_file plots <- image_file(image_path + definition + "/plots.png");
	image_file plants <- image_file(image_path + definition + "/plants.png");
	
	/**
	 * For an exact match with the PNG landscape, the rectangle should be 261.3#m*148.1#m.
	 * But for convenience, values are rounded to ease objects positionning.
	 */
	float env_width <- 255#m;
	float env_height <- 150#m;
	float cell_size <- 15#m;
	geometry shape <- rectangle(env_width, env_height);
}

//grid HelpGrid cell_width:cell_size cell_height:cell_size {
//	rgb color <- rgb(255,255,255, 0.1);
//}

experiment debug_landscape type:gui {	
	output {
		display camisol type:opengl axes: false fullscreen: true {
			image sun position: {0, 0, -0.001} size: {1, 0.33};
			image landscape position: {0, 0, 0} size: {1, 1};
			image plots position: {0, 0, 0} size: {1, 1};
			image plants position: {0, 0, 0.001} size: {1, 1};
//			grid HelpGrid border:#black;
		}
	}
}