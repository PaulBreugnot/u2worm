/**
* Name: landscape
* Based on the internal skeleton template. 
* Author: pbreugno
* Tags: 
*/

model landscape

global {
	file landscape <- file("../../images/landscape.png");
	
	/**
	 * For an exact match with the PNG landscape, the rectangle should be 261.3#m*148.1#m.
	 * But for convenience, values are rounded to ease objects positionning.
	 */
	float env_width <- 260#m;
	float env_height <- 150#m;
	float cell_size <- 15#m;
	geometry shape <- rectangle(env_width, env_height);
}

grid HelpGrid cell_width:cell_size cell_height:cell_size {
	rgb color <- rgb(255,255,255, 0.1);
}

experiment landscape type:gui {	
	output {
		display camisol type:opengl axes: false fullscreen: true {
			image "landscape" position: {0, 0, 0} size: {1, 1} file: landscape;
			// grid HelpGrid border:#black;
		}
	}
}