/**
* Name: CamisolInterface
* Based on the internal empty template. 
* Author: Paul Breugnot
* Tags: 
*/


model CamisolInterface

/* Insert your model definition here */

global {
	file landscape <- file("../images/landscape.png");
	
	/**
	 * For an exact match with the PNG landscape, the rectangle should be 261.3#m*148.1#m.
	 * But for convenience, values are rounded to ease objects positionning.
	 */
	float env_width <- 260#m;
	float env_height <- 150#m;
	float cell_size <- 15#m;
	geometry shape <- rectangle(env_width, env_height);
	
	list<string> buttons;
	int i_button <- 0;
	/**
	 * Current button coordinates.
	 */
	point button_coordinates <- init_crop_button_coordinates();
	/**
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
	
	point init_crop_button_coordinates {
		float num_cell_width <- env_width / cell_size;
		float num_cell_height <- env_height / cell_size;
		
		return {(num_cell_width-2.5)*cell_size, 0.5*cell_size};
	}
	
	action next_button_coordinates {
		button_coordinates <- button_coordinates + button_moves[i_button mod 3];
		i_button <- i_button + 1;
	}
	
	init {
		loop i from: 1 to: 19 {
			buttons <- buttons + ("crop_" + i + ".png");
		}
		loop button over: buttons {
			create Button number:1 with: (location: button_coordinates, image: image_file("../images/crops/" + button));
			do next_button_coordinates;
		}
	}
}

grid HelpGrid cell_width:cell_size cell_height:cell_size {
	rgb color <- rgb(255,255,255, 0.1);
	
}

species Button {
	image_file image;
	
	aspect debug {
		draw circle(4#m) color:#blue;
	}
	aspect button_image {
		draw image size:0.8*cell_size;
	}
}

experiment application type:gui {
	output {
		display camisol type:opengl axes: false fullscreen: true {
			image "landscape" position: {0, 0, 0} size: {1, 1} file: landscape;
			// grid HelpGrid border:#black;
			species Button aspect:button_image;
		}
	}
}