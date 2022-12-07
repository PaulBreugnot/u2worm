/**
* Name: crop
* Based on the internal skeleton template. 
* Author: pbreugno
* Tags: 
*/

model seed

global {
	/** Insert the global definitions, variables and actions here */
}

species Seed {
	image_file image;
	float size;
	
	aspect crop_image {
		draw image size:size;
	}
}
