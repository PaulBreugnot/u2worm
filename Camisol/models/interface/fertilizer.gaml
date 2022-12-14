/**
* Name: fertilizer
* Based on the internal empty template. 
* Author: pbreugno
* Tags: 
*/


model fertilizer

global {
	/** Insert the global definitions, variables and actions here */
	list<image_file> fertilizer_images;
	
	init {
		loop i from:0 to:9 {
			add image_file("../../images/fertilizers/fertilizer_" + (i+1) + ".png") to:fertilizer_images;
		}
	}
}

species Fertilizer {
	int type <- 1 min: 1 max:10;

}

species FertilizerView {
	Fertilizer fertilizer;
	float size;
	image_file image;

	aspect debug {
		draw circle(size/2);
	}
	aspect default {
		draw fertilizer_images[fertilizer.type-1] size:size;
	}
}
