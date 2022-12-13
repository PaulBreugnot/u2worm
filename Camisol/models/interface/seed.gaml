/**
* Name: crop
* Based on the internal skeleton template. 
* Author: pbreugno
* Tags: 
*/

model seed

global {
	/** Insert the global definitions, variables and actions here */
	list<image_file> seed_images;
	
	init {
		loop i from:0 to:18 {
			add image_file("../../images/crops/crop_" + (i+1) + ".png") to:seed_images;
		}
	}
}

species Seed {
	int type <- 1 min: 1 max:19;

}

species SeedView {
	Seed seed;
	float size;
	image_file image;

	aspect debug {
		draw circle(size/2);
	}
	aspect default {
		draw seed_images[seed.type-1] size:size;
	}
}
