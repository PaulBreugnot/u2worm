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
	
	geometry shape <- rectangle(261.3#m, 148.1#m);
}

experiment application type:gui {
	output {
		display camisol type:opengl axes: false fullscreen: true {
			image "landscape" position: {0, 0, 0} size: {1, 1} file: landscape;
		}
	}
}