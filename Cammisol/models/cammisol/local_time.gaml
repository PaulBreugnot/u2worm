/**
* Name: configuration
* Based on the internal empty template. 
* Author: pbreugno
* Tags: 
*/


model local_time

/* Insert your model definition here */

global {
	int local_cycle <- 0;
	float local_step <- 1#h;
	float local_time <- 0.0;

	reflex {
		local_cycle <- local_cycle+1;
		local_time <- local_cycle * local_step;
	}
}