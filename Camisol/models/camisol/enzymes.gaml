/**
* Name: enzymes
* Based on the internal skeleton template. 
* Author: pbreugno
* Tags: 
*/

model enzymes

import "microbes.gaml"

global {
	/** Insert the global definitions, variables and actions here */
}

experiment enzymes type: gui {
	init {
		float carbone_concentration_in_dam <- (729.0#gram * 10^-6)/#gram;
		float azote_concentration_in_dam <- (60.0#gram * 10^-6)/#gram;
		float azote_concentration_in_dim <- (4.74#gram * 10^-6)/#gram;
		float phosphore_concentration_in_dam <- (400.0#gram * 10^-6)/#gram;
		float phosphore_concentration_in_dim <- (1.43#gram * 10^-6)/#gram;
		float model_surface <- 1#cm*1#cm;
		float model_weight <- 1.17#gram * model_surface;
		float C <- 0.02157#gram/(#cm*#cm);
		float N <- 0.00132#gram/(#cm*#cm);
		float P <- 0.00077#gram/(#cm*#cm);
		
		create Dam with: [
				dom: [
					azote_concentration_in_dam * model_weight,
					phosphore_concentration_in_dam * model_weight,
					carbone_concentration_in_dam * model_weight
				],
				dim: [
					azote_concentration_in_dim * model_weight,
					phosphore_concentration_in_dim * model_weight
				]
			] {
			}
		create PoreParticle with: [carrying_capacity::10*total_initial_bacteria_weight] {
			ask Dam {
				myself.dam <- self;
			}
		}
		
		create Copiotrophe_R {
			ask Dam {
				myself.dam <- self;
			}
			ask PoreParticle {
				add myself to: self.populations;
			}
		}
		create OrganicParticle with: [
			N: N*model_surface,
			P: P*model_surface,
			C_recalcitrant: (C/2)*model_surface,
			C_labile: (C/2)*model_surface
		] {
			ask PoreParticle {
				organic_particle <- myself;
				add myself to: accessible_organics;
			}
		}
	}
	
	reflex {
		ask Dam {
			write dom;
			write dim;
		}
	}
	
	/** Insert here the definition of the input and output of the model */
	output {
		display "dam" type: java2D {
			chart "dam" type:series {
				data "N (dom)" value: (sum(Dam collect each.dom[0])) style:spline marker:false thickness:3;
				data "P (dom)" value: (sum(Dam collect each.dom[1])) style:spline marker:false thickness:3;
				data "C (dom)" value: (sum(Dam collect each.dom[2])) style:spline marker:false thickness:3;
				data "N (dim)" value: (sum(Dam collect each.dim[0])) style:spline marker:false thickness:3;
				data "P (dim)" value: (sum(Dam collect each.dim[1])) style:spline marker:false thickness:3;
			}
		}
		
		display "organics" type:java2D {
			chart "Organics composition" type:series {
				data "C_labile" value: (sum(OrganicParticle collect each.C_labile)) style:spline marker:false thickness:3;
				data "C_recalcitrant" value: (sum(OrganicParticle collect each.C_recalcitrant)) style:spline marker:false thickness:3;
				data "N" value: (sum(OrganicParticle collect each.N)) style:spline marker:false thickness:3;
				data "P" value: (sum(OrganicParticle collect each.P)) style:spline marker:false thickness:3;
			}
		}
		
		display "populations" type:java2D {
			chart "Bacteria populations" type:series {
				data "Copiotrophe R" value: (sum(Copiotrophe_R collect each.C)) style:spline color: #red marker:false thickness:3;
				data "Copiotrophe K" value: (sum(Copiotrophe_K collect each.C)) style:spline color: #green marker:false thickness:3;
				data "Oligotrophe K" value: (sum(Oligotrophe_K collect each.C)) style:spline color: #blue marker:false thickness:3;
			}
		}
	}
}
