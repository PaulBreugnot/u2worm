/**
* Name: grid
* Based on the internal skeleton template. 
* Author: pbreugno
* Tags: 
*/

model cammisol_grid

import "pore_particle.gaml"
import "mineral_particle.gaml"

global {
	string ORGANIC <- "organic";
	string MINERAL <- "mineral";
	string PORE <- "pore";

	float organic_rate <- 1/3;
	float mineral_rate <- 1/3;
	float pore_rate <- 1 - organic_rate - mineral_rate;

	int grid_size <- 30;
	
	float soil_size <- 1#cm;
	
	geometry shape <- square(soil_size);
	
	float cell_size <- soil_size/grid_size;
	float cell_area <- cell_size^2;
	
	/**
	 * Solid matter weight by volume of dry soil.
	 */
	 // Razanamalala et al. 2023, supplementary data - https://doi.org/10.3390/microbiolres14040117 
	float bulk_density <- 1.17#gram/(#cm*#cm*#cm);
	
	/**
	 * The total weight of matter in the model.
	 */
	float soil_weight <- bulk_density * world.shape.area * 1#cm;
	
		
	// Default C/N/P concentration in pom for black soils in Razanamalala et al. 2023, supplementary data
	// https://doi.org/10.3390/microbiolres14040117
	
	/**
	 * Initial C from Fresh Organic Matter concentration in soil,
	 * expressed as a weight of C by volume of soil. 
	 * This does not include C in DOM.
	 */
	float C_concentration_in_pom <- 0.02246#gram/(#cm*#cm*#cm);
	
	/**
	 * Initial N from Fresh Organic Matter concentration in soil.
	 * This does not include N in DAM.
	 */
	float N_concentration_in_pom <- 0.00137#gram/(#cm*#cm*#cm);
	
	/**
	 * Initial P from Fresh Organic Matter concentration in soil.
	 * This does not include P in DAM.
	 */
	float P_concentration_in_pom <- 0.00070#gram/(#cm*#cm*#cm);
	
	
	// TODO: Default C/N/P concentrations in dam are currently arbitrary.
	
	/**
	 * Initial C from Dissolved Organic Matter concentration in the soil,
	 * expressed as a weight of C by weight of soil.
	 */
	float C_concentration_in_dom <- (729.0 * 10^-6) #gram/#gram;
	
	/**
	 * Initial N from Dissolved Organic Matter concentration in the soil,
	 * expressed as a weight of N by weight of soil.
	 */
	float N_concentration_in_dom <- (60.0 * 10^-6) #gram/#gram;
	
	/**
	 * Initial P from Dissolved Organic Matter concentration in the soil,
	 * expressed as a weight of P by weight of soil.
	 */
	float P_concentration_in_dom <- (400.0 * 10^-6) #gram/#gram;
	
	// Default C/N/P concentration in dim for black soils in Razanamalala et al. 2023, supplementary data
	// https://doi.org/10.3390/microbiolres14040117
	
	/**
	 * Initial N from Dissolved Inorganic Matter in the soil, expressed as a weight of N by weight of soil.
	 */
	float N_concentration_in_dim <- (10.1 * 10^-6) #gram/#gram;
	/**
	 * Initial P from Dissolved Inorganic Matter in the soil, expressed as a weight of N by weight of soil.
	 */
	float P_concentration_in_dim <- (1.19 * 10^-6) #gram/#gram;

	
	/**
	 * Rate of fresh organic matter considered as labile organic matter.
	 */
	// Arbitrary
	float labile_rate_pom <- 0.2;
	
	
	action init_grid {
		ask Particle {
			type <- rnd_choice([ORGANIC::organic_rate, MINERAL::mineral_rate, PORE::pore_rate]);
			do init_particle;
		}
		// Counts the number of PORES after the initialization of the grid
		int pores_count <- length(PoreParticle);
		
		float C_in_each_dom <- soil_weight * C_concentration_in_dom / pores_count;
		float N_in_each_dom <- soil_weight * N_concentration_in_dom / pores_count;
		float P_in_each_dom <- soil_weight * P_concentration_in_dom / pores_count;
		
		float N_in_each_dim <- soil_weight * N_concentration_in_dim / pores_count;
		float P_in_each_dim <- soil_weight * P_concentration_in_dim / pores_count;
		
		ask PoreParticle {
			/* Init grid structure */
			add self to: pore_neighbors;
			ask Particle[grid_x, grid_y] {
				loop neighbor_cell over: neighbors {
					if(neighbor_cell.type = PORE) {
						add PoreParticle(neighbor_cell.particle) to: myself.pore_neighbors;
					}
				}
			}
			
			add self.organic_particle to: accessible_organics;
			ask Particle[grid_x, grid_y] {
				loop neighbor_cell over: neighbors {
					if(neighbor_cell.type = ORGANIC) { // TODO: and not from neighbor pores?
						add OrganicParticle(neighbor_cell.particle) to: myself.accessible_organics;
					}
				}
			}
			
			/* Init DAM */
			ask dam {
				dom <- [N_in_each_dom, P_in_each_dom, C_in_each_dom];
				dim <- [N_in_each_dim, P_in_each_dim];
			}
		}
			
	 	// Default soil initialization. Can be safely overriden with subsequent init_soil() calls.
		do init_POM(
			C: C_concentration_in_pom,
			N: N_concentration_in_pom,
			P: P_concentration_in_pom
		);
	}
		
	action init_POM(float C, float N, float P) {
		write "Soil C (pom): " + C/(#kg/#m3) + "kg/m3";
		write "Soil N (pom): " + N/(#kg/#m3) + "kg/m3";
		write "Soil P (pom): " + P/(#kg/#m3) + "kg/m3";
//		ask MineralParticle {
//			N <- N*cell_area;
//			P <- P*cell_area;
//		}

		float C_labile_in_each_organic <- labile_rate_pom*C*cell_area*soil_size;
		float N_labile_in_each_organic <- labile_rate_pom*N*cell_area*soil_size;
		float P_labile_in_each_organic <- labile_rate_pom*P*cell_area*soil_size;
		
		float C_recal_in_each_organic <- (1-labile_rate_pom)*C*cell_area*soil_size;
		float N_recal_in_each_organic <- (1-labile_rate_pom)*N*cell_area*soil_size;
		float P_recal_in_each_organic <- (1-labile_rate_pom)*P*cell_area*soil_size;
		
		ask OrganicParticle {
			if(!in_pore) {
				C_labile <- C_labile_in_each_organic;
				N_labile <- N_labile_in_each_organic;
				P_labile <- P_labile_in_each_organic;
				C_recalcitrant <- C_recal_in_each_organic;
				N_recalcitrant <- N_recal_in_each_organic;
				P_recalcitrant <- P_recal_in_each_organic;
			}
		}
	}
}

grid Particle width: grid_size height: grid_size neighbors: 4 {
	string type;
	agent particle;
	
	action init_particle {
		switch type {
			match MINERAL 
			{
				color <- #yellow;
				create MineralParticle with: (
					grid_x: self.grid_x,
					grid_y: self.grid_y,
					location: self.location,
					shape: self.shape
				) {
					myself.particle <- self;
				}
			}
			match ORGANIC {
				color <- #green;
				create OrganicParticle with: (
					grid_x: self.grid_x,
					grid_y: self.grid_y,
					location: self.location,
					shape: self.shape
				) {
					myself.particle <- self;
				}
			}
			match PORE { 
				color <- #black;
				// TODO: autant de nutriment dans la partie organique du pore que dans une vraie particule organique?
				OrganicParticle pore_organic_particle;
				Dam pore_dam;

				create OrganicParticle with: (
					C_labile: 0,
					N_labile: 0,
					P_labile: 0,
					C_recalcitrant: 0,
					N_recalcitrant: 0,
					P_recalcitrant: 0,
					location: self.location,
					grid_x: self.grid_x,
					grid_y: self.grid_y,
					in_pore: true,
					location: self.location,
					shape: self.shape
				) {
					pore_organic_particle <- self;
				}
				
				create Dam {
					pore_dam <- self;
				}
				create PoreParticle with: (
					organic_particle: pore_organic_particle,
					dam: pore_dam,
					grid_x: self.grid_x,
					grid_y: self.grid_y,
					location: self.location,
					shape: self.shape
				) {
					myself.particle <- self;
				}
			}
			
		}
	}
}

experiment cammisol_environment {
	bool show_legend;
	parameter "Legend" var: show_legend <- false;
	
	parameter "Grid size" category: "Grid" var:grid_size;
	parameter "Soil size" category: "Grid" var:soil_size;
	parameter "Organic particle rate" category: "Grid" var:organic_rate;
	parameter "Mineral particle rate" category: "Grid" var:mineral_rate;
	parameter "Bulk density" category: "Grid" var:bulk_density;
	
	parameter "C concentration in POM" category: "POM" var:C_concentration_in_pom;
	parameter "N concentration in POM" category: "POM" var:N_concentration_in_pom;
	parameter "P concentration in POM" category: "POM" var:P_concentration_in_pom;
	parameter "Rate of labile OM in POM" category: "POM" var:labile_rate_pom;
	
	parameter "C concentration in DOM" category: "DAM" var:C_concentration_in_dom;
	parameter "N concentration in DOM" category: "DAM" var:N_concentration_in_dom;
	parameter "P concentration in DOM" category: "DAM" var:P_concentration_in_dom;
	
	parameter "N concentration in DIM" category: "DAM" var:N_concentration_in_dim;
	parameter "P concentration in DIM" category: "DAM" var:P_concentration_in_dim;
	
	init {
		ask simulation {
			do init_grid;
		}
		write cell_size;
		do update_outputs recompute:true;
	}
	
	output {
		display grid type:opengl axes:false {
			grid Particle;
			graphics legend {
				if show_legend {
					draw square(cell_size) at: {soil_size + cell_size, 1.5*cell_size, 0} color: #yellow;
					draw square(cell_size) at: {soil_size + cell_size, 2.8*cell_size, 0} color: #green;
					draw square(cell_size) at: {soil_size + cell_size, 4.1*cell_size, 0} color: #black;
					draw "Mineral particle" font:font("Helvetica", 20 , #bold) at: {soil_size + 2.1*cell_size, 1.7*cell_size, 0} color: #black;
					draw "Organic particle" font:font("Helvetica", 20 , #bold) at: {soil_size + 2.1*cell_size, 3*cell_size, 0} color: #black;
					draw "Pore particle" font:font("Helvetica", 20 , #bold) at: {soil_size + 2.1*cell_size, 4.3*cell_size, 0} color: #black;
				}
			}
		}
		
		display "C compartments" axes:false {
			chart "C compartments (mg/g of soil)" type: pie {
				data "Labile C" value: sum(OrganicParticle collect each.C_labile)/soil_weight / (1e-3#gram/1e3#gram);
				data "Recalcitrant C" value: sum(OrganicParticle collect each.C_recalcitrant)/soil_weight / (1e-3#gram/1e3#gram);
				data "C DOM" value: sum(Dam collect each.dom[2])/soil_weight / (1e-3#gram/1e3#gram);
			}
		}
		
		display "N compartments" axes:false {
			chart "N compartments (mg/g of soil)" type: pie {
				data "Labile N" value: sum(OrganicParticle collect each.N_labile)/soil_weight / (1e-3#gram/1e3#gram);
				data "Recalcitrant N" value: sum(OrganicParticle collect each.N_recalcitrant)/soil_weight / (1e-3#gram/1e3#gram);
				data "N DOM" value: sum(Dam collect each.dom[0])/soil_weight / (1e-3#gram/1e3#gram);
				data "N DIM" value: sum(Dam collect each.dim[0])/soil_weight / (1e-3#gram/1e3#gram);
			}
		}
		
		display "P compartments" axes:false {
			chart "P compartments (mg/g of soil)" type: pie {
				data "Labile P" value: sum(OrganicParticle collect each.P_labile)/soil_weight / (1e-3#gram/1e3#gram);
				data "Recalcitrant P" value: sum(OrganicParticle collect each.P_recalcitrant)/soil_weight / (1e-3#gram/1e3#gram);
				data "P DOM" value: sum(Dam collect each.dom[1])/soil_weight / (1e-3#gram/1e3#gram);
				data "P DIM" value: sum(Dam collect each.dim[1])/soil_weight / (1e-3#gram/1e3#gram);
			}
		}
	}
}
