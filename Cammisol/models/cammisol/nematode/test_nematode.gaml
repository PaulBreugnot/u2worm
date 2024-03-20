/**
* Name: testnematode
* Based on the internal skeleton template. 
* Author: pbreugno
* Tags: 
*/

model testnematode

import "../cammisol_grid.gaml"

global {
	int nematodes_count <- 1;
	float total_C;
	float total_N;
	float total_P;
	
	init {
		mineral_rate <- 0.0;
		do init_grid;
		
		create Nematode number: nematodes_count {
			current_pore <- one_of(PoreParticle); 
			location <- any_location_in(current_pore);
		}
		
		ask PoreParticle {
			float microbe_C <- 10 * nematode_predation_rate * local_step;
			create MicrobePopulation with:[
				C: microbe_C,
				cytosol_C: microbe_C / 5
			] {
				add self to: myself.populations;
				
				total_C <- total_C + C + cytosol_C;
				total_N <- total_N + N + cytosol_N;
				total_P <- total_P + P + cytosol_P;
			}
		}
	}
}

experiment TestNematodeMetabolism type: gui {
	parameter "Grid size" var:grid_size;
	parameter "Seed" var:seed init:10.0;
	parameter "Nematodes count" var:nematodes_count init:1;
	// parameter "Mineral rate" var:mineral_rate init:0.0;

	reflex {
		ask Nematode {
			do life;
		}
	}
	output {
		display grid {
			grid Particle;
			species Nematode aspect: red_dot;
		}
		
		display "C/N/P microbes" {
			chart "C/N/P microbes" {
				data "C (g)" value:sum(PoreParticle collect sum(each.populations collect each.C))/#gram marker:false;
				data "N (g)" value:sum(PoreParticle collect sum(each.populations collect each.N))/#gram marker:false;
				data "P (g)" value:sum(PoreParticle collect sum(each.populations collect each.P))/#gram marker:false;
				data "C (cytosol, g)" value:sum(PoreParticle collect sum(each.populations collect each.cytosol_C))/#gram marker:false;
				data "N (cytosol, g)" value:sum(PoreParticle collect sum(each.populations collect each.cytosol_N))/#gram marker:false;
				data "P (cytosol, g)" value:sum(PoreParticle collect sum(each.populations collect each.cytosol_P))/#gram marker:false;
			}
		}
		display "C/N/P organic/mineral" {
			chart "C/N/P organic" {
				data "C organic (g)" value:sum(PoreParticle collect each.organic_particle.C_labile)/#gram marker:false;
				data "N organic (g)" value:sum(PoreParticle collect each.organic_particle.N_labile)/#gram marker:false;
				data "P organic (g)" value:sum(PoreParticle collect each.organic_particle.P_labile)/#gram marker:false;
				data "N dim (g)" value:sum(Dam collect each.dim[0])/#gram marker:false;
				data "P dim (g)" value:sum(Dam collect each.dim[1])/#gram marker:false;
			}
		}
		display "Nematode awake" {
			chart "Nematode awake" {
				data "Awake (%)" value:100*(Nematode count each.awake)/nematodes_count marker:false;
			}
		}
		display "CO2" {
			chart "CO2" {
				data "CO2 (g)" value:nematode_CO2_emissions/#gram marker:false;
			}
		}
		
		display "C/N/P conservation" {
			chart "C/N/P conservation" {
				data "C" value:
					(sum(OrganicParticle collect each.C_labile) +
					sum(MicrobePopulation collect (each.C + each.cytosol_C)) +
					nematode_CO2_emissions)/total_C marker:false;
				data "N" value:
					(sum(OrganicParticle collect each.N_labile) +
					sum(MicrobePopulation collect (each.N + each.cytosol_N))
					)/total_N marker:false;
				data "P" value:
					(sum(OrganicParticle collect each.P_labile) +
					sum(MicrobePopulation collect (each.P + each.cytosol_P))
					)/total_P marker:false;
			}
		}
	}
}
