/**
* Name: enzymesactivityexploration
* Based on the internal skeleton template. 
* Author: pbreugno
* Tags: 
*/

model enzymatic_activity_exploration

import "../microbes.gaml"

global {
	/* Parameters to explore */
	float C_labile <- 1#h * 0.01#gram/#day;
	float CN_labile <- 10.0;
	float CP_labile <- 17.0;
	float C_recalcitrant <- 1#h * 1#gram/#day;
	float CN_recalcitrant <- 10.0;
	float CP_recalcitrant <- 20.0;
	
	OrganicParticle organic_particle;
	PoreParticle pore_particle;
	/** Insert the global definitions, variables and actions here */
	
	init {
		create OrganicParticle with: [
			C_labile: world.C_labile,
			N_labile: world.C_labile / CN_labile,
			P_labile: world.C_labile / CP_labile,
			C_recalcitrant: world.C_recalcitrant,
			N_recalcitrant: world.C_recalcitrant / CN_recalcitrant,
			P_recalcitrant: world.C_recalcitrant / CP_recalcitrant
			]{
			world.organic_particle <- self;
		}
		create Dam {
			create PoreParticle with: [dam::self] {
				world.pore_particle <- self;
			}
		}
		create Copiotrophe_R with: [
			C: 1#gram,
			C_N: 10.0,
			C_P: 17.0
		] {
		}
		
		create Copiotrophe_K with: [
			C: 1#gram,
			C_N: 10.0,
			C_P: 50.0
		] {
		}
		
		create Oligotrophe_K with: [
			C: 1#gram,
			C_N: 10.0,
			C_P: 17.0
		] {
		}
	}
	
	reflex {
		ask Copiotrophe_R + Copiotrophe_K + Oligotrophe_K {
			do produce_enzymes(world.pore_particle.dam, [world.organic_particle]);
		}
	}
}


experiment ExploreCP type:batch until: cycle>0 repeat: 10 {
	list<string> strategies <- [COPIOTROPHE_R, COPIOTROPHE_K, OLIGOTROPHE_K];

	parameter "CP" var: CP_labile min: 1.0 max: 70.0 step: 1.0;
	method exploration;
	
	permanent {
		display "Enzymatic activity (Y)" {
			chart "Enzymatic activity (Y)" type:xy x_label: "C/P" y_label: "Enzymatic activity (gC/gM/d)" {
				data "T_P" value:{CP_labile, simulations mean_of(sum(Copiotrophe_R collect (each.enzymes.T_P/(#gram/#gram/#d))))};
				data "T_cellulolytic" value:{CP_labile, simulations mean_of(sum(Copiotrophe_R collect (each.enzymes.T_cellulolytic/(#gram/#gram/#d))))};
				data "T_amino" value:{CP_labile, simulations mean_of(sum(Copiotrophe_R collect (each.enzymes.T_amino/(#gram/#gram/#d))))};
				data "T_recal" value:{CP_labile, simulations mean_of(sum(Copiotrophe_R collect (each.enzymes.T_recal/(#gram/#gram/#d))))};
			}
		}
			
		display "Enzymatic activity (A)" {
			chart "Enzymatic activity (A)" type:xy x_label: "C/P" y_label: "Enzymatic activity (gC/gM/d)" {
				data "T_P" value:{CP_labile, simulations mean_of(sum(Copiotrophe_K collect (each.enzymes.T_P/(#gram/#gram/#d))))};
				data "T_cellulolytic" value:{CP_labile, simulations mean_of(sum(Copiotrophe_K collect (each.enzymes.T_cellulolytic/(#gram/#gram/#d))))};
				data "T_amino" value:{CP_labile, simulations mean_of(sum(Copiotrophe_K collect (each.enzymes.T_amino/(#gram/#gram/#d))))};
				data "T_recal" value:{CP_labile, simulations mean_of(sum(Copiotrophe_K collect (each.enzymes.T_recal/(#gram/#gram/#d))))};
			}
		}
			
		display "Enzymatic activity (S)" {
			chart "Enzymatic activity (S)" type:xy x_label: "C/P" y_label: "Enzymatic activity (gC/gM/d)" {
				data "T_P" value:{CP_labile, simulations mean_of(sum(Oligotrophe_K collect (each.enzymes.T_P/(#gram/#gram/#d))))};
				data "T_cellulolytic" value:{CP_labile, simulations mean_of(sum(Oligotrophe_K collect (each.enzymes.T_cellulolytic/(#gram/#gram/#d))))};
				data "T_amino" value:{CP_labile, simulations mean_of(sum(Oligotrophe_K collect (each.enzymes.T_amino/(#gram/#gram/#d))))};
				data "T_recal" value:{CP_labile, simulations mean_of(sum(Oligotrophe_K collect (each.enzymes.T_recal/(#gram/#gram/#d))))};
			}
		}
	}
}

experiment ExploreCP_HighC type:batch parent:ExploreCP until: cycle>0 repeat: 10 {
	init {
		C_labile <- 1#h * 1#gram/#day;
	}
}

experiment ExploreCP_LowC type:batch parent:ExploreCP until: cycle>0 repeat: 10 {
	init {
		C_labile <- 1#h * 0.01#gram/#day;
	}
}

experiment ExploreCN type:batch until: cycle>0 repeat: 10 {
	list<string> strategies <- [COPIOTROPHE_R, COPIOTROPHE_K, OLIGOTROPHE_K];
	
	parameter "CN" var: CN_labile min: 1.0 max: 70.0 step: 1.0;
	method exploration;
	
	permanent {
		display "Enzymatic activity (Y)" {
			chart "Enzymatic activity (Y)" type:xy x_label: "C/N" y_label: "Enzymatic activity (gC/gM/d)" {
				data "T_P" value:{CN_labile, simulations mean_of(sum(Copiotrophe_R collect (each.enzymes.T_P/(#gram/#gram/#d))))};
				data "T_cellulolytic" value:{CN_labile, simulations mean_of(sum(Copiotrophe_R collect (each.enzymes.T_cellulolytic/(#gram/#gram/#d))))};
				data "T_amino" value:{CN_labile, simulations mean_of(sum(Copiotrophe_R collect (each.enzymes.T_amino/(#gram/#gram/#d))))};
				data "T_recal" value:{CN_labile, simulations mean_of(sum(Copiotrophe_R collect (each.enzymes.T_recal/(#gram/#gram/#d))))};
			}
		}
			
		display "Enzymatic activity (A)" {
			chart "Enzymatic activity (A)" type:xy x_label: "C/N" y_label: "Enzymatic activity (gC/gM/d)" {
				data "T_P" value:{CN_labile, simulations mean_of(sum(Copiotrophe_K collect (each.enzymes.T_P/(#gram/#gram/#d))))};
				data "T_cellulolytic" value:{CN_labile, simulations mean_of(sum(Copiotrophe_K collect (each.enzymes.T_cellulolytic/(#gram/#gram/#d))))};
				data "T_amino" value:{CN_labile, simulations mean_of(sum(Copiotrophe_K collect (each.enzymes.T_amino/(#gram/#gram/#d))))};
				data "T_recal" value:{CN_labile, simulations mean_of(sum(Copiotrophe_K collect (each.enzymes.T_recal/(#gram/#gram/#d))))};
			}
		}
			
		display "Enzymatic activity (S)" {
			chart "Enzymatic activity (S)" type:xy x_label: "C/N" y_label: "Enzymatic activity (gC/gM/d)" {
				data "T_P" value:{CN_labile, simulations mean_of(sum(Oligotrophe_K collect (each.enzymes.T_P/(#gram/#gram/#d))))};
				data "T_cellulolytic" value:{CN_labile, simulations mean_of(sum(Oligotrophe_K collect (each.enzymes.T_cellulolytic/(#gram/#gram/#d))))};
				data "T_amino" value:{CN_labile, simulations mean_of(sum(Oligotrophe_K collect (each.enzymes.T_amino/(#gram/#gram/#d))))};
				data "T_recal" value:{CN_labile, simulations mean_of(sum(Oligotrophe_K collect (each.enzymes.T_recal/(#gram/#gram/#d))))};
			}
		}
	}
}

experiment ExploreCN_HighC type:batch parent:ExploreCN until: cycle>0 repeat: 10 {
	init {
		C_labile <- 1#h * 1#gram/#day;
	}
}

experiment ExploreCN_LowC type:batch parent:ExploreCN until: cycle>0 repeat: 10 {
	init {
		C_labile <- 1#h * 0.01#gram/#day;
	}
}
