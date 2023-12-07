/**
* Name: enzymes
* Based on the internal skeleton template. 
* Author: pbreugno
* Tags: 
*/

model enzyme_activity

import "microbes.gaml"

global {
	/** Insert the global definitions, variables and actions here */
	float amino_CN <- 2.0;
	
	ExactCN exact_CN;
	MaxCN max_CN;
	ExactCP exact_CP;
	MaxCP max_CP;
	MaxC max_C;
	
	init {
		create ExactCN {
			exact_CN <- self;
		}
		create MaxCN {
			max_CN <- self;
		}
		create ExactCP {
			exact_CP <- self;
		}
		create MaxCP {
			max_CP <- self;
		}
		create MaxC {
			max_C <- self;
		}
	}
}

species Enzymes {
	// Enzymatic action rates
	float T_cellulolytic <- 0.0#gram/#h;
	float T_amino <- 0.0#gram/#h;
	float T_P <- 0.0#gram/#h;
	float T_recal <- 0.0#gram/#h;
}

species EnzymeProducer {
	// Enzymatic activity budget
	float T_max <- 1.0#gram/#h;
//	
//	float x_cellulolytic <- 0.2/3;
//	float x_amino <- 0.2/3;
//	float x_P <- 0.2/3;
//	float x_recal <- 0.8;
	
	float x_cellulolytic <- 1.0/3;
	float x_amino <- 1.0/3;
	float x_P <- 1.0/3;
	float x_recal <- 0.0;
	
	float T_min <- 0.0;
	
	float beta_cellulolytic_amino <- 0.5;
	float beta_cellulolytic_P <- 0.5;
	float beta_amino_P <- 0.5;
	
	float alpha_P <- 1e-3;
	float alpha_CP <- 0.0;
	float alpha_P_dim <- 0.9;
	
	Enzymes min_enzymes;
	
	// The smallest K is, the more powerful the enzyme is. When K=0, the enzyme instantly imposes the limiting reaction rate, whatever the substrate concentration is.
	// Concentration at which the reaction rate from [C_c] to [C_c_dom] is half of the limiting reaction rate.
	float K_cellulolytic <- 0.0;
	// Concentration at which the reaction rate from [C_om] to [C_dom] is half of the limiting reaction rate.
	float K_amino <- 0.0;
	// Concentration at which the reaction rate from [N_om] to [N_dim] is half of the limiting reaction rate.
	float K_P <- 0.0;
	
	float K_recal <- 0.0;
	
	float L_R_enzyme_rate;
	
	init {
		T_min <- x_cellulolytic * min_enzymes.T_cellulolytic + x_amino * min_enzymes.T_amino + x_P * min_enzymes.T_P + x_recal * min_enzymes.T_recal;	
	}
}

species SimulatedAnnealingState {
	Enzymes enzymes;
	float X_C_cellulolytic <- 0.0;
	float X_C_amino <- 0.0;
	float X_N_amino <- 0.0;
	float X_C_P <- 0.0;
	float X_P_labile_to_dom <- 0.0;
	float X_P_labile_to_dim <- 0.0;
	
	float X_C_recal <- 0.0;
	float X_P_recal_to_dim <- 0.0;

	action _init(
		float total_C_labile, float total_N_labile, float total_P_labile,
		float total_C_recal, float total_N_recal, float total_P_recal,
		EnzymeProducer enzyme_producer,
		float C_microbes,
		float dt
	) {
		// Cellulolytic action
		float E_C_cellulolytic <- total_C_labile > 0.0 ? min(
			(C_microbes * dt * enzymes.T_cellulolytic) * total_C_labile / (enzyme_producer.K_cellulolytic + total_C_labile),
			total_C_labile
		) : 0.0;

		float E_C_P <- min(
			(C_microbes * dt * enzymes.T_P) * total_C_labile / (enzyme_producer.K_P + total_C_labile),
			total_C_labile
		);
		
		float expected_C_amino <- total_C_labile > 0 ? (
			C_microbes * dt * enzymes.T_amino
		) * total_C_labile / (enzyme_producer.K_amino + total_C_labile) : 0.0;
		float expected_N_amino <- expected_C_amino / amino_CN;
		float limiting_amino <- expected_C_amino > 0 ? min(
			1.0, total_C_labile / expected_C_amino, total_N_labile / expected_N_amino 
		) : 0.0;
		// total_C_labile <- total_C_labile - expected_C_amino * limiting_amino;
		// total_N_labile <- total_N_labile - expected_N_amino * limiting_amino;
		float E_C_amino <- expected_C_amino * limiting_amino;
		float E_N_amino <- expected_N_amino * limiting_amino;
		
		X_C_cellulolytic <- E_C_cellulolytic
			- enzyme_producer.beta_cellulolytic_amino * E_C_amino * E_C_cellulolytic / total_C_labile
			- enzyme_producer.beta_cellulolytic_P * E_C_P * E_C_cellulolytic / total_C_labile;
			
		X_C_amino <- E_C_amino
			- (1-enzyme_producer.beta_cellulolytic_amino) * E_C_amino * E_C_cellulolytic / total_C_labile
			- enzyme_producer.beta_amino_P * E_C_amino * E_C_P / total_C_labile;
		X_C_P <- E_C_P
			- (1-enzyme_producer.beta_cellulolytic_P) * E_C_P * E_C_cellulolytic / total_C_labile
			- (1-enzyme_producer.beta_amino_P) * E_C_amino * E_C_P / total_C_labile;
		
		X_N_amino <- X_C_amino / amino_CN;
		float total_X_P <- X_C_P / (total_C_labile / total_P_labile);

		
		X_P_labile_to_dim <- enzyme_producer.alpha_P_dim * total_X_P;
				
		// Only a fraction goes to the dom, the rest is left in the labile section
		X_C_P <- enzyme_producer.alpha_CP * X_C_P;
		X_P_labile_to_dom <- enzyme_producer.alpha_CP * (total_X_P-X_P_labile_to_dim);
		
		float E_P_om_recal <- 0.0;
		float E_N_om_recal <- 0.0;
		float X_C_c_recal <- 0.0;
		float X_C_om_recal <- 0.0;
		
		// Recalcitrant C decomposition
		X_C_recal <- total_C_recal > 0.0 ? min(
			C_microbes * dt * enzymes.T_recal * total_C_recal / (enzyme_producer.K_recal + total_C_recal),
			total_C_recal
		) : 0.0;
		X_P_recal_to_dim <- (total_C_recal > 0 and total_P_recal > 0) ?
			enzyme_producer.alpha_P * X_C_recal / (total_C_recal / total_P_recal) : 0.0;
	}
	
	float C_avail(float C_DOM) {
		return C_DOM + X_C_cellulolytic + X_C_amino + X_C_P;
	}
	
	float N_avail(float N_DOM, float N_DIM) {
		return N_DOM + N_DIM + X_N_amino;
	}
	
	float P_avail(float P_DOM, float P_DIM) {
		return P_DOM + P_DIM + X_P_labile_to_dom + X_P_labile_to_dim;
	}
}

species Objective {
	action value(SimulatedAnnealing simulated_annealing, SimulatedAnnealingState state) virtual: true type: float;
}

species ExactCN parent: Objective {
	action value(SimulatedAnnealing simulated_annealing, SimulatedAnnealingState state) type: float {
		float N_avail <- state.N_avail(simulated_annealing.N_DOM, simulated_annealing.N_DIM);
		float C_avail <- state.C_avail(simulated_annealing.C_DOM);
		if(N_avail > 0 and C_avail > 0) {
			return (1.0 - simulated_annealing.C_N / (C_avail / N_avail)) ^ 2;
		}
		return #infinity;
	}
}

species MaxCN parent: Objective {
	action value(SimulatedAnnealing simulated_annealing, SimulatedAnnealingState state) type: float {
		float N_avail <- state.N_avail(simulated_annealing.N_DOM, simulated_annealing.N_DIM);
		float C_avail <- state.C_avail(simulated_annealing.C_DOM);
		if(N_avail > 0 and C_avail > 0) {
			return (1.0 - simulated_annealing.C_N / max(C_avail / N_avail, simulated_annealing.C_N)) ^ 2;
		}
		return #infinity;
	}
}

species ExactCP parent: Objective {
	action value(SimulatedAnnealing simulated_annealing, SimulatedAnnealingState state) type: float {
		float P_avail <- state.P_avail(simulated_annealing.P_DOM, simulated_annealing.P_DIM);
		float C_avail <- state.C_avail(simulated_annealing.C_DOM);
		if(P_avail > 0 and C_avail > 0) {
			return (1.0 - simulated_annealing.C_P / (C_avail / P_avail)) ^ 2;
		}
		return #infinity;
	}
}

species MaxCP parent: Objective {
	action value(SimulatedAnnealing simulated_annealing, SimulatedAnnealingState state) type: float {
		float P_avail <- state.P_avail(simulated_annealing.P_DOM, simulated_annealing.P_DIM);
		float C_avail <- state.C_avail(simulated_annealing.C_DOM);
		if(P_avail > 0 and C_avail > 0) {
			return (1.0 - simulated_annealing.C_P / max((C_avail / P_avail), simulated_annealing.C_P)) ^ 2;	
		}
		return #infinity;
	}
}

species MaxC parent: Objective {
	action value(SimulatedAnnealing simulated_annealing, SimulatedAnnealingState state) type: float {
		float max_C_decomposition;
		ask simulated_annealing {
			max_C_decomposition <- enzyme_producer.L_R_enzyme_rate * enzyme_producer.T_max /
			min(enzyme_producer.x_cellulolytic, enzyme_producer.x_amino, enzyme_producer.x_P) * dt * C_microbes;
		}
		return (1.0 - (state.X_C_cellulolytic + state.X_C_amino + state.X_C_P)
				/ max_C_decomposition
			) ^ 2;
	}
}

species SimulatedAnnealing {
	float C_N;
	float C_P;
	float C_microbes;
	float dt;
	
	float total_C_labile;
	float total_N_labile;
	float total_P_labile;
	float total_C_recal;
	float total_N_recal;
	float total_P_recal;
	
	float C_DOM;
	float P_DOM;
	float N_DOM;
	float P_DIM;
	float N_DIM;
	
	EnzymeProducer enzyme_producer;
	
	float T_init;
	float T;
	float e;
	SimulatedAnnealingState s;
	list<Objective> objectives;
	
	action init(Enzymes enzymes) {
//		T <- C_N ^ 2 + C_P ^ 2 + total_C_labile ^ 2
//		+ total_N_labile ^ 2 + total_P_labile ^ 2
//		;
		T_init <- 3.0;
		T <- T_init;
		
		create SimulatedAnnealingState with:[enzymes::enzymes] {
			do _init(
				myself.total_C_labile, myself.total_N_labile, myself.total_P_labile,
				myself.total_C_recal, myself.total_N_recal, myself.total_P_recal,
				myself.enzyme_producer, myself.C_microbes, myself.dt
			);
			myself.s <- self;
		}
		e <- E(s);
	}
	
	action step {
		SimulatedAnnealingState s_new <- neighbour();
		float e_new <- E(s_new);
		if e_new < e or rnd(1.0) < exp(-(e_new - e) / T) {
			ask s {
				ask enzymes {
					do die;
				}
				do die;
			}
			s <- s_new;
			e <- e_new;
		} else {
			ask s_new {
				ask enzymes {
					do die;
				}
				do die;
			}
		}
		T <- 0.99 * T;
	}
	
	action optimize(Enzymes enzymes) {
		do init(enzymes);
		loop i from: 0 to: 100 {
			do step;
		}
	}
	
	SimulatedAnnealingState neighbour {
		SimulatedAnnealingState result;
		
		float T_range;
		list<float> _T;
		ask enzyme_producer {
			T_range <- T_max - T_min;
			_T <- [
				x_cellulolytic * (myself.s.enzymes.T_cellulolytic - min_enzymes.T_cellulolytic)/T_range,
				x_amino * (myself.s.enzymes.T_amino - min_enzymes.T_amino)/T_range,
				x_P * (myself.s.enzymes.T_P - min_enzymes.T_P)/T_range,
				x_recal * (myself.s.enzymes.T_recal - min_enzymes.T_recal)/T_range
			];
		}

		list<int> is <- shuffle(enzyme_producer.x_recal > 0 ? [0, 1, 2, 3] : [0, 1, 2]);
		//float delta <- 0.1+0.99 * (T/T_init);
		float delta <- 0.5;
		float range <- 1.0;
		loop i from: 0 to: length(is)-2 {
			_T[is[i]] <- max(0, min(range, gauss(_T[is[i]], 1.0)));
//			_T[i] <- rnd(max(0, _T[i] - delta), min(range, _T[i] + delta));
			range <- range - _T[is[i]];
		}
		_T[is[length(is)-1]] <- range;
		Enzymes enzymes;
		ask enzyme_producer {
			create Enzymes with: [
				T_cellulolytic: min_enzymes.T_cellulolytic + _T[0]*T_range / x_cellulolytic,
				T_amino: min_enzymes.T_amino + _T[1]*T_range / x_amino,
				T_P: min_enzymes.T_P + _T[2]*T_range / x_P,
				T_recal: myself.enzyme_producer.x_recal > 0 ? min_enzymes.T_recal + _T[3]*T_range / x_recal : 0.0
				] {
					enzymes <- self;
			}
		}
		create SimulatedAnnealingState with:[enzymes::enzymes] {
			do _init(
				myself.total_C_labile, myself.total_N_labile, myself.total_P_labile,
				myself.total_C_recal, myself.total_N_recal, myself.total_P_recal,
				myself.enzyme_producer, myself.C_microbes, myself.dt
			);
			result <- self;
		}
		return result;
	}
	
	float E(
		SimulatedAnnealingState state
	) {
		return sum(objectives collect each.value(self, state));
	}
}

experiment TestEnzymes type: gui {
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


experiment TestSimulatedAnnealing type: gui {
	
	Enzymes test_enzymes;
	SimulatedAnnealing simulated_annealing;
	
	init {
		float C_N_microbes <- 5.0;
		float C_N_labile <- 5.0;
		float C_P_microbes <- 20.0;
		float C_P_labile <- 20.0;
		
		TestSimulatedAnnealing current_experiment <- self;
		create Enzymes with: [
			T_cellulolytic::0 #gram / #h,
			T_amino::0 #gram / #h,
			T_P::0 #gram / #h,
			T_recal::0 #gram / #h
		] {
			create EnzymeProducer with: [
				L_R_enzyme_rate::1.0,
				min_enzymes::self
			] {
				create SimulatedAnnealing with:[
					C_N::C_N_microbes,
					C_P::C_P_microbes,
					C_microbes::10#gram,
					dt::1#h,
					total_C_labile::200#gram,
					total_N_labile::(200#gram/C_N_labile),
					total_P_labile::(200#gram/C_P_labile),
					total_C_recal::0.0,
					total_N_recal::0.0,
					total_P_recal::0.0,
					C_DOM::0.0,
					P_DOM::0.0,
					N_DOM::0.0,
					P_DIM::0.0,
					N_DIM::0.0,
					objectives::[exact_CN, exact_CP, max_C],
					enzyme_producer::self
				] {
					current_experiment.simulated_annealing <- self;
				}
			}
		}
		
		EnzymeProducer enzyme_producer <- simulated_annealing.enzyme_producer;
		int n <- enzyme_producer.x_recal > 0 ? 4 : 3;
		float T_range <- enzyme_producer.T_max - enzyme_producer.T_min;
		create Enzymes with: [
//			T_cellulolytic: 1.0#gram/#hour,
//			T_amino: 1.0#gram/#hour,
//			T_P: 1.0#gram/#hour
			T_cellulolytic: enzyme_producer.min_enzymes.T_cellulolytic + 1/n * T_range/enzyme_producer.x_cellulolytic,
			T_amino: enzyme_producer.min_enzymes.T_amino + 1/n * T_range / enzyme_producer.x_amino,
			T_P: enzyme_producer.min_enzymes.T_P + 1/n * T_range / enzyme_producer.x_P,
			T_recal: enzyme_producer.x_recal > 0 ?
				enzyme_producer.min_enzymes.T_recal + 1/n * T_range / enzyme_producer.x_recal : 0.0
			] {
			myself.test_enzymes <- self;
		}
		ask simulated_annealing {
			do init(myself.test_enzymes);
			do step;
		}
	}
	
	reflex {
		ask simulated_annealing {
			do step;
			write "";
			write "s: " + s;
			write "T_cellulolytic: " + s.enzymes.T_cellulolytic / (#gram / #h);
			write "T_amino: " + s.enzymes.T_amino / (#gram / #h);
			write "T_P: " + s.enzymes.T_P / (#gram / #h);
			write "T_recal: " + s.enzymes.T_recal / (#gram/#h);
			write "C/N: " + s.C_avail(C_DOM) / s.N_avail(N_DOM, N_DIM);
			write "C/P: " + s.C_avail(C_DOM) / s.P_avail(P_DOM, P_DIM);
			write "e: " + E(s);

		}
	}
	
	output {
		display Energy type:2d{
			chart "E" type:series {
				data "E" value: sum(SimulatedAnnealing collect (each.e/(#gram/#h)));
			}
		}
	}
}

experiment EnzymaticActivity type: gui {
	string C_N_objective;
	string C_P_objective;
	string C_objective;
	list<Objective> objectives;

	float C_microbes;
	float total_C_labile;
	float C_dom;
	float C_N_microbes;
	float C_P_microbes;
	
	float C_N_labile_init;
	float C_N_labile_final;
	float C_P_labile_init;
	float C_P_labile_final;
	
	float C_N_dom_init;
	float C_N_dom_final;
	float C_P_dom_init;
	float C_P_dom_final;
	
	parameter "C/N" category: "Objectives" var: C_N_objective init: "Exact C/N" among: ["none", "Exact C/N", "Max C/N"];
	parameter "C/P" category: "Objectives" var: C_P_objective init: "Exact C/P" among: ["none", "Exact C/P", "Max C/P"];
	parameter "C" category: "Objectives" var: C_objective init: "none" among: ["none", "Max C"];
	
	parameter "Amino CN" category: "Constants" var: amino_CN init: 2.0;
	parameter "C microbes" category: "Constants" var: C_microbes init: 200#gram;
	parameter "total C labile" category: "Constants" var: total_C_labile init: 10#gram;
	parameter "C dom" category: "Constants" var: C_dom init: 1#gram;
	parameter "Microbe population's C/N" category: "Constants" var:C_N_microbes init:5.0;
	parameter "Microbe population's C/P" category: "Constants" var:C_P_microbes init:20.0;
	
	parameter "Initial C/N (labile)" category: "Labile OM" var: C_N_labile_init init: 20.0 ;
	parameter "Final C/N (labile)" category: "Labile OM" var:C_N_labile_final init:20.0;
	parameter "Initial C/P (labile)" category: "Labile OM" var:C_P_labile_init init:20.0;
	parameter "Final C/P (labile)" category: "Labile OM" var:C_P_labile_final init:20.0;
	
	parameter "Initial C/N (DOM)" category: "DOM" var: C_N_dom_init init: 0.5 ;
	parameter "Final C/N (DOM)" category: "DOM" var:C_N_dom_final init:20.0;
	parameter "Initial C/P (DOM)" category: "DOM" var:C_P_dom_init init:20.0;
	parameter "Final C/P (DOM)" category: "DOM" var:C_P_dom_final init:20.0;
	
	SimulatedAnnealing simulated_annealing;
	int steps <- 1000 parameter: "Count of steps";
	
	init {
		if(C_N_objective = "Exact C/N") {
			add exact_CN to: objectives;
		} else if (C_N_objective = "Max C/N") {
			add max_CN to: objectives;
		}
		
		if(C_P_objective = "Exact C/P") {
			add exact_CP to: objectives;
		} else if (C_P_objective = "Max C/P") {
			add max_CP to: objectives;
		}
		
		if(C_objective = "Max C") {
			add max_C to: objectives;
		}
	}
	
	reflex {
		ask SimulatedAnnealing {
			ask enzyme_producer {
				ask min_enzymes {
					do die;
				}
				do die;
			}
			ask s {
				ask enzymes {
					do die;
				}
				do die;
			}
			do die;
		}
		
		EnzymaticActivity exp <- self;
		create Enzymes with: [
			T_cellulolytic::0 #gram / #h,
			T_amino::0 #gram / #h,
			T_P::0 #gram / #h,
			T_recal::0 #gram / #h
		] {
			create EnzymeProducer with: [
				L_R_enzyme_rate::1.0,
				min_enzymes::self
			] {
				create SimulatedAnnealing with:[
					C_N::exp.C_N_microbes,
					C_P::exp.C_P_microbes,
					C_microbes::exp.C_microbes,
					dt::1#h,
					total_C_labile::exp.total_C_labile,
					total_N_labile::(exp.total_C_labile/(exp.C_N_labile_init + cycle * (exp.C_N_labile_final - exp.C_N_labile_init)/exp.steps)),
					total_P_labile::(exp.total_C_labile/(exp.C_P_labile_init + cycle * (exp.C_P_labile_final - exp.C_P_labile_init)/exp.steps)),
					total_C_recal::0.0,
					total_N_recal::0.0,
					total_P_recal::0.0,
					C_DOM::exp.C_dom,
					N_DOM::(exp.C_dom/(exp.C_N_dom_init + cycle * (exp.C_N_dom_final - exp.C_N_dom_init)/exp.steps)),
					P_DOM::(exp.C_dom/(exp.C_P_dom_init + cycle * (exp.C_P_dom_final - exp.C_P_dom_init)/exp.steps)),
					P_DIM::0.0,
					N_DIM::0.0,
					objectives::exp.objectives,
					enzyme_producer::self
				] {
					Enzymes enzymes;
					int n <- enzyme_producer.x_recal > 0 ? 4 : 3;
					float T_range <- enzyme_producer.T_max - enzyme_producer.T_min;
					create Enzymes with: [
						T_cellulolytic: enzyme_producer.min_enzymes.T_cellulolytic + 1/n * T_range/enzyme_producer.x_cellulolytic,
						T_amino: enzyme_producer.min_enzymes.T_amino + 1/n * T_range / enzyme_producer.x_amino,
						T_P: enzyme_producer.min_enzymes.T_P + 1/n * T_range / enzyme_producer.x_P,
						T_recal: enzyme_producer.x_recal > 0 ?
							enzyme_producer.min_enzymes.T_recal + 1/n * T_range / enzyme_producer.x_recal : 0.0
					] {
						enzymes <- self;
					}
					do init(enzymes);
		
					loop i from: 0 to: 1000 {
						do step;
					}
					write "";
					write "s: " + s;
					write "T_cellylolytic: " + s.enzymes.T_cellulolytic / (#gram / #h);
					write "T_amino: " + s.enzymes.T_amino / (#gram / #h);
					write "T_P: " + s.enzymes.T_P / (#gram / #h);
					write "C/N: " + s.C_avail(C_DOM) / s.N_avail(N_DOM, N_DIM);
					write "C/P: " + s.C_avail(C_DOM) / s.P_avail(P_DOM, P_DIM);
					write "e: " + E(s);
								
					ask s {
						write "X_C_cellulolytic: " + X_C_cellulolytic;
						write "X_C_amino: " + X_C_amino;
						write "X_N_amino: " + X_N_amino;
						write "X_C_P: " + X_C_P;
						write "X_P_labile_to_dom:" + X_P_labile_to_dom;
						write "X_P_labile_to_dim: " + X_P_labile_to_dim;
						write "X_C_recal: " + X_C_recal;
						write "X_P_recal_to_dim:" + X_P_recal_to_dim;
					}
					exp.simulated_annealing <- self;
				}
			}
		}
		
	}
	
	reflex when: cycle=steps {
		ask simulation {
			do pause;
		}
	}
	
	output {
		display "OM" type:2d {
			chart "OM" type:series {
				data "C" value: sum(SimulatedAnnealing collect each.total_C_labile)/#gram;
				data "N" value: sum(SimulatedAnnealing collect each.total_N_labile)/#gram;
				data "P" value: sum(SimulatedAnnealing collect each.total_P_labile)/#gram;
			}
		}
		display "Decomposed OM" type:2d {
			chart "Decomposed OM" type:series {
				data "Decomp C" value: sum(SimulatedAnnealing collect each.s.C_avail(each.C_DOM))/#gram;
				data "Decomp N" value: sum(SimulatedAnnealing collect each.s.N_avail(each.N_DOM, each.N_DIM))/#gram;
				data "Decomp P" value: sum(SimulatedAnnealing collect each.s.P_avail(each.P_DOM, each.P_DIM))/#gram;
			}
		}
		display "Enzymatic rates" type:2d {
			chart "Enzymatic rates" type:series {
				data "T_cellulolytic" value: sum(SimulatedAnnealing collect each.s.enzymes.T_cellulolytic) style:line;
				data "T_amino" value: sum(SimulatedAnnealing collect each.s.enzymes.T_amino) style:line;
				data "T_P" value: sum(SimulatedAnnealing collect each.s.enzymes.T_P) style:line;
				data "T_recal" value: sum(SimulatedAnnealing collect each.s.enzymes.T_recal) style:line;
			}
		}
		display "C/N" type:2d {
			chart "C/N" type:series {
				data "organic C/N" value: sum(SimulatedAnnealing collect (each.total_C_labile/each.total_N_labile));
				data "dom C/N" value: sum(SimulatedAnnealing collect (each.N_DOM > 0 ? each.C_DOM/each.N_DOM : 0.0));
				data "available C/N" value: sum(SimulatedAnnealing collect (each.s.C_avail(each.C_DOM)/each.s.N_avail(each.N_DOM, each.N_DIM)));
				data "C/N" value: sum(SimulatedAnnealing collect each.C_N);
			}
		}
		display "C/P" type:2d {
			chart "C/P" type:series {
				data "organic C/P" value: sum(SimulatedAnnealing collect (each.total_C_labile/each.total_P_labile));
				data "dom C/P" value: sum(SimulatedAnnealing collect (each.P_DOM > 0 ? each.C_DOM/each.P_DOM : 0.0));
				data "available C/P" value: sum(SimulatedAnnealing collect (each.s.C_avail(each.C_DOM)/each.s.P_avail(each.P_DOM, each.P_DIM)));
				data "C/P" value: sum(SimulatedAnnealing collect each.C_P);
			}
		}
	}
}
