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
}
species Enzymes {
	// Enzymatic action rates
	float T_C_c <- 0.0#gram/#h;
	float T_C_om <- 0.0#gram/#h;
	float T_N <- 0.0#gram/#h;
	float T_P <- 0.0#gram/#h;
	// list<float> t <- [0.25, 0.5, 0.75];
}

species EnzymeProducer {
	// Enzymatic activity budget
	float T_max <- 1.0#gram/#h;
	
	float x_C_c <- 1.0/4;
	float x_C_om <- 1.0/4;
	float x_P <- 1.0/4;
	float x_N <- 1.0/4;
	
	Enzymes min_enzymes;
	
	float beta_C <- 0.5;
	float beta_N <- 0.5;
	float beta_P <- 0.5;
	
	float gama_C_c <- 1.0;
	float gama_C_om <- 1.0;
	float gama_P <- 1.0;
	float gama_N <- 1.0;
	
	// The smallest K is, the more powerful the enzyme is. When K=0, the enzyme instantly imposes the limiting reaction rate, whatever the substrate concentration is.
	// Concentration at which the reaction rate from [C_c] to [C_c_dom] is half of the limiting reaction rate.
	float K_C_c <- 0.0;
	// Concentration at which the reaction rate from [C_om] to [C_dom] is half of the limiting reaction rate.
	float K_C_om <- 0.0;
	// Concentration at which the reaction rate from [N_om] to [N_dim] is half of the limiting reaction rate.
	float K_N <- 0.0;
	// Concentration at which the reaction rate from [P_om] to [P_dim] is half of the limiting reaction rate.
	float K_P <- 0.0;
	
	float L_R_enzyme_rate;
}

species SimulatedAnnealingState {
	Enzymes enzymes;
	float X_C_c_labile_to_dom <- 0.0;
	float X_C_om_labile_to_dom <- 0.0;
	float X_P_labile_to_dom <- 0.0;
	float X_P_labile_to_dim <- 0.0;
	float X_N_labile_to_dom <- 0.0;
	float X_N_labile_to_dim <- 0.0;
	
	float X_C_c_recal_to_dom <- 0.0;
	float X_C_c_recal_to_labile <- 0.0;
	float X_C_om_recal_to_dom <- 0.0;
	float X_C_om_recal_to_labile <- 0.0;
	float X_P_om_recal_to_dom <- 0.0;
	float X_P_om_recal_to_labile <- 0.0;
	float X_N_om_recal_to_dom <- 0.0;
	float X_N_om_recal_to_labile <- 0.0;
	
	float X_P_recal_to_dim <- 0.0;
	float X_P_recal_to_labile <- 0.0;
		
	float X_N_recal_to_dim <- 0.0;
	float X_N_recal_to_labile <- 0.0;

	action _init(
		float total_C_labile, float total_N_labile, float total_P_labile,
		float total_C_recal, float total_N_recal, float total_P_recal,
		EnzymeProducer enzyme_producer,
		float C_microbes,
		float dt
	) {
		float E_P_om_labile <- 0.0;
		float E_N_om_labile <- 0.0;
		if(total_C_labile > 0.0) {
			// Labile C decomposition
			// Pure C decomposition
			float E_C_c_labile <- min(
				(enzyme_producer.L_R_enzyme_rate * C_microbes * dt * enzymes.T_C_c) * total_C_labile / (enzyme_producer.K_C_c + total_C_labile),
				total_C_labile
			);
			// Organic C/N/P decomposition
			float E_C_om_labile <- min(
				(enzyme_producer.L_R_enzyme_rate * C_microbes * dt * enzymes.T_C_om) * total_C_labile / (enzyme_producer.K_C_om + total_C_labile),
				total_C_labile
			);
			
			X_C_c_labile_to_dom <- E_C_c_labile - enzyme_producer.beta_C * E_C_c_labile * E_C_om_labile / total_C_labile;
			X_C_om_labile_to_dom <- E_C_om_labile - (1-enzyme_producer.beta_C) * E_C_c_labile * E_C_om_labile / total_C_labile;
			
			// Organic P decomposition, according to the total C/P rate of the decomposed matter
			E_P_om_labile <- E_C_om_labile * total_P_labile / total_C_labile;

			// Organic N decomposition, according to the total C/N rate of the decomposed matter
			E_N_om_labile <- E_C_om_labile * total_N_labile / total_C_labile;
		}
		
		if(total_P_labile > 0) {
			// Pure P decomposition
			float E_P_labile <- min(
				(enzyme_producer.L_R_enzyme_rate * C_microbes * dt * enzymes.T_P) * total_P_labile / (enzyme_producer.K_P + total_P_labile),
				total_P_labile
			);
		
			X_P_labile_to_dom <- E_P_om_labile - enzyme_producer.beta_P * (E_P_om_labile * E_P_labile) / total_P_labile;
			X_P_labile_to_dim <- E_P_labile - (1-enzyme_producer.beta_P) * (E_P_om_labile * E_P_labile) / total_P_labile;
		}
		
		if(total_N_labile > 0.0) {
			// Pure N decomposition
			float E_N_labile <- min(
				(enzyme_producer.L_R_enzyme_rate * C_microbes * dt * enzymes.T_N) * total_N_labile / (enzyme_producer.K_N + total_N_labile),
				total_N_labile
			);
		
			X_N_labile_to_dom <- E_N_om_labile - enzyme_producer.beta_N * (E_N_om_labile * E_N_labile) / total_N_labile;
			X_N_labile_to_dim <- E_N_labile - (1-enzyme_producer.beta_N) * (E_N_om_labile * E_N_labile) / total_N_labile;
		}
		
		float E_P_om_recal <- 0.0;
		float E_N_om_recal <- 0.0;
		float X_C_c_recal <- 0.0;
		float X_C_om_recal <- 0.0;
		if (total_C_recal > 0.0) {
			// Recalcitrant C decomposition
			float E_C_c_recal <- min(
				((1-enzyme_producer.L_R_enzyme_rate) * C_microbes * dt * enzymes.T_C_c) * total_C_recal / (enzyme_producer.K_C_c + total_C_recal),
				total_C_recal
			);
			float E_C_om_recal <- min(
				((1-enzyme_producer.L_R_enzyme_rate) * C_microbes * dt * enzymes.T_C_om) * total_C_recal / (enzyme_producer.K_C_om + total_C_recal),
				total_C_recal
			);
			
			X_C_c_recal <- E_C_c_recal - enzyme_producer.beta_C * E_C_c_recal * E_C_om_recal / total_C_recal;
			X_C_om_recal <- E_C_c_recal - (1-enzyme_producer.beta_C) * E_C_c_recal * E_C_om_recal / total_C_recal;
			
			// Recalcitrant P decomposition
			E_P_om_recal <- E_C_om_recal * total_P_recal / total_C_recal;
			
			// Recalcitrant N decomposition
			E_N_om_recal <- E_C_om_recal * total_N_recal / total_C_recal;
		}

		float X_P_om_recal <- 0.0;
		float X_P_recal <- 0.0;
		if(total_P_recal > 0.0) {
			float E_P_recal <- min(
				((1-enzyme_producer.L_R_enzyme_rate) * C_microbes * dt * enzymes.T_P) * total_P_recal / (enzyme_producer.K_P + total_P_recal),
				total_P_recal
			);
			
			X_P_om_recal <- E_P_om_recal - enzyme_producer.beta_P * (E_P_om_recal * E_P_recal) / total_P_recal;
			X_P_recal <- E_P_recal - (1-enzyme_producer.beta_P) * (E_P_om_recal * E_P_recal) / total_P_recal;
		}
		
		float X_N_om_recal <- 0.0;
		float X_N_recal <- 0.0;
		if(total_N_recal > 0.0) {
			float E_N_recal <- min(
				((1-enzyme_producer.L_R_enzyme_rate) * C_microbes * dt * enzymes.T_N) * total_N_recal / (enzyme_producer.K_N + total_N_recal), total_N_recal
			);
			
			X_N_om_recal <- E_N_om_recal - enzyme_producer.beta_N * (E_N_om_recal * E_N_recal) / total_N_recal;
			X_N_recal <- E_N_recal - (1-enzyme_producer.beta_N) * (E_N_om_recal * E_N_recal) / total_N_recal;
		}
		
		X_C_c_recal_to_dom <- enzyme_producer.gama_C_c * X_C_c_recal;
		X_C_c_recal_to_labile <- X_C_c_recal - X_C_c_recal_to_dom;
		
		X_C_om_recal_to_dom <- enzyme_producer.gama_C_om * X_C_om_recal;
		X_C_om_recal_to_labile <- X_C_om_recal - X_C_om_recal_to_dom;
		X_P_om_recal_to_dom <- enzyme_producer.gama_C_om * X_P_om_recal;
		X_P_om_recal_to_labile <- X_P_om_recal - X_P_om_recal_to_dom;
		X_N_om_recal_to_dom <- enzyme_producer.gama_C_om * X_N_om_recal;
		X_N_om_recal_to_labile <- X_N_om_recal - X_N_om_recal_to_dom;
		
		X_P_recal_to_dim <- enzyme_producer.gama_P * X_P_recal;
		X_P_recal_to_labile <- X_P_recal - X_P_recal_to_dim;
		
		X_N_recal_to_dim <- enzyme_producer.gama_N * X_N_recal;
		X_N_recal_to_labile <- X_N_recal - X_N_recal_to_dim;

	}
	
	float C_avail(float C_DOM_gl, float C_DOM) {
		return C_DOM_gl + C_DOM + X_C_c_labile_to_dom + X_C_om_labile_to_dom + X_C_c_recal_to_dom + X_C_om_recal_to_dom;
	}
	
	float N_avail(float N_DOM, float N_DIM) {
		return N_DOM + N_DIM + X_N_labile_to_dom + X_N_labile_to_dim + X_N_om_recal_to_dom + X_N_recal_to_dim;
	}
	
	float P_avail(float P_DOM, float P_DIM) {
		return P_DOM + P_DIM + X_P_labile_to_dom + X_P_labile_to_dim + X_P_om_recal_to_dom + X_P_recal_to_dim;
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
	
	float C_DOM_gl;
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
			float T_min <- x_C_c * min_enzymes.T_C_c + x_C_om * min_enzymes.T_C_om + x_P * min_enzymes.T_P + x_N * min_enzymes.T_N;
			T_range <- T_max - T_min;
			_T <- [
				x_C_c * (myself.s.enzymes.T_C_c - min_enzymes.T_C_c)/T_range,
				x_C_om * (myself.s.enzymes.T_C_om - min_enzymes.T_C_om)/T_range,
				x_N * (myself.s.enzymes.T_N - min_enzymes.T_N)/T_range,
				x_P * (myself.s.enzymes.T_P - min_enzymes.T_P)/T_range
			];
		}

		list<int> is <- shuffle([0, 1, 2, 3]);
		//float delta <- 0.1+0.99 * (T/T_init);
		float delta <- 0.5;
		float range <- 1.0;
		loop i from: 0 to: 2 {
			_T[is[i]] <- max(0, min(range, gauss(_T[is[i]], 1.0)));
//			_T[i] <- rnd(max(0, _T[i] - delta), min(range, _T[i] + delta));
			range <- range - _T[is[i]];
		}
		_T[is[3]] <- range;
		Enzymes enzymes;
		ask enzyme_producer {
			create Enzymes with: [
				T_C_c: min_enzymes.T_C_c + _T[0]*T_range / x_C_c,
				T_C_om: min_enzymes.T_C_om + _T[1]*T_range / x_C_om,
				T_P: min_enzymes.T_P + _T[2]*T_range / x_P,
				T_N: min_enzymes.T_N + _T[3]*T_range / x_N
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

//		float delta <- 0.1;
//		list<float> t <- [s.t[0], s.t[1], s.t[2]];
//		loop i from: 0 to: 2 {
//			t[i] <- rnd(max(0, t[i] - delta), min(1.0, t[i] + delta));
//		}
//
//		t <- t sort_by (each);
		
//		ask enzyme_producer {
//			create Enzymes with: [
//				T_C_c: min_enzymes.T_C_c + t[0]*T_range / x_C_c,
//				T_C_om: min_enzymes.T_C_om + (t[1] - t[0])*T_range / x_C_om,
//				T_P: min_enzymes.T_P + (t[2] - t[1])*T_range / x_P,
//				T_N: min_enzymes.T_N + (1 - t[2])*T_range / x_N,
//				t: t
//				] {
//					result <- self;
//			}
//		}
		return result;
	}
	
	float E(
		SimulatedAnnealingState state
	) {
		float N_avail <- state.N_avail(N_DOM, N_DIM);
		float P_avail <- state.P_avail(P_DOM, P_DIM);
		
		float result;
		if (N_avail > 0 and P_avail > 0) {
			float C_avail <- state.C_avail(C_DOM_gl, C_DOM);
			float max_C_decomposition <- enzyme_producer.T_max / min(enzyme_producer.x_C_c, enzyme_producer.x_C_om) * dt * C_microbes;
			
			ask state {
				result <-
				(1.0 - myself.C_N / max(C_avail / N_avail, myself.C_N)) ^ 2 +
				(1.0 - myself.C_P / max(C_avail / P_avail, myself.C_P)) ^ 2
				+
				(1.0 - (X_C_c_labile_to_dom + X_C_om_labile_to_dom + X_C_c_recal_to_dom + X_C_c_recal_to_labile + X_C_om_recal_to_dom + X_C_om_recal_to_labile)
					/ max_C_decomposition
				) ^ 2
				;
			}
		} else {
			result <- #infinity;
		}
		return result;
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
	SimulatedAnnealing test_simulated_annealing;
	
	init {
		float C_N_microbes <- 5.0;
		float C_N_labile <- 5.0;
		float C_P_microbes <- 20.0;
		float C_P_labile <- 20.0;
		
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
			C_DOM_gl::0.0,
			C_DOM::0.0,
			P_DOM::0.0,
			N_DOM::0.0,
			P_DIM::0.0,
			N_DIM::0.0
		] {
			create EnzymeProducer with: [
				L_R_enzyme_rate::1.0
			] {
				create Enzymes with: [
					T_C_c::0 #gram / #h,
					T_C_om::0 #gram / #h,
					T_P::0 #gram / #h,
					T_N::0 #gram / #h
				] {
					myself.min_enzymes <- self;
				}
				myself.enzyme_producer <- self;
			}
			myself.test_simulated_annealing <- self;
		}
		
		create Enzymes with: [
			T_C_c: 1.0#gram/#hour,
			T_C_om: 1.0#gram/#hour,
			T_N: 1.0#gram/#hour,
			T_P: 1.0#gram/#hour
			] {
			myself.test_enzymes <- self;
		}
		ask test_simulated_annealing {
			do init(myself.test_enzymes);
			do step;
		}
	}
	
	reflex {
		ask test_simulated_annealing {
			do step;
			write "";
			write "s: " + s;
			write "T_C_c: " + s.enzymes.T_C_c / (#gram / #h);
			write "T_C_om: " + s.enzymes.T_C_om / (#gram / #h);
			write "T_P: " + s.enzymes.T_P / (#gram / #h);
			write "T_N: " + s.enzymes.T_N / (#gram / #h);
			write "C/N: " + s.C_avail(C_DOM_gl, C_DOM) / s.N_avail(N_DOM, N_DIM);
			write "C/P: " + s.C_avail(C_DOM_gl, C_DOM) / s.P_avail(P_DOM, P_DIM);
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

experiment SimulatedAnnealingVarC_N type: gui {
	float C_N_labile_init;
	float C_N_labile_final;
	float C_N_microbes;
	float C_P_microbes;
	float C_P_labile_init;
	float C_P_labile_final;
	
	parameter "Initial C/N (labile)" var: C_N_labile_init init: 0.5 ;
	parameter "Final C/N (labile)" var:C_N_labile_final init:20.0;
	parameter "Microbe population's C/N" var:C_N_microbes init:5.0;
	parameter "Microbe population's C/P" var:C_P_microbes init:20.0;
	parameter "Initial C/P (labile)" var:C_P_labile_init init:20.0;
	parameter "Final C/P (labile)" var:C_P_labile_final init:20.0;
	
	SimulatedAnnealing simulated_annealing;
	int steps <- 1000 parameter: "Count of steps";
	
	init {
	
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
		create SimulatedAnnealing with:[
			C_N::C_N_microbes,
			C_P::C_P_microbes,
			C_microbes::200#gram,
			dt::1#h,
			total_C_labile::200#gram,
			total_N_labile::(200#gram/(C_N_labile_init + cycle * (C_N_labile_final - C_N_labile_init)/steps)),
			total_P_labile::(200#gram/(C_P_labile_init + cycle * (C_P_labile_final - C_P_labile_init)/steps)),
			total_C_recal::0.0,
			total_N_recal::0.0,
			total_P_recal::0.0,
			C_DOM_gl::0.0,
			C_DOM::0.0,
			P_DOM::0.0,
			N_DOM::0.0,
			P_DIM::0.0,
			N_DIM::0.0
		] {
			create EnzymeProducer with: [
				L_R_enzyme_rate::1.0
			] {
				create Enzymes with: [
					T_C_c::0 #gram / #h,
					T_C_om::0 #gram / #h,
					T_P::0 #gram / #h,
					T_N::0 #gram / #h
				] {
					myself.min_enzymes <- self;
				}
				myself.enzyme_producer <- self;
			}
			Enzymes enzymes;
			create Enzymes with: [
			T_C_c: 1.0#gram/#hour,
			T_C_om: 1.0#gram/#hour,
			T_N: 1.0#gram/#hour,
			T_P: 1.0#gram/#hour
			] {
				enzymes <- self;
			}
			do init(enzymes);

			loop i from: 0 to: 1000 {
				do step;
			}
			write "";
			write "s: " + s;
			write "T_C_c: " + s.enzymes.T_C_c / (#gram / #h);
			write "T_C_om: " + s.enzymes.T_C_om / (#gram / #h);
			write "T_P: " + s.enzymes.T_P / (#gram / #h);
			write "T_N: " + s.enzymes.T_N / (#gram / #h);
			write "C/N: " + s.C_avail(C_DOM_gl, C_DOM) / s.N_avail(N_DOM, N_DIM);
			write "C/P: " + s.C_avail(C_DOM_gl, C_DOM) / s.P_avail(P_DOM, P_DIM);
			write "e: " + E(s);
						
			ask s {
				write "X_C_c_labile_to_dom: " + X_C_c_labile_to_dom;
				write "X_C_om_labile_to_dom: " + X_C_om_labile_to_dom;
				write "X_P_labile_to_dom:" + X_P_labile_to_dom;
				write "X_P_labile_to_dim: " + X_P_labile_to_dim;
				write "X_N_labile_to_dom: " + X_N_labile_to_dom;
				write "X_N_labile_to_dim: " + X_N_labile_to_dim;
				write "X_C_c_recal_to_dom: " + X_C_c_recal_to_dom;
				write "X_C_c_recal_to_labile: " + X_C_c_recal_to_labile;
				write "X_C_om_recal_to_dom: " + X_C_om_recal_to_dom;
				write "X_C_om_recal_to_labile: " + X_C_om_recal_to_labile;
				write "X_P_om_recal_to_dom: " + X_P_om_recal_to_dom;
				write "X_P_om_recal_to_labile: " + X_P_om_recal_to_labile;
				write "X_N_om_recal_to_dom: " + X_N_om_recal_to_dom;
				write "X_N_om_recal_to_labile: " + X_N_om_recal_to_labile;
				write "X_P_recal_to_dim:" + X_P_recal_to_dim;
				write "X_P_recal_to_labile:" + X_P_recal_to_labile;
				write "X_N_recal_to_dim:" + X_N_recal_to_dim;
				write "X_N_recal_to_labile:" + X_N_recal_to_labile;
			}
			myself.simulated_annealing <- self;
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
				data "Decomp C" value: sum(SimulatedAnnealing collect each.s.C_avail(each.C_DOM_gl, each.C_DOM))/#gram;
				data "Decomp N" value: sum(SimulatedAnnealing collect each.s.N_avail(each.N_DOM, each.N_DIM))/#gram;
				data "Decomp P" value: sum(SimulatedAnnealing collect each.s.P_avail(each.P_DOM, each.P_DIM))/#gram;
			}
		}
		display "Enzymatic rates" type:2d {
			chart "Enzymatic rates" type:series {
				data "T_C_c" value: sum(SimulatedAnnealing collect each.s.enzymes.T_C_c) style:line;
				data "T_C_om" value: sum(SimulatedAnnealing collect each.s.enzymes.T_C_om) style:line;
				data "T_P" value: sum(SimulatedAnnealing collect each.s.enzymes.T_P) style:line;
				data "T_N" value: sum(SimulatedAnnealing collect each.s.enzymes.T_N) style:line;
			}
		}
		display "C/N" type:2d {
			chart "C/N" type:series {
				data "organic C/N" value: sum(SimulatedAnnealing collect (each.total_C_labile/each.total_N_labile));
				data "available C/N" value: sum(SimulatedAnnealing collect (each.s.C_avail(each.C_DOM_gl, each.C_DOM)/each.s.N_avail(each.N_DOM, each.N_DIM)));
				data "C/N" value: sum(SimulatedAnnealing collect each.C_N);
			}
		}
		display "C/P" type:2d {
			chart "C/P" type:series {
				data "organic C/N" value: sum(SimulatedAnnealing collect (each.total_C_labile/each.total_P_labile));
				data "available C/P" value: sum(SimulatedAnnealing collect (each.s.C_avail(each.C_DOM_gl, each.C_DOM)/each.s.P_avail(each.P_DOM, each.P_DIM)));
				data "C/P" value: sum(SimulatedAnnealing collect each.C_P);
			}
		}
	}
}
