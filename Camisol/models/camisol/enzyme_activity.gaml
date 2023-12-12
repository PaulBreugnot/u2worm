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
	MaxN max_N;
	MaxP max_P;
	MaxRecalC max_recal_C;
	MaxTotalC max_total_C;
	
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
		create MaxN {
			max_N <- self;
		}
		create MaxP {
			max_P <- self;
		}
		create MaxRecalC {
			max_recal_C <- self;
		}
		create MaxTotalC {
			max_total_C <- self;
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
	float beta_cellulolytic_amino <- 0.5;
	float beta_cellulolytic_P <- 0.5;
	float beta_amino_P <- 0.5;
	
	float alpha_P <- 1e-3;
	float alpha_CP <- 0.0;
	float alpha_P_dim <- 0.9;
	
	Enzymes min_enzymes;
	Enzymes max_enzymes;
	
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
	}
}

species SimulatedAnnealingState {
	EnzymeProducer enzyme_producer;
	Enzymes enzymes;
	list<float> budget;
	
	float X_C_cellulolytic <- 0.0;
	float X_C_amino <- 0.0;
	float X_N_amino <- 0.0;
	float X_C_P <- 0.0;
	float X_P_labile_to_dom <- 0.0;
	float X_P_labile_to_dim <- 0.0;
	
	float X_C_recal <- 0.0;
	float X_P_recal_to_dim <- 0.0;

	init {
		enzymes <- build_enzyme();
	}
	
	Enzymes build_enzyme {
		Enzymes new_enzymes;
		ask enzyme_producer {
			create Enzymes with: [
				T_cellulolytic: min_enzymes.T_cellulolytic + myself.budget[0]*(max_enzymes.T_cellulolytic - min_enzymes.T_cellulolytic),
				T_amino: min_enzymes.T_amino + myself.budget[1]*(max_enzymes.T_amino - min_enzymes.T_amino),
				T_P: min_enzymes.T_P + myself.budget[2]*(max_enzymes.T_P - min_enzymes.T_P),
				T_recal: min_enzymes.T_recal + myself.budget[3]*(max_enzymes.T_recal - min_enzymes.T_recal)
			] {
				new_enzymes <- self;
			}
		}
		return new_enzymes;
	}
	
	action compute(
		float total_C_labile, float total_N_labile, float total_P_labile,
		float total_C_recal, float total_N_recal, float total_P_recal,
		float C_microbes,
		float dt
	) {
		if(total_C_labile > 0) {
			// Cellulolytic action
			float E_C_cellulolytic <- min(
				(C_microbes * dt * enzymes.T_cellulolytic) * total_C_labile / (enzyme_producer.K_cellulolytic + total_C_labile),
				total_C_labile
			);
	
			float E_C_P <- min(
				(C_microbes * dt * enzymes.T_P) * total_C_labile / (enzyme_producer.K_P + total_C_labile),
				total_C_labile
			);
			
			float expected_C_amino <- (
				C_microbes * dt * enzymes.T_amino
			) * total_C_labile / (enzyme_producer.K_amino + total_C_labile);
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
				
			if(total_P_labile > 0) {
				X_N_amino <- X_C_amino / amino_CN;
				float total_X_P <- X_C_P / (total_C_labile / total_P_labile);
		
				
				X_P_labile_to_dim <- enzyme_producer.alpha_P_dim * total_X_P;
						
				// Only a fraction goes to the dom, the rest is left in the labile section
				X_C_P <- enzyme_producer.alpha_CP * X_C_P;
				X_P_labile_to_dom <- enzyme_producer.alpha_CP * (total_X_P-X_P_labile_to_dim);
			}
		}

		if(total_C_recal > 0) {
			float E_P_om_recal <- 0.0;
			float E_N_om_recal <- 0.0;
			float X_C_c_recal <- 0.0;
			float X_C_om_recal <- 0.0;
			
			// Recalcitrant C decomposition
			X_C_recal <- min(
				C_microbes * dt * enzymes.T_recal * total_C_recal / (enzyme_producer.K_recal + total_C_recal),
				total_C_recal
			);
			if(total_P_recal > 0) {
				X_P_recal_to_dim <-
					enzyme_producer.alpha_P * X_C_recal / (total_C_recal / total_P_recal);
			}
		}
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
	float weight <- 1.0;
	
	action value(SimulatedAnnealing simulated_annealing, SimulatedAnnealingState state) virtual: true type: float;
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
	
	list<float> init_budget <- [1.0/4, 1.0/4, 1.0/4, 1.0/4];
	
	init {
//		T <- C_N ^ 2 + C_P ^ 2 + total_C_labile ^ 2
//		+ total_N_labile ^ 2 + total_P_labile ^ 2
//		;
		T_init <- float(length(objectives));
		T <- T_init;
		
		create SimulatedAnnealingState with:[enzyme_producer::enzyme_producer, budget::init_budget] {
			do compute(
				myself.total_C_labile, myself.total_N_labile, myself.total_P_labile,
				myself.total_C_recal, myself.total_N_recal, myself.total_P_recal,
				myself.C_microbes, myself.dt
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
	
	action optimize {
		loop i from: 0 to: 100 {
			do step;
		}
	}
	
	SimulatedAnnealingState neighbour {
		SimulatedAnnealingState result;
		
		list<float> _T;
		loop item over: s.budget {
			add item to: _T;
		}

		// The indexes list is shuffled to determine in which order each component will vary.
		// This the total budget (1.0) is fixed, this allows a fair access to the budget for
		// each component (i.e. the last component does not always get the last part left after
		// all other get their own part).
		list<int> indexes;
		loop i from: 0 to: length(_T)-1 {
			add i to: indexes;
		}
		indexes <- shuffle(indexes);
		float delta <- 0.1;
		float range <- 1.0;
		loop i from: 0 to: length(indexes)-2 {
			_T[indexes[i]] <- max(0, min(range, gauss(_T[indexes[i]], delta)));
			range <- range - _T[indexes[i]];
		}
		_T[indexes[length(indexes)-1]] <- range;
		
		create SimulatedAnnealingState with:[enzyme_producer::enzyme_producer, budget::_T] {
			do compute(
				myself.total_C_labile, myself.total_N_labile, myself.total_P_labile,
				myself.total_C_recal, myself.total_N_recal, myself.total_P_recal,
				myself.C_microbes, myself.dt
			);
			result <- self;
		}
		return result;
	}
	
	float E(
		SimulatedAnnealingState state
	) {
		return sum(objectives collect (each.weight * each.value(self, state)))/sum(objectives collect each.weight);
	}
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
			max_C_decomposition <- min(
				max(
					enzyme_producer.max_enzymes.T_cellulolytic,
					enzyme_producer.max_enzymes.T_amino, // TODO: limiting factor
					enzyme_producer.alpha_CP * enzyme_producer.max_enzymes.T_P
				) * dt * C_microbes,
				total_C_labile
			);
		}
		return (1.0 - (state.X_C_cellulolytic + state.X_C_amino + state.X_C_P)
				/ max_C_decomposition
			) ^ 2;
	}
}

species MaxP parent: Objective {
	action value(SimulatedAnnealing simulated_annealing, SimulatedAnnealingState state) type: float {
		float max_P_decomposition;
		ask simulated_annealing {
			max_P_decomposition <- max(
				min(
					total_C_labile > 0 and total_P_labile > 0 ?
					dt * C_microbes * enzyme_producer.max_enzymes.T_P / (total_C_labile / total_P_labile) : 0.0,
					total_P_labile
				),
				min(
					total_C_recal > 0 and total_P_recal > 0 ?
					enzyme_producer.alpha_P * dt * C_microbes * enzyme_producer.max_enzymes.T_recal / (total_C_recal / total_P_recal) : 0.0,
					total_P_recal
				)
			);
		}
		return max_P_decomposition > 0 ?
				(1.0 - (state.X_P_labile_to_dom + state.X_P_labile_to_dim + state.X_P_recal_to_dim) / max_P_decomposition) ^ 2
			: #infinity;
	}
}

species MaxN parent: Objective {
	action value(SimulatedAnnealing simulated_annealing, SimulatedAnnealingState state) type: float {
		float max_N_decomposition;
		ask simulated_annealing {
			max_N_decomposition <- total_C_labile > 0 and total_N_labile > 0 ?
				dt * C_microbes * enzyme_producer.max_enzymes.T_amino / (total_C_labile / total_N_labile) : 0.0;
		}
		return max_N_decomposition > 0 ?
			(1.0 - state.X_N_amino / max_N_decomposition) ^ 2
			: #infinity;
	}
}

species MaxRecalC parent: Objective {
	action value(SimulatedAnnealing simulated_annealing, SimulatedAnnealingState state) type: float {
		float max_recal_C_decomposition;
		ask simulated_annealing {
			max_recal_C_decomposition <- min(
				enzyme_producer.max_enzymes.T_recal * dt * C_microbes,
				total_C_recal
				);
		}
		return max_recal_C_decomposition > 0 ? 
			(1.0 - state.X_C_recal / max_recal_C_decomposition) ^ 2
			: (simulated_annealing.enzyme_producer.max_enzymes.T_recal - state.enzymes.T_recal > 0 ?
			exp(10000 * state.enzymes.T_recal / (simulated_annealing.enzyme_producer.max_enzymes.T_recal - state.enzymes.T_recal)) - 1.0
			: #infinity);
	}
}
	
species MaxTotalC parent: Objective {
	action value(SimulatedAnnealing simulated_annealing, SimulatedAnnealingState state) type: float {
		float result;
		ask simulated_annealing {
			float max_labile_C_enzyme <- max(
							enzyme_producer.max_enzymes.T_cellulolytic,
							enzyme_producer.max_enzymes.T_amino, // TODO: limiting factor
							enzyme_producer.alpha_CP * enzyme_producer.max_enzymes.T_P
						);
			if(enzyme_producer.max_enzymes.T_recal > 0 and max_labile_C_enzyme > 0) {
				float max_C_decomposition <- min(
					max_labile_C_enzyme * dt * C_microbes,
					total_C_labile
				);
				float max_C_recal_decomposition <- min(
					enzyme_producer.max_enzymes.T_recal * dt * C_microbes,
					total_C_recal
				);
				// It is possible to produced recalcitrant enzymes
				if(total_C_recal > 0 and total_C_labile > 0) {
					result <- (1.0 - (state.X_C_recal + state.X_C_cellulolytic + state.X_C_amino + state.X_C_P) / (max_C_decomposition + max_C_recal_decomposition)) ^ 2;
				} else if (total_C_labile > 0) {
					// Maximize labile C decomposition while prohibiting recalcitrant enzyme production (since total_C_recal = 0).
					result <- (simulated_annealing.enzyme_producer.max_enzymes.T_recal - state.enzymes.T_recal) > 0 ?
						(1.0 - (state.X_C_cellulolytic + state.X_C_amino + state.X_C_P) / (max_C_decomposition)) ^ 2
							* exp(state.enzymes.T_recal / (simulated_annealing.enzyme_producer.max_enzymes.T_recal - state.enzymes.T_recal))
						: #infinity;
				} else if (total_C_recal > 0) {
					// Maximize recalcitrant C decomposition while prohibiting labile enzyme production (since total_C_labile = 0).
					float T_labile <- state.enzymes.T_cellulolytic + state.enzymes.T_amino + enzyme_producer.alpha_CP * enzyme_producer.max_enzymes.T_P;
					float T_max_labile <- enzyme_producer.max_enzymes.T_cellulolytic + enzyme_producer.max_enzymes.T_amino
						+ enzyme_producer.alpha_CP * enzyme_producer.max_enzymes.T_P;
					result <- (T_max_labile - T_labile) > 0 ?
							(1.0 - (state.X_C_cellulolytic + state.X_C_amino + state.X_C_P) / (max_C_decomposition)) ^ 2
								* exp(10000 * T_labile / (T_max_labile - T_labile))
						: #infinity;
				} else {
					result <- 0.0;
				}
			} else if (max_labile_C_enzyme > 0) {
				// It is not possible to produce recalcitrant enzyme, so act like if they do not exist
				result <- max_C.value(simulated_annealing, state);
			} else if(enzyme_producer.max_enzymes.T_recal > 0) {
				// It is not possible to produce labile enzyme, so act like if they do not exist
				result <- max_recal_C.value(simulated_annealing, state);
			} else {
				// No C enzyme can be produced.
				result <- 0.0;
			}
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
	SimulatedAnnealing simulated_annealing;
	
	init {
		float C_N_microbes <- 5.0;
		float C_N_labile <- 5.0;
		float C_P_microbes <- 20.0;
		float C_P_labile <- 20.0;
		
		TestSimulatedAnnealing current_experiment <- self;
		
		Enzymes min_enzymes;
		create Enzymes with: [
			T_cellulolytic::0 #gram / #h,
			T_amino::0 #gram / #h,
			T_P::0 #gram / #h,
			T_recal::0 #gram / #h
		] {
			min_enzymes <- self;
		}
		
		Enzymes max_enzymes;
		create Enzymes with: [
			T_cellulolytic::1 #gram / #h,
			T_amino::1 #gram / #h,
			T_P::1 #gram / #h,
			T_recal::0.2 #gram / #h
		] {
			max_enzymes <- self;
		}
		
		create EnzymeProducer with: [
			L_R_enzyme_rate::1.0,
			min_enzymes::min_enzymes,
			max_enzymes::max_enzymes
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
				objectives::[max_total_C],
				enzyme_producer::self
			] {
				current_experiment.simulated_annealing <- self;
			}
		}
		
		ask simulated_annealing {
			do step;
		}
	}
	
	reflex {
		ask simulated_annealing {
			do step;
			float C_avail <- s.C_avail(C_DOM);
			float N_avail <- s.N_avail(N_DOM, N_DIM);
			float P_avail <- s.P_avail(P_DOM, P_DIM);
			write "";
			write "s: " + s;
			write "T_cellulolytic: " + s.enzymes.T_cellulolytic / (#gram / #h);
			write "T_amino: " + s.enzymes.T_amino / (#gram / #h);
			write "T_P: " + s.enzymes.T_P / (#gram / #h);
			write "T_recal: " + s.enzymes.T_recal / (#gram/#h);
			write "C/N: " + (N_avail > 0 ? C_avail / N_avail : 0.0);
			write "C/P: " + (P_avail > 0 ? C_avail / P_avail : 0.0);
			write "e: " + E(s);

		}
	}
	
	output {
		display Energy type:2d{
			chart "E" type:series {
				data "E" value: sum(SimulatedAnnealing collect (each.e/(#gram/#h))) marker:false;
			}
		}
	}
}

experiment EnzymaticActivity type: gui {
	string C_N_objective;
	float C_N_objective_weight;
	string C_P_objective;
	float C_P_objective_weight;
	string C_objective;
	float C_objective_weight;
	string N_objective;
	float N_objective_weight;
	string P_objective;
	float P_objective_weight;
	string C_recal_objective;
	float C_recal_objective_weight;
	list<Objective> objectives;

	float T_cellulolytic;
	float T_amino;
	float T_P;
	float T_recal;
	float C_microbes;
	float C_dom;
	float C_N_microbes;
	float C_P_microbes;
	
	float C_labile_init;
	float C_labile_final;
	float C_N_labile_init;
	float C_N_labile_final;
	float C_P_labile_init;
	float C_P_labile_final;
	
	float C_N_dom_init;
	float C_N_dom_final;
	float C_P_dom_init;
	float C_P_dom_final;
	
	float C_recal_init;
	float C_recal_final;
	
	int steps;
	
	parameter "C/N" category: "Objectives" var: C_N_objective init: "Max C/N" among: ["none", "Exact C/N", "Max C/N"];
	parameter "C/N weight" category: "Objectives" var: C_N_objective_weight init: 100.0;
	parameter "C/P" category: "Objectives" var: C_P_objective init: "Max C/P" among: ["none", "Exact C/P", "Max C/P"];
	parameter "C/P weight" category: "Objectives" var: C_P_objective_weight init: 100.0;
	parameter "C" category: "Objectives" var: C_objective init: "Max total C" among: ["none", "Max C", "Max total C"];
	parameter "C weight" category: "Objectives" var: C_objective_weight init: 1.0;
	parameter "N" category: "Objectives" var: N_objective init: "none" among: ["none", "Max N"];
	parameter "N weight" category: "Objectives" var: N_objective_weight init: 1.0;
	parameter "P" category: "Objectives" var: P_objective init: "none" among: ["none", "Max P"];
	parameter "P weight" category: "Objectives" var: P_objective_weight init: 1.0;
	parameter "C recalcitrant" category: "Objectives" var: C_recal_objective init: "none" among: ["none", "Max recalcitrant C"];
	parameter "C recalcitrant weight" category: "Objectives" var: C_recal_objective_weight init: 1.0;
	
	parameter "Max cellulolytic (g/h)" category: "Constants" var: T_cellulolytic init: 1.0;
	parameter "Max amino (g/h)" category: "Constants" var: T_amino init: 1.0;
	parameter "Max P (g/h)" category: "Constants" var: T_P init: 1.0;
	parameter "Max recalcitrant (g/h)" category: "Constants" var: T_recal init: 0.0;
	parameter "Amino CN" category: "Constants" var: amino_CN init: 6.0;
	parameter "C microbes" category: "Constants" var: C_microbes init: 200#gram;
	parameter "C dom" category: "Constants" var: C_dom init: 1#gram;
	parameter "Microbe population's C/N" category: "Constants" var:C_N_microbes init:5.0;
	parameter "Microbe population's C/P" category: "Constants" var:C_P_microbes init:20.0;
	
	parameter "Initial C (labile)" category: "Labile OM" var: C_labile_init init: 10#gram;
	parameter "Final C (labile)" category: "Labile OM" var: C_labile_final init: 10#gram;
	parameter "Initial C/N (labile)" category: "Labile OM" var: C_N_labile_init init: 20.0 ;
	parameter "Final C/N (labile)" category: "Labile OM" var:C_N_labile_final init:20.0;
	parameter "Initial C/P (labile)" category: "Labile OM" var:C_P_labile_init init:20.0;
	parameter "Final C/P (labile)" category: "Labile OM" var:C_P_labile_final init:20.0;
	
	parameter "Initial C/N (DOM)" category: "DOM" var: C_N_dom_init init: 0.5 ;
	parameter "Final C/N (DOM)" category: "DOM" var:C_N_dom_final init:20.0;
	parameter "Initial C/P (DOM)" category: "DOM" var:C_P_dom_init init:20.0;
	parameter "Final C/P (DOM)" category: "DOM" var:C_P_dom_final init:20.0;
	
	parameter "Initial C (recalcitrant)" category: "Recalcitrant OM" var: C_recal_init init: 0.0;
	parameter "Final C (recalcitrant)" category: "Recalcitrant OM" var: C_recal_final init: 0.0;
	
	SimulatedAnnealing simulated_annealing;
	parameter "Count of steps" category: "Experiment" var: steps init: 1000;
	
	init {
//		write "Labile C threshold: " + C_microbes * T_cellulolytic #gram/#h * 1#h / #gram;
//		write "Recalcitrant C threshold: " + C_microbes * T_recal #gram/#h * 1#h / #gram;
		Objective _C_N_objective;
		if(C_N_objective = "Exact C/N") {
			_C_N_objective <- exact_CN;
		} else if (C_N_objective = "Max C/N") {
			_C_N_objective <- max_CN;
		}
		if _C_N_objective != nil {
			ask _C_N_objective {
				weight <- myself.C_N_objective_weight;
				add self to: myself.objectives;
			}	
		}
		
		Objective _C_P_objective;
		if(C_P_objective = "Exact C/P") {
			_C_P_objective <- exact_CP;
		} else if (C_P_objective = "Max C/P") {
			_C_P_objective <- max_CP;
		}
		if _C_P_objective != nil {
			ask _C_P_objective {
				weight <- myself.C_P_objective_weight;
				add self to: myself.objectives;
			}	
		}
		
		Objective _C_objective;
		if(C_objective = "Max C") {
			_C_objective <- max_C;
		} else if (C_objective = "Max total C") {
			_C_objective <- max_total_C;
		}
		if(_C_objective != nil) {
			ask _C_objective {
				weight <- myself.C_objective_weight;
				add self to: myself.objectives;
			}	
		}
				
		if(N_objective = "Max N") {
			add max_N to: objectives;
			ask max_N {
				weight <- myself.N_objective_weight;
			}
		}
				
		if(P_objective = "Max P") {
			add max_P to: objectives;
			ask max_P {
				weight <- myself.P_objective_weight;
			}
		}
		
		if(C_recal_objective = "Max recalcitrant C") {
			add max_recal_C to: objectives;
			ask max_recal_C {
				weight <- myself.C_recal_objective_weight;
			}
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
		
		Enzymes min_enzymes;
		create Enzymes with: [
			T_cellulolytic::0 #gram / #h,
			T_amino::0 #gram / #h,
			T_P::0 #gram / #h,
			T_recal::0 #gram / #h
		] {
			min_enzymes <- self;
		}
		Enzymes max_enzymes;
		create Enzymes with: [
			T_cellulolytic::T_cellulolytic #gram/#h,
			T_amino::T_amino #gram/#h,
			T_P::T_P #gram/#h,
			T_recal::T_recal #gram/#h
		] {
			max_enzymes <- self;
		}
		EnzymaticActivity exp <- self;
		
		create EnzymeProducer with: [
			L_R_enzyme_rate::1.0,
			min_enzymes::min_enzymes,
			max_enzymes::max_enzymes
		] {
			float total_C_labile <- exp.C_labile_init + cycle * (exp.C_labile_final - exp.C_labile_init)/exp.steps;
			float total_N_labile <- total_C_labile/(exp.C_N_labile_init + cycle * (exp.C_N_labile_final - exp.C_N_labile_init)/exp.steps);
			float total_P_labile <- total_C_labile/(exp.C_P_labile_init + cycle * (exp.C_P_labile_final - exp.C_P_labile_init)/exp.steps);
			
			float total_C_recal <- exp.C_recal_init + cycle * (exp.C_recal_final - exp.C_recal_init)/exp.steps;
			
			create SimulatedAnnealing with:[
				C_N::exp.C_N_microbes,
				C_P::exp.C_P_microbes,
				C_microbes::exp.C_microbes,
				dt::1#h,
				total_C_labile::total_C_labile,
				total_N_labile::total_N_labile,
				total_P_labile::total_P_labile,
				total_C_recal::total_C_recal,
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
				loop i from: 0 to: 1000 {
					do step;
				}
				float C_avail <- s.C_avail(C_DOM);
				float N_avail <- s.N_avail(N_DOM, N_DIM);
				float P_avail <- s.P_avail(P_DOM, P_DIM);
				write "";
				write "s: " + s;
				write "T_cellylolytic: " + s.enzymes.T_cellulolytic / (#gram / #h);
				write "T_amino: " + s.enzymes.T_amino / (#gram / #h);
				write "T_P: " + s.enzymes.T_P / (#gram / #h);
				write "C/N: " + (N_avail > 0 ? C_avail / N_avail : 0.0);
				write "C/P: " + (P_avail > 0 ? C_avail / P_avail : 0.0);
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
	
	reflex when: cycle=steps {
		ask simulation {
			do pause;
		}
	}
	
	output {
		display "OM" type:2d {
			chart "OM" type:series style:line {
				data "C (labile)" value: sum(SimulatedAnnealing collect each.total_C_labile)/#gram marker:false;
				data "N (labile)" value: sum(SimulatedAnnealing collect each.total_N_labile)/#gram marker:false;
				data "P (labile)" value: sum(SimulatedAnnealing collect each.total_P_labile)/#gram marker:false;
				data "C (recal)" value: sum(SimulatedAnnealing collect each.total_C_recal)/#gram marker:false;
				data "N (recal)" value: sum(SimulatedAnnealing collect each.total_N_recal)/#gram marker:false;
				data "P (recal)" value: sum(SimulatedAnnealing collect each.total_P_recal)/#gram marker:false;
			}
		}
		display "Enzymatic rates" type:2d {
			chart "Enzymatic rates" type:series style:line {
				data "T_cellulolytic" value: sum(SimulatedAnnealing collect (each.s.enzymes.T_cellulolytic)) marker:false;
				data "T_amino" value: sum(SimulatedAnnealing collect (each.s.enzymes.T_amino)) marker:false;
				data "T_P" value: sum(SimulatedAnnealing collect (each.s.enzymes.T_P)) marker:false;
				data "T_recal" value: sum(SimulatedAnnealing collect (each.s.enzymes.T_recal)) marker:false;
				data "T_cellulolytic (max)" value: sum(SimulatedAnnealing collect (each.enzyme_producer.max_enzymes.T_cellulolytic)) marker:false;
				data "T_amino (max)" value: sum(SimulatedAnnealing collect (each.enzyme_producer.max_enzymes.T_amino)) marker:false;
				data "T_P (max)" value: sum(SimulatedAnnealing collect (each.enzyme_producer.max_enzymes.T_P)) marker:false;
				data "T_recal (max)" value: sum(SimulatedAnnealing collect (each.enzyme_producer.max_enzymes.T_recal)) marker:false;
			}
		}
		display "C/N" type:2d {
			chart "C/N" type:series style:line {
				data "organic C/N" value: sum(SimulatedAnnealing collect (each.total_N_labile > 0 ? each.total_C_labile/each.total_N_labile : 0.0)) marker:false;
				data "dom C/N" value: sum(SimulatedAnnealing collect (each.N_DOM > 0 ? each.C_DOM/each.N_DOM : 0.0)) marker:false;
				data "available C/N" value: sum(SimulatedAnnealing collect (
					each.s.N_avail(each.N_DOM, each.N_DIM) > 0 ?
					each.s.C_avail(each.C_DOM)/each.s.N_avail(each.N_DOM, each.N_DIM) : 0.0
				)) marker:false;
				data "C/N" value: sum(SimulatedAnnealing collect each.C_N) marker:false;
			}
		}
		display "C/P" type:2d {
			chart "C/P" type:series style:line {
				data "organic C/P" value: sum(SimulatedAnnealing collect (each.total_P_labile > 0 ? each.total_C_labile/each.total_P_labile : 0.0)) marker:false;
				data "dom C/P" value: sum(SimulatedAnnealing collect (each.P_DOM > 0 ? each.C_DOM/each.P_DOM : 0.0)) marker:false;
				data "available C/P" value: sum(SimulatedAnnealing collect (
					each.s.P_avail(each.P_DOM, each.P_DIM) > 0 ?
					each.s.C_avail(each.C_DOM) / each.s.P_avail(each.P_DOM, each.P_DIM) : 0.0
				)) marker:false;
				data "C/P" value: sum(SimulatedAnnealing collect each.C_P) marker:false;
			}
		}
	}
}
