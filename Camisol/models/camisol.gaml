/***
* Name: camisole
* Author: Laetitia BERNARD, Lucas GROSJEAN, Nicolas MARILLEAU
***/

model camisol

import "csv_function.gaml"

//todo calcul facteur efficacité enzymatique récalcitrant
//todo facteur mortalité virale (besoin des chiffres)
//todo absorption mineral (besoin des chiffres) + desorption 
//todo fixation d'un partie des enzyme sur les minéraux + relachement pendant la désorption
//todo améliorer mouvement nématode (besoin des chiffres)
//todo passage en 3D + fractalité

//todo ajout des echelles + repartition des cases selons les chiffres des échelles
//todo sortie sol (production + porosité) -> déja fait mais besoin de faire les échelles
//todo entrée organique -> déja fait mais besoin de faire les échelles
//flux hydrique -> à faire mais besoin de faire les échelles

//todo intégration de u2worm dans le modèle
//todo ajout des popultation dans les nouveaux pores (anciens organics)

global
{
	string ORGANIC <- "organic";
	string MINERAL <- "mineral";
	string PORE <- "pore";
	
	int grid_width <- 30;
	int grid_height <- 30;
	
	geometry shape <- square(1#cm);
	
	//Charactéristiques des décomposeurs :
	float temps_division_opportuniste <- 1#h;
	float temps_division_decomposeur <- 2#h;
	
	float dividing_time_copiotrophe_R <- 1#h;
	float dividing_time_copiotrophe_K <- 24#h;
	float dividing_time_oligotrophe <- 368#h;
	
	float growth_decomposeur <- 0.5/#mn;
	
	float CO2_produit; 
	
	float carbone_concentration_in_dam <- (729.0#gram * 10^-6)/#gram;
	float azote_concentration_in_dam <- (60.0#gram * 10^-6)/#gram;
	float phosphore_concentration_in_dam <- (400.0#gram * 10^-6)/#gram;
	
	int nb_copioK <- 0;
	int nb_copioR <- 0;
	int nb_oligo <- 0;
	int nb_nematode <- 30;
	
	float min_heatmap <- 0.0;
	float max_heatmap <- 0.0;
	
	float facteur_enzyme_cr <- 1.0;
	float facteur_enzyme_cl <- 1.0;
	
	float total_model_weight <- bulk_density * world.shape.area * 1#cm;
	float nb_bacterie_per_cm <- 5 * 10^9; 
	
	int nb_pores -> {particle count ( each.my_type = PORE) };
	int nb_organics -> {particle count ( each.my_type = ORGANIC) };
	int nb_minerals -> {particle count ( each.my_type = MINERAL) };
	
	float poids_sec_bacterie <-  0.0000000000003#gram;
	float weight_population_copioK <- (poids_sec_bacterie * nb_bacterie_per_cm * world.shape.area) / 3;
	float weight_population_copioR <- (poids_sec_bacterie * nb_bacterie_per_cm * world.shape.area) / 3;
	float weight_population_oligo <- (poids_sec_bacterie * nb_bacterie_per_cm * world.shape.area) / 3;
	
	//  C_init de pop => weight_init
	// nbC * weight_init / C_init => weight_at_T  
	float weight_bacterie_in_one_pore -> { (weight_population_copioK + weight_population_copioR + weight_population_oligo) / nb_pores};
	
	//  0.22gC => 0.135day-1 / 24 (qty de C par heure)
	// weight_at_T * (activite_enzymatique_per_hour) / 0.22 => qty de Carbone à produire
	
	/* Julien Sainte-marie */
	float masse_microbienne_etude <- 0.22#gram;
	float activite_enzymatique <- 0.135#day-1;
	float activite_enzymatique_per_hour <- activite_enzymatique / 24;
	/* Julien Sainte-marie */
	
	
	float total_masseC <- 0.0;
	float total_masseC_init <- 0.0;
	
	float total_co2_exp;
	
	float total_particle_C_labile -> { sum(particle collect(each.C_labile)) } ;
	float total_particle_C_recalcitrant -> { sum(particle collect(each.C_recalcitrant)) };
	
	float total_pores_C_available_C -> { sum(pores_particles collect(each.dam.available_C)) };
	
	float C_labile -> {(sum(particle where( each.my_type = PORE) collect(each.C_labile)))};
	float C_recalcitrant -> {(sum(particle where( each.my_type = PORE) collect(each.C_recalcitrant)))};
	
	float total_enzyme_cl -> {sum(pores_particles collect(each.dam.enzyme_Cl))};
	float total_enzyme_cr -> {sum(pores_particles collect(each.dam.enzyme_Cr))};
	
	// #g
	float pop_nematode -> {sum(nematode collect(each.C))};
	float total_nematode_stomackC -> {sum(nematode collect(each.stomack_C))};
	float pop_copio_R -> {sum (copiotrophe_R collect(each.C))};
	float pop_copio_K -> {sum (copiotrophe_K collect(each.C))};
	float pop_oligo -> {sum (oligotrophe_K collect(each.C))};
	
	float total_bacteries_cytosol_C -> { sum( particle collect( sum(each.populations collect(each.cytosol_C)))) };
	float total_incomeC -> { sum( particle collect( sum(each.populations collect(each.total_perception_C_lost)))) };
	
	float pop_copio_R_perception_C -> {sum (copiotrophe_R collect(each.perception_C))};
	float pop_copio_K_perception_C -> {sum (copiotrophe_K collect(each.perception_C))};
	float pop_oligo_perception_C -> {sum (oligotrophe_K collect(each.perception_C))};
	float total_perception -> { pop_copio_R_perception_C + pop_copio_K_perception_C + pop_oligo_perception_C};
	// #g
		 
	float total_co2_exp_init;
	float total_particle_C_labile_init ;
	float total_particle_C_recalcitrant_init;
	float total_pores_C_available_C_init;
	float total_bacteries_cytosol_C_init;
	float total_incomeC_init;
	float total_nematode_stomackC_init;
	float C_labile_init;
	float C_recalcitrant_init;
	float pop_copio_R_init;
	float pop_copio_K_init;
	float pop_oligo_init;
	float pop_nematode_init;
	float total_enzyme_cl_init;
	float total_enzyme_cr_init;
	
	list<particle> pores_particles -> {particle where(each.my_type = PORE)};
	list<particle> organics_particles -> {particle where(each.my_type = ORGANIC)};
	list<particle> minerals_particles -> {particle where(each.my_type = MINERAL)};
	
	
	float masse_atome_N <- 17.0067;
	float masse_atome_P <- 30.97;
	
	float masse_atomique_N_in_kg <- 2.3258671 * 10^-26#kg;
	float masse_atomique_P_in_kg <- 5.14332 * 10^-26#kg;
		
	float masse_total_N -> {sum(pores_particles collect (each.dam.dim[0])) * masse_atomique_N_in_kg};
	float masse_total_P -> {sum(pores_particles collect (each.dam.dim[1])) * masse_atomique_P_in_kg};
	
	float pop_init_Cop_R;
	float pop_init_Cop_K;
	float pop_init_Oligo;
	
	string species_to_produce <- "riz";
	string type_apport <- "Mada Compost";
	float kilo_of_production <- 0.0;
	
	
	// 5E8 -> 5E9 bacterie / gramme de sol
	/*
	 * A modifier................................................ 
	 */
	float pop_init <-  0.05*1.5#gram/(#cm*#cm)*world.shape.area;
	
	float bulk_density <- 1.17#gram/(#cm*#cm*#cm);
	
	/*
	 * 
	 * 4% de la surface dans les sols sont colonisés
	 * Calculer la carring capacity (une bactérie mesure 1micron) et comparé au 4% de la litterature.
	 */
	float K_init <- 10*pop_init;
	
	map<string,float> soil_characteristics <-   map(["C"::0.02157#gram/(#cm*#cm),"N"::0.00132#gram/(#cm*#cm),"P"::0.00077#gram/(#cm*#cm)]);	
	float step_value <- 1#hour;


 	init
	{
		
		step <- step_value;
		pop_init <- 0.05 * 1.5#gram/(#cm*#cm)*world.shape.area;
		pop_init <- (pop_init / 2);
	
		K_init <- 10*pop_init;
			
		
		ask particle
		{
			my_neighbors_organics <- neighbors where(each.my_type = ORGANIC);
		}
		
		/*Jean*/
		float carbon_in_all_pore <- total_model_weight * carbone_concentration_in_dam; 
		float carbon_in_pore <- carbon_in_all_pore / nb_pores;
		
		float azote_in_all_pore <- total_model_weight * azote_concentration_in_dam;
		float azote_in_pore <- azote_in_all_pore / nb_pores;
		
		float phosphore_in_all_pore <- total_model_weight * phosphore_concentration_in_dam;
		float phosphore_in_pore <- phosphore_in_all_pore / nb_pores;
		/*Jean*/
		
		/*Bernard Et Al 2021 -> Bouzac 2015*/
		float phospore_in_dim <- (1.43#gram * 10^-6) * total_model_weight;
		float phospore_in_dim_in_pore <- phospore_in_dim / nb_pores;
		
		float azote_in_dim <- (4.74#gram * 10^-6) * total_model_weight;
		float azote_in_dim_in_pore <- phospore_in_dim / nb_pores;
		/*Bernard Et Al 2021 -> Bouzac 2015*/
		
		ask pores_particles
		{
			K_C <- K_init / nb_pores;
			particle pp <- self;
			
			C_labile <- C_labile + (pop_init / 2) / nb_pores;
			
			create Dam number:1
			{
				associated_particle <- pp;
				associated_particle.dam <- self;
				dom <- [azote_in_pore, phosphore_in_pore, carbon_in_pore];
				dim <- [azote_in_dim_in_pore, phospore_in_dim_in_pore];
			}
			/**
			 * voir fichier de paramétrisation excel
			 * 		Espece	C	N	P
			 *	Copiotrophe R	1,42857E-05	2,85714E-06	1,42857E-06
			 *	Copiotrophe  K	3,33333E-05	6,66667E-06	3,33333E-06
			 *	Oligotrohphe	0,00002	0,000004	0,000002
			 * 
			 */			
						
			create copiotrophe_R number:1
			{
				porous_cell <- pp;
				location <- any_location_in(pp.shape);
				add self to:pp.populations;
				nb_copioR <- nb_copioR + 1;
				//N/min	P/min	C/min	C/N	C/P
				//0,000001	0,0000005	0,00001	10	20	
			}
			create copiotrophe_K number:1
			{
				porous_cell <- pp;
				location <- any_location_in(pp.shape);
				add self to:pp.populations;
				nb_copioK <- nb_copioK + 1;
			}
			create oligotrophe_K number:1
			{
				porous_cell <- pp;
				location <- any_location_in(pp.shape);
				add self to:pp.populations;
				nb_oligo <- nb_oligo + 1;
			} 
		} 
		
		pop_init_Cop_R <- sum (copiotrophe_R collect(each.C));
		pop_init_Cop_K <- sum (copiotrophe_K collect(each.C));
		pop_init_Oligo <- sum (oligotrophe_K collect(each.C));
		
		create nematode number: nb_nematode
		{
			my_pore <- one_of(pores_particles); 
			location <- my_pore.location;
			nb_nematode <- nb_nematode + 1;
		}
	}
	
	reflex parameter_production
	{
			write("species_to_produce "+species_to_produce);
			write("type_apport "+type_apport);
	}
	reflex scheduleSimulation 
	{
		
		CO2_produit <- 0.0;
		ask Dam
		{
			do decompe;
		}
		ask pores_particles
		{
			do microbe_eat;
		}
		ask organics_particles
		{
			 if(self.C_labile < self.C_labile_init){
			 	
				//do organic_to_pore;
			 }
		}
		if(cycle mod 10 = 0 and cycle > 10){
			do production(species_to_produce); // besoin de faire les échelles avant d'utiliser
		}
		
		if(cycle mod 100 = 0 and cycle > 50 and type_apport != nil)
		{
			do apport_MO(type_apport); // besoin de faire les échelles avant d'utiliser
		}
	}
	
	action production(string name_species)
	{	
		switch(lower_case(name_species))
		{
			// n et p dans le dim
			match "riz"
			{
				do produce_riz_or_mais("riz");
			} 
			match "mais" 
			{
				do produce_riz_or_mais("mais");
			}
			default {
				
				int index <- search_index_species(name_species);
				if(index != -1)
				{
					float N_graine <- get_N(index);
					float P_graine <- get_P(index);
					float N_plante <- N_graine;
					float P_plante <- P_graine;	
					
					float harvest_index <- get_HI(index);
					
					float N_for_graine <- masse_total_N * harvest_index;
					float N_for_plante <- masse_total_N * (1 - harvest_index);
					
					float P_for_graine <- masse_total_P * harvest_index;
					float P_for_plante <- masse_total_P * (1 - harvest_index);
					
					float ratio_dispo_N_graine <- N_for_graine / N_graine;
					float ratio_dispo_P_graine <- P_for_graine / P_graine;
					
					float ratio_dispo_N_plante <- N_for_plante / N_plante;
					float ratio_dispo_P_plante <- N_for_plante / N_plante;
					
					if(lower_case(name_species) = "vesce" or lower_case(name_species) = "arachide " or lower_case(name_species) = "haricot" or lower_case(name_species) = "soja" or lower_case(name_species) = "guianensis" )
					{
						// all N is taken from ambient air
						float min_ratio <- min([ratio_dispo_P_graine,ratio_dispo_P_plante]);
						float qty_of_P_to_retrieve <- min_ratio * (ratio_dispo_P_plante+ratio_dispo_P_graine);
			
						kilo_of_production <- kilo_of_production + min_ratio * 1000#kg;
						
						ask Dam
						{
							dim <- [dim[0], dim[1]-qty_of_P_to_retrieve];
						}
			
						return;
					}
					
					float min_ratio_N <- min([ratio_dispo_N_graine,ratio_dispo_N_plante]);
					float min_ratio_P <- min([ratio_dispo_P_graine,ratio_dispo_P_plante]);
					
					float min_ratio <- min([min_ratio_N,min_ratio_P]);
				
					float qty_of_N_to_retrieve <- min_ratio * (ratio_dispo_N_plante+ratio_dispo_N_graine);
					float qty_of_P_to_retrieve <- min_ratio * (ratio_dispo_P_plante+ratio_dispo_P_graine);
					
			
					kilo_of_production <- kilo_of_production + min_ratio * 1000#kg;
					
					ask Dam
					{
						dim <- [dim[0] - qty_of_N_to_retrieve, dim[0] - qty_of_P_to_retrieve];
					}
					
				}else
				{
					write("species : "+name_species+" not found in file");	
				}
			}
		}
	}
	
	action produce_riz_or_mais(string name_species){
		
			int index_graine;
			int index_plante;
			if(lower_case(name_species) = "riz")
			{
				index_graine <- search_index_species("riz_g");
				index_plante <- search_index_species("riz_p");
				
			}else{
				index_graine <- search_index_species("mais_g");
				index_plante <- search_index_species("mais_p");
			}
			
			float N_graine <- get_N(index_graine);
			float P_graine <- get_P(index_graine);
			
			float N_plante <- get_N(index_plante);
			float P_plante <- get_P(index_plante);
			
			float harvest_index <- get_HI(index_graine);
			
			float N_for_graine <- masse_total_N * harvest_index;
			float N_for_plante <- masse_total_N * (1 - harvest_index);
			
			float P_for_graine <- masse_total_P * harvest_index;
			float P_for_plante <- masse_total_P * (1 - harvest_index);
			
			float ratio_dispo_N_graine <- N_for_graine / N_graine;
			float ratio_dispo_P_graine <- P_for_graine / P_graine;
			
			float ratio_dispo_N_plante <- N_for_plante / N_plante;
			float ratio_dispo_P_plante <- N_for_plante / N_plante;
			
			float min_ratio_N <- min([ratio_dispo_N_graine,ratio_dispo_N_plante]);
			float min_ratio_P <- min([ratio_dispo_P_graine,ratio_dispo_P_plante]);
			
			float min_ratio <- min([min_ratio_N,min_ratio_P]);
			
			float qty_of_N_to_retrieve <- min_ratio * (N_for_graine + N_for_plante);
			float qty_of_P_to_retrieve <- min_ratio * (P_for_graine + P_for_plante);
			
			kilo_of_production <- kilo_of_production + min_ratio * 1000#kg;
			
			ask Dam
			{
				dim <- [dim[0] - qty_of_N_to_retrieve, dim[1] - qty_of_P_to_retrieve];
			}
		
	}
	action apport_MO(string name_MO)
	{
		int index_MO <- search_OM(name_MO);
		if(index_MO != -1)
		{
			
			// 80% soluble -> C_labile 
			// 20% soluble -> DOM
			
			// hemi + celluose -> C_labile_organics
			// lignine -> C_recalcitrant_organics
			
			float taux_labile_soluble <- 0.8;
			
			//10% de carbone
			
			float hectare <- 10000#meter ^2;
			float aire_echelle_1_in_cm2 <- 20#cm ^2; // todo selon l'echelle
			float aire_echelle_1_in_m2 <- aire_echelle_1_in_cm2 / 10000;
			
			float ratio_echelle_modele_realite <- aire_echelle_1_in_m2 / hectare;
			
			write("ratio aire "+ratio_echelle_modele_realite);
			
			float qty_MO_per_hectare <- 12000.0#kg; // a affiner
			float qty_MO_for_model <- qty_MO_per_hectare * ratio_echelle_modele_realite;
			
			write("qty_MO_for_model = "+qty_MO_for_model);
			
			float soluble <- get_soluble_with_index(index_MO);
			float hemi <- get_hemicellulose_with_index(index_MO);
			float celluose <- get_celluose_with_index(index_MO);
			float lignine <- get_lignine_with_index(index_MO);
			float CN <- get_CN_with_index(index_MO);
			float CP <- get_CP_with_index(index_MO);
			
			float soluble_Cl <- ((qty_MO_for_model * soluble / 100)) * taux_labile_soluble;
			float N_pore <- soluble_Cl * CN;
			float P_pore <- soluble_Cl * CP;
			
			float soluble_dom <- ((qty_MO_for_model * soluble / 100)) * (1 - taux_labile_soluble);
			float N_dom <- soluble_dom * CN;
			float P_dom <- soluble_dom * CN;
			
			float labile_organic <- ((qty_MO_for_model * hemi/ 100)) + ((qty_MO_for_model * celluose/ 100));
			float recalcitrant_organic <- ((qty_MO_for_model * lignine/ 100));
			float N_organic <- (labile_organic+recalcitrant_organic) * CN;
			float P_organic <- (labile_organic+recalcitrant_organic) * CN;
			
			ask pores_particles{
				self.C_labile <- self.C_labile + soluble_Cl;
				self.N <- self.N + N_pore;
				self.P <- self.P + P_pore;
				
				self.dam.dom <- [self.dam.dom[0] + (N_dom/nb_pores),self.dam.dom[1] + (P_dom/nb_pores),self.dam.dom[2] + (soluble_dom/nb_pores)];
			}
			ask organics_particles{
				self.C_labile <- self.C_labile + labile_organic;
				self.C_recalcitrant <- self.C_recalcitrant + recalcitrant_organic;
				self.N <- self.N + N_organic;
				self.P <- self.P + P_organic;
			}
			
		}	
	}
	
	reflex heatmap 
	{
		list<float> co2_value;
		
		ask particle
		{
			add sum(populations collect(each.CO2_producted)) to: co2_value;
		}
		
		min_heatmap <- min(co2_value);
		
		if(max(co2_value) > max_heatmap)
		{
			max_heatmap <- max(co2_value);
		}
	}
	
	float process_masseC
	{
		float masseC <- 0.0;
		
		masseC <- masseC + total_pores_C_available_C;
		masseC <- masseC + total_particle_C_labile;
		masseC <- masseC + total_particle_C_recalcitrant;
		masseC <- masseC + total_bacteries_cytosol_C;
		masseC <- masseC + pop_copio_R;
		masseC <- masseC + pop_copio_K;
		masseC <- masseC + pop_oligo;
		
		masseC <- masseC + pop_nematode;
		masseC <- masseC + total_nematode_stomackC;
		
		masseC <- masseC + total_enzyme_cl;
		masseC <- masseC + total_enzyme_cr;
		
		masseC <- masseC + total_co2_exp;
		
		masseC<- masseC + total_perception;
		
		write("---------------------step "+cycle +"----------------------------------");
		
		return masseC;
	}
	
	reflex total_masseC_init_process when: cycle = 1
	{
		 total_masseC_init <- process_masseC();
		 
		 total_co2_exp_init <-  total_co2_exp;
		 total_particle_C_labile_init <- total_particle_C_labile;
		 total_particle_C_recalcitrant_init <- total_particle_C_recalcitrant;
		 total_pores_C_available_C_init <- total_pores_C_available_C;
		 total_bacteries_cytosol_C_init <- total_bacteries_cytosol_C;
		 total_incomeC_init <- total_incomeC;
		 total_nematode_stomackC_init <- total_nematode_stomackC;
		 C_labile_init <- C_labile_init;
		 C_recalcitrant_init <- C_recalcitrant;
		 pop_copio_R_init <- pop_copio_R;
		 pop_copio_K_init <- pop_copio_K;
		 pop_oligo_init <- pop_oligo;
		 pop_nematode_init <- pop_nematode;
		 total_enzyme_cl_init <- total_enzyme_cl;
		 total_enzyme_cr_init <- total_enzyme_cr;	
	}
	
	reflex total_masseC_process when: cycle >= 1
	{
		total_masseC <- process_masseC();
	}
	
}

species soil  schedules:[]//skills:[camisole]
{
	init
	{
		//location <- point(10#cm,10#cm,10#cm);
	}
}


// Dam solution de matière accessible
// dissolved available matter
species Dam 
{	
	//solution organique du sol
	//					N	P	C
	list<float> dom <- [0.0,0.0,0.0];
	
	//solution inorganique du sol
	//					N	 P
	list<float> dim <- [0.0,0.0];
	particle associated_particle;
	float C_N -> {(dom at 2) / ((dom at 0) + (dim at 0))};
	float C_P -> {(dom at 2) / ((dom at 1) + (dim at 1))};
	
	float available_C -> {dom at 2};
	float available_N -> {((dom at 0) + (dim at 0))};
	float available_P -> {((dom at 1) + (dim at 1))};
	
	//---------------------- DEM ------------------------ //
	
	//les Enzymes
	//enzyme carbone labile
	float enzyme_Cl <- 0.0;
	//enzyme carbone récalcitrant
	float enzyme_Cr <- 0.0;
	// enzyme Phosphore
	float enzyme_P <- 0.0;
	// enzyme Azone
	float enzyme_N <- 0.0;
	
	
	//les Fixateurs
	
	// enzyme Acide Carboxilique (CA+)
	float acide_C;
	// H+ pour le ph
	float ph_h;
	
	float process_masseC_dam 
	{
		
		/*float masseC <- 0.0;
		
		masseC <- masseC + total_particle_C_labile;
		masseC <- masseC + total_particle_C_recalcitrant;
		masseC <- masseC + total_pores_C_available_C;
		masseC <- masseC + total_bacteries_cytosol_C;
		 
		masseC <- masseC + pop_copio_R;
		masseC <- masseC + pop_copio_K;
		masseC <- masseC + pop_oligo;
		
		masseC <- masseC + pop_nematode;
		masseC <- masseC + total_nematode_stomackC;
		
		masseC <- masseC + total_enzyme_cl;
		masseC <- masseC + total_enzyme_cr;
		
		masseC <- masseC + total_co2_exp;
		masseC<- masseC + total_perception;
		
		return masseC;*/
		return 0.0;
	}
	
	action inject_enzim(float e_Cl,float e_Cr,float e_P,float e_N)
	{
		enzyme_Cl <- enzyme_Cl + e_Cl;
		enzyme_Cr <- enzyme_Cr + e_Cr;
		
		enzyme_P <- enzyme_P + e_P;
		enzyme_N <- enzyme_N + e_N;
	}
		
	action decompe
	{
		float init <-  process_masseC_dam();
		
		float rate_oligoK <- 0.0;
		float rate_copioK <- 0.0;
		float rate_copioR <- 0.0;
		
		list<particle> organics <- shuffle(associated_particle.my_neighbors_organics + self.associated_particle);
		float add_C <- 0.0;
		float add_P <- 0.0;
		float add_N <- 0.0;
		float taille <- length(organics) + 0.0;
		
		float facteur_cl <- facteur_enzyme_cl; // todo
		float facteur_cr <- facteur_enzyme_cr;
		float facteur_p_n <- 1.0;
		
		//write("recal "+ sum(organics collect (each.C_recalcitrant)));
		//write("labile "+ sum(organics collect (each.C_labile)));
		//write(" length(organics) "+  length(organics));
		
		
		ask organics 
		{	
			float qty_enzyme_Cr_basic <- (myself.enzyme_Cr/taille);
			float carbone_expected_cr <- (qty_enzyme_Cr_basic * facteur_cr);
			
			//write("poids  "+(weight_bacterie_in_pore * ( 0.135 / 24 ) / 0.22#gram) );
			
			//write("diff C_recalcitrant enzyme "+ (C_recalcitrant - qty_enzyme_Cr));
			float dec_Cr <- min([C_recalcitrant, carbone_expected_cr]);
			C_recalcitrant  <- C_recalcitrant - dec_Cr;
			
			float qty_enzyme_Cl_basic <- (myself.enzyme_Cl/taille);
			float carbone_expected_cl <- (qty_enzyme_Cl_basic * facteur_cl);
			
			//write("facteur cl " + carbone_expected_cl );
			//write("diff C_labile enzyme "+ (C_recalcitrant - qty_enzyme_Cl));
			float dec_Cl <- min([C_labile, carbone_expected_cl]);
			C_labile  <- C_labile - dec_Cl;
			
			//add_C <- add_C + dec_Cl+ dec_Cr;
			add_C <- add_C + dec_Cl + dec_Cr + qty_enzyme_Cl_basic + qty_enzyme_Cr_basic;
			//add_C <- add_C + dec_Cl +  dec_Cr + qty_enzyme_Cr + qty_enzyme_Cl;
			
			float qty_enzyme_P <- ((myself.enzyme_P/taille) * facteur_p_n);
			float dec_P <- max([0,min([P,qty_enzyme_P])]);
			P <- P - dec_P;
			add_P <- add_P + dec_P;
			
			float qty_enzyme_N <- ((myself.enzyme_N/taille) * facteur_p_n);
			float dec_N <- max([0,min([N,qty_enzyme_N])]);
			N <- N - dec_N;
			add_N <- add_N + dec_N;
		}
		
		
		dom <- [ dom[0] + add_N, dom[1] + add_P, dom[2] + add_C];
		
		enzyme_Cl <- 0.0;
		enzyme_Cr <- 0.0;
		enzyme_P <- 0.0;
		enzyme_N <- 0.0;
		
		
		float end <- process_masseC_dam();
		
		if(end != init){
			write("decompe = "+ (init - end));
		}
		//assert(end = init);
	}
	
}

grid particle width: grid_width height: grid_height neighbors: 4
{
	string my_type <- one_of([ORGANIC, MINERAL, PORE]); // todo pas dependre de l'aléa mais des chiffres pour la répartition des types
	float C_labile;
	float C_labile_init;
	float C_recalcitrant;
	float C_recalcitrant_init;
	
	float N;
	float P;
	Dam dam;

	float ratio_organics_total <- 0.0; // selon l'echlle changer la facon de répartir les type de cases
	float ratio_pores_total <- 0.0;
	float ratio_minerals_total <- 0.0;

	list<particle> my_neighbors <- neighbors;
	list<particle> my_neighbors_organics <- nil; //ps where(each.my_type = "organic")
	list<particle> my_neighbors_local <- neighbors;
	
	float total_bacteria_C -> {sum(populations collect(each.C))};
	float K_C;
	
	list<microbes_P> populations<- [];
	
	float pop_init ;
	
	init
	{
		//qte carbone labile et recalcitrant par gram volume de sol
		N <- gauss(12, 3)#gram;
		P <- gauss(3, 0.3)#gram;
		
		switch(my_type)
		{
			match MINERAL 
			{ 
					C_labile <- 0.0;
					C_recalcitrant <- 0.0; 
					N <- (soil_characteristics at "N")*self.shape.area; 
					P <- (soil_characteristics at "P")*self.shape.area;  
			}
			match ORGANIC { 
					C_labile_init <- (soil_characteristics at "C"/2)*self.shape.area;
					C_labile <- C_labile_init;
					C_recalcitrant_init <- (soil_characteristics at "C"/2)*self.shape.area;
					C_recalcitrant <- (soil_characteristics at "C"/2)*self.shape.area;
					N <- (soil_characteristics at "N")*self.shape.area; 
					P <- (soil_characteristics at "P")*self.shape.area; 
			}
			default 
			{ 
					C_labile_init <- (soil_characteristics at "C"/2)*self.shape.area;
					C_labile <- C_labile_init;
					C_recalcitrant_init <- (soil_characteristics at "C"/2)*self.shape.area;
					C_recalcitrant <- (soil_characteristics at "C"/2)*self.shape.area;
					N <- (soil_characteristics at "N")*self.shape.area; 
					P <- (soil_characteristics at "P")*self.shape.area;  
			}
			
		}
	//write "start "+ sum(particle collect(each.N))+" "+ sum(particle collect(each.P));
	//write "test "+	(soil_characteristics at "N")*self.shape.area;
	}
		
	/*
	 *  [Croissance, Ratio Carbone/respiration, Qte de matiere ingérée,
	 *  4,5 10-13g masse d'une bactérie (avec eau) (2micrometre) (Masse volumique 0.9 Kg/l)  et 10-14g de carbone 
	 *  K masse Carbone biomasse
	 */
	
	action microbe_eat
	{
		
		float init <- process_masseC_particle();
		
		float total_C_wanted <- sum(populations collect(each.C_wanted));
		float total_N_wanted <- sum(populations collect(each.N_wanted));
		float total_P_wanted <- sum(populations collect(each.P_wanted));
		
		float C_Soil <- min([dam.available_C, total_C_wanted]);
		float N_Soil <- min([dam.available_N, total_N_wanted]);
		float P_Soil <- min([dam.available_P, total_P_wanted]);
		
		
		ask populations 
		{
			float C_rate  <- total_C_wanted > 0 ? (C_wanted / total_C_wanted) : 0;
			float N_rate  <- total_N_wanted > 0 ? (N_wanted / total_N_wanted) : 0;
			float P_rate  <- total_P_wanted > 0 ? (P_wanted / total_P_wanted) : 0;
			
			float C_consum <- C_Soil * C_rate;
			float N_consum <- N_Soil * N_rate;
			float P_consum <- P_Soil * P_rate;
			
			perception_C <- perception_C + C_consum;
			perception_N <- perception_N + N_consum;
			perception_P <- perception_P + P_consum;
			
			C_consum <- min([myself.dam.dom[2],C_consum]);
			P_consum <- min([myself.dam.dom[1],P_consum]);
			N_consum <- min([myself.dam.dom[0],N_consum]);

			myself.dam.dom[2]<- myself.dam.dom[2] - C_consum;			
			myself.dam.dom[1]<- myself.dam.dom[1] - P_consum;			
			myself.dam.dom[0]<- myself.dam.dom[0] - N_consum;
		}
		
		float end <- process_masseC_particle();
		if(end != init)
		{
			write("microbe_eat ="+(init - end));
		}
		//assert(end = init);	
	}
	
	action organic_to_pore
	{
		write("organics to pore : "+nb_pores);
				
		my_type <- PORE;
		K_C <- K_init / nb_pores;
		add self to: pores_particles;
		remove self from: organics_particles;
		remove self from: my_neighbors_organics;
		ask particle
		{
			my_neighbors_organics <- neighbors where(each.my_type = ORGANIC);
		}
		
		particle pp <- self;
		
		create Dam number:1
		{
			associated_particle <- pp;
			associated_particle.dam <- self;
			dom <- [myself.N, myself.P, myself.C_recalcitrant + myself.C_labile]; // todo qty of carbone restant + n + p
			dim <- [0.0, 0.0];
			myself.C_recalcitrant <- 0.0;
			myself.C_labile <- 0.0;
		}
		
		list<particle> my_neighbors_pores <- neighbors where(each.my_type = PORE);
		int nb_adjacent_pore <- length(my_neighbors_pores);
		
		float qty_C_neighbors <- 0.0;
		float qty_N_neighbors <- 0.0; 
		float qty_P_neighbors <- 0.0;
		
		/*ask my_neighbors_pores{ // todo fix C conservation
			ask populations {
				float tmp_C <- (C / (nb_adjacent_pore+1)) / 3;
				C <- C - tmp_C;	
				qty_C_neighbors <- qty_C_neighbors + tmp_C;
				
				float tmp_N <- (N / (nb_adjacent_pore+1)) / 3;
				N <- N - tmp_N;	
				qty_N_neighbors <- qty_C_neighbors + tmp_N;
				
				float tmp_P <- (P / (nb_adjacent_pore+1)) / 3;
				P <- P - tmp_P;	
				qty_P_neighbors <- qty_P_neighbors + tmp_P;
			}
		}
			
			
		create copiotrophe_R number:1
		{
			C <- qty_C_neighbors / 3;
			N <- qty_N_neighbors / 3;
			P <- qty_P_neighbors / 3;
			
			porous_cell <- myself;
			location <- any_location_in(pp.shape);
			add self to: myself.populations;
			nb_copioR <- nb_copioR + 1;
		}
		create copiotrophe_K number:1
		{
			C <- qty_C_neighbors / 3;
			N <- qty_N_neighbors / 3;
			P <- qty_P_neighbors / 3;
			
			porous_cell <- myself;
			location <- any_location_in(pp.shape);
			add self to: myself.populations;
			nb_copioK <- nb_copioK + 1;
		}
		create oligotrophe_K number:1
		{
			C <- qty_C_neighbors / 3;
			N <- qty_N_neighbors / 3;
			P <- qty_P_neighbors / 3;
			
			porous_cell <- myself;
			location <- any_location_in(pp.shape);
			add self to: myself.populations;
			nb_oligo <- nb_oligo + 1;
		}*/
	}
	
	float process_masseC_particle
	{
		
		/*float masseC <- 0.0;
		
		masseC <- masseC + total_pores_C_available_C;
		masseC <- masseC + total_particle_C_labile;
		masseC <- masseC + total_particle_C_recalcitrant;
		masseC <- masseC + total_bacteries_cytosol_C;
		masseC <- masseC + pop_copio_R;
		masseC <- masseC + pop_copio_K;
		masseC <- masseC + pop_oligo;
		
		masseC <- masseC + pop_nematode;
		masseC <- masseC + total_nematode_stomackC;
		
		masseC <- masseC + total_enzyme_cl;
		masseC <- masseC + total_enzyme_cr;
		
		masseC <- masseC + total_co2_exp;
		masseC<- masseC + total_perception;
		
		return masseC;	*/
		return 0.0;
	}
	
	rgb selectColor
	{
		switch(my_type)
		{
			match "pore" { return #black;}
			match "mineral" { return #yellow;}
			match "organic" { return #green;}
			default { return #blue;}
		}
	}
	aspect heatmap
	{
		rgb heat <- #red;
		
		switch(my_type)
		{
			match MINERAL {heat <- #yellow;}
			match ORGANIC {heat <- #green;}
			match PORE 
			{
				float current <- sum(populations collect(each.CO2_producted));
				float normalized_value <- (255-0)* ((current-min_heatmap)/(max_heatmap-min_heatmap)) +0;
				heat <- rgb(normalized_value as int, 0, 0);
			}	
		}
		
		draw shape color: heat;
	}
	
	aspect default
	{
		rgb mcolor <- #blue;
		switch(my_type)
		{
			match "pore" { mcolor <- #black;}
			match "mineral" { mcolor <- #yellow;}
			match "organic" { mcolor <- #green;}
		}
		draw shape color:mcolor;
	}
	
}


species microbes_P 
{
	float threshold_DOC <- 0.5; // seuil de dormance
	float taux_respiration <- 0.5;
	float dividing_time <- 0.0;
	//taux labile récalcitrant
	float L_R_rate; // labile_recalcitrante_rate
	float C_N;
	float C_P;
	bool awake <- false;
	
	// quantity of C N P in the colony
	float C <- 60.0;
	float N <- 7.0;
	float P <- 1.0;
	
	float facteur_dormance <- 0.5;
	float wakeup_factor <- 0.50;
	float C_actif -> {C * facteur_dormance};
	
	float C_init <- C;
	// Mon estomac le cytosol 
	float cytosol_C <- 0#gram; 
	float cytosol_N <- 0#gram; 
	float cytosol_P <- 0#gram; 
	
	// ce que je percois
	float perception_C <- 0#gram; 
	float perception_N <- 0#gram;
	float perception_P <- 0#gram;
	
	float total_perception_C_lost <- 0.0;


	string name_species;
	/* 
	 R : C * (step/1)  /  (1-0.2)
	 K : C * (step/24)  / (1-0.4)
	 O : C * (step/368) / (1-0.7)s
	 */
	float C_assimilation_speed -> {C_actif * (step / dividing_time) / (1 - taux_respiration)};
	
	
	float CO2_producted <- 0.0;
	
	particle porous_cell;
	
	//float decomposing_speed <- 0.00001#gram/#mn;
	float C_wanted -> { C_assimilation_speed };
	float N_wanted -> { (C+cytosol_C+C_wanted) / C_N - (cytosol_N + N)};
	float P_wanted -> { (C+cytosol_C+C_wanted) / C_P - (cytosol_P + P)};
	float DOC_assimilation_friction;
	float death_rate <- 0.0;
	
	float division_enzyme_rate <- 0.0;
	float cytosol_enzyme_facteur <- 0.0;
	
	float process_masseC_microbe 
	{
		
		/*float masseC <- 0.0;
		
		masseC <- masseC + total_particle_C_labile;
		masseC <- masseC + total_particle_C_recalcitrant;
		masseC <- masseC + total_pores_C_available_C;
		masseC <- masseC + total_bacteries_cytosol_C;
		 
		masseC <- masseC + pop_copio_R;
		masseC <- masseC + pop_copio_K;
		masseC <- masseC + pop_oligo;
		
		masseC <- masseC + pop_nematode;
		masseC <- masseC + total_nematode_stomackC;
		
		masseC <- masseC + total_enzyme_cl;
		masseC <- masseC + total_enzyme_cr;
		
		masseC <- masseC + total_co2_exp;
		masseC<- masseC + total_perception;
		
		return masseC;*/
		return 0.0;
	}
	
	action respirate(float cytosol_respiration) // pas de perte c
	{
		float init <- process_masseC_microbe();
		cytosol_C <- cytosol_C  - cytosol_respiration; 
		
		CO2_producted <- CO2_producted + cytosol_respiration;
		total_co2_exp <- total_co2_exp + cytosol_respiration;
		
		float end <- process_masseC_microbe();
		
		assert (init = end);
		if(end != init)
		{
			write("respirate = "+(init-end));
		}
	}
	
	action growth(float cytosol_division) // pas de perte c
	{
		
		float init <- process_masseC_microbe();
		
		// quantity de carbone cytosol -> structurel
		float new_individual <- C * ( step / dividing_time) * (1- self.porous_cell.total_bacteria_C / self.porous_cell.K_C);
		
		//new_individual <- min([new_individual , cytosol_C]);
		new_individual <- min([new_individual , cytosol_division]);
		
		
		float assi_C_N <- cytosol_N * C_N;  //(assimilated_N=0) ? 0:(assimilated_C/assimilated_N);
		float assi_C_P <- cytosol_P * C_P; //(assimilated_P=0) ? 0:(assimilated_C/assimilated_P);
		
		new_individual <- min([new_individual, (assi_C_N) ]);
		new_individual <- min([new_individual, (assi_C_P) ]);
		
		float transfert_C <- new_individual;
		float transfert_N <- transfert_C / C_N;
		float transfert_P <- transfert_C / C_P;
		
		cytosol_C <- cytosol_C - transfert_C;
		cytosol_N <- cytosol_N - transfert_N;
		cytosol_P <- cytosol_P - transfert_P;
		
		
		C <- C + transfert_C;
		N <- N + transfert_N;
		P <- P + transfert_P;
		
		
		float end <- process_masseC_microbe();
		
		assert (init = end);
		if(init != end)
		{
			write("growth = "+(init-end));
		}
		
	}
	
	action death
	{
		// A partir d'un seuil d'augmentation de population, le virus se déclanche
		
		// plus la population crois vite, plus le virus est actif
		
		// un taux de mortalité par espèce....
	}
	
	float decompose(float cytosol_enzyme) // perte c
	{
		float init <- process_masseC_microbe();
		
		cytosol_C <- cytosol_C - cytosol_enzyme;
		
		float wanted_P <- cytosol_enzyme / C_P;
		float wanted_N <- cytosol_enzyme / C_N;
		
		ask porous_cell.dam
		{
			do inject_enzim(myself.L_R_rate*cytosol_enzyme, (1-myself.L_R_rate)*cytosol_enzyme, wanted_P, wanted_N);
		}
			
		float end <- process_masseC_microbe();
		assert(init = end);
		
		if(init != end)
		{
			write("decompose = "+(init-end));
		}
		
		return 0.0;
	}
	
	reflex life
	{
		
		float init <- process_masseC_microbe();
		CO2_producted <- 0.0;
		
		if(C_wanted = 0)
		{
			facteur_dormance <- 0.0;
			
		}else
		{
			facteur_dormance <- min([(perception_C / C_wanted), 1]);
		}
		
		facteur_dormance <- max([facteur_dormance,wakeup_factor]);
		
	
		cytosol_C <- cytosol_C + perception_C * facteur_dormance;
		cytosol_N <- cytosol_N + perception_N * facteur_dormance;
		cytosol_P <- cytosol_P + perception_P * facteur_dormance;
		
		float init_cytosol_C <- cytosol_C;
		
		float cytosol_respiration <- (taux_respiration) * cytosol_C ;
		float cytosol_division <- ((1-taux_respiration) * cytosol_C) * division_enzyme_rate;
		float cytosol_enzyme <- ((1-taux_respiration) * cytosol_C) * (1 - division_enzyme_rate);
		
		float weight_at_T <- (C_actif * weight_bacterie_in_one_pore) / C_init;
		float qty_C_for_population <- weight_at_T * (activite_enzymatique_per_hour) / masse_microbienne_etude; // qty C / per hour
		
		//write("masse micro "+qty_C_for_population);
		
		cytosol_enzyme_facteur <- cytosol_enzyme;
		
		//assert(init_cytosol_C = (cytosol_respiration + cytosol_division + cytosol_enzyme) );
		if( init_cytosol_C != (cytosol_respiration + cytosol_division + cytosol_enzyme )){
			//write(" "+init_cytosol_C+" ?? "+ (cytosol_respiration + cytosol_division + cytosol_enzyme) );
		}
		
		
		do respirate(cytosol_respiration);
		do growth(cytosol_division);
		do decompose(cytosol_enzyme);
		
	
		ask porous_cell
		{
			dam.dom[2]<- dam.dom[2] + myself.perception_C * (1 - myself.facteur_dormance);
			dam.dom[1]<- dam.dom[1] + myself.perception_P * (1 - myself.facteur_dormance);
			dam.dom[0]<- dam.dom[0] + myself.perception_N * (1 - myself.facteur_dormance);
			
		}
		
		
		total_perception_C_lost <- total_perception_C_lost + perception_C;
		
		perception_C <- 0.0;
		perception_N <- 0.0;
		perception_P <- 0.0;
		
		do death;
		
		float end <- process_masseC_microbe();
		assert(end = init);
		
		if(end-init != 0.0)
		{
			write("life = "+( init - end));
		}
	}
}

species copiotrophe_R parent:microbes_P
{
	init
	{
		dividing_time <- dividing_time_copiotrophe_R;
		L_R_rate <- 1.0;
		C_N <- 5.0;
		C_P <- 10.0;
		C <- pop_init/nb_pores/3;
		C_init <- C;
		cytosol_C <- C;
		cytosol_N <- C / C_N;
		cytosol_P <- C / C_P;
		
		wakeup_factor <- 0.002;
		
		N <- C / C_N;
		P <- C / C_P;
		//threshold_DOC_act <- 0.2;
		threshold_DOC <- 0.18;
		taux_respiration <- 0.2;
		death_rate <- 0.0001;
		name_species <- "copiotrophe_R";
		division_enzyme_rate <- 0.714;	
	}	
	aspect base
	{
		draw circle(0.01#cm) color:#blue;
	}
}
species copiotrophe_K parent:microbes_P
{
	init{
		dividing_time <- dividing_time_copiotrophe_K;
		L_R_rate <- 0.2;
		C_N <- 5.0;
		C_P <- 10.0;
		C <- pop_init/nb_pores/3; //*0.2;//0.000001#gram;
		C_init <- C;
		cytosol_C <- C;
		cytosol_N <- C / C_N;
		cytosol_P <- C / C_P;
		
		wakeup_factor <- 0.018;
		 
		N <- C/ C_N;
		P <- C/ C_P;
		threshold_DOC <- 0.05;
		taux_respiration <- 0.4;
		death_rate <- 0.00001;
		name_species <- "copiotrophe_K";
		division_enzyme_rate <- 0.5;		
	}	
	aspect base
	{
		draw circle(0.01#cm) color:#red;
	}
}

species oligotrophe_K parent:microbes_P
{
	init
	{
		dividing_time <- dividing_time_oligotrophe;
		L_R_rate <- 0.8;
		C_N <- 5.0;
		C_P <- 10.0;
		C <- pop_init/nb_pores/3;   //*0.7;//0.000001#gram;
		
		C_init <- C;
		cytosol_C <- C;
		cytosol_N <- C / C_N;
		cytosol_P <- C / C_P;
		
		wakeup_factor <- 1.0;
		
		N <- C/ C_N;
		P <- C/ C_P;
		threshold_DOC <- 0.0;
		taux_respiration <- 0.7;
		death_rate <- 0.0;
		name_species <- "oligotrophe_K";
		division_enzyme_rate <- 0.5;	
	}
	aspect base
	{
		draw circle(0.01#cm) color:#pink;
	}
}


species nematode
{
	float threshold_DOC <- 0.5; // à vérifier
	
	float C_N <- 6.0;
	float C_P <- 47.0;
	float CO2_efficiency <- 0.84;
	float t1 <- 0.5; // à vérifier
	
	// quantity of C N P
	float C <- (470#gram * 10^-6);
	float N <- (470#gram * 10^-6)/C_N;
	float P <- (470#gram * 10^-6)/C_P;
	
	float stomack_C;
	float stomack_N;
	float stomack_P;
	
	float CO2_production <- 0.0 ;
	particle my_pore <- nil;
	
	bool awake <- true;
	float predation_rate <- 0.00001#gram/#day; // 0.5E-9#gram/#mn; 
	
	//un nematode mange entre 10K et 100k bacteries par jour -> moins de 10k il dort (ref article de colman)
	float wanted_C <- predation_rate * step;
	
	//perception of the nematode
	float available_C <- 0.0;
	float available_N <- 0.0; 
	float available_P <- 0.0; 
	
	reflex life
	{
		do perceive;
		awake <- (available_C / wanted_C) >= threshold_DOC;
		
		if(awake)
		{
			do move;
			do predating;
			do respirate;	
			do anabolize;
			
		}else
		{
			//write("sleeping");
		}
	}
	
	action perceive
	{
		available_C <- sum(my_pore.populations collect(each.C + each.cytosol_C));
		available_N <- sum(my_pore.populations collect(each.N + each.cytosol_N));
		available_P <- sum(my_pore.populations collect(each.P + each.cytosol_P));
	}
	
	action move
	{
		list<particle> neighbors_pores <- (my_pore.neighbors where (each.my_type = PORE)) + my_pore;
		
		float max_perceive <- sum(neighbors_pores[0].populations collect(each.C + each.cytosol_C));
		particle most_attractive <- neighbors_pores[0];
		ask neighbors_pores
		{
			float tmp <- sum(self.populations collect(each.C + each.cytosol_C));
			if(tmp >= max_perceive)
			{
				max_perceive <- tmp;
				most_attractive <- self;
			}
		}
		my_pore <- most_attractive;
		location <- my_pore.location;
	}
	
	action respirate
	{
		float carbone_respiration <- CO2_efficiency * stomack_C;
		stomack_C <- stomack_C  - carbone_respiration; 
		
		CO2_production <- CO2_production + carbone_respiration;
		total_co2_exp <- total_co2_exp + carbone_respiration;
	}
	
	action anabolize
	{
		float carbone_assimilation <- (1 - CO2_efficiency) * stomack_C;
		stomack_C <- stomack_C  - carbone_assimilation; 
		
		float qty_N_necromasse <- carbone_assimilation/C_N;
		float qty_P_necromasse <- carbone_assimilation/C_P;
		
		// necromasse nematode
		my_pore.C_labile <- my_pore.C_labile + carbone_assimilation;
		my_pore.N <- my_pore.N + qty_N_necromasse;
		my_pore.P <- my_pore.P + qty_P_necromasse;
		
		stomack_N <- stomack_N - qty_N_necromasse;
		stomack_P <- stomack_P - qty_P_necromasse;
		
		// le reste de l'estomac va dans la dim
		my_pore.dam.dim[0] <- my_pore.dam.dim[0] + stomack_N;
		my_pore.dam.dim[1] <- my_pore.dam.dim[1] + stomack_P;
		
		stomack_N <- 0.0;
		stomack_P <- 0.0;		
	}
	
	action recracher(float carbone, float azote, float phosphore)
	{
		my_pore.C_labile <- my_pore.C_labile + carbone;
		my_pore.N <- my_pore.N + azote;
		my_pore.P <- my_pore.P + phosphore;
	}
	
 	action predating
 	{
		float C_to_catch <- min([available_C, wanted_C]);
		float C_acc <-  C_to_catch;
		
		float total_C <- 0.0;
		float total_N <- 0.0;
		float total_P <- 0.0;
		
		ask shuffle(my_pore.populations)
		{	
			float localCC <- (C + cytosol_C) / myself.available_C *C_to_catch;
			float local_C_to_catch <- min([localCC, C_acc]);
			
			float localC;
			//élément dans la partie structure
			if(C+cytosol_C = 0){
				localC <- 0.0;
			}else{
				localC <- C/(C+cytosol_C) *localCC;
			}
			float localN <- min([N, localC / C_N]);
			float localP <- min([P, localC / C_P]);
			
			//élément dans la partie cytosol
			float localCytosol <- localCC - localC ;
			float local_cyt_N <- min([(localCytosol / C_N), cytosol_N ]);
			float local_cyt_P <- min([(localCytosol / C_P), cytosol_P]);

			C <- C - localC ;
			N <- N - localN;
			P <- P - localP;
			
			cytosol_C <- cytosol_C  - localCytosol;
			cytosol_N <- cytosol_N  - local_cyt_N;
			cytosol_P <- cytosol_P  - local_cyt_P;
			
			total_C <- total_C + (localCytosol + localC);
			total_N <- total_N + (local_cyt_N +  localN);
			total_P <- total_P + (local_cyt_P +  localP);
			
		}
		
		do recracher(total_C*t1,total_N*t1,total_P*t1);
		
		stomack_C <- stomack_C + ((1-t1) * total_C);
		stomack_N <- stomack_N + ((1-t1) * total_N);
		stomack_P <- stomack_P + ((1-t1) * total_P);
 	}
 	
 	
	aspect base
	{
		draw square(0.05#cm) color:#white;
	}
}

experiment ex type:gui
{
	// todo parameter
	
	/*
	parameter "dividing_time_copiotrophe_R" var: dividing_time_copiotrophe_R <- 1#h;
	parameter "dividing_time_copiotrophe_K" var: dividing_time_copiotrophe_K <- 24#h;
	parameter "dividing_time_oligotrophe" var: dividing_time_oligotrophe <- 368#h;
	
	parameter "growth_decomposeur" var: growth_decomposeur <- 0.5/#mn;
	
	parameter "carbone_concentration_in_dam" var:carbone_concentration_in_dam <- (729.0#gram * 10^-6)/#gram;
	parameter "azote_concentration_in_dam" var:azote_concentration_in_dam <- (60.0#gram * 10^-6)/#gram;
	parameter "phosphore_concentration_in_dam" var:phosphore_concentration_in_dam <- (400.0#gram * 10^-6)/#gram;
	
	parameter "nb_nematode" var: nb_nematode <- 30;
	
	parameter "facteur_enzyme_cr" var: facteur_enzyme_cr <- 1.0;
	parameter "facteur_enzyme_cl" var: facteur_enzyme_cl <- 200.0;
	
	//parameter "pop_init" var: pop_init <-  0.05*1.5#gram/(#cm*#cm)*world.shape.area;
	
	parameter "bulk_density" var: bulk_density <- 1.17#gram/(#cm*#cm*#cm);
	parameter "K_init" var: K_init <- 10*pop_init;
	
	parameter "plantation" var: species_to_produce <- "riz";
	parameter "compost" var: type_apport <- "Mada Compost";
	
	parameter "soil_characteristics" var: soil_characteristics <-   map(["C"::0.02157#gram/(#cm*#cm),"N"::0.00132#gram/(#cm*#cm),"P"::0.00077#gram/(#cm*#cm)]);	
	parameter "step" var: step <- 1#hour;
	*/
	
	
	
	//Qualité des apports et azote et phosphore
	output
	{
		display map  background:#black
		{
			species particle ;
		/* 	species copiotrophe_R aspect:base;
			species copiotrophe_K aspect:base;
			species oligotrophe_K aspect:base;*/
		}
		
		/*display "Population" type: java2D 
		{
			chart "Population" type: series 
			{
				data "Copiotrophe R" value: pop_copio_R   - pop_init_Cop_R color: #blue ;// style: spline;
				data "Copiotrophe K" value: pop_copio_K - pop_init_Cop_K color: #red;// style: spline;
				data "oligotrophe" value: pop_oligo - pop_init_Oligo color: #green ;//style: spline;
			}
		}
		display "N and P in" type: java2D 
		{
			chart "N and P" type: series
			{
				data "N" value: sum(pores collect(each.dam.dim[0])) color: #blue style: spline;
				data "P" value: sum(pores collect(each.dam.dim[1])) color: #red style: spline;
			}
	
		}
		display "N and P org dis" type: java2D 
		{
			chart "N and P org dis" type: series
			{
				data "N" value: sum(pores collect(each.dam.dom[0])) color: #blue style: spline;
				data "P" value: sum(pores collect(each.dam.dom[1])) color: #red style: spline;
			}
	
		}
		display "N and P org" type: java2D 
		{
			chart "N and P in" type: series
			{
				data "N" value: sum(particle collect(each.N)) color: #blue style: spline;
				data "P" value: sum(particle collect(each.P)) color: #red style: spline;
			}
	
		}*/
		display "Carboon in soil" type: java2D 
		{
			chart "Carboon in soil" type: series
			{
				data "C Labile" value: sum(particle collect(each.C_labile)) color: #blue style: spline;
				data "C Recalcitrant" value: sum(particle collect(each.C_recalcitrant)) color: #orange style: spline;
				data "C in DOM" value: sum(particle where (each.my_type = PORE) collect(each.dam.dom[2])) color: #red style: spline;
			}
		}
		display "C in pore" type: java2D 
		{
			chart "C in pore" type: series
			{
				data "C_labile" value: sum(pores_particles collect(each.C_labile)) color: #blue style: spline;
				data "C_recalcitrant" value: sum(pores_particles collect(each.C_recalcitrant)) color: #red style: spline;
			}
	
		}
		display "C in organics" type: java2D 
		{
			chart "C in organics" type: series
			{
				data "C_labile" value: sum(organics_particles collect(each.C_labile)) color: #blue style: spline;
				data "C_recalcitrant" value: sum(organics_particles collect(each.C_recalcitrant)) color: #red style: spline;
			}
	
		}
		display "CO2" type: java2D 
		{
			chart "C02" type: series
			{
				data "CO2 copiotrophe_R" value: sum(copiotrophe_R collect(each.CO2_producted)) color: #blue style: spline;
				data "CO2 copiotrophe_K" value: sum(copiotrophe_K collect(each.CO2_producted)) color: #red style: spline;
				data "CO2 oligotophe" value: sum(oligotrophe_K collect(each.CO2_producted)) color: #green style: spline;
				data "CO2 nematode" value: sum(nematode collect(each.CO2_production)) color: #purple style: spline;
			}
		
		}
		display "awake population" type: java2D 
		{
			chart "awake population" type: series
			{
				data "CopioR awake %" value: (sum(copiotrophe_R collect (each.facteur_dormance)) / nb_copioR) * 100;
				data "CopioK awake %" value: (sum(copiotrophe_K collect (each.facteur_dormance)) / nb_copioK) * 100;
				data "Oligo awake %" value: (sum(oligotrophe_K collect (each.facteur_dormance)) / nb_oligo) * 100;
				data "Nematode awake %" value: (sum(nematode collect (each.awake as int)) / nb_nematode) * 100;
			}
		}
		display "microbes colony" type: java2D
		{
			chart "microbes colony" type: series
			{
				data "C colony" value: sum(copiotrophe_R collect (each.C)) + sum(copiotrophe_K collect (each.C)) + sum(oligotrophe_K collect (each.C));
				data "N colony" value: sum(copiotrophe_R collect (each.N)) + sum(copiotrophe_K collect (each.N)) + sum(oligotrophe_K collect (each.N));
				data "P colony" value: sum(copiotrophe_R collect (each.P)) + sum(copiotrophe_K collect (each.P)) + sum(oligotrophe_K collect (each.P));	
			}
		}
		display "Mass C conservation "
		{
			chart "Mass C conservation"
			{	
				data "diff mass initial" value: total_masseC_init - total_masseC color: #red;	
			}
		}
		display "heat map co2"
		{
			species particle aspect: heatmap;
	
		}
		
		monitor "diff mass initial" value: total_masseC_init - total_masseC color: #red;
			
	}
}



experiment calcul_facteur_enzyme type:gui
{
	parameter "nb_nematode" var: nb_nematode <- 22;
	parameter "facteur_enzyme_cr" var: facteur_enzyme_cr <- 1.0;
	parameter "facteur_enzyme_cl" var: facteur_enzyme_cl <- 200.0;
	parameter "step_value" var: step_value <- 1#hour; 
	parameter "grid_width" var: grid_width <- 12; 
	parameter "grid_height" var: grid_height <- 12; 
	
	output
	{
		display map  background:#red
		{
			species particle;
			species copiotrophe_R aspect:base;
			species copiotrophe_K aspect:base;
			species oligotrophe_K aspect:base;
			species nematode aspect:base;
		}
		display "CO2" type: java2D 
		{
			chart "C02" type: series
			{
				data "CO2 copiotrophe_R" value: sum(copiotrophe_R collect(each.CO2_producted)) color: #blue style: spline;
				data "CO2 copiotrophe_K" value: sum(copiotrophe_K collect(each.CO2_producted)) color: #red style: spline;
				data "CO2 oligotophe" value: sum(oligotrophe_K collect(each.CO2_producted)) color: #green style: spline;
			}
		}
		display "C in pore" type: java2D 
		{
			chart "C in pore" type: series
			{
				data "C_labile" value: sum(pores_particles collect(each.C_labile)) color: #blue style: spline;
				data "C_recalcitrant" value: sum(pores_particles collect(each.C_recalcitrant)) color: #red style: spline;
			}
	
		}
		display "C in organics" type: java2D 
		{
			chart "C in organics" type: series
			{
				data "C_labile" value: sum(organics_particles collect(each.C_labile)) color: #blue style: spline;
				data "C_recalcitrant" value: sum(organics_particles collect(each.C_recalcitrant)) color: #red style: spline;
			}
	
		}
		display "Carboon in soil" type: java2D 
		{
			chart "Carboon in soil" type: series
			{
				data "C Labile" value: sum(particle collect(each.C_labile)) color: #blue style: spline;
				data "C Recalcitrant" value: sum(particle collect(each.C_recalcitrant)) color: #orange style: spline;
				data "C in DOM" value: sum(particle where (each.my_type = PORE) collect(each.dam.dom[2])) color: #red style: spline;
			}
		}
		display "awake population" type: java2D 
		{
			chart "awake population" type: series
			{
				data "CopioR awake %" value: (sum(copiotrophe_R collect (each.facteur_dormance)) / nb_copioR) * 100;
				data "CopioK awake %" value: (sum(copiotrophe_K collect (each.facteur_dormance)) / nb_copioK) * 100;
				data "Oligo awake %" value: (sum(oligotrophe_K collect (each.facteur_dormance)) / nb_oligo) * 100;
				data "Nematode awake %" value: (sum(nematode collect (each.awake as int)) / nb_nematode) * 100;
			}
	
		}
		display "microbes colony" type: java2D
		{
			chart "microbes colony" type: series
			{
				
				data "C colony" value: sum(copiotrophe_R collect (each.C + each.cytosol_C)) + sum(copiotrophe_K collect (each.C+ each.cytosol_C)) + sum(oligotrophe_K collect (each.C+ each.cytosol_C));
				data "N colony" value: sum(copiotrophe_R collect (each.N + each.cytosol_N)) + sum(copiotrophe_K collect (each.N+ each.cytosol_N)) + sum(oligotrophe_K collect (each.N+ each.cytosol_N));
				data "P colony" value: sum(copiotrophe_R collect (each.P + each.cytosol_P)) + sum(copiotrophe_K collect (each.P+ each.cytosol_P)) + sum(oligotrophe_K collect (each.P+ each.cytosol_P));
					
			}
		}
		/*display "nematode colony" type: java2D
		{
			chart "nematode colony" type: series
			{
				
				data "C" value: sum(nematode collect (each.C + each.stomack_C));
				data "N" value: sum(nematode collect (each.N + each.stomack_N));
				data "P" value: sum(nematode collect (each.P + each.stomack_P));		
			}
		}*/
		display "Mass C conservation "
		{
			chart "Mass C conservation"
			{	
				data "diff mass initial" value: total_masseC_init - total_masseC color: #red;	
			}
		}
		/*display "qualité Sol"
		{
			chart "porosité"
			{
				data "porosité" value: nb_pores/(nb_pores+nb_organics+nb_minerals);
			}
		}
		
		display "sortie Sol"
		{
			chart "production"
			{
				data "production "+species_to_produce value: kilo_of_production;
			}
		}*/
		monitor "diff mass initial" value: total_masseC_init - total_masseC color: #red;
		//monitor "production in kg" value: kilo_of_production color: #blue;
	}
}
