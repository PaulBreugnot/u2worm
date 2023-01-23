/**
* Name: csvfunction
*  
* Author: Lucas GROSJEAN
* 
* Library function for csv file used in camisole model
*/


model csv_loader

global
{
	// Fertilizer fields headers
	string CSV_QUANTITE_SOLUBLE <- "Solubles";
	string CSV_HEMICELLULOSE <- "Hemicellulose";
	string CSV_CELLULOSE <- "Cellulose";
	string CSV_LIGNINE <- "Lignine";
	string CSV_C_N <- "C:N";
	string CSV_C_P <- "C:P";
	string CSV_DOSE_ESSAI <- "Dose d'essai (kg/ha)";
	
	// Crop fields headers
	string CSV_N <- "N(kg/tonne)";
	string CSV_P <- "P(kg/tonne)";
	string CSV_HARVEST_INDEX <- "Harvest index";
	string CSV_N_FROM_SOIL <- "N from soil";
	
	
	map<string, map<string, float>> fertilizers_data;
	map<string, map<string, float>> crops_data;
	
	action load_csv_into_map(csv_file csv, map<string, map<string, float>> data) {
		list<string> headers <- csv.contents row_at 0;
		loop i from: 1 to: int(csv.contents.dimension.y)-1 {
			list<unknown> row <- csv.contents row_at i;
			data[string(row[0])] <- [];
			loop j from: 1 to: int(csv.contents.dimension.x)-1 {
				data[string(row[0])][headers[j]] <- float(row[j]);
			}
		}
	}
	
	action load_csv_data {
		write "Loading fertilizers CSV data...";
		csv_file fertilizer_csv <- csv_file("../includes/qualite_amendements_organiques.csv", false); // Do not handle header
		do load_csv_into_map(fertilizer_csv, fertilizers_data);
		write fertilizers_data;
		
		write "Loading crops CSV data...";
		csv_file crops_csv <- csv_file("../includes/exportation_NPK_cultures.csv", false);
		do load_csv_into_map(crops_csv, crops_data);
		write crops_data;
	}
}

experiment debug_csv_loader type:gui {
	action _init_ {
		create simulation {
			ask world {
				do load_csv_data;
			}
		}
	}
	output {
		display test {
			
		}
	}
}