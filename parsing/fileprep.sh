#Remove malaspina from abundance
grep -v Malaspina Normalized_Viral_Abundances_ALL_5kb.txt > Normalized_Viral_Abundances_ALL_5kb_noMalaspina.txt

#Jens nitrate_levels file with nitrates
perl -pe 's/TARA_0+/Station/g' tara_levels.csv | perl -pe 's/TARA_(\d{3})/Station$1/g' | perl -pe 's/SRF/SUR/g' > nitrate_levels2.csv;

#Daniels station environmental file with depths
cut -f3,8,9 -d "," station_environmental_data.csv | perl -pe 's/TARA_0+/Station/g' | perl -pe 's/TARA_(\d{3})/Station$1/g' | perl -pe 's/\[(.*)\].*(,.*)/$1$2/g' | perl -pe 's/SRF/SUR/g' > station_environmental_data2.csv 

#Run station_zone_means.R #to get lat, lon, depth, and nitrate means per station and zone
#Run all_info2.R #To get a single file with all info
