require(reshape2)
#retrieve info
contig_info = read.csv("genome_covariates.csv")
abund = read.table("Normalized_Viral_Abundances_ALL_5kb_noMalaspina.txt",header = T)
station_info = read.csv("station_zone_means.csv")

#Check that contigs and abundance have same info
seq_contigs = unique(sort(contig_info$ID))
abund_contigs = unique(sort(abund$Contig))
setdiff(abund_contigs,seq_contigs) #they are the same

#order by contig name
contig_info2 = contig_info[order(contig_info$ID),]
abund2 = abund[order(abund$Contig),1:132] #Remove MPS

#melt abundance matrix
abund2_colnames = colnames(abund2)[2:132]
all_table = data.frame(abund2,
                       contig_info2[,c(2:7,10:12)])

all_table_melted = melt(all_table, measure.vars = abund2_colnames,
                  variable.name = "station_abundance",
                  value.name="abundance")

all_table_melted = data.frame(all_table_melted,
                               "station_id"=NA,
                               "level_code"=NA,
                               "latitude"=NA,
                               "longitude"=NA,
                               "depth"=NA,
                               "tara_nitrate"=NA,
                               "ebi_nitrate"=NA,
                              "darwin_nitrate"=NA,
                               "nitrate"=NA,
                              "tara_temperature" = NA,
                              "ebi_temperature" =NA,
                              "temperature"=NA,
                              "tara_salinity_level"=NA,
                              "ebi_salinity_level"=NA,
                              "salinity"=NA,
                              "tara_oxygen_level"=NA,
                              "ebi_oxygen_level"=NA,
                              "darwin_oxygen_level"=NA,
                              "oxygen"=NA,
                              "darwin_alkalinity"=NA)

#Remove 0 valued abundances for memory sake
all_table_melted = all_table_melted[all_table_melted$abundance>0,]

#add environmental station information
uniq_stations=abund2_colnames

for(j in 1:length(uniq_stations)){
  tmp = station_info[station_info$station_zone==uniq_stations[j],2:21]
  tmp$station_id = as.character(tmp$station_id)
  tmp$level_code = as.character(tmp$level_code)
  if(nrow(tmp)>0)
    all_table_melted[all_table_melted$station_abundance==uniq_stations[j],13:32] = tmp
  print(j)
}
ids = which(as.numeric(all_table_melted$nitrate)<0)
all_table_melted$nitrate[ids] = all_table_melted$darwin_nitrate[ids]

##Taking out those missing some station data
#no nitrate
nonitrate_ids = which(is.na(all_table_melted$nitrate))
nooxygen_ids = which(is.na(all_table_melted$oxygen))
notemperature_ids = which(is.na(all_table_melted$temperature))
nosalinity_ids = which(is.na(all_table_melted$salinity))
ids = unique(sort(c(nonitrate_ids,nooxygen_ids,notemperature_ids,nosalinity_ids)))
all_table_melted = all_table_melted[-ids,]

#save
write.csv(all_table_melted,"alldata.csv",quote = F,row.names = F)

