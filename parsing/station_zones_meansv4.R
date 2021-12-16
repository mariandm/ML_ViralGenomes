#Jens table for nitrate
ns = read.csv("nitrate_levels2.csv")
station_zone = paste(ns$station_id,ns$level_code,sep="_")
uniq_station_zone = unique(sort(station_zone))

#Removing empty data from tara and ebi
ns2 = cbind(station_zone,ns)
ns2$tara_nitrate_level[ns2$tara_nitrate_level==0]=NA
ns2$ebi_nitrate_level[ns2$ebi_nitrate_level==99999.000000]=NA
ns2$ebi_nitrate_level[ns2$ebi_nitrate_level==9999.000000]=NA
ns2$tara_temperature[ns2$tara_temperature==0]=NA
ns2$ebi_temperature[ns2$ebi_temperature==99999.000000]=NA
ns2$ebi_temperature[ns2$ebi_temperature==9999.000000]=NA
ns2$tara_salinity_level[ns2$tara_salinity_level==0]=NA
ns2$ebi_salinity_level[ns2$ebi_salinity_level==99999.000000]=NA
ns2$ebi_salinity_level[ns2$ebi_salinity_level==9999.000000]=NA
ns2$tara_oxygen_level[ns2$tara_oxygen_level==0]=NA
ns2$ebi_oxygen_level[ns2$ebi_oxygen_level==99999.000000]=NA
ns2$ebi_oxygen_level[ns2$ebi_oxygen_level==9999.000000]=NA

#Use Daniel-s table to get depth
depths = read.csv("station_environmental_data2.csv")
colnames(depths) = c("station_id","level_code","depth")
station_zone = paste(depths$station_id,depths$level_code,sep="_")
depths2 = cbind(station_zone,depths)

#Map everything
ns3 = data.frame("station_zone" = uniq_station_zone,
                 "station_id" =NA,
                 "level_code"=NA,
                 "latitude" = NA,
                 "longitude" = NA,
                 "depth" = NA,
                 "tara_nitrate" = NA,
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

for(i in 1:length(uniq_station_zone)){
  ns3[i,"latitude"] = mean(ns2[ns2$station_zone==uniq_station_zone[i],"latitude"],na.rm=TRUE)
  ns3[i,"longitude"] = mean(ns2[ns2$station_zone==uniq_station_zone[i],"longitude"],na.rm=TRUE)
  ns3[i,"tara_nitrate"] = mean(ns2[ns2$station_zone==uniq_station_zone[i],"tara_nitrate_level"],na.rm=TRUE)
  ns3[i,"ebi_nitrate"] = mean(ns2[ns2$station_zone==uniq_station_zone[i],"ebi_nitrate_level"],na.rm=TRUE)
  ns3[i,"darwin_nitrate"] = mean(ns2[ns2$station_zone==uniq_station_zone[i],"darwin_nitrate_level"],na.rm=TRUE)
  ns3[i,"nitrate"] = mean(ns2[ns2$station_zone==uniq_station_zone[i],"nitrate_level"],na.rm=TRUE)
  ns3[i,"tara_temperature"] = mean(ns2[ns2$station_zone==uniq_station_zone[i],"tara_temperature"],na.rm=TRUE)
  ns3[i,"ebi_temperature"] = mean(ns2[ns2$station_zone==uniq_station_zone[i],"ebi_temperature"],na.rm=TRUE)
  ns3[i,"temperature"] = mean(ns2[ns2$station_zone==uniq_station_zone[i],"temperature"],na.rm=TRUE)
  ns3[i,"tara_salinity_level"] = mean(ns2[ns2$station_zone==uniq_station_zone[i],"tara_salinity_level"],na.rm=TRUE)
  ns3[i,"ebi_salinity_level"] = mean(ns2[ns2$station_zone==uniq_station_zone[i],"ebi_salinity_level"],na.rm=TRUE)
  ns3[i,"salinity"] = mean(ns2[ns2$station_zone==uniq_station_zone[i],"salinity_level"],na.rm=TRUE)
  ns3[i,"tara_oxygen_level"] = mean(ns2[ns2$station_zone==uniq_station_zone[i],"tara_oxygen_level"],na.rm=TRUE)
  ns3[i,"ebi_oxygen_level"] = mean(ns2[ns2$station_zone==uniq_station_zone[i],"ebi_oxygen_level"],na.rm=TRUE)
  ns3[i,"darwin_oxygen_level"] = mean(ns2[ns2$station_zone==uniq_station_zone[i],"darwin_oxygen_level"],na.rm=TRUE)
  ns3[i,"oxygen"] = mean(ns2[ns2$station_zone==uniq_station_zone[i],"oxygen_level"],na.rm=TRUE)
  ns3[i,"darwin_alkalinity"] = mean(ns2[ns2$station_zone==uniq_station_zone[i],"darwin_alkalinity"],na.rm=TRUE)
  
  ns3[i,"depth"] = mean(depths2[depths2$station_zone==uniq_station_zone[i],"depth"],na.rm =TRUE)
  
  ns3[i,"station_id"] = strsplit(uniq_station_zone[i],"_")[[1]][1]
  ns3[i,"level_code"] = strsplit(uniq_station_zone[i],"_")[[1]][2]
}

write.csv(ns3,"station_zone_means.csv",quote = F,row.names = F)

