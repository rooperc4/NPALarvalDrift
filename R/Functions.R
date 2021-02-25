#' Function to extract start-end positions and distance
#'
#' This function reads a net cdf file output from Ocean Parcels and compiles the final position of the particle, as well as the start position,
#' beginning and ending dates, the cumulative distance travelled and the straight-line distance travelled. It also attaches the year of release
#' and the seamount of release to the output.
#' @param nc_file Output ncdf file from Ocean Parcels
#' @param drift_days number of days for which to calculate the drift
#' @keywords globcurrent larval drift duration
#' @export
#' @examples
#' Particle_drift(kkjkh, 120)


Particle_drift<-function(nc_file,drift_days){
  #nc_file<-"C:/Users/rooperc/Desktop/OceanParcels/particle/y00_d0.nc"
  #drift_days<-120
  require(ncdf4)
  nc_file<-nc_open(nc_file)
  Lon<-c(ncvar_get(nc_file,varid="lon"))
  Lat<-c(ncvar_get(nc_file,varid="lat"))
  Time<-c(ncvar_get(nc_file,varid="time"))
  Time<-as.Date(as.POSIXct(Time, origin = "1995-01-01 12:00:00"))
  Depth<-c(ncvar_get(nc_file,varid="z"))
  NPA<-c(ncvar_get(nc_file,varid="trajectory"))
  nc_close(nc_file)

  data1<-data.frame(Lon=Lon,Lat=Lat,Time=Time,Depth=Depth,NPA=NPA)
  data1<-data1[is.na(data1$Time)==FALSE,]
  yc<-paste(unique(format(Time[1],"%Y"))[1],unique(format(Time,"%Y"))[2],sep="-")

  seamounts<-data.frame(Seamount=c("Suiko","Showa","Youmei","Nintoku","Jingu","Ojin","Koko","Kinmei","Yuryaku","Kammu","Colahan","C-H","NW Hancock","SE Hancock"),
                        Longitude=c(170.3, 170.4, 170.4, 170.6, 171.2, 170.5, 171.6, 171.5, 172.3, 173, 176, 177.6, 178.7, 179.1),
                        Latitude=c(44.6, 43, 42.3, 41.1, 38.8, 38, 35.3, 33.7, 32.7, 32.2, 31.3, 30.4, 30.3, 29.8),stringsAsFactors = FALSE)


   outdata<-data.frame(particle_id=numeric(),seamount=character(),year=character(),start_date=character(),end_date=character(),start_lon=numeric(),start_lat=numeric(),end_lon=numeric(),end_lat=numeric(),distance=numeric(),cumulative_distance=numeric())
   particles<-unique(data1$NPA)

  for(i in 1:length(particles)){
    temp<-subset(data1,data1$NPA==particles[i])
    p1<-particles[i]
    s1<-seamounts$Seamount[which(round(seamounts$Latitude,2)==round(temp$Lat[1],2)&round(seamounts$Longitude,2)==round(temp$Lon[1],2))]
    t3<-yc
    t1<-temp$Time[which.min(temp$Time)]
    t2<-temp$Time[which.min(temp$Time)+drift_days]
    l1<-temp$Lon[which.min(temp$Time)]
    l2<-temp$Lat[which.min(temp$Time)]
    l3<-temp$Lon[which.min(temp$Time)+drift_days]
    l4<-temp$Lat[which.min(temp$Time)+drift_days]
    d1<-dist_xy(l1,l2,l3,l4,"km")
    d2<-sum(dist_xy(temp$Lon[1:(drift_days-1)],temp$Lat[1:(drift_days-1)],temp$Lon[2:drift_days],temp$Lat[2:drift_days],"km"))

    temp1<-data.frame(particle_id=p1,seamount=s1,year=t3,start_date=t1,end_date=t2,start_lon=l1,start_lat=l2,end_lon=l3,end_lat=l4,distance=d1,cumulative_distance=d2)
    outdata<-rbind(outdata,temp1)
  }

  return(outdata)

}


#' Function to extract start-end positions and distance
#'
#' This function reads a net cdf file output from Ocean Parcels and outputs the trajectory of the particals on an interval in days. The function
#' also attaches the year of release and the seamount of release to the output.
#' @param nc_file Output ncdf file from Ocean Parcels
#' @param drift_days number of days for which to calculate the drift
#' @param interval interval in number of days over which to output the larval position (e.g. daily (the default), weekly = 7)
#' @keywords globcurrent larval drift duration
#' @export
#' @examples
#' Particle_drift(kkjkh, 120, 7)


Particle_trajectory<-function(nc_file,drift_days,interval=1){
  #nc_file<-"C:/Users/rooperc/Desktop/OceanParcels/particle/y00_d0.nc"
  #drift_days<-120
  require(ncdf4)

   int1<-seq(1,drift_days,interval)
    if(int1[length(int1)]!=drift_days){int1<-c(int1,drift_days)}

  nc_file<-nc_open(nc_file)
  Lon<-c(ncvar_get(nc_file,varid="lon"))
  Lat<-c(ncvar_get(nc_file,varid="lat"))
  date<-c(ncvar_get(nc_file,varid="time"))
  date<-as.Date(as.POSIXct(date, origin = "1995-01-01 12:00:00"))
  Depth<-c(ncvar_get(nc_file,varid="z"))
  particle_id<-c(ncvar_get(nc_file,varid="trajectory"))
  nc_close(nc_file)

  data1<-data.frame(Lon=Lon,Lat=Lat,date=date,Depth=Depth,particle_id=particle_id)
  data1<-data1[is.na(data1$date)==FALSE,]
  yc<-paste(unique(format(date[1],"%Y"))[1],unique(format(date,"%Y"))[2],sep="-")

  seamounts<-data.frame(Seamount=c("Suiko","Showa","Youmei","Nintoku","Jingu","Ojin","Koko","Kinmei","Yuryaku","Kammu","Colahan","C-H","NW Hancock","SE Hancock"),
                        Longitude=c(170.3, 170.4, 170.4, 170.6, 171.2, 170.5, 171.6, 171.5, 172.3, 173, 176, 177.6, 178.7, 179.1),
                        Latitude=c(44.6, 43, 42.3, 41.1, 38.8, 38, 35.3, 33.7, 32.7, 32.2, 31.3, 30.4, 30.3, 29.8),stringsAsFactors = FALSE)


  outdata<-data.frame(particle_id=numeric(),seamount=character(),year=character(),date=character(),lon=numeric(),lat=numeric())
  particles<-unique(data1$particle_id)

  for(i in 1:length(particles)){
    temp<-subset(data1,data1$particle_id==particles[i])
    temp<-temp[int1,]
    temp$seamount<-seamounts$Seamount[which(round(seamounts$Latitude,2)==round(temp$Lat[1],2)&round(seamounts$Longitude,2)==round(temp$Lon[1],2))]
    temp$year<-yc

   outdata<-rbind(outdata,temp)
  }

  return(outdata)

}



#' Function to calculate distance
#'
#' This function calculates the distance between two points on a round earth.
#' @param Longitude_start Start longitude
#' @param Latitude_start Start latitude
#' @param Longitude_end  End longitude
#' @param Latitude_end End latitude
#' @param unit units to use for calculation (either m or km)
#' @keywords distance calculation
#' @export
#' @examples
#' dist_xy(-120,58,-121,59,"km")



dist_xy<-function(Longitude_start,Latitude_start,Longitude_end,Latitude_end,unit){
  #lat1<-56
  #lat2<-57
  #long1<-(-174)
  #long2<-(-171)
  #unit="m"


  r<-ifelse(unit == "m", 6371200,6371.2)

  la1 = Latitude_start * pi / 180
  la2 = Latitude_end *pi / 180
  lo1 = Longitude_start * pi / 180
  lo2 = Longitude_end * pi/ 180


  rads <- acos((sin(la1) * sin(la2)) + (cos(la1) * cos(la2) * cos(lo2 - lo1)))
  return(r * rads)
}
