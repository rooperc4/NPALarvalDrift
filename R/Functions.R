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



#' Function to extract the sea surface temperature at a date and location
#'
#' This function downloads daily Optimum Interpolation Sea Surface Temperatures from the NCDC website (https://www.ncdc.noaa.gov/oisst)
#' and stores them to a local destination. It will then extract the sea surface temperature for a set of dates at a set of locations.
#' THe SST is gridded on 1/4 degrees, so the output will be on the same grid size.
#' @param Longitude Longitude of interest
#' @param Latitude Latitude of interest
#' @param date  date of interest
#' @param destination_path Local path for location of AVHRR daily files
#' @param do_download TRUE (default) or FALSE should the AVHRR daily files be downloaded
#' @param Area_of_interest Either a polygon or x and y min and max defining the area of interest (defauilt is the North Pacific)
#' @keywords Larval drift, sea surface temperature, AVHRR
#' @export
#' @examples
#' temps<-AVHRR_daily_extract_function(alltraj$Lon[1:4],alltraj$Lat[1:4],alltraj$date[1:4],"D:/AVHRR Data/",do_download=FALSE)
#'
#'
#'
AVHRR_daily_extract_function<-function(Longitude,Latitude,date,destination_path,do_download=TRUE,Area_of_interest=c(-230,-120,15,60)){
  require(stringr)
  require(RCurl)
  require(ncdf4)
  require(httr)
  require(sp)

  #  file1<-nc_open("D:/AVHRR Data/oisst-avhrr-v02r01.20001115.nc")
  # Longitude<-alltraj$Lon[1:4]
  # Latitude<-alltraj$Lat[1:4]
  # date<-alltraj$date[1:4]
  # do_download<-FALSE
  # Area_of_interest<-c(-230,-120,15,60)
  # destination_path<-"D:/AVHRR Data/"

  Latitude_trim<-ceiling(abs(Latitude)*4)/4*Latitude/abs(Latitude)
  Longitude_trim<-ceiling(abs(Longitude)*4)/4*Longitude/abs(Longitude)



  if(is.vector(Area_of_interest)){
    #  Area_of_interest<-c(130,-120,15,60)
    x_coords<-c(Area_of_interest[1],Area_of_interest[1],Area_of_interest[2],Area_of_interest[2],Area_of_interest[1])
    y_coords<-c(Area_of_interest[3],Area_of_interest[4],Area_of_interest[4],Area_of_interest[3],Area_of_interest[3])
    Area_of_interest <- sp::Polygon(cbind(x_coords,y_coords))
    Area_of_interest <- sp::Polygons(list(Area_of_interest), ID = "Area")
    Area_of_interest <- sp::SpatialPolygons(list(Area_of_interest),proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))
  }

  year<-as.numeric(unique(format(date,"%Y")))
  #year<-seq(Year_min,Year_max,1)
  destpath<-destination_path
  urlpath<-"https://www.ncei.noaa.gov/data/sea-surface-temperature-optimum-interpolation/v2.1/access/avhrr/"

  month<-as.numeric(unique(format(date,"%m")))
  #month<-seq(1,12,1)
  #if(Season=="Spring"){month<-seq(3,5,1)}
  #if(Season=="Winter"){month<-c(12,1,2)}
  #if(Season=="Summer"){month<-seq(6,8,1)}
  #if(Season=="Fall"){month<-seq(9,11,1)}

  month<-formatC(month,width=2,format="d",flag="0")
  day<-as.numeric(unique(format(date,"%d")))
  #day<-seq(1,31,1)
  day<-formatC(day,width=2,format="d",flag="0")

  if(do_download==TRUE){
    ##DOWNLOAD RELEVANT MONTHS AND YEARS OF DATA
    for(i in 1:length(year)){
      for(j in 1:length(month)){
        for(k in 1:length(day)){
          urlname<-paste(urlpath,year[i],month[j],"/oisst-avhrr-v02r01.",year[i],month[j],day[k],".nc",sep="")
          if(!http_error(urlname)==TRUE){
            if(file.exists(paste(destpath,"oisst-avhrr-v02r01.",year[i],month[j],day[k],".nc",sep=""))==FALSE){
              download.file(urlname, destfile=paste(destpath,"oisst-avhrr-v02r01.",year[i],month[j],day[k],".nc",sep=""),mode="wb", quiet = FALSE)
            }}}}}}


  file1<-nc_open(paste(destpath,"oisst-avhrr-v02r01.",year[1],month[1],day[1],".nc",sep=""))
  Lat<-as.vector(ncvar_get(file1,varid="lat"))
  Lon<-as.vector(ncvar_get(file1,varid="lon"))
  Lon<-rep(Lon,720)-360
  Lat<-rep(Lat,each=1440)

  spots<-data.frame(Lon=Lon,Lat=Lat)
  meuse_sf = SpatialPoints(spots,proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))
  goodpts<-sp::over(meuse_sf,Area_of_interest)
  Lon<-Lon[is.na(goodpts)==FALSE]
  Lat<-Lat[is.na(goodpts)==FALSE]

  AVHRR1<-data.frame(lon=numeric(),lat=numeric(),value=numeric(),year=numeric(),month=numeric(),day=numeric())
  for(i in 1:length(year)){
    for(j in 1:length(month)){
      for(k in 1:length(day)){
        if(file.exists(paste(destpath,"oisst-avhrr-v02r01.",year[i],month[j],day[k],".nc",sep=""))){
          file1<-nc_open(paste(destpath,"oisst-avhrr-v02r01.",year[i],month[j],day[k],".nc",sep=""))
          AVHRR_sst<-ncvar_get(file1,varid="sst")
          AVHRR_sst<-as.vector(AVHRR_sst)
          AVHRR_sst<-AVHRR_sst[is.na(goodpts)==FALSE]
          data1<-data.frame(Longitude=ceiling(abs(Lon)*4)/4*Lon/abs(Lon),Latitude=ceiling(abs(Lat)*4)/4*Lat/abs(Lat),SST=AVHRR_sst,date=as.Date(paste(year[i],month=month[j],day=day[k],sep="-")),stringsAsFactors=FALSE)
          nc_close(file1)}
        AVHRR1<-rbind(AVHRR1,data1)}}}

  AVHRR2<-data.frame(Longitude,Latitude,Longitude_trim,Latitude_trim,date)
  AVHRR1$Longitude[AVHRR1$Longitude<=(-180)]<-AVHRR1$Longitude[AVHRR1$Longitude<=(-180)]+360

  AVHRR3<-merge(AVHRR2,AVHRR1,by.x=c("Longitude_trim","Latitude_trim","date"),by.y=c("Longitude","Latitude","date"),all.x=TRUE)
  AVHRR3<-AVHRR3[,3:5]
  return(AVHRR3)
}



#' Function to extract the Modis data at a date and location
#'
#' This function downloads MODIS satellite ocean productivity from the Oregon State University Ocean Productivity website
#' (http://sites.science.oregonstate.edu/ocean.productivity/standard.product.php) and stores them to a local destination.
#' It will then extract the chlorophyll data for a set of dates at a set of locations and merge to the original data set.
#' The satellite data is only available on 8 day intervals, so the output will contain NA where no data were available.
#' THe Ocean Color is gridded on 1/4 degrees in this output.
#' @param Longitude Longitude of interest
#' @param Latitude Latitude of interest
#' @param date  date of interest
#' @param destination_path Local path for location of MODIS daily files
#' @param do_download TRUE (default) or FALSE should the MODIS daily files be downloaded
#' @param Area_of_interest Either a polygon or x and y min and max defining the area of interest (defauilt is the North Pacific)
#' @keywords Larval drift, sea surface temperature, AVHRR
#' @export
#' @examples
#' temps<-MODIS_daily_extract_function(alltraj$Lon[1:4],alltraj$Lat[1:4],alltraj$date[1:4],"D:/MODIS data/",do_download=FALSE)
#'
#'
#'


MODIS_extract_function<-function(Longitude,Latitude,date,destination_path,do_download=TRUE,Area_of_interest=c(-230,-120,15,60)){
  require(stringr)
  require(RCurl)
  #
  #    Longitude<-alltraj$Lon[2200000:220500]
  #    Latitude<-alltraj$Lat[2000000:2000500]
  #    date<-alltraj$date[2000000:2000500]
  #    do_download<-FALSE
  #    Area_of_interest<-c(-230,-120,15,60)
  #    destination_path<-"D:/MODIS data/"
  #


  ##DOWNLOAD RELEVANT YEARS OF DATA
  datematrix<-paste(format(date,"%Y"),formatC(format(date,"%j"),width=3,format="d",flag="0"),sep="")
  datematrix<-unique(datematrix)
  year<-as.numeric(unique(format(date,"%Y")))

  urlpath<-"http://orca.science.oregonstate.edu/data/1x2/8day/vgpm.r2018.m.chl.m.sst/xyz/"
  destpath<-destination_path

  if(do_download==TRUE){
    for(i in 1:length(year)){
      urlname<-paste(urlpath,"vgpm.m.",year[i],".xyz.tar",sep="")
      if(url.exists(urlname)==TRUE){
        if(file.exists(paste(destpath,"vgpm.m.",year[i],".xyz.tar",sep=""))==FALSE){
          download.file(urlname, destfile=paste(destpath,"vgpm.m.",year[i],".xyz.tar",sep=""),mode="wb", quiet = FALSE)
          untar(paste(destpath,"vgpm.m.",year[i],".xyz.tar",sep=""),exdir=paste0(destpath,"GZ/"))}}}}
  #vgpm = net primary production (units of mg C / m**2 / day) based on the standard vgpm algorithm
  #yyyy = year
  #ddd = day of year of the start of each 8day file
  #all = all pixels output, including those with no data (nodata = -9999)
  #xyz = text file, where x = longitude, y = latitude and z = npp value
  #gz = compressed with gzip (uncompress with gunzip)


  setwd(paste0(destpath,"GZ/"))

  f1<-list.files(full.names=TRUE)
  f2<-which(as.numeric(str_sub(f1, -18, -12))%in%as.numeric(datematrix))
  f1<-f1[f2]


  if(is.vector(Area_of_interest)){
    #  Area_of_interest<-c(130,-120,15,60)
    x_coords<-c(Area_of_interest[1],Area_of_interest[1],Area_of_interest[2],Area_of_interest[2],Area_of_interest[1])
    y_coords<-c(Area_of_interest[3],Area_of_interest[4],Area_of_interest[4],Area_of_interest[3],Area_of_interest[3])
    Area_of_interest <- sp::Polygon(cbind(x_coords,y_coords))
    Area_of_interest <- sp::Polygons(list(Area_of_interest), ID = "Area")
    Area_of_interest <- sp::SpatialPolygons(list(Area_of_interest),proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))
  }

  data1<-read.table(f1[1],skip=1)
  spots<-data.frame(Lon=data1[,1],Lat=data1[,2])
  spots$Lon[spots$Lon>0]<-spots$Lon[spots$Lon>0]-360
  meuse_sf = SpatialPoints(spots,proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs"))
  goodpts<-which(sp::over(meuse_sf,Area_of_interest)>0)

  modis1<-array(dim=c(0,4))
  colnames(modis1)<-c("lon","lat","color","date")
  for(i in 1:length(f1)){

    data1<-read.table(f1[i],skip=1)
    data1<-data1[goodpts,]
    data1<-subset(data1,data1[,3]>-9999)
    if(length(data1[,1])>0){
      colnames(data1)<-c("lon","lat","color")
      data1$date<-as.Date(paste(as.numeric(str_sub(f1[i], -18, -15)),as.numeric(str_sub(f1[i], -14, -12)),sep="-"),"%Y-%j")
      modis1<-rbind(modis1,data1)}
  }


  modis1$lon[modis1$lon<0]<-modis1$lon[modis1$lon<0]+360

  Latitude_trim<-ceiling(abs(Latitude)*4)/4*Latitude/abs(Latitude)
  Longitude_trim<-ceiling(abs(Longitude)*4)/4*Longitude/abs(Longitude)

  modis1$lat<-ceiling(abs(modis1$lat)*4)/4*modis1$lat/abs(modis1$lat)
  modis1$lon<-ceiling(abs(modis1$lon)*4)/4*modis1$lon/abs(modis1$lon)
  modis1<-aggregate(color~lat+lon+date,data=modis1,FUN="mean")


  modis2<-data.frame(Longitude,Latitude,Longitude_trim,Latitude_trim,date)


  modis1<-merge(modis2,modis1,by.x=c("Longitude_trim","Latitude_trim","date"),by.y=c("lon","lat","date"),all.x=TRUE)
  modis1<-modis1[,3:5]

  return(modis1)
}
