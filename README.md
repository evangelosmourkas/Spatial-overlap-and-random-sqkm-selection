# Spatial Overlap and Random sqkm selection in a geographical space
1. Randomly select squares in the same geographical space (e.g. country)
2. Estimating spatial overlap between data points from two different populations occupying the same geographical space.

## Preprint
The following code has been developed in collaboration with **Jose Valdebenito** and functions of code used in by Marcel Cardillo & Dan L. Warren in _"Analysing patterns of spatial and niche overlap among species at multiple resolutions"_ paper published in 2016. The methodology underlying the use of the scripts has been detailed in the preprint _"Urbanization spreads antimicrobial resistant enteric pathogens in wild bird microbiomes"_. Please refer to citation information at the bottom of this document.

## Files Summary
The files include data for two populations, each strain representing a different condition. Each directory has the following:
* **_Population_data_1_**: data for population 1 (latitude, longitude, population density score)
* **_Population_data_2_**: data for population 2 (latitude, longitude, population density score)

# 1. Randomly select 100 squares (1x1 degree of latitude and longitude) within the territory of a country
## R Packages for Selecting Random Squares in the same geographical space (e.g. country)
* ggplot2
* sf
* sp
* terra
* spatialEco
* rnaturalearth
* dplyr

## Generate 100 random points within the territory of a country (e.g. US)
```
set.seed(2022)
shp.us <- ne_countries(country = "United States of America",scale = "large", returnclass = "sp") #it has to be class sp
plot(shp.us) # country silhouette
```
![image](https://github.com/evangelosmourkas/Spatial-overlap-in-the-same-geographical-area/assets/73548463/45f36c48-2be1-47f7-8c54-1c70b088dbfa)

## Select 100 points within the country
```
random_poly_us <- spsample(shp.us, n=100, "random") 
plot(random_poly_us)
```
![image](https://github.com/evangelosmourkas/Spatial-overlap-in-the-same-geographical-area/assets/73548463/22abc355-1180-4d4a-9eaa-44c4fdf28f68)

## Saving the coordinates of the 100 points and adding +/- 0.5 degrees to points to create squaremeters of 110square km
```
df.us<-as.data.frame(random_poly_us@coords)
df.us$x.min <- df.us$x-0.5
df.us$x.max <- df.us$x+0.5
df.us$y.min <- df.us$y-0.5
df.us$y.max <- df.us$y+0.5
```

## Loading population 1 file. In this case this is human density data in the US. This file is too big so we split it into two
```
file1 = read.csv("Population_1_1.csv")
file2 = read.csv("Population_2_1.csv")
us.pop.r <- rbind(file1, file2)
us.pop.r <- us.pop.r[us.pop.r$Population>0,]
```
## Filtering loop for population 1 (Human density data in the US)
```
us.pop.f <- data.frame()
for(i in 1:100){
  us1.1 <- us.pop.r %>% filter(Lat > df.us[i,5], Lat < df.us[i,6])
  us1.2 <- us1.1 %>% filter(Lon > df.us[i,3], Lon < df.us[i,4])
  us.pop.f <- rbind(us.pop.f, us1.2)
}
```

## Plotting map, squares and human population
```
sf.us <- ne_countries(country = "United States of America",scale = "large", returnclass = "sf") #to plot in ggplot it needs to be class sf
ggplot() +
  geom_sf(data= sf.us, color = "dark green", fill = "white") + theme_classic()+
  geom_point(data=us.pop.r, aes(y=Lat, x=Lon), size=0.01, color="black", alpha = 0.01) + #adding population raw
  #geom_point(data=us.pop.f, aes(y=Lat, x=Lon), size=0.01, color="black", alpha = 0.01) + #adding population filtered
  #geom_point(data=swe.ancre.f, aes(y=Lat, x=Lon), size=0.1, color="yellow",  alpha = 0.4)+ #adding birds
  #geom_rect(data= df.us, aes(xmin = x.min, xmax = x.max, ymin = y.min, ymax = y.max), color = "red", fill = NA, size=0.3)+#adding squares
  coord_sf(crs = "+proj=lonlat +lat_0=45 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs")+
  coord_sf(xlim = c(-167, -66.98292), ylim = c(19.03319, 69.64097), expand = TRUE)
```
![image](https://github.com/evangelosmourkas/Spatial-overlap-in-the-same-geographical-area/assets/73548463/03a922f0-3c63-406a-bdfc-8733968760f2)

## Adding new coordinate system, which will replace the existing one
```
ggplot() +
  geom_sf(data= sf.us, color = "dark green", fill = "white") + theme_classic()+
  #geom_point(data=us.pop.r, aes(y=Lat, x=Lon), size=0.01, color="black", alpha = 0.01) + #adding population raw
  geom_point(data=us.pop.f, aes(y=Lat, x=Lon), size=0.01, color="black", alpha = 0.01) + #adding population filtered
  #geom_point(data=swe.ancre.f, aes(y=Lat, x=Lon), size=0.1, color="yellow",  alpha = 0.4)+ #adding birds
  geom_rect(data= df.us, aes(xmin = x.min, xmax = x.max, ymin = y.min, ymax = y.max), color = "red", fill = NA, size=0.3)+#adding squares
  coord_sf(crs = "+proj=lonlat +lat_0=45 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs")+
  coord_sf(xlim = c(-167, -66.98292), ylim = c(19.03319, 69.64097), expand = TRUE)
```
![image](https://github.com/evangelosmourkas/Spatial-overlap-in-the-same-geographical-area/assets/73548463/be418af5-0a61-44cd-b7bc-776077fa7724)

## Adding population 2 values. In our case this will be wild bird observational data for _Anas acuta_
```
us.anasacu.r <- read.csv("/Volumes/EXT\ 1T\ JOSE/AMR/sp\ by\ countries/US-AnaAcu\ Galaxy46-Select_random_lines_on_data_26.csv")
head(us.anasacu.r)
```
## Filtering loop for population 2 (wild bird)
```
us.anasacu.f <- data.frame() #filterred population
for(i in 1:100){
  us1.1 <- us.anasacu.r %>% filter(Lat > df.us[i,5], Lat < df.us[i,6])
  us1.2 <- us1.1 %>% filter(Lon > df.us[i,3], Lon < df.us[i,4])
  us.anasacu.f <- rbind(us.anasacu.f, us1.2)
}
```

## Zoom in figure to show details
```
ggplot() +
  geom_sf(data= sf.us, color = "dark green", fill = "white") + theme_classic()+
  geom_point(data=us.pop.f, aes(y=Lat, x=Lon), size=0.01, color="black", alpha = 0.01) + #adding population filtered
  geom_point(data=us.anasacu.f, aes(y=Lat, x=Lon), size=0.1, color="yellow",  alpha = 0.4)+ #adding birds
  geom_rect(data= df.us, aes(xmin = x.min, xmax = x.max, ymin = y.min, ymax = y.max), color = "red", fill = NA, size=0.3)+#adding squares
  coord_sf(crs = "+proj=lonlat +lat_0=45 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs")+
  coord_sf(xlim = c(-100, -75), ylim = c(22, 40), expand = TRUE)
```
![image](https://github.com/evangelosmourkas/Spatial-overlap-in-the-same-geographical-area/assets/73548463/a23249b7-aefd-4698-83bb-d0228fd7d28c)


# 2. Calculating Spatial Overlap
## R Packages for calculating overlap
* fields
Note: for function rdist()
* ENMTools
Note: for function nncluster().
If ENMTools cannot be installed, paste nncluster() directly from source (https://github.com/cran/nnclust/blob/master/R/nncluster.R)

```
nncluster <- function(xy, species){
    species.names <- unique(species)
  df <- as.data.frame(matrix (nrow=length(species.names), ncol=length(species.names)))
  rownames(df) <- species.names
  colnames(df) <- species.names
  
  for(i in 1:length(species.names)){
    for(j in i:length(species.names)){
      if(i != j){
        
        if(nrow(xy[species == species.names[i],]) > 1 && nrow(xy[species == species.names[j],]) > 1){
          
          print(paste("Calculating", species.names[i], "vs", species.names[j]))
          within1 <- rdist(xy[species == species.names[i],],xy[species == species.names[i],]) # eucl dist within sp1 datapoints
          within2 <- rdist(xy[species == species.names[j],],xy[species == species.names[j],]) # eucl dist within sp2 datapoints
          between <- rdist(xy[species == species.names[i],],xy[species == species.names[j],]) # eucl dist betwen sp1 and sp2 datapoints
          score1 <- rep(NA,length(within1[1,])) #the length equals the number of datapoints initially included
          score2 <- rep(NA,length(within1[2,])) #the length equals the number of datapoints initially included
          
          for(k in 1:length(within1[1,])){
            thisscore <- min(within1[k,-(k)])/min(between[k,]) # O(xi) = w(xi) / b(xi)
            score1[k] <- thisscore
          }
          
          for(k in 1:length(within2[1,])){
            thisscore <- min(within2[k,-(k)])/min(between[,k]) # O(xi) = w(xi) / b(xi)
            score2[k] <- thisscore
          }
          
          score1 <- length(which(score1>1))/length(score1) # p = n(O > 1) / n(O)
          score2 <- length(which(score2>1))/length(score2) # p = n(O > 1) / n(O)
          w <- mean(c(score1, score2)) # O = p(x) + p(y) / 2
          
        }
        else{w <- NA}
        df[i,j] <- w
      }
    }
    
  }
  return(df)
}
```
Note: nncluster() function works best if provided equal data points to determine the overlap from. Thus we randomly select 5000 data points from the human and bird population data. Remember that our currently our human and bird data corresponds to coordinates within the 100 squares we selected earlier. 5000 point coordinates should be an even representation of both datasets.

## Calculating spatial overlap of population 1 (human) and population 2 (wild bird) distributions. The calculated value in the end is referred to the paper as the _proximity to urbanisation_ score
```
set.seed(2022)
spess <- c(NA)
prox <- c(NA)
country <- c(NA)
prox.result <- as.data.frame(cbind(spess,prox,country))
asd <- matrix(1:130, nrow=130,ncol = 1)
for(i in 1:10){
us <- rbind(us.pop.f[sample(nrow(us.pop.f), 5000, replace = T),1:2],
            us.anasacu.f[sample(nrow(us.anasacu.f), 5000, replace = T),3:4])
us.sp <- c(rep("pop", 5000), rep("anacuta", 5000))
us <- cbind(us,us.sp)
asd[i] <- nncluster(xy=us[us$us.sp==c("anacuta","pop"),1:2],species=us[us$us.sp==c("anacuta","pop"),3])[1,2]
}
spess <- c(rep("anacuta", 10))
asd <- as.data.frame(cbind(asd,spess))
asd$V1 <- as.numeric(asd$V1)
asd$spess <- as.factor(asd$spess)
over <- aggregate(V1~spess,data=asd,mean)
over

colnames(over)[2] <- "prox"
over$country <- rep("USA", 13)
prox.result <- rbind(prox.result, over)
```
# How to cite
Mourkas E, Valdebenito JO, Marsh H, Hitchings MD, Cooper KK, Parker CT, Székely T, Johansson H, Ellström P, Pascoe B, Waldenström J, Sheppard SK (2023) **Urbanisation spreads antimicrobial resistant enteric pathogens in wild bird microbiomes**.
bioRxiv doi: 10.1101/2023.07.11.548564
