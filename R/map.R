library(rgl)
library(raster)
library(viridis)
library(tiff)
library(rworldmap)

rotateccw = function(x) t(apply(t(x),2,rev));
rotatecw = function(x) t(apply(x,2,rev));

rendermap = function(datatype="", infected) {
	filename = "img/ethiopia";
	if(datatype != "") {
		filename = paste(filename, "_", datatype, sep="")

	}

	filename = paste(filename, ".tif", sep="")

	readimg = function(raster, filename) {
		col_rgb = readTIFF(filename);

		col = matrix(nrow=raster@nrows,ncol=raster@ncols);

		for (x in 1:raster@nrows) {
			for (y in 1:raster@ncols) {
				col[x,y] = rgb(col_rgb[x,y,1], col_rgb[x,y,2], col_rgb[x,y,3]);
			}
		}

		return(rotatecw(col));
	}

	africa_map_altitude <- raster("img/ethiopia.tif")
	ethiopia<-africa_map_altitude/max(getValues(africa_map_altitude),na.rm=TRUE) #rescale between 0 and 1 (otherwise scale is in minutes)

	scale = 5;

	heights = rotateccw(scale * matrix(ethiopia@data@values,nrow=ethiopia@nrows));

	col = readimg(ethiopia, filename);


	rgl.clear();

	str(heights);

	surface3d(1:nrow(heights),1:ncol(heights),heights,color=col)

	data("countriesLow") #load map border data
	countriesLow = subset(countriesLow,LAT>0 & LAT<20 & LON>30 & LON<50) #Filter to only local countries

	for(j in 1:length(countriesLow)) {
		points = countriesLow@polygons[[j]]@Polygons[[1]]@coords;

		x = ((points[,1] - 30) / 20 * 240);
		y = ((1 - (points[,2] - 0) / 20) * 240);

		filter = (x > 0 & x < 240 & y > 0 & y < 240);
		x = x[filter];
		y = y[filter];

		z = c();
		for(i in 1:length(x)) {
			z = c(z, 1.1 * heights[round(x[i]), 240 - round(y[i])]);
		}

		lines3d(x,240 - y,z,col="black",lwd=4);

	}

	rgl.viewpoint(phi=5, zoom=.75);
	return(heights);
}

points = function(infected, heights) {
	visited = rep(FALSE, length(infected[[1]]));

	for(i in 1:length(infected)) {
		for(j in 1:length(infected[[i]])){
			count = 0;
			if(infected[[i]][[j]] & !visited[[j]]) {
				x = ((cities$lon[[j]] - 30) / 20 * 240);
				y = ((1 - (cities$lat[[j]]- 0) / 20) * 240);

				z = 1.1 * heights[round(x), 240 - round(y)];

				spheres3d(x,240 - y,z,col="red",radius=1);
				count = count + 1;

				visited[[j]] = TRUE;
			}
		}



		filename=paste("/mnt/ram/frame",i,".png",sep="");
		rgl.snapshot(filename=filename,fmt="png");
	}
}

INAPlot = function() {
	cities = read.csv("farms.csv");
	count = 1259;
	coords = t(matrix(nrow=2,c(cities$lon[1:count], cities$lat[1:count]), byrow=T));

	initinfo = round(runif(count, 0, 1));
	initbio = round(runif(count, 0, 1));

	#res = src();
	#rgl.open();
	rgl.clear();

	heights=rendermap();
	points(infected=res$multdetails[[1]]$multout[[1]]$estabvecL, heights=heights);

}

ll = function() {
	source("map.R");
	INAPlot();
}
