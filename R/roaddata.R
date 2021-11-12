#' Evaluates scenarios in an impact network analysis (INA)
#'
#' This function implements and summarizes multiple simulations across a designated range of parameter values
#'
#' Updated 2021-10-22

#' @param geocoords5 matrix of x,y coordinates of nodes
#' @param roaddistfilepath5n filepath to a R Data (.RData extension) file containing an adjacency matrix with the distances via the road network between nodes in geocoords - mutually exclusive with roadtimefilepath5n
#' @param roadtimefilepath5n filepath to a R Data (.RData extension) file containing an adjacency matrix with the travel times via the road network between nodes in geocoords - mutually exclusive with roaddistfilepath5n
#' @keywords road
#' @export 
#' @import osrm

#' @examples

#NEED TO WRITE EXAMPLES



roaddata = function(geocoords5, roaddistfilepath5=NA, roadtimefilepath5=NA) {
	filepath = "";
	isDist = F;
    if(!is.na(roaddistfilepath5)) {
	filepath = roaddistfilepath5;
	isDist = T;
    } else if(!is.na(roadtimefilepath5))  {
	filepath = roadtimefilepath5;
	isDist = F;
    } else {
	stop("Please specify either roaddistfilepath or roadtimefilepath");
    }

	if(file.exists(filepath)) {
		load(filepath);

		if(isDist) {
			return(distMat); # Loaded from file
		} else {
			return(timeMat); # Loaded from file
		}

	} else {

		message("The data file containing the adjacency matrix with distances between locations has not been found.");
		message("It will be downloaded and generated using the OSRM API.");
		message("If you have this file, ensure that the correct filepath has been entered.");
		message("This download may take a while.");

    continue = readline(prompt="Are you sure you would like to continue (Y/n): ");

    if(!(continue == "Y" || continue == "y")) {
      stop("Either do not specify roaddatafilepath, specify a valid roaddatafilepath, or allow download to occur for simulation to proceed");
    }

		# Load osrm library if not already loaded
		# Fail if osrm is not installed
		allLibs = library()$results[,1];

		if(!is.element("osrm", allLibs)) {
			message("osrm library not installed, exiting");
      message("Please install osrm package to use road distance functionality");
			return;
		}

		loadedLibs = (.packages());

		if(!is.element("osrm", loadedLibs)) {
			library(osrm);
		}

		# With osrm loaded, proceed to get data


    #BIG ASSUMPTION: geocoords (x,y) data coorresponds to coordinate (lon, lat) values
    #Convert coordinates into list, needed for function
    #Also need to add id for OSRM package
    lon = geocoords5[,1]
    lat = geocoords5[,2]
    id = 1:length(lon)
    #List varible order is important, OSRM breaks if out of order
    locations = data.frame(id=id, lon=lon, lat=lat)

	count = length(locations$id);
	distMat = matrix(nrow=count, ncol=count, dimnames=list(locations$id,locations$id));
	timeMat = matrix(nrow=count, ncol=count, dimnames=list(locations$id,locations$id));

	#The OSRM API's default server does not allow the generation of a matrix with greater than or equal to
	#10000 entries (a 100x100 square matrix, 50x200 rectangular matrix, etc.)
	#Thus, for lists with 100 points or more, the matrix needs to be loaded in 99x99 submatrices, and stitched together
	#into the final matrix

	repCount = ceiling(count / 99);

	loadedCount = 0;
	for(y in 1:repCount) {
		for(x in 1:repCount) {
			#Declare lower and upper bounds for locations to be downloaded
			xl = 1 + 99 * (x - 1);
			yl = 1 + 99 * (y - 1);
			xu = min(99 * x, count); #Don't exceed maximum matrix dimension
			yu = min(99 * y, count); #Don't exceed maximum matrix dimension

			#Get result from OSRM
			result = osrmTable(
				src=locations[yl:yu, c("id", "lon", "lat")],
				dst=locations[xl:xu, c("id", "lon", "lat")],
				measure=c("distance", "duration"), osrm.profile="car"
			);

			distMat[yl:yu,xl:xu] = result$distances #Copy from result into large matrix
			timeMat[yl:yu,xl:xu] = result$durations #Copy from result into large matrix
			loadedCount = loadedCount + (xu-xl+1) * (yu-yl+1); #Update progress counter

			#Output progress
			message(paste("Downloaded ", loadedCount, " of ", count^2, " distances, ", trunc(100 * loadedCount / count^2, digits=5), "% done", sep=""));
		}
	}

    #Check for NA values, should only need to do for either distMat or timeMat, but not both (would be redundant)
    temp = c(distMat);
    #idx = 1:length(temp)
    if(length(temp[is.na(temp)]) > 0) {
      message("There are some NA entries in the adjacency matrix.");
      message("This is because a road network route could not be found between at least one pair of locations.");
      message("This may occur if, for instance, one location is on an island.")
      message("These can either be replaced with Inf (infinity), these locations can be removed from the dataset.");
      continue = readline(prompt="Substitute for infinity, or terminate and manually filter dataset (I/r): ");

      if(!(continue == "I" || continue == "i")) {
        stop("Please rerun with an approppriate dataset");
      }

      #Replace all NA with Inf
      for(i in 1:count) {
        for(j in 1:count) {
          if(is.na(distMat[i,j])) {
            distMat[i,j] = Inf;
          }
        }
      }


      for(i in 1:count) {
        for(j in 1:count) {
          if(is.na(distMat[i,j])) {
            timeMat[i,j] = Inf;
          }
        }
      }

      #Probably more efficient way
      #temp[is.na(temp)] = Inf;
      #asdf=matrix(nrow=nrow(distMat), temp);

    }

		save( distMat, timeMat, file=filepath);

		if(isDist) {
			return(distMat); # Loaded from file
		} else {
			return(timeMat); # Loaded from file
		}
	}
}

validate = function(locations, distMat) {
	x = round(runif(1, 1, 228));
	y = round(runif(1, 1, 228));
	#x=228;
	#y=228;
	width = round(runif(1, 1, 98));
	height = round(runif(1, 1, 98));

	print(x)
	print(y)
	print(width)
	print(height)


	#waste = file.remove("roaddist.RData");

	#distMat = genDist(locations, "roaddist.RData");

	download = osrmTable(
			src=locations[y:(y+98), c("id", "lon", "lat")],
			dst=locations[x:(x+98), c("id", "lon", "lat")],
			measure="distance", osrm.profile="car"
		);

	download2 = osrmTable(
			src=locations[y:(y+height - 1), c("id", "lon", "lat")],
			dst=locations[x:(x+width - 1), c("id", "lon", "lat")],
			measure="distance", osrm.profile="car"
		);

	#str(download$distances);
	#str(distMat[y:(y+98),x:(x+98)]);

	#print(download$distances - distMat[y:(y+98),x:(x+98)]);

	#print(sum(abs(distMat[y:(y+98),x:(x+98)] - download$distances)));
	#print(sum((distMat[y:(y+98),x:(x+98)] - download$distances)));

	print(sum(abs(download$distances[1:(height),1:(width)] - download2$distances)));
	print(sum((download$distances[1:(height),1:(width)] - download2$distances)));

	#print(download$distances);
	#print(download2$distances);

	write.table(download$distances, file="download.txt");
	write.table(download2$distances, file="download2.txt");
}

#test = function() {
	#cities = read.csv("cities.csv");
	#filepath = "roaddist.RData";
	#waste = genDist(cities,filepath);
#}
