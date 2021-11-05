src = function() {
  source("INAscene.R");
  source("estab.R");
  source("estinfo.R");
  source("genlocs.R");
  source("genmovnet.R");
  source("initvals.R");
  source("layerplot.R");
  source("load.R");
  source("makedec.R");
  source("multistart.R");
  source("multsame2.R");
  source("ntsteps2.R");
  source("onestart.R");
  source("reduce_nodes_to_focus.R");
  source("roaddata.R");
  source("setup2.R");
  source("smartsurv.R");
  source("smartsurv.weight.R");
  source("spreadstep.R");
  source("squish.R");

  test();
}

test = function() {
  #cities = read.csv("cities.csv");
  #count = 326;
  #coords = t(matrix(nrow=2,c(cities$lon[1:count], cities$lat[1:count]), byrow=T));

  #Use Ethiopia data
  cities = read.csv("farms.csv");
  count = 1259;
  coords = t(matrix(nrow=2,c(cities$lon[1:count], cities$lat[1:count]), byrow=T));

  initinfo = round(runif(count, 0, 1));
  initbio = round(runif(count, 0, 1));

j25.readgeocoords <- INAscene(nreals=3,
  ntimesteps=3,
  doplot=F,
  readgeocoords=T,
  #geocoords=matrix(c(1, 1, 1, 2, 1, 3, 2, 1, 2, 2, 2, 3), byrow=T, ncol=2),
  geocoords=coords,
  roaddatafilepath="roaddata.RData",
  numnodes=NA,
  xrange=NA,
  yrange=NA,
  randgeo=F,
  readinitinfo=T,
  initinfo=initinfo,
  initinfo.norp=NA,
  initinfo.n=NA,
  initinfo.p=NA,
  initinfo.dist=NA,
  readinitbio=T,
  initbio=initbio,
  initbio.norp=NA,
  initbio.n=NA,
  initbio.p=NA,
  initbio.dist=NA,
  readseam=F,
  seam=NA,
  seamdist='powerlaw',
  seamrandp=c(0.01, 0.05, 0.1, 0.5),
  seampla=1,
  seamplb=1,
  readbpam=F,
  bpam=NA,
  bpamdist='random',
  bpamrandp=0.1,
  bpampla=NA,
  bpamplb=NA,
  readprobadoptvec=F,
  probadoptvec=NA,
  probadoptmean=0.1,
  probadoptsd=0.1,
  readprobestabvec=F,
  probestabvec=NA,
  probestabmean=0.1,
  probestabsd=0.1,
  maneffdir='decrease_estab',
  maneffmean=0.5,
  maneffsd=0.1,
  usethreshman=F,
  maneffthresh=NA,
  sampeffort=NA
) 
}
