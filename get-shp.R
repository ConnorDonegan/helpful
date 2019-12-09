# provide url to download a shapefile and unzip it into the working directory
get.shp <- function(url, folder = "shape") {
	tmp.dir <- tempfile(fileext=".zip")
	download.file(url, destfile = tmp.dir)
	unzip(tmp.dir, exdir = folder)
	list.files(folder, full.names = TRUE)
	}

