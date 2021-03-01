# Functions for importing (vectors) figures to Powerpoint slides.
# R Wu. Oct 2018
# revised Feb 2021
# add support for png and tiff


# install.packages("officer")
# install.packages("magrittr") #the package for function "%>%"
# install.packages("rvg") # for function "ph_with_vg"
# install.packages("ggplot2")
#

# Sept 2020 new dml object

create_pptx <- function(plot, path, width = 8, height = 6, left=0.05, top=0.05, png=TRUE, tiff=FALSE, res=100){
  library("officer")
  library("magrittr")
  library("rvg")
  library("stringr")
  if(!file.exists(path)) {
    out <- read_pptx()
  } else {
    out <- read_pptx(path)
  }

  out %>%
    add_slide(layout = "Title and Content", master = "Office Theme") %>%
    ph_with(dml(code = print(plot)), # rvg:::ph_with.dml(dml(code = print(plot)),
                      # ph_location_fullsize()
            ph_location(lef=left, top=top, width=width, height = height)
                      ) %>%
    print(target = path)
  
  ### create png ### 
  # check series numbers.
  if (png){
  n <- 1
  while (n){
    file <- paste0(substr(path,1, nchar(path)-4), str_pad(n, 3, "left", "0"), ".png")
    if (! file.exists(file))
      break
    n = n + 1
  }
  png(file, width=width, height = height, units = "in", res = res)
  print(plot)
  dev.off()
  }
  
  ### create tiff
  if (tiff){
    n <- 1
    while (n){
      file <- paste0(substr(path,1, nchar(path)-4), str_pad(n, 3, "left", "0"), ".tif")
      if (! file.exists(file))
        break
      n = n + 1
    }
    tiff(file, width=width, height = height, units = "in", res = res)
    print(plot)
    dev.off()
  }
}
