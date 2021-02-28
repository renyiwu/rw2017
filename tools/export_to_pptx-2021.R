# Functions for importing (vectors) figures to Powerpoint slides.
# R Wu. Oct 2018
# revised Nov 2020

# install.packages("officer")
# install.packages("magrittr") #the package for function "%>%"
# install.packages("rvg") # for function "ph_with_vg"
# install.packages("ggplot2")
#

# Sept 2020 new dml object

create_pptx <- function(plot, path, width = 8, height = 6, left=0.05, top=0.05){
  library("officer")
  library("magrittr")
  library("rvg")
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
}
