# Functions for importing (vectors) figures to Powerpoint slides.
# R Wu. Oct 2018
# revised Nov 2020

# install.packages("officer")
# install.packages("magrittr") #the package for function "%>%"
# install.packages("rvg") # for function "ph_with_vg"
# install.packages("ggplot2")
#

# Sept 2020 new dml object

create_pptx <- function(plot, path, w=10, h = 7.5, l = 0.05, t = 0.05 ){
  library("officer")
  library("magrittr")
  library("rvg")
  
 if(file.exists(path)) {
   if (!file.exists(paste0("~$", path))){
    out <- read_pptx(path)
   }
  } 
  else {
    path <- paste0("temp_", path)
     out <- read_pptx(path)
   }

   out %>%
     add_slide(layout = "Title and Content", master = "Office Theme") %>%
     ph_with(dml(code = print(plot)), 
             # rvg:::ph_with.dml(dml(code = print(plot)),
             # ph_location_fullsize()
             ph_location(left = l, top = t, width = w,  height = h, newlabel = "", bg = NULL, rotation = NULL)
             ) %>%
     print(target = path)

}


## to use in a loop
#1. first create an object outside the loop
# pptx_out <- read_pptx()

#2. then put this in the loop
pptx_out %>%
add_slide(layout = "Title and Content", master = "Office Theme") %>%
ph_with(dml(code = print(plot)), # rvg:::ph_with.dml(dml(code = print(plot)),
        ph_location_fullsize()
)

#3. finally write to file
# print(pptx_out, path)