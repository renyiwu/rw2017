# Funtions for importing (vectors) figures to Powerpoint slides.
# R Wu. Oct 2018
# revised May 2019

# install.packages("officer")
# install.packages("magrittr") #the package for function "%>%"
# install.packages("rvg") # for function "ph_with_vg"
# install.packages("ggplot2")
#

# function 1, for plots saved in variables. eg by pheatmap,
## deprecated as of Sept 2020
create_pptx <- function(plot, path, width = 6, height = 6, pointsize = 12){
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
    ph_with_vg(code = print(plot),
                  type = "body",
                  width = width, # Width and heigh define the size of the plot on powerpoint
                  height = height,
                  pointsize = pointsize,
               offx = 0,
               offy = 0
                  ) %>%
    print(target = path)
}

# Sept 2020 new dml object

# create_pptx <- function(plot, path, width = 6, height = 6, pointsize = 12){
#   library("officer")
#   library("magrittr")
#   library("rvg")
#   if(!file.exists(path)) {
#     out <- read_pptx()
#   } else {
#     out <- read_pptx(path)
#   }
#   
#   out %>%
#     add_slide(layout = "Title and Content", master = "Office Theme") %>%
#     ph_with(dml(code = print(plot)), # rvg:::ph_with.dml(dml(code = print(plot)),
#                       ph_location_fullsize()
#                       ) %>%
#     print(target = path)
# }




#### below methods not recommended

# Function 2. insert ggplot2 objects as bitmap

  create_pptx_gg <- function(plot, path, width = 6, height = 6, pointsize = 12){
    library(officer)
    library(magrittr)
    library(rvg)
    library(ggplot2)
    if(!file.exists(path)) {
      out <- read_pptx()
    } else {
      out <- read_pptx(path)
    }

    out %>%
      add_slide(layout = "Title and Content", master = "Office Theme") %>%
      # ph_with_vg(code = print(plot),
      #            type = "body",
      #            width = width, # Width and heigh define the size of the plot on powerpoint
      #            height = height,
      #            pointsize = pointsize,
      #            offx = 0,
      #            offy = 0
      # ) %>%
      ph_with_gg(value = plot) %>% ### Embeded a bitmap rater than a vector.
      print(target = path)
  }

# Dunction 3. Manually put code within "{}".
#
read_pptx() %>%
  add_slide(layout = "Title and Content", master = "Office Theme") %>%
  ph_with_vg(code = {},
             # type = "body",
             # width = width, # Width and heigh define the size of the plot on powerpoint
             # height = height,
             # pointsize = pointsize,
             offx = 0,
             offy = 0
  ) %>%
  print(target = "Rplots.svg2.pptx")

## Function 3 example #####
read_pptx() %>%
  add_slide(layout = "Title and Content", master = "Office Theme") %>%
  ph_with_vg(code = {plot(dtp1$diff ~ dtp1$mu,
                          pch = ".",
                          xlab = "Mean",
                          ylab = "Difference",
                          main = "Proportion of Methylated Reads in Control\nat Week 15 vs. Week 2, FDR < 0.1")
    points(dtp1$diff[dtp1$fdrs < 0.1 & dtp1$diff > 0] ~ dtp1$mu[dtp1$fdrs < 0.1 & dtp1$diff > 0] ,
           pch = "x",
           col = "green")
    points(dtp1$diff[dtp1$fdrs < 0.1 & dtp1$diff < 0] ~ dtp1$mu[dtp1$fdrs < 0.1 & dtp1$diff < 0] ,
           pch = "x",
           col = "red")
    abline(h = c(-0.2, 0.2),
           lty = 2)},
             # type = "body",
             # width = width, # Width and heigh define the size of the plot on powerpoint
             # height = height,
             # pointsize = pointsize,
             offx = 0,
             offy = 0
  ) %>%
  print(target = "Rplots.svg2.pptx")


#











# https://www.datacamp.com/community/tutorials/15-questions-about-r-plots





create_pptx(p1, "data/tmp333.pptx")



doc <- read_pptx()
doc <- add_slide(doc, "Title and Content", "Office Theme")
doc <- ph_with_vg(doc, code = barplot(1:5, col = 2:6), type = "body")
doc <- add_slide(doc, "Title and Content", "Office Theme")
doc <- ph_with_vg_at(doc, code = barplot(1:5, col = 2:6),
                     left = 1, top = 2, width = 6, height = 4)
print(doc, target = "data/tmp/vg.pptx")
