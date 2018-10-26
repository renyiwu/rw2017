

create_pptx <- function(plot, path, width = 6, height = 6, pointsize = 12){
  library(officer)
  library(magrittr)
  library(rvg)
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



















# https://www.datacamp.com/community/tutorials/15-questions-about-r-plots





create_pptx(p1, "data/tmp333.pptx")



doc <- read_pptx()
doc <- add_slide(doc, "Title and Content", "Office Theme")
doc <- ph_with_vg(doc, code = barplot(1:5, col = 2:6), type = "body")
doc <- add_slide(doc, "Title and Content", "Office Theme")
doc <- ph_with_vg_at(doc, code = barplot(1:5, col = 2:6),
                     left = 1, top = 2, width = 6, height = 4)
print(doc, target = "data/tmp/vg.pptx")