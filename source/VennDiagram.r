#Drew Venn Diagram with all values 
require(VennDiagram)
venn.plot <- draw.triple.venn(
  area1 = 49,
  area2 = 33,
  area3 = 70,
  n12 = 10,
  n23 = 10,
  n13 = 7,
  n123 = 3,
  category = c("First", "Second", "Third"),
  fill = c("blue", "red", "green"),
 lty = "blank",
 #lty = rep("solid", 3),
  cex = 2,
  cat.cex = 2,
 euler.d=TRUE,
 scaled=TRUE)
tiff(filename = "test3.tiff", compression = "none",type = "windows",antialias = "none")
grid.draw(venn.plot)
dev.off()
