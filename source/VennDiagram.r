#Drew Venn Diagram with all values 
#install.packages("VennDiagram")
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

#Draw a Venn Diagram with 4 sets:
venn.plot <- draw.quad.venn(
  area1 = 72,
  area2 = 86,
  area3 = 50,
  area4 = 52,
  n12 = 44,
  n13 = 27,
  n14 = 32,
  n23 = 38,
  n24 = 32,
  n34 = 20,
  n123 = 18,
  n124 = 17,
  n134 = 11,
  n234 = 13,
  n1234 = 6,
  category = c("First", "Second", "Third", "Fourth"),
  fill = c("orange", "red", "green", "blue"),
  lty = "dashed",
  cex = 2,
  cat.cex = 2,
  cat.col = c("orange", "red", "green", "blue")
);
# Writing to file
tiff(filename = "Quad_Venn_diagram.tiff", compression = "lzw");
grid.draw(venn.plot);
dev.off()

#And more, see http://www.stats.bris.ac.uk/R/web/packages/VennDiagram/VennDiagram.pdf

