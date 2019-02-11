return_overlap <- function(x1=8, x2=9, y1=0, y2=8.3) {
  # Function takes input minimum and maximum and target min and max,
  # then returns the amount of overlap between the sections.
  return(max(min(x2,y2) - max(x1,y1), 0))
}