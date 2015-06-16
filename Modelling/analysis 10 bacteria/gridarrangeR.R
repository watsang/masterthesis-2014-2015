arrangegrids<-function(listggplots)
{
  #this function takes as an argument a list object of ggplot objects
  #TODO: implement sanity check for listggplots object: are all constiuents ggplot objects
  #TODO: minimal working example
  plots <- listggplots
  grobs <- list()
  widths <- list()
  # collect the widths for each grob of each plot
  
  for (i in 1:length(plots)){
      grobs[[i]] <- ggplotGrob(plots[[i]])
      widths[[i]] <- grobs[[i]]$widths[2:5]
  }
  # use do.call to get the max width
  
  maxwidth <- do.call(grid::unit.pmax, widths)
  # asign the max width to each grob
  
  for (i in 1:length(grobs)){
       grobs[[i]]$widths[2:5] <- as.list(maxwidth)
  }
  plot
  
  do.call("grid.arrange", c(grobs, ncol = 1))
}