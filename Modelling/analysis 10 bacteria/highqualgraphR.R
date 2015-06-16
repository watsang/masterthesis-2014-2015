highqualgraphR<-function(x,filename,res=1200,pointsize=12,extension = "pdf")
{
  #here x is a ggplot2 object
  #the filename is the filename without an extension
  #@extension: string "pdf", "png", "tiff", "svg" or "postscript"
  #TODO: minimal working example
  #TODO: sanity checks for x and filename
  
  
  #######################################################################
  # If you receive Error: 
  #     unknown family 'Arial' 
  # Then run the following code:
  # 
  # library("extrafont")
  # font_import()
  # loadfonts()
  # loadfonts(device = "postscript")
  # 
  #######################################################################
  
  # Create subfolder figure to dump your files
  working_directory <- getwd()
  subdirectory <- "figure/"
  dir.create(file.path(working_directory, subdirectory), showWarnings = FALSE)
  
  filename <- paste(subdirectory, filename, sep = "")
  
  # define function to print your graphs
  print_graph <- function(format){
    if (format == "pdf")
    {
      pdf(file=paste(filename,".pdf",sep=""), height = 6, width = 9,
          family = "Arial", paper = "special", onefile = FALSE,
          colormodel="cmyk")
    }
    else if (format == "png")
    {
      png(filename=paste(filename,".png",sep=""),width=12000,height=7000,units="px",res=res,pointsize=pointsize,family="Arial",type="cairo")
    }
    else if (format == "tiff")
    {
      tiff(filename=paste(filename,".tiff",sep=""),width=12000,height=7000,res=res,pointsize=pointsize,family="Arial",type="cairo",compression="lzw")  
    }
    else if (format == "svg")
    {
      svg(filename=paste(filename,".svg",sep=""),width=9,height=6,pointsize=12,family="Arial")  
    }
    else if (format == "postscript")
    {
      postscript(file=paste(filename,".eps",sep=""), height = 6, width = 9,  
                 family = "Arial", paper = "special", onefile = FALSE,
                 horizontal = FALSE,colormodel="cmyk")
    }
    else
    {
      stop("Not a valid extension")
    }
    # Print your files
    print(x)
    dev.off()
  }
  
  # Check whether to print one or multiple plots
  if (is.vector(extension))
  {
    for (format in extension)
    {
      print_graph(format)
    }
  }
  else
  {
    print_graph(extension)
  }
}