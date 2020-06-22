#' Function to plot time-course label enrichment both as average_labeling and mole_equivalents_labeling
#' in a non-stationary isotopic labeling experiment
#'
#' This functions allows plotting of the average_labeling and mol_equivalent_labeling objects for user interpretation of label enrichment over the course of the experiment
#'
#' @import data.table
#' @import purrr
#' @import stringr
#' @import devtools
#' @import ggfortify
#' @import ggplot2
#' @import ggpubr
#' @import gridExtra
#' @import grid
#' @import lattice
#' @import RColorBrewer
#' @param mydata object (1) The 'average_labeling' or the 'mol_equivalent_labeling' object created using get_table_objects() function
#' @param Category string (1) name of the category by which you are binning, for the example data set - Lipidomics has two categories (i.e. tissue types),
#' "Cotyledon" and "EA" and the DualLabel dataset has two catergories (i.e. genotypes), "WT" and "GAT"
#' @param ylim numeric (1) sets the limit of y-axis, default is set to maximum
#' @param xlim numeric (1) sets the limit of the x-axis, default is set to maximum
#' @param axisTitle string (1) The label to be used for the Y-Axis. Typically it is "% labeling" for average_labeling and "mol equivalents of label" for mol_equivalent_labeling
#' @param plotTitle string (1) The title to be used for the plot - has to be one of the columns within the object.
#' Typically, it is "Compound", however, if the string is too long in the 9x9 output generated, "Bin" can be used as a substitute
#' @param plotTitle2 string (1) Additional title to be used for the plot, default is set to NULL
#' @param outputName string (1) The name to be appended to the the output pdf
#' @keywords untargeted metabolomics, stable isotopes, non-stationary isotopic labeling, dual labels, MS1
#' @seealso PCA_and_heatmap(), MIDplot(), getClustersAndPlots(), get_table_objects()
#' @return The function returns a pdf file with plots for label enrichment for all the compounds identified within the XCMS_data
#' @author Shrikaar Kambhampati, Allen Hubbard
#' @export
#' @examples
#' Average_Labeling <- label_enrichment_plot(test_13C15NGlutamine$average_labeling, Category = "WT", yLim=NULL,xLim=NULL, axisTitle=" % Labeling", plotTitle2= "Bin", outputName = "Average")
#' Label_Enrichment <- label_enrichment_plot(test_13C15NGlutamine$mol_equivalent_labeling, Category = "WT", yLim=NULL,xLim=NULL, axisTitle="mol equivalents of Label", plotTitle="Compound", plotTitle2= "Bin", outputName = "mol Equi")

##########do the plotting
label_enrichment_plot <- function(mydata, Category, yLim=NULL,xLim=NULL, axisTitle="Labeling", plotTitle="Bin", plotTitle2=NULL, outputName = "_")
{

  #print("this is allen's june 19th version")
  #pull out the first field which correpsonds to timepoint of non-stationary labeling experiment
  data_clean <- function(x) sapply (strsplit(x , '[_]' ), `[` , 1)

  #read in the data table or use the object from previous function
  mydata = mydata
  mydata = na.omit(mydata)


  #set the axes limit
  yLim = yLim
  xLim = xLim

  #get the outputName
  outputName = outputName

  #set the axisTitle
  axisTitle = axisTitle
  plotTitle = plotTitle
  plotTitle2 = plotTitle2

  #workaround to the yLim resetting issue
  #we reset yLim to the max( avg + sd) + 5
  #and then it is no longer null, so we need to know
  #if it was originally null
  YoriginallyNull = FALSE
  if(is.null(yLim) == TRUE)
  {
    YoriginallyNull = TRUE
  }

  XoriginallyNull = FALSE
  if(is.null(xLim) == TRUE)
  {
    XoriginallyNull = TRUE
  }

  #combine the multiple columns into a single column with unique names
  names = make.names(paste(mydata$Bin,mydata$rt,mydata$mz,sep = "_"), unique = TRUE)
  rownames(mydata) = names

  #what is the identifier for this group (i.e. the category, a condition or tissue group)
  Category = Category

  #pull out the condition/species specific information
  #by only those columns which match the Category
  mydata =mydata[, !(colnames(mydata) %in% colnames(mydata)[(colnames(mydata) %like% Category)])]

  #generate the pdf
  name = paste(Category,outputName,"label_enrichment.pdf", sep = "_")
  pdf(name)

  par(mfrow=c(3,3))
  for(i in 1:nrow(mydata))
  {

    #reset the plotting frame every 9 plots
    #use modular division to determine if we're at the 9th
    if(i %% 9 == 0)
    {
      par(mfrow=c(3,3))
    }
    theVecOfInfo = mydata[i,]


    #if using id's for title
    #set title before we make it numeric
    title = theVecOfInfo[colnames(theVecOfInfo) == plotTitle][1,]


    if(!is.null(plotTitle2))
    {
      title = paste(theVecOfInfo[colnames(theVecOfInfo) == plotTitle][1,], theVecOfInfo[colnames(theVecOfInfo) == plotTitle2][1,], sep = "  ")
    }

    ##calculate each time point and each std dev

    #get the set of times from the column names
    myTimes = data_clean(as.character(colnames(mydata)))
    #pull out any non-numeric elements of the timepoints
    myTimes =  as.numeric(gsub("X","",as.character(myTimes)))
    #remove any non-valid timepoints
    myTimes = myTimes[is.na(myTimes) == FALSE]
    #get a unique vector of the timepoints
    myTimesUnique = unique(myTimes)

    #hold all the Avgs and standard deviations
    allMyInfoAvg = vector()
    allMyInfoSD = vector()

    ###get the average and the sd for each of the metabolites
    #for each time
    for(k in 1:length(myTimesUnique))
    {

      #get the time
      myTime = myTimesUnique[k]

      #pull out first field,i.e time
      myTimepointMatch = gsub("(*.*)_.*_.*", "\\1", names(theVecOfInfo))

      #pull out non-numeric characters
      myTimepointMatch = gsub("[^0-9.-]", "", myTimepointMatch)

      #rename the colnames with just the time information
      #I think that this is redundant, actually
      names(theVecOfInfo) = gsub("[^0-9.-]", "", names(theVecOfInfo))

      #make sure that it's numeric now
      myTimepointMatch = as.numeric(myTimepointMatch)

      #set the names to just the timepoint now, it's all we need
      names(theVecOfInfo) = as.numeric(myTimepointMatch)

      #calculate the mean for this timepoint
      myTimeMean = mean(as.numeric(theVecOfInfo[,as.numeric(colnames(theVecOfInfo)) %in% c(myTime)]))

      #calculate the SD for this timepoint
      myTimeSD = sd(as.numeric(theVecOfInfo[,as.numeric(colnames(theVecOfInfo)) %in% c(myTime)]))

      allMyInfoAvg = c(allMyInfoAvg, myTimeMean)
      allMyInfoSD = c(allMyInfoSD, myTimeSD)
    }

    if(is.null(yLim) == TRUE)
    {
      yLim = max(allMyInfoAvg) + max(allMyInfoSD) + 5
      #print(yLim)
      #print("is yLim")
      #print(max(allMyInfoAvg+allMyInfoSD))
      #print("other potential max")
    }


    #if the user doesn't specify the x-limit, we're going
    #to use the max as the highest timepoint
    if(is.null(xLim) == TRUE)
    {
      xLim = max(myTimesUnique)
    }

    #set the y-limit to at least the value of the highest error bar
    if(yLim == "default")
    {
      yLim = max(allMyInfoAvg[!is.na(allMyInfoAvg)]) + allMyInfoAvg+allMyInfoSD

    }



    #declare the times which will be the x-axis of the plots
    x = myTimesUnique

    #make the plots (time on the x-axis)
    #the label enrichment on the y-axis

    #if we have an xLim but not an yLim
    if(is.null(xLim) == FALSE & is.null(yLim) == TRUE)
    {
      plot(x, allMyInfoAvg,
           xlim=range(c(0, xLim)),
           pch=19, xlab="Time", ylab=axisTitle,
           main=title
      )
    }

    #if we have a yLim but not an xLim
    if(is.null(yLim) == FALSE & is.null(xLim) == TRUE)
    {
      plot(x, allMyInfoAvg,
           ylim=range(c(0, yLim)),
           pch=19, xlab="Time", ylab=axisTitle,
           main=title
      )
    }

    #if we have a yLim and an xLim
    if(is.null(yLim) == FALSE & is.null(xLim) == FALSE)
    {
      plot(x, allMyInfoAvg,
           ylim=range(c(0, yLim)), xlim=range(c(0, xLim)),
           pch=19, xlab="Time", ylab=axisTitle,
           main=title
      )
    }

    #add in the error bars as specified by the clever
    #stack overflow post
    # hack: we draw arrows but with very special "arrowheads"
    arrows(x, allMyInfoAvg-allMyInfoSD, x, allMyInfoAvg+allMyInfoSD, length=0.05, angle=90, code=3)

    #reset yLlim to null if needed
    if(YoriginallyNull == TRUE)
    {
      yLim = NULL
    }

    #reset xLlim to null if needed
    if(XoriginallyNull == TRUE)
    {
      xLim = NULL
    }

  }
  dev.off()

}
