#' Function to plot time-course label enrichment both as average_labeling and mole_equivalents_labeling
#'in a non-stationary isotopic labeling experiment
#'   XCMS data and start getting isotopic information from the annotation file
#' and add information columns to be used as input for next steps
#'
#' @param mydata you can use average_labeling or mole_equivalents_labeling as an input for this function
#' @param Category will be the tissue type, i.e. WT/Mut, or treatment, i.e treated vs untreated
#' @param ylim is the limit of y-axis
#' @param xlim is the limit of the x-axis
#' @param axisTitle This will be name for the Y-Axis typically it is "% label  enrichment" for average_labeling and "mole_equivalents" for mole_equivalents_labeling)

#' @param plotTitle This will be name for the Y-Axis typically it is "% label  enrichment" for average_labeling and "mole_equivalents" for mole_equivalents_labeling)
#' @param plotTitle2 This will be name for the Y-Axis typically it is "% label  enrichment" for average_labeling and "mole_equivalents" for mole_equivalents_labeling)


#' @keywords untargeted metabolomics, stable isotopes, non-stationary isotopic labeling, dual labels, MS1
#' @export
#' @examples
#' makeAvgsplot()

##########do the plotting
label_enrichment_plot <- function(mydata, Category, yLim=NULL,xLim=NULL, axisTitle1="Labeling", plotTitle=Bin, plotTitle2=NULL)
{



  #read in the data table
  mydata = mydata
  #parameter about to which percent we want to plot
  #for example, if yLim= 65 we will only plot up
  #until 65 percent
  yLim = yLim
  xLim = xLim

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

  #what is the identifier for this group (i.e. is it a condition or tissue group)
  Category = Category

  #pull out the condition/species specific information
  #by only those columns which match the Category
  mydata =mydata[, !(colnames(mydata) %in% colnames(mydata)[(colnames(mydata) %like% Category)])]

  #read all the data in
  name = paste(Category, "allen_home_averages_plotting_june11th_AvgsII.pdf", sep = "_")
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
    title = paste(theVecOfInfo[colnames(theVecOfInfo) == "ID"][1,], theVecOfInfo[colnames(theVecOfInfo) == "Bin"][1,])

    ##calculate each time point and each std dev

    #get the set of times from the column names
    myTimes = data_clean(as.character(colnames(mydata)))
    #pull out any non-numeric elements of the timepoints
    myTimes =  as.numeric(gsub("X","",as.character(myTimes)))
    #remove any non-valid timepoints
    myTimes = myTimes[is.na(myTimes) == FALSE]
    #get a unique vector of the timepoints
    myTimesUnique = unique(myTimes)
    #isotopologueName = as.character(myVecOfAllInfo$Isotopologue)

    #if the user doesn't specify the x-limit, we're going
    #to use the max as the highest timepoint
    if(is.null(xLim) == TRUE)
    {
      xLim = max(myTimesUnique)
    }

    #hold all the Avgs and standard deviations
    allMyInfoAvg = vector()
    allMyInfoSD = vector()

    ###get the average and the sd for each of the metabolites

    #for each time
    for(k in 1:length(myTimesUnique))
    {

      #get the time
      myTime = myTimesUnique[k]



      #pull out first field
      myTimepointMatch = gsub("(*.*)_.*_.*", "\\1", names(theVecOfInfo))

      #pull out non-numeric stuff
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
      #df = rbind(df, data.frame(Isotopologue=isotopologueName, Time=myTime, Mean=myTimeMean, stdDev=myTimeSD))
    }


    print(allMyInfoAvg)
    print("averages")
    print(allMyInfoSD)
    print("SD")
    print(yLim)
    print("is yLim")
    if(is.null(yLim) == TRUE)
    {
      yLim = max(allMyInfoAvg) + max(allMyInfoSD) + 5
      print(yLim)
      print("is yLim")
      print(max(allMyInfoAvg+allMyInfoSD))
      print("other potential max")
    }

    #declare the times which will be the x-axis of the plots
    x = myTimesUnique

    #make the plots (time on the x-axis)
    #the averages on the y-axis

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
