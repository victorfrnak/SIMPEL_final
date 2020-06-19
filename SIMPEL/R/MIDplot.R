#' Function to plot Mass Isotopologue Distribution, i.e. MIDs, in a non stationary isotopic
#' labeling experiment
#'
#' @param myData you can use MIDs or scaled_MIDs as an input for this function
#' @param Category will be the tissue type, i.e. WT/Mut, or treatment, i.e treated vs untreated
#' @param ylim is the limit of y-axis
#' @param xlim is the limit of the x-axis
#' @param axisTitle This will be name for the Y-Axis typically it is "MIDs" for MIDs and "mole equivalents of signals" for scaled_MIDs)
#' @param plotTitle This will be name for the Y-Axis typically it is "% label  enrichment" for average_labeling and "mole_equivalents" for mole_equivalents_labeling)
#' @param plotTitle2 This will be name for the Y-Axis typically it is "% label  enrichment" for average_labeling and "mole_equivalents" for mole_equivalents_labeling)
#' @keywords untargeted metabolomics, stable isotopes, non-stationary isotopic labeling, dual labels, MS1
#' @export
#' @examples
#' MIDplot()

#MIDplot(tableObjectToUse$scaled_MIDs, "GAT", "C0N0",axisTitle = "Mean", plotTitle = "customize")
#MIDplot(tableObjectToUse$scaled_MIDs, "GAT", "C0N0",axisTitle = "Mean", plotTitle = "customize",plotTitle2 = NULL,c(5,10))

#######Do all the plottings for the MIDs
#MIDplot(tableObjectToUse$scaled_MIDs, "GAT", "default","C0N0","default")


MIDplot <- function(myData, Category, splitIsotopologue = "C0N0", axisTitle="MID", plotTitle="Bin", plotTitle2=NULL, yLimit="default", xLimit="default")
{


  #pull out the first field which correpsonds to timepoint of non-stationary labeling experiment
  data_clean <- function(x) sapply (strsplit(x , '[_]' ), `[` , 1)


  #default options
  splitIsotopologue = splitIsotopologue
  yLimit = yLimit
  myData = myData
  plotTitle = plotTitle
  plotTitle2 = plotTitle2
  axisTitle = axisTitle
  Category = Category
  #print(yLimit)
  #print(" is yLimit")

  #focus on the just the columns with the category/tissue of interest
  myData = myData[, !(colnames(myData) %in% colnames(myData)[(colnames(myData) %like% Category)])]

  #make the rownames all the relevant identifying information (bin, retention time, m.z)
  #but make sure these are unique (just in case) using the make.names function
  names = make.names(paste(myData$rt,myData$mz,myData$Bin, sep = "_"), unique = TRUE)
  rownames(myData) = names
  #list of the plots
  p<-list()
  #plots with all of them together
  pAll <- list()
  #loop through each bin
  for(i in 1:length(unique(myData$Bin)))
  {

    theBin = unique(myData$Bin)[i]
    df <- data.frame(Isotopologue=character(),
                     Time=numeric(),
                     Mean=numeric(),
                     stdDev=numeric())
    subsetBin = subset(myData, Bin  == theBin)
    subsetBin$Isotopologue =as.factor(make.unique(as.character(subsetBin$Isotopologue), sep = "."))
    ###this is where we get the mass and retention time
    #to format
    theMass = subsetBin[subsetBin$Isotopologue == splitIsotopologue,]$m.z
    retentTime = subsetBin[subsetBin$Isotopologue == splitIsotopologue,]$rt

    #CompName = subsetBin[subsetBin$Isotopologue == splitIsotopologue,]$Bin
    #I've updated this on June 11th 2020 to make sure that the user can customize the column
    #they want to use in order to plot
    CompName = subsetBin[subsetBin$Isotopologue == splitIsotopologue,][plotTitle]

    title = paste(CompName)

    #if the user wants to include another column we'll include that in the title
    if(is.null(plotTitle2) == FALSE)
    {
      title = paste(title, subsetBin[subsetBin$Isotopologue == splitIsotopologue,][plotTitle2], sep = " ")
    }

    #take all the data and put into a new table
    #we well need to have this data formatted in order
    #to be compatible with the ggplot2 function
    for(j in 1:nrow(subsetBin))
    {

      #get the row with all of the information
      myVecOfAllInfo = subsetBin[j,]

      #determine how many timepoint

      ##calculate each time point mean and each std dev
      myTimes = data_clean(as.character(colnames(myData)))
      myTimes =  as.numeric(gsub("X","",as.character(myTimes)))
      myTimes = myTimes[is.na(myTimes) == FALSE]
      myTimesUnique = unique(myTimes)

      if(xLimit == "default")
      {
        xLim = max(myTimesUnique)
      }

      isotopologueName = as.character(myVecOfAllInfo$Isotopologue)

      #make sure that we get the information for each timepoint
      for(k in 1:length(myTimesUnique))
      {
        myTime = myTimesUnique[k]

        #pull out first field
        myTimepointMatch = gsub("(*.*)_.*_.*", "\\1", names(myVecOfAllInfo))

        #pull out non-numeric stuff
        myTimepointMatch = gsub("[^0-9.-]", "", myTimepointMatch)
        names(myVecOfAllInfo) = gsub("[^0-9.-]", "", names(myVecOfAllInfo))

        #make sure that it's numeric now
        myTimepointMatch = as.numeric(myTimepointMatch)

        #set the names to just the timepoint now, it's all we need
        names(myVecOfAllInfo) = as.numeric(myTimepointMatch)

        myTimeMean = mean(as.numeric(myVecOfAllInfo[,as.numeric(colnames(myVecOfAllInfo)) %in% c(myTime)]))
        myTimeSD = sd(as.numeric(myVecOfAllInfo[,as.numeric(colnames(myVecOfAllInfo)) %in% c(myTime)]))
        df = rbind(df, data.frame(Isotopologue=isotopologueName, Time=myTime, Mean=myTimeMean, stdDev=myTimeSD))
      }
    }

    #just M0
    smallDf = subset(df, Isotopologue  == splitIsotopologue)
    #everything except M0
    largerDF = subset(df, Isotopologue  != splitIsotopologue)

    #store all of these in the list of images
    #we will have to plot by default the entire range
    #as well as using the limits provided by the user if they want to use the option
    if(nrow(largerDF) > 0)
    {
      if(length(yLimit) > 1)
      {

        p[[i]] = ggarrange(qplot(Time, Mean, data=smallDf, colour=Isotopologue, geom=c("line","point")) +
                             geom_errorbar(aes(ymin=Mean-stdDev, ymax=Mean+stdDev), width = .3) + ggtitle(title) + ylim(yLimit) + ylab(axisTitle) ,qplot(Time, Mean, data=largerDF, colour=Isotopologue, geom=c("line","point")) +
                             geom_errorbar(aes(ymin=Mean-stdDev, ymax=Mean+stdDev, width = .3)) + ggtitle(title) + ylim(yLimit) + ylab(axisTitle))

      }


      if(length(yLimit) == 1)
      {

        p[[i]] = ggarrange(qplot(Time, Mean, data=smallDf, colour=Isotopologue, geom=c("line","point")) +
                             geom_errorbar(aes(ymin=Mean-stdDev, ymax=Mean+stdDev), width = .3) + ggtitle(title) ,qplot(Time, Mean, data=largerDF, colour=Isotopologue, geom=c("line","point")) +
                             geom_errorbar(aes(ymin=Mean-stdDev, ymax=Mean+stdDev, width = .3)) + ggtitle(title))
      }
    }
    if(nrow(largerDF) == 0)
    {

      if(length(yLimit) == 1)
      {

        #print("we've got null for ylimit")
        p[[i]] = qplot(Time, Mean, data=df, colour=Isotopologue, geom=c("line","point")) + geom_errorbar(aes(ymin=Mean-stdDev, ymax=Mean+stdDev), width = .3) + ggtitle(title)


        #p[[i]] = qplot(Time, Mean, data=df, colour=Isotopologue, geom=c("line","point")) + geom_errorbar(aes(ymin=Mean-stdDev, ymax=Mean+stdDev), width = .3) + ggtitle(title)
      }

      if(length(yLimit) > 1)
      {
        print("we have a window")

        p[[i]] = qplot(Time, Mean, data=df, colour=Isotopologue, geom=c("line","point")) + geom_errorbar(aes(ymin=Mean-stdDev, ymax=Mean+stdDev), width = .3) + ggtitle(title) + ylim(yLimit) + ylab(axisTitle)

        #p[[i]] = qplot(Time, Mean, data=df, colour=Isotopologue, geom=c("line","point")) + geom_errorbar(aes(ymin=Mean-stdDev, ymax=Mean+stdDev), width = .3) + ggtitle(title) + ylab(axisTitle) + xlim(xmin,xmax) + ylim(ymin, ymax)
      }


    }

    #pAll[[i]] = qplot(Time, Mean, data=df, colour=Isotopologue, geom=c("line","point")) + geom_errorbar(aes(ymin=Mean-stdDev, ymax=Mean+stdDev), width = .3) + ggtitle(title)

    if(length(yLimit) == 1)
    {
      pAll[[i]] = qplot(Time, Mean, data=df, colour=Isotopologue, geom=c("line","point")) + geom_errorbar(aes(ymin=Mean-stdDev, ymax=Mean+stdDev), width = .3) + ggtitle(title) + ylab(axisTitle)

    }

    if(length(yLimit) > 1)
    {
      #print("we've got null for the no split")
      pAll[[i]] = qplot(Time, Mean, data=df, colour=Isotopologue, geom=c("line","point")) + geom_errorbar(aes(ymin=Mean-stdDev, ymax=Mean+stdDev), width = .3) + ggtitle(title) + ylab(axisTitle) + ylim(yLimit)

    }

  }


  #print out the side by side plots
  Split_MIDs <- marrangeGrob(p, nrow=2, ncol=1)

  splitMIDSname = paste(Category, "split_MIDs.pdf",sep = "_")
  MIDSnonSplitname = paste(Category, "MIDs.pdf",sep = "_")
  ggsave(splitMIDSname, Split_MIDs, width = 11, height = 11)
  #print out the all by all plots
  MIDs <- marrangeGrob(pAll, nrow=3, ncol=3)
  ggsave(MIDSnonSplitname, MIDs, width = 12, height = 12)
}


