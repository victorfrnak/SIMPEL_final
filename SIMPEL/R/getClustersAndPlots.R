#' Function to do statistical analyses and create plots from non-stationary labeling experiments
#'
#' This function allows you to performs kmeans clustering analyses using label enrichment data
#' @param mydata1 will be your average labeling OR your mole equivalent average labeling
#' @param mydata2 will be your MIDs or scaled MIDs
#' @param Category name of the category by which we'll be binning
#' @param  nClust how many kmeans clusters
#' @param  labels label datapoints by either the "Bin" or "Compound" column
#' @param doMIDs do we want to do the plots for the MIDs as well?
#' @param outputName This will be the name that will be appended to the the output pdf
#' @keywords untargeted metabolomics, stable isotopes, non-stationary isotopic labeling, dual labels, MS1
#' @export
#' @examples
#' getClustersAndPlots()


#rownames(tableObjectToUse$mol_equivalent_labeling) = tableObjectToUse$mol_equivalent_labeling[labels]
#rownames(tableObjectToUse$mol_equivalent_labeling) = tableObjectToUse$mol_equivalent_labeling[labels]
#getClusterAndPlots(tableObjectToUse$mol_equivalent_labeling, tableObjectToUse$scaled_MIDs,"GAT", 5,5)
#getClusterAndPlots(testWholePackage$mol_equivalent_labeling, testWholePackage$scaled_MIDs,"GAT", 5,5)
getClusterAndPlots <- function(mydata1, mydata2, Category, nClust=0, labels="Bin", doMIDs=FALSE, outputName = "_")
{

  ##function that will create all the dataframes
  createDataframe = function(columnnames)
  {
    columnnames = columnnames
    lengthColNames = length(columnnames)
    df =  data.frame(matrix(nrow = 0, ncol = lengthColNames))
    colnames(df) = columnnames
    return(df)
  }


  #pull out the first field which correpsonds to timepoint of non-stationary labeling experiment
  data_clean <- function(x) sapply (strsplit(x , '[_]' ), `[` , 1)


  listOfAllOutputs = list()
  #first things first, make sure that we've used the right rownames
  #for the input tables

  #how many iterations for the kmeans - 25 should be more than enough for convergence
  howmanyN = 25
  mydata1 = mydata1
  mydata2 = mydata2
  labels = labels
  outputName = outputName

  #labels = "Bin"
  rownames(mydata1) = mydata1[,colnames(mydata1) %in% c(labels)]
  rownames(mydata2) = make.unique( mydata2[,colnames(mydata2) %in% c(labels)])

  #get the category name, pare down these variables, note well
  Category = Category
  #how many clusters for the kmeans

  #if nClust is 0, we will implement the sqrt
  nClust = nClust

  #subset the table by just the columns matching our category of interest
  catSubset = mydata1[,colnames(mydata1) %like% Category]

  #scale the data from 0 to 1 for each row of the table.
  #This will enable comparison between different clusters
  for(row in 1:nrow(catSubset))
  {
    #get the max
    myMax = max(catSubset[row,])

    #process each entry
    for(column in 1:ncol(catSubset))
    {
      #scale every entry by this max
      toUpdate = catSubset[row,column] / myMax
      #replace with the scaled value
      catSubset[row,column] = toUpdate
    }
  }

  #there'll be some na's so set them to zero
  catSubset[is.na(catSubset)] = 0

  #if the user is using the default kmeans setting
  #we'll using the sqrt formula to determine the number of clusters
  if(nClust == 0)
  {
    nClust = as.integer(sqrt(nrow(catSubset) / 2))
  }

  #go ahead and do the kmeans right away on the normalized table
  kmeans_object = kmeans(catSubset,nClust,nstart = howmanyN)

  #bring the kmeans and the PCA in here
  #have one object to store all of the plots
  plotsList <- list()


  #this data frame is going to have the global information for a complete k-mean cluster
  #across all timepoints
  dfAllClusters = createDataframe(c("Time","mean","sd","Cluster"))
  dfAllClusters$Time = as.numeric(dfAllClusters$Time)
  dfAllClusters$mean = as.numeric(dfAllClusters$mean)
  dfAllClusters$sd = as.numeric(dfAllClusters$sd)
  dfAllClusters$Cluster = as.numeric(dfAllClusters$Cluster)



  #once we've done the PCA's
  #progress to analyzing the kmeans analysis

  #for each of the clusters, we will want to know all of the data
  #to plot

  #for each kmeans cluster in the normalized label enrichment table
  #we're going to subset by just the compounds in that cluster
  #then we're going to get all the MIDs for these compounds
  #and we will further cluster them

  #we will collect information for all the metabolites + MID's in each cluster
  #through plots that focus on all of the compounds in a single cluster
  #as well as all the clusters to reflect the metabolism of pooled, similar compounds
  for(i in 1:nClust)
  {
    #cluster specific information
    theCluster = i

    #we reuse i before making it up here, so have a placeholder for the cluster number
    #since we will need it after the resetting
    theClusterAvgNum = i

    #plot out everything associated with each cluster
    #this is the kmeans results on the label enrichment table
    toSubsetTable = names(kmeans_object$cluster[which(kmeans_object$cluster == i)])

    #for each cluster, pull out just those compounds which fall in that cluster
    for_plottingFurther = subset(catSubset, rownames(catSubset) %in% toSubsetTable)

    #get their MIDs
    #make sure that the naming is consistent with previous table


    #moved this from 201
    #this is softcoded to get all of the timepoints from the column names
    vecOfExpTimes = unique(data_clean(colnames(for_plottingFurther)))


    ####get all the timepoint specific information about each cluster
    for(timepoint in 1:length(vecOfExpTimes))
    {
      whichTimepoint = vecOfExpTimes[timepoint]
      #this is going to be the avg's and standard deviation for all of the compounds associated with a cluster
      myMeanToAdd = mean(rowMeans(for_plottingFurther[,colnames(for_plottingFurther) %like% whichTimepoint]))
      mySdToAdd = sd(rowMeans(for_plottingFurther[,colnames(for_plottingFurther) %like% whichTimepoint]))
      dfAllClusters = rbind(dfAllClusters, data.frame(Time=whichTimepoint, mean=myMeanToAdd, sd=mySdToAdd,Cluster=theClusterAvgNum))
    }


    #if we want to do the MIDs as well
    if( doMIDs == TRUE)
    {
      #we have to match up the label enrichment and the MIDs
      #only pull out those MIDs  of compounds within the cluster
      for_plottingMIDs <- mydata2[ which(mydata2[,colnames(mydata2) %in% c(labels)] %in% rownames(for_plottingFurther)), ]

      #get the bin info
      justBinInfo = for_plottingMIDs[labels]

    #make sure that it's unique
      justBinInfo = make.names(justBinInfo, unique = TRUE)

      for_plottingMIDs[labels] = NULL


    ####get all the timepoint specific information about each cluster
      #for(timepoint in 1:length(vecOfExpTimes))
      #{
      #  whichTimepoint = vecOfExpTimes[timepoint]

        #this is going to be the avg's and standard deviation for all of the compounds associated with a cluster
      #  myMeanToAdd = mean(rowMeans(for_plottingFurther[,colnames(for_plottingFurther) %like% whichTimepoint]))
      #  mySdToAdd = sd(rowMeans(for_plottingFurther[,colnames(for_plottingFurther) %like% whichTimepoint]))
      #  dfAllClusters = rbind(dfAllClusters, data.frame(Time=whichTimepoint, mean=myMeanToAdd, sd=mySdToAdd,Cluster=theClusterAvgNum))
      #}


    #get the number of clusters for the MIDs kmeans
      numberOfClusters = as.integer(sqrt(nrow(for_plottingMIDs) / 2))

    #make sure that our MIDs are just our Category of interest
      for_plottingMIDs = for_plottingMIDs[,colnames(for_plottingMIDs) %like% Category]

    #if we have any Na's set them to zero
      for_plottingMIDs[is.na(for_plottingMIDs)] = 0

    #now, do kmeans on the many MIDs associated
    #with the Avg concentrations that cluster together
      kmeans_MIDs_objects = kmeans(for_plottingMIDs, numberOfClusters,nstart = howmanyN)

    #as we iterate through clusters, this list will store our plots
      allAverageAndMIDSCluster = list()

    #store each MID cluster's information

    #for each of the clusters of MIDs
    #pool together all of the information
    #of the MIDs in that cluster
      dfAllClustersMIDs <- data.frame(
      Time=numeric(),
      mean=numeric(),
      sd=numeric(),
      Cluster=factor()
      )

    #within each cluster of MIDs
    #keep track of all of the abundance data
      dfMIDs <- data.frame(
      Time=numeric(),
      mean=numeric(),
      sd=numeric(),
      Compound=factor(),
      cluster=factor()
      )

    #build up a list of everything for MIDs clusters
      MidsPlotList = list()

    ##for each MIDs cluster
    #we're going to go through and collect the information to fill in the data frames
      for(j in 1:length(unique(kmeans_MIDs_objects$cluster)))
      {

        toSubsetTableII = names(kmeans_MIDs_objects$cluster[which(kmeans_MIDs_objects$cluster == j)])
        for_plottingSub = subset(mydata2, rownames(mydata2) %in% toSubsetTableII)
        for_plottingSub= for_plottingSub[,colnames(for_plottingSub) %like% Category]

      #theMIDs cluster is here
        theCluster = j

      ####get all the specific information about each compounds cluster
        for(timepoint in 1:length(vecOfExpTimes))
        {
          whichTimepoint = vecOfExpTimes[timepoint]
          myMeanToAdd =  mean(rowMeans(for_plottingSub[,colnames(for_plottingSub) %like% whichTimepoint]))
          mySdToAdd =  sd(as.vector(as.matrix(for_plottingSub[,colnames(for_plottingSub) %like% whichTimepoint])))

        #get all the MIDs associated with that
          dfAllClustersMIDs = rbind(dfAllClustersMIDs, data.frame(Time=whichTimepoint, mean=myMeanToAdd, sd=mySdToAdd,Cluster=theCluster))
        }


      #Get the compound specific information in each MIDs cluster
        for(k in 1:length(rownames(for_plottingSub)))
        {
          theCompound = rownames(for_plottingSub)[k]
          for(timepoint in 1:length(vecOfExpTimes))
          {
            whichTimepoint = vecOfExpTimes[timepoint]

            myMeanToAdd =  mean(as.numeric(for_plottingSub[rownames(for_plottingSub) == theCompound,colnames(for_plottingSub) %like% whichTimepoint]))
            mySdToAdd =  sd(as.numeric(for_plottingSub[rownames(for_plottingSub) == theCompound,colnames(for_plottingSub) %like% whichTimepoint]))
          #we're looking at all of the MIDs in a single cluster
            dfMIDs = rbind(dfMIDs, data.frame(Time=whichTimepoint, mean=myMeanToAdd, sd=mySdToAdd, Compound=theCompound, cluster = j))
          }

        }

        theTotalName = paste(Category, "compounds cluster", theClusterAvgNum, "MIDs cluster", j)
        whichMIDsCluster = j
        toSubsetPlotMIDs = dfMIDs[dfMIDs$cluster == j,]
        toSubsetPlotMIDs$Time = as.numeric(gsub("X","",as.character(toSubsetPlotMIDs$Time)))

      #we will re-order after sorting by the slope

      ###do the slope and sorting now
      #at the mids and averages level
        for(l in 1:length(unique(toSubsetPlotMIDs$Compound)))
        {


        #extract all of the MIDs
          myCompound = unique(toSubsetPlotMIDs$Compound)[l]

        #get the slope through the midpoint time
          midTime = unique(toSubsetPlotMIDs)$Time[length(unique(toSubsetPlotMIDs$Time)) / 2]
          startTime = unique(toSubsetPlotMIDs$Time)[1]
          p1 = subset(toSubsetPlotMIDs, subset=(Time=="0" & Compound== myCompound))$mean
          p2 = subset(toSubsetPlotMIDs, subset=(Time==midTime & Compound== myCompound))$mean


          theSlope = p2 - p1

          toSubsetPlotMIDs[toSubsetPlotMIDs$Compound ==  myCompound, "ClusterSlope"]  = rep(theSlope,length(toSubsetPlotMIDs$Compound[toSubsetPlotMIDs$Compound == myCompound]))

        }

      ##sort by the slope
        toSubsetPlotMIDs = toSubsetPlotMIDs[order(toSubsetPlotMIDs$ClusterSlope,decreasing=T),]

        labeledSortVec = vector()
        for(i in 1:length(unique(toSubsetPlotMIDs$Compound)))
        {
          labeledSortVec=  c(labeledSortVec,rep(i,length(unique(vecOfExpTimes))))
        }

        toSubsetPlotMIDs$slopeSorted = labeledSortVec

      #set the factor levels based off of the reordered by slope, now
        toSubsetPlotMIDs$SortNames = toSubsetPlotMIDs$Compound %>% factor(levels = unique(toSubsetPlotMIDs$Compound))

      #plot all of the MIDs from this cluster of MIDs associated with the Avgs Cluster
        MidsPlotList[[j]] =  ggplot(toSubsetPlotMIDs, aes(x=Time,y=mean,colour=SortNames,group=SortNames)) + geom_line() + ggtitle(theTotalName) + geom_ribbon(aes(ymax=mean + sd, ymin=mean - sd, linetype=NA, fill = factor(SortNames)), show.legend = F, alpha=.1)
      }


      dfAllClustersMIDs$Cluster = as.factor(dfAllClustersMIDs$Cluster)
      MIDsAllClustersPlotName = paste(" Global MIDs cluster for average cluster", theClusterAvgNum)

    #plot everything for the MIDs cluster, now
      dfAllClustersMIDs$Time = as.numeric(gsub("X","",as.character(dfAllClustersMIDs$Time)))


    #for each cluster of MIDs, calculate the collective slope
      for(j in 1:length(unique(dfAllClustersMIDs$Cluster)))
      {

        p1 = subset(dfAllClustersMIDs, subset=(Time==unique(dfAllClustersMIDs$Time)[1] & Cluster==j))$mean
        p2 = subset(dfAllClustersMIDs, subset=(Time==unique(dfAllClustersMIDs$Time)[(length(unique(dfAllClustersMIDs$Time)) / 2)] & Cluster==j))$mean

        theSlope = p2 - p1
        dfAllClustersMIDs[dfAllClustersMIDs$Cluster == j, "ClusterSlope"]  = rep(theSlope,length(dfAllClustersMIDs$Cluster[dfAllClustersMIDs$Cluster == j]))

      #get the slop of the first to numerbs
      #lm(dfAllClusters[dfAllClusters$Cluster == j,]$mean ~ dfAllClusters[dfAllClusters$Cluster == j,]$Time)$coeff[[2]]


      }


      dfAllClustersMIDs = dfAllClustersMIDs[order(dfAllClustersMIDs$ClusterSlope,decreasing=T),]


      labeledSortVec = vector()

    #each cluster will contains the information
    #of many compounds
      for(i in 1:length(unique(dfAllClustersMIDs$Cluster)))
      {
        labeledSortVec=  c(labeledSortVec,rep(i,  length(unique(vecOfExpTimes))))
      }

      dfAllClustersMIDs$slopeSorted = labeledSortVec

    #make sure that we have the legends sorted here, now
      dfAllClustersMIDs$SortNames = dfAllClustersMIDs$Cluster %>% factor(levels = unique(dfAllClustersMIDs$Cluster))

      listToAdd1 = list()
      listToAdd1[[1]] = ggplot(dfAllClustersMIDs, aes(x=Time,y=mean,colour=SortNames,group=SortNames)) + geom_line() + ggtitle(MIDsAllClustersPlotName) + geom_ribbon(aes(ymax=mean + sd, ymin=mean - sd, linetype=NA, fill = factor(SortNames)), show.legend = F, alpha=.1)

    #store all of the MIDs plots in a list
      allAverageAndMIDSCluster = prepend(MidsPlotList,listToAdd1)
      plotsList = prepend(allAverageAndMIDSCluster,plotsList)
    }
  }

  ######End the MIDs here !!!! #####



  #determine the cluster by slope
  #working with the averages
  dfAllClusters$ClusterSlope = dfAllClusters$Cluster
  dfAllClusters$Cluster = as.factor(dfAllClusters$Cluster)

  #determine which cluster has the best slope
  #subset by each cluster
  for(j in 1:length(unique(dfAllClusters$Cluster)))
  {

    p1 = 0
    p2 = 0

    p1 = subset(dfAllClusters, subset=(Time==unique(dfAllClusters$Time)[1] & Cluster==j))$mean
    p2 = subset(dfAllClusters, subset=(Time==unique(dfAllClusters$Time)[(length(unique(dfAllClusters$Time)) / 2)] & Cluster==j))$mean



    theSlope = 0
    theSlope = as.numeric(p2) - as.numeric(p1)

    dfAllClusters[dfAllClusters$Cluster == j, "ClusterSlope"]  =  theSlope
  }

  forAllAvgsCluster = paste(Category, ": K-means Clustering")


  ##for now, we want the ordering to be consistent with the kmeans biplot
  dfAllClusters = dfAllClusters[order(dfAllClusters$Cluster,decreasing=F),]


  #we've re-ordered after sorting by the slope
  allTimesVector = vector()
  for(time in 1:nClust)
  {
    allTimesVector  = c(allTimesVector,rep(time,  length(unique(vecOfExpTimes))))
  }
  dfAllClusters$slopeSorted = allTimesVector

  dfAllClusters$slopeSorted = as.factor(dfAllClusters$slopeSorted)
  dfAllClusters$Cluster = as.factor(dfAllClusters$Cluster)
  dfAllClusters$Time = as.numeric(gsub("X","",as.character(dfAllClusters$Time)))



  ##add the overall kmeans plot of the averages
  #allTotClusters = prepend(allTotClusters,allClusterForBeginning)
  listToAdd = list()
  listToAdd[[1]] = autoplot(kmeans_object, data = catSubset, label = TRUE, label.size = 3, frame = TRUE, frame.type = 'norm', main = forAllAvgsCluster)
  listToAdd2 = list()
  listToAdd2[[1]] = ggplot(dfAllClusters, aes(x=Time,y=mean,colour=Cluster,group=Cluster)) + geom_line() + ggtitle(forAllAvgsCluster) + geom_ribbon(aes(ymax=mean + sd, ymin=mean - sd, linetype=NA, fill = factor(Cluster)), show.legend = F, alpha=.1)


  autoplotAndGlobalClusters = prepend(listToAdd2,listToAdd)


  #allTotClusters = prepend(allTotClusters,autoplotAndGlobalClusters)

  #list to store global kmeans plot
  #allClusterForBeginning = list()

  #global k means plot
  #allClusterForBeginning[[1]] = ggplot(dfAllClusters, aes(x=Time,y=mean,colour=Cluster,group=Cluster)) + geom_line() + ggtitle(forAllAvgsCluster) + geom_ribbon(aes(ymax=mean + sd, ymin=mean - sd, linetype=NA, fill = factor(Cluster)), show.legend = F, alpha=.1)

  ##prepend our global kmeans plots to the full list of plot
  allTotClusters = prepend(plotsList, autoplotAndGlobalClusters)



  pdf(paste(Category, outputName, "kmeans_plots.pdf", sep = "_"))
  for(i in allTotClusters){print(i)}
  dev.off()



  return(allTotClusters)
}



