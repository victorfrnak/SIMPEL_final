#' Function to process XCMS data for isotopic enrichment and create MID and average labeling objects
#'from non-stationary isotopic labeling experiments
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
#' @param XCMS_data this will be your edited XCMS output
#' @param ppm will be ppm for theupper and lower bounds
#' @param rt_tolerance will be the retention time window for the upper and lower bounds
#' @param compounds_data This will be your annotation file.
#' @param output This is the name which will be appended to the files if the user wishes to print out the analysis tables
#' @keywords untargeted metabolomics, stable isotopes, non-stationary isotopic labeling, dual labels, MS1
#' @export
#' @examples
#' get_table_objects()

#XCMS_data = read.table(file = "/Users/ahubbard/Downloads/xcms_data_DualLabel_SingleCategory.csv", sep = ",", header = TRUE)
#compounds_data = read.table(file = "/Users/ahubbard/Downloads/Compound_data_DualLabel_REA_june5thL.txt", sep = "\t", header = TRUE)
#compounds_data = head(compounds_data, n = 1)
#testFullSmallFile = get_table_objects(XCMS_data, compounds_data)

#XCMS_data = read.table(file = "/Users/ahubbard/Downloads/xcms_data_DualLabel_REAL_june10th.csv", sep = ",", header = TRUE)

#compounds_data = read.table(file = "/Users/ahubbard/Downloads/Compound_data_DualLabel_REAL_june10th.txt", sep = "\t", header = TRUE)

#XCMS_data = read.table(file = "/Users/ahubbard/Downloads/xcms_data_Lipidomics_REAL_june15th.csv", sep = ",", header = TRUE)


#XCMS_data = read.table(file = "/Users/ahubbard/Downloads/xcms_data_Lipidomics_REAL_june15th.csv", sep = ",", header = TRUE)

#compounds_data = read.table(file = "/Users/ahubbard/Downloads/Compound_data_Lipidomics_REAL_june15th.txt", sep = "\t", header = TRUE)



#compounds_data = head(compounds_data, n = 1)
#testFullSmallFile = get_table_objects(XCMS_data, compounds_data)

#testWholePackage = get_table_objects(XCMS_data, compounds_data)

get_table_objects <- function(XCMS_data, compounds_data, ppm=10, rt_tolerance=.1, output=NULL){

  #set the variables from the input, i.e. ppm
  ppm = ppm
  rt_tolerance = rt_tolerance
  XCMS_data = XCMS_data
  compounds_data = compounds_data

  #add the columns to the table that we want all the outputs to have
  XCMS_data <- XCMS_data %>%
    mutate(comp_result = NA)

  XCMS_data <- XCMS_data %>%
    mutate(carbon = NA)

  XCMS_data <- XCMS_data %>%
    mutate(nitrogen = NA)

  XCMS_data <- XCMS_data %>%
    mutate(total_isotopes = NA)

  XCMS_data <- XCMS_data %>%
    mutate(Bin = NA)

  XCMS_data <- XCMS_data %>%
    mutate(Compound = NA)

  XCMS_data <- XCMS_data %>%
    mutate(Formula = NA)

  #These will be the functions to separate the fields from the file labels
  #pull out the first field which correpsonds to timepoint of non-stationary labeling experiment
  data_clean <- function(x) sapply (strsplit(x , '[_]' ), `[` , 1)

  #pull out the second, or middle, field which correpsonds to category [treatment/mutant, etc]
  data_cleanMiddle <- function(x) sapply (strsplit(x , '[_]' ), `[` , 2)

  #pull out the ending field which correpsonds to replicate
  data_cleanEnd <- function(x) sapply (strsplit(x , '[_]' ), `[` , 3)


  #looping through the list of mass features in the XCMS data
  #to identfiy isotopologs and bins
  for(i in 1:nrow(compounds_data))
  {
    #for each compound we're gonna pull out formula and retention time
    #from the annotation file
    compounds_data[i, ] %>%
      unlist() %>%
      as.vector() %>%
      print()

    ##for each compound in the annotation file with rt and formula
    #determine its m/z
    comp_formula <- as.character(compounds_data[i, "formula"])

    #retention time and polarity
    r_time  <- as.numeric(compounds_data[i, "rt"])
    polarity = compounds_data$polarity[i]

    #send that m/z to look up all of the C and N isotopes
    #and get the upper and lower bounds on each m/z based on the specified ppm error

    #filter the compounds_data by polarity
    subset_comp_lookup_table = compounds_data[compounds_data$polarity == polarity,]
    comp_lookup_table <- get_comp_mz_lookup(subset_comp_lookup_table,  comp_formula, r_time, ppm, polarity)

    #information for the bin as found in the prefix column of the annotation
    #file
    myBin = as.character(compounds_data[i, "prefix"])
    myFormula = as.character(compounds_data[i, "formula"])
    myCompound = as.character(compounds_data[i, "Compound"])

    #looping through the xcms data

    #for each feature in the XCMS data, id'd by RT and m/z
    #determine if it is one of the C or N or a combination of isotopes by whether or not
    #it falls within the error in ppm that are stored in comp_lookup_table
    for(j in 1:nrow(XCMS_data)){
      if(is.na(XCMS_data[j, "comp_result"]))
      {
        #only make sure to search the correct polarity
        if(XCMS_data[j, "polarity"] == polarity)
        {
          #print(j)
          x = as.numeric(XCMS_data[j, "mz"])
          y = as.numeric(XCMS_data[j, "rt"])

          #determine if this feature, id'd by m/z and RT is one of the C and N isotopes
          #whose ppm error is stored in comp_lookup_table
          #print("in the second loop")
          #add the number of N's and the number of C's
          #add the isotopolog columns
          val <- get_comp_stage(x, y, comp_lookup_table, r_time, rt_tolerance)

          if(is.null(val) == FALSE)
          {
            print("val is false")

            XCMS_data[j, "comp_result"] <- val$compound
            XCMS_data[j, "carbon"] <- val$carbon
            XCMS_data[j, "nitrogen"] <- val$nitrogen

            isotopeNumbers = val$carbon$carbon + val$nitrogen$nitrogen

            XCMS_data[j, "total_isotopes"] <- isotopeNumbers
            XCMS_data[j, "Bin"] <- myBin
            XCMS_data[j, "Compound"] <- myCompound
            XCMS_data[j, "Formula"] <- myFormula


          }


          if(is.null(val) == TRUE)
          {
            XCMS_data[j, "comp_result"] <- NA
            XCMS_data[j, "carbon"] <- 0
            XCMS_data[j, "nitrogen"] <- 0
            XCMS_data[j, "total_isotopes"] <- 0
            XCMS_data[j, "Bin"] <- myBin
            #XCMS_data[j, "Formula"] <- myFormula
          }



        }
      }
    }
  }


  #create a subset table that only contains binned rows with
  #labeled data
  XCMS_data = XCMS_data[is.na(XCMS_data$comp_result) == FALSE,]


  #calculate the median normalized data
  forMedianNormalization = XCMS_data


  #we're going to exclude the column labels whose rows do not include the numeric data
  vecToExclude = c("mz","polarity", "rt", "comp_result","Formula", "carbon", "nitrogen" , "total_isotopes", "Bin", "Compound" )

  #this vector of the medians is going to be critical because we are going to use
  #and reuse it
  columnsToSearch= setdiff(colnames(forMedianNormalization), vecToExclude)

  #get the median for each column
  for(column in 1:length(columnsToSearch))
  {

    myColumn = columnsToSearch[column]
    myMedian = median(forMedianNormalization[,colnames(forMedianNormalization) ==  myColumn])

    #divide the intensity of each entry by the median of the column
    forMedianNormalization[,colnames(forMedianNormalization) ==  myColumn] = forMedianNormalization[,colnames(forMedianNormalization) ==  myColumn] / myMedian

  }



  #declare the proxyPoolTable variable.
  #Since we don't have pool size in metabolomics data we obtain proxy pools
  #by median normalization followed by summing all the isotopologues within the bin
  names <- c("Bin",columnsToSearch)
  proxyPoolTable <- data.frame(matrix(ncol = length(names), nrow = length(unique(forMedianNormalization$Bin))))
  x <-  names
  colnames(proxyPoolTable) <- x

  proxyPoolTable$Bin = unique(forMedianNormalization$Bin)

  #for each sample, sum up the measurements of the bin
  for(theBin in 1:length(unique(forMedianNormalization$Bin)))
  {


    theBinToSum = unique(forMedianNormalization$Bin)[theBin]

    theBinSumsToBind = c(theBinToSum)
    for(column in 1:length(columnsToSearch))
    {
      myColumn = columnsToSearch[column]
      sumTheIntensities = forMedianNormalization[forMedianNormalization$Bin == theBinToSum, colnames(forMedianNormalization) == myColumn]
      print(sumTheIntensities)


      mySumToAdd = sum(sumTheIntensities)
      proxyPoolTable[proxyPoolTable$Bin == theBinToSum, colnames(proxyPoolTable) == myColumn] = mySumToAdd


    }

  }

  #we will use this proxyPoolTable to adjust the MIDs


  #get the xcms data as input
  MIDs_table = XCMS_data

  #get the bins
  allBins = MIDs_table$Bin
  allBins = unique(allBins)

  #have a vector to keep track of all
  #the compounds without an MID.
  #We don't use it here, but we may in the future
  vecOfNoMIDs = vector()

  #calculate the average labeling description
  average_labeling = XCMS_data

  #we're going to have a column by which we bin
  #and we don't need a lot of the information in here
  average_labeling = subset(average_labeling,  select = -c(carbon, nitrogen, total_isotopes, comp_result))

  #only one average per bin
  average_labeling = average_labeling[!duplicated(average_labeling$Bin),]


  #subset for each bin and creae the MIDs table
  for(i in 1:length(allBins))
  {

    myBin = allBins[i]

    #focus on just the xcms data associated with a bin
    #we will subsequently subset this by category (treatments/mutants) and then replicates
    AllMIDSubsBeforeCategories = subset(MIDs_table,  Bin == myBin)

    #determine the formula of whatever is in the bin
    myFormula = unique(AllMIDSubsBeforeCategories$Formula)

    #get the carbons and nitrogens in the formula
    carbonsFormula = get_element_count( myFormula )[['C']]
    nitrogenFormula = get_element_count( myFormula )[['N']]

    #get all of the isotopologues
    isotopologueList = AllMIDSubsBeforeCategories$comp_result


    #get the M0 for the bin
    M0 = subset(AllMIDSubsBeforeCategories,  carbon == 0 & nitrogen == 0, select = -c(Bin,Formula, mz, rt, carbon, nitrogen, total_isotopes, comp_result,polarity))

    #get the total number of nitrogens
    nitrogens = AllMIDSubsBeforeCategories$nitrogen
    #get the total number of carbons
    carbons = AllMIDSubsBeforeCategories$carbon

    #now we've got everything that we need
    AllMIDSubsBeforeCategories = subset(MIDs_table,  Bin == myBin, select = -c(Bin,Formula, mz, rt, carbon, nitrogen, total_isotopes, comp_result,polarity) )

    #if there's no M0 keep track of the bin so that we can remove it.  If there is no M0, the bin is likely noise so keep track of the bin so we can remove it
    if(nrow(M0) == 0)
    {
      vecOfNoMIDs = c(vecOfNoMIDs, myBin)
    }

    #if we have mutiple MIDs, we're gonna start binning them
    if(nrow(M0) > 0)
    {

      #get the replicates

      #see if we have replicates for multiple categories
      reps = unique(data_cleanEnd(colnames(AllMIDSubsBeforeCategories)))

      #get these replicates
      reps = unique(reps)

      reps = reps[is.na(reps) == FALSE]
      #and get the categories as well, which may be treatment coniditions
      categories = unique(data_cleanMiddle(colnames(AllMIDSubsBeforeCategories)))
      categories = unique(categories)
      categories = categories[is.na(categories) == FALSE]

      #subset each row of MIDs by its replicate
      for(rep in 1:length(reps))
      {
        #in each category
        for(category in 1:length(categories))
        {


          justReps = colnames(AllMIDSubsBeforeCategories)[colnames(AllMIDSubsBeforeCategories) %like% paste0("_",reps[rep], sep = NULL)]

          #subset by the caetegorical variable (i.e. WT or MUT, treatment1 or treatment 2) and their associated reps
          repsAndCategories = justReps[justReps %like% categories[category]]


          whichRep = reps[rep]
          whichCategory = categories[category]
          AllMIDSubs = subset(AllMIDSubsBeforeCategories, select = repsAndCategories)

          #this will take us across each timepoint for the given
          #category and replicate so that we can track how the MIDs change across time

          for(j in 1:ncol(AllMIDSubs))
          {
            myToGetPercent = AllMIDSubs[,j]
            justToSum = colSums(AllMIDSubs)[j]

            justPercent =  (myToGetPercent / justToSum) * 100

            #effectively "refill" the table
            #and also calculate the avgs information

            #store all the avgs information

            #will store all of the information at a single timepoint for all the MIDs
            vectorOfAvgsInfo = vector()
            allMIDsNumber = sum(AllMIDSubsBeforeCategories$total_isotopes)

            #this will take us across all the isotoplogues for a timepoint
            for(k in 1:length(justPercent))
            {
              nameOfColumn = colnames(AllMIDSubs)[j]
              toReplace = justPercent[k]

              myIsoBin = isotopologueList[k]

              #fill in all the elements in the jth row
              #by each of their columns
              MIDs_table[MIDs_table$comp_result == myIsoBin, colnames(MIDs_table) == nameOfColumn] = toReplace


              #multiply by the number of labeled carbons and nitrogens of that isotopologue
              AvgsInfo = toReplace * (carbons[k] + nitrogens[k])

              #sum up all of the isotopologues to create the mass distribution vector
              vectorOfAvgsInfo = c(vectorOfAvgsInfo, AvgsInfo)

            }

            #in order to calculate average labeling, first multiply each isotologue by the
            #number of labeled carbons and nitrogens in that isotopologue, then sum all the isotopologues
            #and divide the sum by total number of carbons and nitrogens in the compound

            #calculate the total number of carbons and nitrogens in the compound
            CandNnumbers = carbonsFormula + nitrogenFormula

            #collapse all the isotologue info into one datapoint
            vectorOfAvgs = sum(vectorOfAvgsInfo) / CandNnumbers



            if(sum(CandNnumbers) == 0)
            {
              vectorOfAvgs = 0

            }


            #repopulate the average_labeling table with all of the updated data
            average_labeling[average_labeling$Bin == myBin, colnames(average_labeling) == nameOfColumn] = vectorOfAvgs

          }
        }
      }
    }
  }

  #hold the MIDs_table before we adjust by proxy pool size
  MIDs_tableBeforeScaling = MIDs_table


  #initialize the isotopologue column for the MIDs_tableBeforeScalin, which we'll update in the later loop
  MIDs_tableBeforeScaling$Isotopologue =  MIDs_tableBeforeScaling$total_isotopes


  for(MIDsBin in 1:length(unique(MIDs_table$Bin)))
  {
    myMIDsBinInfo = unique(MIDs_table$Bin)[MIDsBin]

    for(myMIDsColumn in 1:length(columnsToSearch))
    {
      MIDsColumnLookUp = columnsToSearch[myMIDsColumn]

      myMIDtoMultiply = MIDs_table[MIDs_table$Bin == myMIDsBinInfo, colnames(MIDs_table) == MIDsColumnLookUp]

      #now that we have the MID, let's multiply by the pool size
      poolSizeToMultiply = proxyPoolTable[proxyPoolTable$Bin == myMIDsBinInfo, colnames(proxyPoolTable) == MIDsColumnLookUp]
      #print(poolSizeToMultiply)

      #update the MIDs table with the MID scaled by pool size
      MIDs_table[MIDs_table$Bin == myMIDsBinInfo, colnames(MIDs_table) == MIDsColumnLookUp] = myMIDtoMultiply * poolSizeToMultiply

    }

  }

  #remove the bins the have no M0s from the MID table
  MIDs_table = MIDs_table[!MIDs_table$Bin %in% vecOfNoMIDs,]


  #remove the bins the have no M0s from the average labeling table
  average_labeling = average_labeling[!average_labeling$Bin %in% vecOfNoMIDs,]

  #initialize the isotopologue column
  MIDs_table$Isotopologue = MIDs_table$total_isotopes

  #for each MID, go through the MID table and fill in the number of labeled carbon and nitrogen
  for(i in 1:nrow(MIDs_table))
  {

    carbonNum = MIDs_table[i,]$carbon
    nitrogenNum = MIDs_table[i,]$nitrogen


    nameToAdd = paste0("C", carbonNum, "N", nitrogenNum)
    MIDs_table[i, colnames(MIDs_table) == "Isotopologue"] = nameToAdd



  }



  #for each MID, go through the unscaled MID table and fill in the number of labeled carbon and nitrogen
  for(i in 1:nrow(MIDs_tableBeforeScaling))
  {

    carbonNum = MIDs_tableBeforeScaling[i,]$carbon
    nitrogenNum = MIDs_tableBeforeScaling[i,]$nitrogen


    nameToAdd = paste0("C", carbonNum, "N", nitrogenNum)
    MIDs_tableBeforeScaling[i, colnames(MIDs_tableBeforeScaling) == "Isotopologue"] = nameToAdd



  }




  #Next we calculate label enrichment which is the nmole equivalent of labeled compound at a given timepoint
  #Label enrichment is the sum of labeled isotopologues of each bin within the pool size scaled MID table
  #note that this does not include M0 [unlabeled pool is removed]

  #declare the label enrichment
  names <- c("Bin",columnsToSearch)
  label_enrichment <- data.frame(matrix(ncol = length(names), nrow = length(unique(forMedianNormalization$Bin))))
  x <-  names
  colnames(label_enrichment) <- x
  label_enrichment$Bin = unique(forMedianNormalization$Bin)

  #remove the M0's from the MID tables
  MIDs_tableNoM0 = MIDs_table[MIDs_table$total_isotopes > 0,]


  #for each bin, sum up all of the MIDs
  for(column in 1:length(columnsToSearch))
  {
    columnSearch = columnsToSearch[column]

    for(noMIDBins in 1:length(label_enrichment$Bin))
    {
      myBinSearch = label_enrichment$Bin[noMIDBins]
      mySubsetToSum = MIDs_tableNoM0[MIDs_tableNoM0$Bin == myBinSearch, colnames(MIDs_tableNoM0) == columnSearch]

      MIDsSumToAdd = sum(mySubsetToSum)
      label_enrichment[label_enrichment$Bin == myBinSearch,colnames(label_enrichment) == columnSearch ] = MIDsSumToAdd
    }
  }

  average_labeling = average_labeling[!average_labeling$Bin %in% vecOfNoMIDs,]


  if (!is.null(output))
  {
    write.table(MIDs_tableBeforeScaling, file = paste(output, "MIDs.txt", sep = "_"), sep = "\t")
    write.table(MIDs_table, file = paste(output, "scaled_MIDs.txt", sep = "_"), sep = "\t")
    write.table(average_labeling, file = paste(output, "average_labeling.txt", sep = "_"), sep = "\t")
    write.table(label_enrichment, file = paste(output, "nanomol_equivalents_of_label.txt", sep = "_"), sep = "\t")
  }

  #return all of outputs as a get_table_objects() object
  listReturn = list(MIDs = MIDs_tableBeforeScaling, scaled_MIDs = MIDs_table,  average_labeling = average_labeling, mol_equivalent_labeling = label_enrichment)
}

