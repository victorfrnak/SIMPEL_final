#' Create natural abundance corrected post-processing data tables for pre-processed, stable isotope labeling based metabolomics datasets
#'
#' Function to process XCMS data for stable isotope enriched metabolomics experiments. This function first calculates all the
#' isotopologue m/z for a list of chemical formulae provided by the user as "compound_data" and identifies the matched m/z features within
#' the "xcms_data". In addition, natural abundance correction is achieved using the package IsoCorrectoR.
#' Four objects are created using this function, an MIDs table containing mass isotopologue distribution matrix for all the
#' compounds, a scaled_MIDs table where the MIDs are adjusted by a proxy for pool size to account for differences in Pool sizes of compounds,
#' an average_labeling table where the label enrichment for each compound is represented as a percentage and finally a mol_equivalent_labeling
#' table where the lable enrichment is represented as mol equivalents of labeled compound at each time point.  Four additional objects similar to the above
#' but corrected for natural abundance are also created.
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
#' @param XCMS_data XCMS pre-processed data provided as user input in .csv format
#' @param ppm numeric (1) define the maximal tolerated deviation in m/z (as parts per million) from calculate labeled isotopologue within XCMS_data
#' @param rt_tolerance  numeric (1) the maximal tolerate retention time deviation in minutes (from specified value provided in the annotation file) when searching for labeled isotopologues within XCMS_data.
#' @param compounds_data User annotation file (in .txt format) containing a list of compounds for which isotopologues are to be identified.
#' @param output A title to append to the files if the user wishes to print out the analyses tables/objects.
#' @keywords untargeted metabolomics, stable isotopes, non-stationary isotopic labeling, dual labels, MS1, natural abundance correction
#' @seealso get_table_objects(),PCA_and_heatmap(),label_enrichment_plot(),MIDplot(),getClustersAndPlots()
#' @return The function returns a class instance with eight objects 1. MIDs 2. scaled_MIDs 3. average_labeling 4. mol_equivalent_labeling as well as natural
#' abundance corrected versions of 1-4.
#' in addition, the function also exports the eight objects as .txt files containing the specific data tables for user verification into two
#' separate directories.
#' @author Shrikaar Kambhampati, Allen Hubbard
#' @export
#' @examples
#' test_13C15NGlutamine <- get_table_objectsII(XCMS_data, compounds_data, ppm = 5, rt_tolerance = 0.1, output = "13C15N_Glutamine")

get_table_objects_NA_corrected <- function(XCMS_data, compounds_data, ppm=10, rt_tolerance=.1, output=NULL){

  #we're going to have an output directory for SIMPEL
  #in the installed package location
  setWorkingDir = file.path(.libPaths(), "SIMPEL/Output_directory")

  #create the output directory for the user which will hold all their data
  dir.create(file.path(".", setWorkingDir), showWarnings = FALSE)


  #make sure that we've created the output directory as well
  dir.create(file.path(".", setWorkingDir), showWarnings = FALSE)

  #we are now going to be working in the directory we have created to hold all
  #of the outputs
  dir.create(file.path(setWorkingDir, output), showWarnings = FALSE)


  setwd(file.path(setWorkingDir, output))

  #make the directories for the not NA_adjusted table
  dir.create(file.path(".","get_table_objects"), showWarnings = FALSE)

  #make the directories for the NA corrected table
  dir.create(file.path(".", "NA_corrected"), showWarnings = FALSE)


  #set the variables from the input, i.e. ppm
  ppm = ppm
  rt_tolerance = rt_tolerance
  XCMS_data = XCMS_data
  compounds_data = compounds_data
  output = output

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
            #print("val is false")

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

  #move this up before removing only the things that we've matched
  #calculate the median normalized data
  forMedianNormalization = XCMS_data

  #create a subset table that only contains binned rows with
  #labeled data
  XCMS_data = XCMS_data[is.na(XCMS_data$comp_result) == FALSE,]


  #we're going to exclude the column labels whose rows do not include the numeric data
  vecToExclude = c("mz","polarity", "rt", "comp_result","Formula", "carbon", "nitrogen" , "total_isotopes", "Bin", "Compound" )

  #this vector of the medians is going to be critical because we are going to use
  #and reuse it
  columnsToSearch= setdiff(colnames(forMedianNormalization), vecToExclude)

  #also save the columns to add back
  subsetOfTableJustAnnotationData = forMedianNormalization[,colnames(forMedianNormalization) %in% vecToExclude]

  #get the median for each column
  #make sure to include all the data just for the median normalizaton
  for(column in 1:length(columnsToSearch))
  {
    myColumn = columnsToSearch[column]
    myMedian = median(forMedianNormalization[,colnames(forMedianNormalization) ==  myColumn])

    #divide the intensity of each entry by the median of the column
    forMedianNormalization[,colnames(forMedianNormalization) ==  myColumn] = forMedianNormalization[,colnames(forMedianNormalization) ==  myColumn] / myMedian

  }

  ##now ... remove all the compounds we don't have hits for
  forMedianNormalization = forMedianNormalization[is.na(forMedianNormalization$comp_result) == FALSE,]

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
      #print(sumTheIntensities)


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

  #run a master function to get the MIDs table, the averages table and also the list of
  #bins with incomplete sets of MIDs
  listOfAllOutputs = MIDsTableNoScale(MIDs_table, XCMS_data)

  #This is going to get us the 1st object of the list which will be the unscaled MIDs table
  MIDs_table = listOfAllOutputs[[1]]
  #This is going to get us the second object of the list, which will be the averages table
  average_labeling = listOfAllOutputs[[2]]
  #This is going to get us the third object of the list, which will be any of the MIDs without
  #sufficient species present
  vecOfNoMIDs = listOfAllOutputs[[3]]

  #hold the MIDs_table before we adjust by proxy pool size
  MIDs_tableBeforeScaling = MIDs_table

  #initialize the isotopologue column for the MIDs_tableBeforeScalin, which we'll update in the later loop
  MIDs_tableBeforeScaling$Isotopologue =  MIDs_tableBeforeScaling$total_isotopes


  #scale the MIDs table
  MIDs_tableScaled = prepareScaledMIDsTable(MIDs_tableBeforeScaling,proxyPoolTable, columnsToSearch,vecOfNoMIDs)


  #Next we calculate label enrichment which is the nmole equivalent of labeled compound at a given timepoint
  #Label enrichment is the sum of labeled isotopologues of each bin within the pool size scaled MID table
  #note that this does not include M0 [unlabeled pool is removed]


  #call the label enrichment fxn
  label_enrichment = getLabelEnrichment(MIDs_table,forMedianNormalization,columnsToSearch,subsetOfTableJustAnnotationData,compounds_data)

  if (!is.null(output))
  {
    write.table(MIDs_tableBeforeScaling, file = file.path("./get_table_objects", paste(output, "MIDs.txt", sep = "_")), sep = "\t")
    write.table(MIDs_table, file = file.path("./get_table_objects", paste(output, "scaled_MIDs.txt", sep = "_")), sep = "\t")
    write.table(average_labeling, file = file.path("./get_table_objects", paste(output, "average_labeling.txt", sep = "_")), sep = "\t")
    write.table(label_enrichment, file = file.path("./get_table_objects", paste(output, "mole_equivalents_of_label.txt", sep = "_")), sep = "\t")
  }



  #NA correct the MIDs now!
  MIDS_NACorrected = NACorrectionFxn(MIDs_table)


  #get average and MIDs table
  listOfAllOutputsNAcorrected = MIDsTableNoScale(MIDS_NACorrected, XCMS_data)



  MIDs_tableNAcorrected = listOfAllOutputsNAcorrected[[1]]
  average_labelingNAcorrected =listOfAllOutputsNAcorrected[[2]]
  vecOfNoMIDsNAcorrected = listOfAllOutputsNAcorrected[[3]]


  #scale the MIDs table by proxy pool
  scaledMIDsTableNAcorrected = prepareScaledMIDsTable(MIDS_NACorrected,proxyPoolTable, columnsToSearch,vecOfNoMIDs)

  #get the enrichment adjusted MIDs
  labelEnrichmentMIDsNAcorrected = getLabelEnrichment(MIDS_NACorrected,forMedianNormalization,columnsToSearch,subsetOfTableJustAnnotationData,compounds_data)


  if (!is.null(output))
  {
    write.table(MIDS_NACorrected, file =  file.path("NA_corrected", paste(output, "MIDsNAcorrected.txt", sep = "_")), sep = "\t")
    write.table(scaledMIDsTableNAcorrected, file = file.path("NA_corrected", paste(output, "scaled_MIDsNAcorrected.txt", sep = "_")), sep = "\t")
    write.table(average_labelingNAcorrected, file = file.path("NA_corrected", paste(output, "average_labelingNAcorrected.txt", sep = "_")), sep = "\t")
    write.table(labelEnrichmentMIDsNAcorrected, file = file.path("NA_corrected", paste(output, "mole_equivalents_of_labelNAcorrected.txt", sep = "_")), sep = "\t")
  }

  #return all of outputs as a get_table_objects() object
  listReturn = list(MIDs = MIDs_table, scaled_MIDs = MIDs_tableScaled,  average_labeling = average_labeling, mol_equivalent_labeling = label_enrichment,MIDS_Corrected = MIDS_NACorrected,scaledMIDsCorrected = scaledMIDsTableNAcorrected, averageLabeling_corrected = average_labelingNAcorrected, molEquivalent_corrected = labelEnrichmentMIDsNAcorrected)
}




##install devtools if you don't already have it
install.packages("devtools")##load the library devtools
library("devtools")##install SIMPEL using the install_github function within devtools
install_github("victorfrnak/SIMPEL_final/SIMPEL")



Compound_data<-read.table(system.file("extdata","Compound_data_DualLabel_SIMPEL.txt",package="SIMPEL"),sep ="\t",header =TRUE)
Compound_data = head(Compound_data)

xcms_data <-read.table(system.file("extdata","xcms_data_DualLabel_SIMPEL.csv",package="SIMPEL"),sep =",",header =TRUE)
