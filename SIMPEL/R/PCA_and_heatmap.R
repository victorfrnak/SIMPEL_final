#' Function to perform global analysis of label enrichment and generate PCA plots and heatmaps
#'
#' This function performs Principle Component Analysis and plots the average_labeling or mol_equivalent_labeling data in a two dimensional space
#' to allow comparison of time couse labeling enrichment data obtained from a non-stationary isotopic labeling experiment. In addition, this function
#' also plots the label enrichment for each of the compounds as a heat map to easily identify compounds that have significant labeling using global view
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
#' @param mydata1 object (1) The 'average_labeling' or the 'mol_equivalent_labeling' object created using get_table_objects() function
#' @param heatMapCategories string (1) or (2) c("Category1") or c("Category1", "Category2") The category to be used for plotting a heatmap
#' @param PCMax numeric (1) maximum numbers of PC's to plot
#' @param labels string (1) label to be used for heat maps i.e. "Bin" or "Compound" column
#' @param outputName string (1) The name to be appended to the the output pdf
#' @keywords untargeted metabolomics, stable isotopes, non-stationary isotopic labeling, dual labels, MS1
#' @seealso label_enrichment_plot(), MIDplot(), getClustersAndPlots(), get_table_objects()
#' @return The function returns a pdf file with PCA plots comparing all PCs upto the user specified number and heat maps with label enrichment
#' for all compounds using the user specified Categories
#' @author Shrikaar Kambhampati, Allen Hubbard
#' @export
#' @examples
#' PCA_Heatmap_Avg <- PCA_and_heatmap(test_13C15NGlutamine$average_labeling, PCMax = 3, heatMapCategories = c("WT", "GAT"), labels="Bin", outputName = "Average")
#' PCA_Heatmap_molEqui <- PCA_and_heatmap(test_13C15NGlutamine$mol_equivalent_labeling, PCMax = 3, heatMapCategories = c("WT", "GAT"), labels="Bin", outputName = "mol_equivalents")

PCA_and_heatmap <- function(mydata1, PCMax=3, heatMapCategories, labels="Bin", outputName = "_")
{

  #print("june 19th")
  #function to extract the middle field of the column labels
  #which corresponds to category
  data_cleanII <- function(x) sapply (strsplit(x , '[_]' ), `[` , 2)


  PCMax = PCMax
  allTotFigures = list()
  heatMapCategories = heatMapCategories
  mydata1 = mydata1
  outputName = outputName

  #subset the table by just the columns matching our category of interest
  #catSubset = mydata1[,colnames(mydata1) %like% Category]


  #remove bin column from table to do PCA on
  rownames(mydata1) =  mydata1[,colnames(mydata1) %in% labels]
  mydata1Backup = mydata1
  mydata1[labels] = NULL



  #make sure that there are now all zeroes columns or na containing columns
  #we're going to exclude the column labels whose rows do not include the numeric data
  vecToExclude = c("mz","polarity", "rt", "comp_result","Formula", "carbon", "nitrogen" , "total_isotopes", "Bin", "Compound", "CompoundPlaceholder", "Isotopologue" )

  #we'll only want to use the columns with numeric values
  columnsToUse = setdiff(colnames(mydata1), vecToExclude)

  #exclude the other columns
  mydata1 = mydata1[,colnames(mydata1) %in% columnsToUse]


  #exclude all of the nun-numeric columns
  mydata1 = mydata1[complete.cases(mydata1), ]
  mydata1 = mydata1[rowSums(mydata1) != 0 ,]


  #make a metadata table (DataFrameLabel) in order to do the PCA
  DataFrameLabel = colnames(mydata1)
  DataFrameLabel = as.data.frame(DataFrameLabel)
  colnames(DataFrameLabel) = c("Samples")

  #make sure unintended characters are removed with this function
  #and that we're including a column with the category + time as well as a separate
  #one with just the time
  DataFrameLabel$CategoryTime = paste(data_cleanII(as.character(DataFrameLabel$Samples)),data_clean(as.character(DataFrameLabel$Samples)), sep = "_")
  DataFrameLabel$Time = data_clean(as.character(DataFrameLabel$Samples))

  #just make sure that we're always having category come first
  #originally this was Category
  DataFrameLabel$Category = gsub("*.*_(.*)_.*", "\\1", DataFrameLabel$Samples)

  #get PCA object to determine maximum number of PC's
  prcompObject = prcomp(t(mydata1))

  if(PCMax > 0)
  {
    #now make the PCA autoplots
    for(x in 1:PCMax)
    {
      for(y in 1:PCMax)
      {
        if(x > y)
        {
          if( x <= ncol(prcompObject$x) & y <=  ncol(prcompObject$x))
          {
          #second plot, coloring by time and then Category as well
            autoplotList = list()
            autoplotList[[1]] = print(autoplot(prcomp(t(mydata1)), data =  DataFrameLabel, colour = 'CategoryTime', shape = 'Category', frame.colour = 'CategoryTime', size = 3, x = x, y = y))
            allTotFigures = prepend(allTotFigures,autoplotList)
          }
        }
      }
    }
  }

  #the heatmaps may be for one or both conditions


  #list to hold heatmap to prepend to the main list
  heatMapListAll = list()

  heatMapList = list()
  toPrintFile = paste(paste0(heatMapCategories),outputName,"heatmap_and_PCA.pdf")
  pdf(toPrintFile)
  for(i in 1:length(heatMapCategories))
  {
    #print(colnames(mydata1))
    #print("is colnames(mydata1)")
    #print(heatMapCategories[i])
    #print("is heatMapCategories[i]")
    catSubset = mydata1[,colnames(mydata1) %like% heatMapCategories[i]]
    heatmapName = paste("Heatmap for", heatMapCategories[i], sep = " ")
    print(heatmap(as.matrix(catSubset), Colv = NA, scale = c("none"), col = cm.colors(300), main = heatmapName))
  }

  for(i in  allTotFigures)
  {
    print(i)
  }

  dev.off()

  return(allTotFigures)
}

