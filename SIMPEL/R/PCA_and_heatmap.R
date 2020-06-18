#' Function to plot Mass Isotopologue Distribution, i.e. MIDs, in a non stationary isotopic
#' labeling experiment
#'
#' @param mydata you can use average_labeling or mole_equivalents_labeling as an input for this function
#' @param heatMapCategories will be the tissue type(s), i.e. WT/Mut, or treatment, i.e treated vs untreated or both.  A heatmap will be conducted for the data associated with each category
#' @param  PCMax maximum numbers of PC's to plot
#' @param  labels label datapoints by either the "Bin" or "Compound" column
#' @keywords untargeted metabolomics, stable isotopes, non-stationary isotopic labeling, dual labels, MS1
#' @export
#' @examples
#' PCA_and_heatmap()

#testPCAHeatmap = PCA_and_heatmap(tableObjectToUse$mol_equivalent_labeling,"GAT", 5)
#testPCAHeatmap = PCA_and_heatmap(tableObjectToUse$mol_equivalent_labeling,"GAT", 3)
#testPCAHeatmap = PCA_and_heatmap(testWholePackage$mol_equivalent_labeling,"GAT", 3)

PCA_and_heatmap <- function(mydata1, PCMax=3, heatMapCategories, labels="Bin")
{

  #function to extract the middle field of the column labels
  #which corresponds to category
  data_cleanII <- function(x) sapply (strsplit(x , '[_]' ), `[` , 2)


  PCMax = PCMax
  allTotFigures = list()
  heatMapCategories = heatMapCategories
  mydata1 = mydata1

  #subset the table by just the columns matching our category of interest
  #catSubset = mydata1[,colnames(mydata1) %like% Category]


  #remove bin column from table to do PCA on
  rownames(mydata1) =  mydata1[,colnames(mydata1) %in% labels]
  mydata1Backup = mydata1
  mydata1[labels] = NULL

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

  if(PCMax > 0)
  {
    #now make the PCA autoplots
    for(x in 1:PCMax)
    {
      for(y in 1:PCMax)
      {
        if(x > y)
        {
          #second plot, coloring by time and then Category as well
          autoplotList = list()
          autoplotList[[1]] = print(autoplot(prcomp(t(mydata1)), data =  DataFrameLabel, colour = 'CategoryTime', shape = 'Category', frame.colour = 'CategoryTime', size = 3, x = x, y = y))
          allTotFigures = prepend(allTotFigures,autoplotList)
        }
      }
    }
  }


  #the heatmaps may be for one or both conditions


  #list to hold heatmap to prepend to the main list
  heatMapListAll = list()

  heatMapList = list()
  toPrintFile = paste(Category,"heatmap_and_PCA.pdf")
  pdf(toPrintFile)
  for(i in 1:length(heatMapCategories))
  {
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

