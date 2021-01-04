#calculate a scaled MIDs table from the proxyPoolTable and a pre-existing
#MIDs table
prepareScaledMIDsTable = function(MIDs_table, proxyPoolTable, columnsToSearch, vecOfNoMIDs)
{
  #hold the MIDs_table before we adjust by proxy pool size
  #MIDs_tableBeforeScaling = MIDs_table

  MIDs_table = MIDs_table
  proxyPoolTable = proxyPoolTable
  columnsToSearch = columnsToSearch
  vecOfNoMIDs = vecOfNoMIDs
  #initialize the isotopologue column for the MIDs_tableBeforeScalin, which we'll update in the later loop
  #MIDs_tableBeforeScaling$Isotopologue =  MIDs_tableBeforeScaling$total_isotopes
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
      #print(MIDs_table[MIDs_table$Bin == myMIDsBinInfo, colnames(MIDs_table) == MIDsColumnLookUp])
      #print(myMIDtoMultiply * poolSizeToMultiply)
    }
  }

  #remove the bins the have no M0s from the MID table
  MIDs_table = MIDs_table[!MIDs_table$Bin %in% vecOfNoMIDs,]

  #remove the bins the have no M0s from the average labeling table
  #average_labeling = average_labeling[!average_labeling$Bin %in% vecOfNoMIDs,]

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

  #return the scaled MIDs_table
  return(MIDs_table)
}

