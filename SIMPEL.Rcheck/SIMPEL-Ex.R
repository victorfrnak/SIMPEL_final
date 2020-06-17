pkgname <- "SIMPEL"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "SIMPEL-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('SIMPEL')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("getClusterAndPlots")
### * getClusterAndPlots

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: getClusterAndPlots
### Title: Function to do statistical analyses and create plots from
###   non-stationary labeling experiments
### Aliases: getClusterAndPlots
### Keywords: MS1 dual isotopes, isotopic labeling, labels, metabolomics,
###   non-stationary stable untargeted

### ** Examples

getClustersAndPlots()



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("getClusterAndPlots", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("get_comp_mass")
### * get_comp_mass

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: get_comp_mass
### Title: This function is calculating the m/z for the chemical formulae
###   provided
### Aliases: get_comp_mass
### Keywords: MS1 dual isotopes, isotopic labeling, labels, metabolomics,
###   non-stationary stable untargeted

### ** Examples

get_comp_mass()



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("get_comp_mass", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("get_comp_mz_lookup")
### * get_comp_mz_lookup

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: get_comp_mz_lookup
### Title: This function is going to evaluate all of the XCMS_data to
###   identify isotopologues
### Aliases: get_comp_mz_lookup
### Keywords: MS1 dual isotopes, isotopic labeling, labels, metabolomics,
###   non-stationary stable untargeted

### ** Examples

get_comp_mz_lookup()



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("get_comp_mz_lookup", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("get_comp_stage")
### * get_comp_stage

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: get_comp_stage
### Title: This function is going to evaluate all of the XCMS_data to
###   identify isotopologues
### Aliases: get_comp_stage
### Keywords: MS1 dual isotopes, isotopic labeling, labels, metabolomics,
###   non-stationary stable untargeted

### ** Examples

get_comp_stage()



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("get_comp_stage", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("get_element_count")
### * get_element_count

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: get_element_count
### Title: This function is calculating the number of each of the elements
###   present in the formula
### Aliases: get_element_count
### Keywords: MS1 dual isotopes, isotopic labeling, labels, metabolomics,
###   non-stationary stable untargeted

### ** Examples

get_element_count()



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("get_element_count", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("get_table_objects")
### * get_table_objects

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: get_table_objects
### Title: Function to process XCMS data for isotopic enrichment and create
###   MID and average labeling objects from non-stationary isotopic
###   labeling experiments
### Aliases: get_table_objects
### Keywords: MS1 dual isotopes, isotopic labeling, labels, metabolomics,
###   non-stationary stable untargeted

### ** Examples

get_table_objects()



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("get_table_objects", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("label_enrichment_plot")
### * label_enrichment_plot

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: label_enrichment_plot
### Title: Function to plot time-course label enrichment both as
###   average_labeling and mole_equivalents_labeling in a non-stationary
###   isotopic labeling experiment XCMS data and start getting isotopic
###   information from the annotation file and add information columns to
###   be used as input for next steps
### Aliases: label_enrichment_plot
### Keywords: MS1 dual isotopes, isotopic labeling, labels, metabolomics,
###   non-stationary stable untargeted

### ** Examples

label_enrichment_plot()
makeAvgsplot()



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("label_enrichment_plot", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
