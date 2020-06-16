#'This function is going to evaluate all of the XCMS_data to identify isotopologues
#'
#' @import tidyverse
#' @param compound_data is the annotation file
#' @param comp_formula is the formula
#' @param polarity will be the polarity from both Compund_Data and XCMS_data
#' @param comp_lookup_table This will be your m/z table that was created based on the annotation file
#' @param r_time This will be the rt from your annotation file
#' @keywords untargeted metabolomics, stable isotopes, non-stationary isotopic labeling, dual labels, MS1
#' @export
#' @examples
#' get_comp_mz_lookup()

get_comp_mz_lookup <- function(compound_data, comp_formula, r_time, ppm,polarity){

   #### INNITIALIZATION ####
  w_c = 1.003355
  w_n = 0.997035

  d = list()
  gln_mass = get_comp_mass(comp_formula, polarity)
  c = get_element_count(comp_formula)[['C']]
  n = get_element_count(comp_formula)[['N']]
  comp_prefix <- compound_data %>%
    filter(formula == comp_formula & rt == r_time) %>%
    pull(prefix)
  ###two loops -
  ##one with compunds with just C and no N
  #one for compounds with C + N
  if(length(c) > 0)
  {
    #still calculate the mass of the compound
    #even if there are no Nitrogens and just look at all of the possible
    #carbon isotopes
    if(length(n) == 0)
    {
      #calculate all the possible carbon isotopes
      #and
      for(i in 0:c)
      {
        #number of Carbons and numbere of Nitrogens
        no_of_c = i
        no_of_n = 0
        #name the compound according to the 13C and 15N molecules
        compound = paste0(comp_prefix, '_', no_of_n, 'N', no_of_c, 'C')
        weight_compound = gln_mass + (i * w_c)
        #print(weight_compound)
        #print(" is weight_compound ")
        deviance <- ppm*0.000001*weight_compound
        #set upper and lower bounds on the masses
        #of the C and N isotopes
        upper_bound = weight_compound + deviance
        lower_bound = weight_compound - deviance
        #store, for each compound, the upper and lower bound masses
        #in a dictionary

        #if n is null, set to 0
        n = ifelse(is.null(n), 0, n)
        isotopeNumbers = n + c
        d[[compound]] = list('lb' = lower_bound,
                             'wc' = weight_compound,
                             'ub' = upper_bound,
                             'carbon' = no_of_c,
                             'nitrogen'= 0,
                             'isotope_numbers' =  isotopeNumbers)
      }
    }
    #if there are nitrogrens, no go ahead and calculate all masses
    #of all the possible C and N isotopes
    if(length(n) > 0)
    {
      #all possible combinations of 13C incorporation
      for(i in 0:c)
      {
        #possible combinations of 15N incorporation
        for(j in 0:n)
        {
          #number of Carbons and numbere of Nitrogens
          no_of_c = i
          no_of_n = j
          #name the compound according to the 13C and 15N molecules
          compound = paste0(comp_prefix, '_', no_of_n, 'N', no_of_c, 'C')
          weight_compound = gln_mass + (i * w_c) + (j * w_n)
          #print(weight_compound)
          #print("is the weight of compound")
          deviance <- ppm*0.000001*weight_compound
          #set upper and lower bounds on the masses
          #of the C and N isotopes
          upper_bound = weight_compound + deviance
          lower_bound = weight_compound - deviance
          #store, for each compound, the upper and lower bound masses
          #in a dictionary
          #if n is null, set to 0
          n = ifelse(is.null(n), 0, n)
          isotopeNumbers = n + c
          d[[compound]] = list('lb' = lower_bound,
                               'wc' = weight_compound,
                               'ub' = upper_bound,
                               'carbon' = no_of_c,
                               'nitrogen'= no_of_n,
                               'isotope_numbers' =  isotopeNumbers)

        }
      }
    }
  }
  return(d)
}
