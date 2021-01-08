##################################################
# Read and merge SEER incidence and incidence-based
# mortality data for women diagnosed with 1 of 6
# cancers at ages 50-54, 60-64, and 70-74 in the
# calendar period 2000-2002.
##################################################
library(tidyverse)
library(here)

datestamp <- '2021-01-07'

##################################################
# Specify filenames for SEER data
##################################################
incfile <- 'seer_incidence_2000-2002.csv'
ibmfile <- 'seer_ibm_2000-2002.csv'

##################################################
# Read and format SEER data
##################################################
read_data <- function(filename, radix=1e5){
    dset <- read_csv(here('data', filename), col_types='ccdddddd')
    dset <- dset %>% select(-'Standard Error')
    dset <- dset %>% rename(Age='Age recode (50-54, 60-64, 70-74)',
                            Site='Site recode (lung, breast, colorectum, ovary, pancreas, liver)',
                            Rate='Crude Rate',
                            Lower='Lower Confidence Interval',
                            Upper='Upper Confidence Interval')
    dset <- dset %>% mutate(Age=sub(' years$', '', Age),
                            Rate=Rate/radix,
                            Lower=Lower/radix,
                            Upper=Upper/radix)
    if(grepl('incidence', filename)){
        dset <- dset %>% rename(Diagnosis.Count=Count,
                                Diagnosis.Rate=Rate,
                                Diagnosis.Lower=Lower,
                                Diagnosis.Upper=Upper)
    } else {
        dset <- dset %>% rename(Death.Count=Count,
                                Death.Rate=Rate,
                                Death.Lower=Lower,
                                Death.Upper=Upper)
        dset <- dset %>% select(-Population)
    }
    return(dset)
}

##################################################
# Control analysis
##################################################
control <- function(incfile, ibmfile, saveit=FALSE){
    iset <- read_data(incfile)
    mset <- read_data(ibmfile)
    dset <- full_join(iset, mset, by=c('Age', 'Site'))
    dset <- dset %>% select(Age,
                            Site,
                            Population,
                            Diagnosis.Count,
                            Diagnosis.Rate,
                            Diagnosis.Lower,
                            Diagnosis.Upper,
                            Death.Count,
                            Death.Rate,
                            Death.Lower,
                            Death.Upper)
    if(saveit){
        outfile <- paste('seer', 'merged', datestamp, sep='_')
        outfile <- paste(outfile, 'csv', sep='.')
        write_csv(dset, here('data', outfile))
    }
}
control(incfile, ibmfile, saveit=TRUE)

