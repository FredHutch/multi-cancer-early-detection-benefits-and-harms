##################################################
# Read and merge SEER incidence and incidence-based
# mortality data by sex, race, and primary site
# diagnosed at ages 50-54, 60-64, or 70-74 in the
# calendar period 2000-2002 with 10- or 15-year
# follow-up for cancer death
##################################################
library(tidyverse)
library(here)

datestamp <- '2021-03-26'

##################################################
# Specify filenames for SEER data and target cancers
##################################################
fset <- tibble(incfile='seer_incidence_2000-2002_group=5_extended.csv',
               ibmfile=str_glue('seer_ibm_2000-2002_followup={c(10, 15)}_extended.csv'),
               cancersfile='cancers_2021-03-26.csv')

##################################################
# Read and format SEER data
##################################################
read_data <- function(filename, radix=1e5){
    dset <- read_csv(here('data', filename), col_types='ccccdddddd')
    dset <- dset %>% select(-'Standard Error')
    dset <- dset %>% rename(Sex='Sex (Male, Female)',
                            Race='Race recode (All, Black)',
                            Age='Age recode (50-54, 60-64, 70-74)',
                            Site='Site recode ICD-O-3/WHO 2008 (Liu et al groupings)',
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
        dset <- dset %>% mutate(Death.Rate=Death.Rate*Population,
                                Death.Lower=Death.Lower*Population,
                                Death.Upper=Death.Upper*Population)
        dset <- dset %>% select(-Population)
    }
    return(dset)
}

##################################################
# Read and format target cancers
##################################################
read_cancers <- function(filename)
    read_csv(here('data', filename), comment='#', col_types='c')

##################################################
# Control analysis
##################################################
control <- function(fset, saveit=FALSE){
    iset <- read_data(fset$incfile)
    mset <- read_data(fset$ibmfile)
    dset <- full_join(iset, mset, by=c('Sex', 'Race', 'Age', 'Site'))
    dset <- dset %>% mutate(Death.Rate=Death.Count/Population,
                            Death.Lower=Death.Lower/Population,
                            Death.Upper=Death.Upper/Population)
    dset <- dset %>% select(Sex,
                            Race,
                            Age,
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
    cset <- read_cancers(fset$cancersfile)
    dset <- left_join(cset, dset, by='Site')
    if(saveit){
        followup <- str_extract(fset$ibmfile, 'followup=(10|15)')
        outfile <- str_glue('seer_merged_2000-2002_{followup}_extended_{datestamp}.csv')
        write_csv(dset, here('data', outfile))
    }
    return(dset)
}
fset <- fset %>% group_by(incfile, ibmfile)
fset <- fset %>% do(control(., saveit=TRUE))

