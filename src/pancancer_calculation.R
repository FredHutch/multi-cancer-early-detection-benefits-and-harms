##################################################
# Project expected outcomes of pan-cancer testing
# at age-specific screening occasions:
# 1. expected unnecessary confirmation tests (UCTs)
# 2. expected cancers detected (CDs)
# 3. expected lives saved (LSs)
##################################################
library(tidyverse)
library(here)
library(scales)
library(viridis)

#datestamp <- '2021-03-03'
#datestamp <- '2021-03-08'
#datestamp <- '2021-03-12'
datestamp <- '2021-03-13'

##################################################
# Project outcomes for k-cancer test
#
# Notation in paper         Variable name
# ---------------------------------------
# P_A                       prevalence (of A)
# P_A(T+)                   sensitivity (of A)
# L_A(A│T+)                 localization (of A)
# (unspecified)             mortality (of A)
# Sp                        specificity
##################################################
k_cancer_test <- function(dset, specificity=0.99, size=1000){
    # group results by cancer site
    gset <- dset %>% group_by(site)
    # expected unnecessary confirmation tests
    tset <- gset %>% summarize(prevs=value[feature == 'prevalence'],
                               tests=value[feature == 'prevalence']*
                                     value[feature == 'sensitivity'] *
                                     (1-value[feature == 'localization']))
    tests <- with(tset, size*(sum(tests)+(1-specificity)*(1-sum(prevs))))
    # expected cancers detected
    iset <- gset %>% summarize(cancers=value[feature == 'prevalence']*
                                       value[feature == 'sensitivity']*
                                       value[feature == 'localization'])
    cancers <- with(iset, size*sum(cancers))
    # expected lives saved
    lset <- gset %>% summarize(lives=value[feature == 'effect']*
                                     value[feature == 'mortality'])
    lives <- with(lset, size*sum(lives))
    tibble(UCT=tests, CD=cancers, LS=lives)
}

##################################################
# Analysis 1: pan-cancer test for two hypothetical
# cancers across selected scenarios:
# 1. varying sensitivity of both cancers=0.5 to 0.9
# 2. fixed prevalence of cancer A=0.001
# 3. varying prevalence of cancer B=0 to 0.01
# 4. varying localization of both cancers=0.5 to 0.8
# 5. varying overall specificity=0.95 to 0.99
##################################################
hypothetical_test <- function(specificity){
    dset <- expand_grid(prevalence.a=0.001,
                        prevalence.b=c(0, 0.0005, 0.001, 0.005, 0.01),
                        sensitivity.a=seq(0.5, 0.9, by=0.1),
                        localization.a=c(0.5, 0.8),
                        mortality.a=0)
    dset <- dset %>% mutate(scenario=seq(nrow(dset)),
                            sensitivity=sensitivity.a,
                            localization=localization.a,
                            prevalence=prevalence.b,
                            sensitivity.b=sensitivity.a,
                            localization.b=localization.a,
                            mortality.b=mortality.a)
    dset <- dset %>% pivot_longer(-c(scenario,
                                     sensitivity,
                                     localization,
                                     prevalence),
                                  names_to=c('feature', 'site'),
                                  names_pattern='(.*)[.]([ab])',
                                  values_to='value')
    dset <- dset %>% mutate(specificity=specificity)
    dset <- dset %>% group_by(scenario,
                              sensitivity,
                              localization,
                              specificity,
                              prevalence)
    dset <- dset %>% do(k_cancer_test(., specificity=specificity))
    dset <- dset %>% ungroup()
}

##################################################
# Read merged SEER incidence and IBM data
##################################################
read_data <- function(filename){
    dset <- read_csv(here('data', filename), col_types='ccddddddddd')
    dset <- dset %>% rename(age='Age', site='Site')
    dset <- dset %>% select(-Population,
                            -Diagnosis.Count,
                            -Death.Count,
                            -ends_with('Lower'),
                            -ends_with('Upper'))
    dset <- dset %>% rename(prevalence='Diagnosis.Rate',
                            mortality='Death.Rate')
    dset <- dset %>% mutate(site=factor(site,
                                        levels=c('Breast',
                                                 'Colon and Rectum',
                                                 'Liver and Intrahepatic Bile Duct',
                                                 'Lung and Bronchus',
                                                 'Ovary',
                                                 'Pancreas'),
                                        labels=c('Breast',
                                                 'Colon',
                                                 'Liver',
                                                 'Lung',
                                                 'Ovary',
                                                 'Pancreas')))
}

##################################################
# Age-specific outcomes from SEER data
##################################################
age_analysis <- function(bset, aset){
    dset <- bind_rows(bset, aset)
    dset <- dset %>% mutate(candidate=unique(bset$site))
    dset <- dset %>% pivot_longer(-c(age, site, candidate),
                                  names_to='feature',
                                  values_to='value')
    dset <- dset %>% group_by(age)
    dset <- dset %>% do(k_cancer_test(.))
    dset <- dset %>% ungroup()
}

##################################################
# Age-specific incremental impact of candidate cancer
##################################################
age_analysis_incremental <- function(dset, existing){
    aset <- dset %>% filter(site %in% existing)
    bset <- dset %>% filter(!site %in% existing)
    bset <- bset %>% group_by(site)
    rset <- bset %>% do(age_analysis(., aset))
    rset <- rset %>% ungroup()
}

##################################################
# Control plot aesthetics
##################################################
gg_theme <- function(...){
    theme_set(theme_classic())
    theme_update(axis.text=element_text(size=10, colour='black'),
                 axis.line=element_line(colour='black'),
                 axis.title=element_text(size=14),
                 axis.ticks.length=unit(0.2, 'cm'),
                 strip.background=element_rect(colour=NA, fill=NA),
                 strip.text=element_text(size=14))
    theme_update(...)
}

##################################################
# Visualize UCTs in hypothetical analysis
##################################################
hypothetical_uct_plot <- function(ext='png', saveit=FALSE){
    dset <- bind_rows(hypothetical_test(specificity=0.97),
                      hypothetical_test(specificity=0.98),
                      hypothetical_test(specificity=0.99),
                      hypothetical_test(specificity=1.0))
    dset <- dset %>% filter(prevalence == 0.01,
                            localization == 0.8)
    dset <- dset %>% select(-prevalence, -localization, -CD, -LS)
    dset <- dset %>% mutate(sensitivity=factor(sensitivity,
                                               levels=unique(sensitivity),
                                               labels=sprintf('%2.0f%%', 100*unique(sensitivity))))
    gg_theme()
    gg <- ggplot(data=dset)
    gg <- gg+geom_line(aes(x=specificity,
                           y=UCT,
                           alpha=sensitivity),
                       size=0.6)
    gg <- gg+scale_x_continuous('\nSpecificity for cancers A and B',
                                labels=percent_format(accuracy=1),
                                breaks=seq(0.97, 1, by=0.01),
                                limits=c(0.969, 1.001))
    gg <- gg+scale_y_continuous('Unnecessary confirmation tests per 1,000 women\n',
                                limits=c(0, 33),
                                expand=c(0, 0))
    gg <- gg+guides(alpha=guide_legend(title='Sensitivity for\ncancers A and B',
                                       keywidth=unit(1, 'cm'),
                                       title.theme=element_text(size=12),
                                       label.theme=element_text(size=12, angle=0)))
    print(gg)
    if(saveit){
        filename <- paste('hypothetical_uct', datestamp, sep='_')
        filename <- paste(filename, ext, sep='.')
        ggsave(here('plots', filename),
               plot=gg,
               width=6,
               height=6)
    }
}

##################################################
# Visualize CDs in hypothetical analysis
##################################################
hypothetical_cd_plot <- function(ext='png', saveit=FALSE){
    dset <- hypothetical_test(specificity=0.99)
    dset <- dset %>% ungroup()
    dset <- dset %>% filter(localization == 0.8)
    dset <- dset %>% select(-specificity, -localization, -UCT, -LS)
    dset <- dset %>% mutate(sensitivity=factor(sensitivity,
                                               levels=unique(sensitivity),
                                               labels=sprintf('%2.0f%%', 100*unique(sensitivity))))
    gg_theme()
    gg <- ggplot(data=dset)
    gg <- gg+geom_line(aes(x=prevalence,
                           y=CD,
                           alpha=sensitivity),
                       size=0.6)
    gg <- gg+scale_x_continuous('\nPrevalence of cancer B',
                                labels=percent_format(accuracy=0.1),
                                breaks=seq(0, 0.01, by=0.002))
    gg <- gg+scale_y_continuous('Cancers detected per 1,000 women\n',
                                limits=c(0, 9),
                                breaks=seq(0, 8, by=2),
                                expand=c(0, 0))
    gg <- gg+guides(alpha=guide_legend(title='Sensitivity for\ncancers A and B',
                                       keywidth=unit(1, 'cm'),
                                       title.theme=element_text(size=12),
                                       label.theme=element_text(size=12, angle=0)))
    print(gg)
    if(saveit){
        filename <- paste('hypothetical_cd', datestamp, sep='_')
        filename <- paste(filename, ext, sep='.')
        ggsave(here('plots', filename),
               plot=gg,
               width=6,
               height=6)
    }
}

##################################################
# Visualize outcomes in empirical analysis
##################################################
empirical_age_plot <- function(dset, figureno, ext='png', sensitivity=FALSE, saveit=FALSE){
    dset <- dset %>% mutate(UCT.CD=UCT/CD, UCT.LS=UCT/LS)
    dset <- dset %>% select(-UCT, -CD, -LS)
    dset <- dset %>% pivot_longer(cols=c(UCT.CD, UCT.LS),
                                  names_to='outcome',
                                  values_to='value')
    dset <- dset %>% ungroup()
    dset <- dset %>% arrange(value)
    dset <- dset %>% mutate(age=sub('-[567]4', ' y', age),
                            outcome=factor(outcome,
                                           levels=c('UCT.CD', 'UCT.LS'),
                                           labels=c('Unnecessary\nconfirmation\ntests per cancer\ndetected',
                                                    'Unnecessary\nconfirmation\ntests per life\nsaved')),
                            site=factor(site, levels=unique(site)))
    if(sensitivity){
        dset <- dset %>% filter(outcome == 'Unnecessary\nconfirmation\ntests per life\nsaved')
        height <- 4
    } else {
        height <- 6
    }
    ymax <- switch(as.character(figureno), '2'=80, '3'=50, 'S1'=150, 'S2'=80)
    gg_theme(axis.text.x=element_text(size=10, angle=90, vjust=0.5, hjust=1),
             axis.ticks.x=element_blank(),
             panel.spacing=unit(0.02, 'npc'),
             strip.text.y=element_text(size=10, angle=0))
    gg <- ggplot(data=dset)
    gg <- gg+geom_bar(aes(x=site, y=value),
                      stat='identity',
                      position='dodge')
    gg <- gg+geom_hline(aes(yintercept=0), colour='black')
    if(!sensitivity)
        gg <- gg+geom_blank(data=dset %>% filter(outcome == 'Unnecessary\nconfirmation\ntests per cancer\ndetected'), aes(y=8))
    gg <- gg+geom_blank(data=dset %>% filter(outcome == 'Unnecessary\nconfirmation\ntests per life\nsaved'), aes(y=ymax))
    gg <- gg+facet_grid(outcome~age, scales='free_y')
    gg <- gg+scale_x_discrete(name='')
    gg <- gg+scale_y_continuous(name='',
                                expand=c(0, 0))
    print(gg)
    if(saveit){
        filename <- str_glue('figure{figureno}_{datestamp}')
        filename <- paste(filename, ext, sep='.')
        ggsave(here('plots', filename),
               plot=gg,
               width=8,
               height=height)
    }
}

##################################################
# Format table for hypothetical two-cancer test
##################################################
format_hypothetical <- function(dset, saveit=FALSE){
    dset <- dset %>% select(-scenario, -LS)
    dset <- dset %>% filter(prevalence > 0, sensitivity %in% c(0.5, 0.8))
    dset <- dset %>% mutate(sensitivity=sprintf('%2.0f', 100*sensitivity),
                            specificity=sprintf('%2.0f', 100*specificity),
                            localization=sprintf('%2.0f', 100*localization),
                            prevalence=sprintf('%5.2f', 100*prevalence),
                            UCT=sprintf('%4.1f', UCT),
                            CD=sprintf('%3.1f', CD))
    dset <- dset %>% rename('Sensitivity, %'='sensitivity',
                            'Localization, %'='localization',
                            'Specificity, %'='specificity',
                            'Prevalence, %'='prevalence',
                            'Unnecessary confirmation tests, n'='UCT',
                            'Cancers detected, n'='CD')
    if(saveit){
        filename <- str_glue('supplemental_table1_{datestamp}.csv')
        write_csv(dset, here('data', filename))
    }
}

##################################################
# Format table for realistic six-cancer test
##################################################
format_empirical <- function(dset, saveit=FALSE){
    dset <- dset %>% select(-site)
    dset <- dset %>% mutate(age=sub('-[567]4', '', age),
                            UCT.CD=UCT/CD,
                            UCT.LS=UCT/LS)
    dset <- dset %>% mutate(UCT=sprintf('%4.1f', UCT),
                            CD=sprintf('%3.1f', CD),
                            LS=sprintf('%3.1f', LS),
                            UCT.CD=sprintf('%3.1f', UCT.CD),
                            UCT.LS=sprintf('%3.1f', UCT.LS))
    dset <- dset %>% arrange(age)
    dset <- dset %>% rename('Screening age, y'='age',
                            'Unnecessary confirmation tests, n'='UCT',
                            'Cancers detected, n'='CD',
                            'Lives saved, n'='LS',
                            'UCT/CD'='UCT.CD',
                            'UCT/LS'='UCT.LS')
    if(saveit){
        filename <- str_glue('table3_{datestamp}.csv')
        write_csv(dset, here('data', filename))
    }
}

##################################################
# Format table for sensitivity analyses
##################################################
format_supplemental <- function(dset, tableno, saveit=FALSE){
    dset <- dset %>% select(age, site, UCT, CD, LS)
    dset <- dset %>% mutate(age=sub('-[567]4', '', age))
    dset <- dset %>% mutate(UCT=sprintf('%4.1f', UCT),
                            CD=sprintf('%3.1f', CD),
                            LS=sprintf('%3.1f', LS))
    dset <- dset %>% arrange(age, site)
    dset <- dset %>% rename('Tissue of origin'='site',
                            'Screening age, y'='age',
                            'Unnecessary confirmation tests, n'='UCT',
                            'Cancers detected, n'='CD',
                            'Lives saved, n'='LS')
    if(saveit){
        filename <- str_glue('supplemental_table{tableno}_{datestamp}.csv')
        write_csv(dset, here('data', filename))
    }
}

##################################################
# Table 1
##################################################
#sset <- read_data('seer_merged_2000-2002_followup=15_2021-03-03.csv')

##################################################
# Table 2
##################################################
#pset <- tribble(~site, ~sensitivity, ~localization,
#                'Breast',   0.64, 0.96,
#                'Colon',    0.74, 0.97,
#                'Lung',     0.59, 0.92,
#                'Ovary',    0.67, 0.96,
#                'Pancreas', 0.78, 0.79,
#                'Liver',    0.68, 0.72)
#pset %>% mutate(marginal=sensitivity*localization)

##################################################
# Table 3
##################################################
#pset <- pset %>% mutate(effect=0.2)
#sset <- full_join(sset, pset, by='site')
#iset6 <- age_analysis_incremental(sset, setdiff(pset$site, 'Breast'))
#format_empirical(iset6, saveit=TRUE)

##################################################
# Figure 1
##################################################
#hypothetical_uct_plot(saveit=TRUE)
#hypothetical_cd_plot(saveit=TRUE)

##################################################
# Figure 2
##################################################
#iset1 <- age_analysis_incremental(sset, 'Breast')
#empirical_age_plot(iset1, figureno=2, sensitivity=FALSE, saveit=TRUE)

##################################################
# Figure 3
##################################################
#iset2 <- age_analysis_incremental(sset, c('Breast', 'Lung'))
#empirical_age_plot(iset2, figureno=3, sensitivity=FALSE, saveit=TRUE)

##################################################
# Supplemental Figure 1
##################################################
#sset1s <- sset %>% mutate(effect=ifelse(site %in% c('Breast', 'Colon', 'Lung'), 0.1, 0.5))
#iset1s <- age_analysis_incremental(sset1s, 'Breast')
#empirical_age_plot(iset1s, figureno='S1', sensitivity=TRUE, saveit=TRUE)

##################################################
# Supplemental Figure 2
##################################################
#sset10 <- read_data('seer_merged_2000-2002_followup=10_2021-03-03.csv')
#iset10 <- age_analysis_incremental(sset, 'Breast')
#empirical_age_plot(iset10, figureno='S2', sensitivity=TRUE, saveit=TRUE)

##################################################
# Supplemental Table 1
##################################################
#hset <- bind_rows(hypothetical_test(specificity=0.95),
#                  hypothetical_test(specificity=0.99))
#format_hypothetical(hset, saveit=TRUE)

##################################################
# Supplemental Table 2
##################################################
#format_supplemental(iset1, tableno=2, saveit=TRUE)

##################################################
# Supplemental Table 3
##################################################
#format_supplemental(iset2, tableno=3, saveit=TRUE)

