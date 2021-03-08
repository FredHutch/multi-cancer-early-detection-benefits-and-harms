##################################################
# Project expected outcomes of pan-cancer testing
# at age-specific screening occasions:
# 1. expected unnecessary confirmation tests (UCTs)
# 2. expected cancers detected (CDs)
# 3. expected cancer deaths prevented (CDPs)
##################################################
library(tidyverse)
library(here)
library(scales)
library(viridis)

#datestamp <- '2021-03-03'
datestamp <- '2021-03-08'

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
k_cancer_test <- function(dset, specificity=0.99, effect=0.2, size=1000){
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
    # expected cancer deaths prevented
    saved <- with(dset, size*effect*sum(value[feature == 'mortality']))
    return(tibble(UCT=tests, CD=cancers, CDP=saved))
}

##################################################
# Analysis 1: pan-cancer test for two hypothetical
# cancers across selected scenarios:
# 1. fixed sensitivity of both cancers=0.7
# 2. fixed prevalence of cancer A=0.001
# 3. varying prevalence of cancer B=0.0005 to 0.01
# 4. varying localization of both cancers=0.5 to 0.8
# 5. varying overall specificity=0.95 to 0.99
##################################################
analysis1 <- function(specificity){
    dset <- expand.grid(prevalence.a=0.001,
                        prevalence.b=c(0, 0.0005, 0.001, 0.005, 0.01),
                        #sensitivity.a=c(0.5, 0.7),
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
}
hset <- bind_rows(analysis1(specificity=0.95),
                  analysis1(specificity=0.99))

##################################################
# Pan-cancer test performance characteristics
##################################################
pset <- tribble(~site, ~sensitivity, ~localization,
                'Breast', 0.64, 0.96,
                'Colon and Rectum', 0.74, 0.97,
                'Lung and Bronchus', 0.59, 0.92,
                'Ovary', 0.67, 0.96,
                'Pancreas', 0.78, 0.79,
                'Liver and Intrahepatic Bile Duct', 0.68, 0.72)
#pset <- pset %>% mutate(Marginal=Sensitivity*Localization)

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
    dset
}
#sset <- read_data(str_glue('seer_merged_2000-2002_followup=10_2021-03-03.csv'))
sset <- read_data(str_glue('seer_merged_2000-2002_followup=15_2021-03-03.csv'))
sset <- full_join(sset, pset, by='site')

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
}

##################################################
# Age-specific incremental impact of candidate cancer
##################################################
age_analysis_incremental <- function(dset, existing){
    aset <- dset %>% filter(site %in% existing)
    bset <- dset %>% filter(!site %in% existing)
    bset <- bset %>% group_by(site)
    rset <- bset %>% do(age_analysis(., aset))
}

##################################################
# Selected incremental analyses
##################################################
iset1 <- age_analysis_incremental(sset, 'Breast')
iset2 <- age_analysis_incremental(sset, c('Breast', 'Lung and Bronchus'))
iset6 <- age_analysis_incremental(sset, setdiff(pset$site, 'Breast'))

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
# Visualize patterns in hypothetical analysis
##################################################
hypothetical_plot <- function(ext='png', saveit=FALSE){
    dset <- bind_rows(analysis1(specificity=0.97),
                      analysis1(specificity=0.98),
                      analysis1(specificity=0.99))
    dset <- dset %>% filter(prevalence != 0.0005,
                            localization != 0.5)
    dset <- dset %>% select(-CDP)
    dset <- dset %>% pivot_longer(cols=c(UCT, CD),
                                  names_to='outcome',
                                  values_to='value')
    dset <- dset %>% mutate(specificity=factor(specificity),
                            localization=factor(localization),
                            outcome=factor(outcome,
                                           levels=c('UCT', 'CD'),
                                           labels=c('Unnecessary confirmation tests',
                                                    'Expected cancers detected')))
    gg_theme(legend.position='bottom')
    gg <- ggplot(data=dset)
    gg <- gg+geom_line(aes(x=prevalence,
                           y=value,
                           colour=specificity,
                           alpha=sensitivity,
                           group=interaction(sensitivity,
                                             specificity)),
                       show.legend=c(colour=TRUE, alpha=FALSE),
                       size=0.6)
    gg <- gg+geom_blank(data=dset %>% filter(outcome == 'Unnecessary confirmation tests'),
                        aes(y=40))
    gg <- gg+geom_blank(data=dset %>% filter(outcome == 'Expected cancers detected'),
                        aes(y=8))
    gg <- gg+geom_blank(aes(y=0))
    gg <- gg+geom_hline(aes(yintercept=0), colour='black')
    gg <- gg+facet_wrap(~outcome, ncol=1, scales='free_y')
    gg <- gg+scale_x_continuous('\nPrevalence of cancer B',
                                labels=percent_format(accuracy=0.1),
                                breaks=seq(0, 0.01, by=0.002),
                                limits=c(0, 0.011))
    gg <- gg+scale_y_continuous('Events per 1,000 women\n',
                                expand=c(0, 0))
    gg <- gg+scale_colour_viridis(discrete=TRUE,
                                  begin=0.2,
                                  end=0.8,
                                  guide=guide_legend(title='Specificity',
                                                     keywidth=unit(1, 'cm'),
                                                     title.theme=element_text(size=12, face='bold'),
                                                     label.theme=element_text(size=12, angle=0)))
    print(gg)
    if(saveit){
        filename <- paste('hypothetical', datestamp, sep='_')
        filename <- paste(filename, ext, sep='.')
        ggsave(here('plots', filename),
               plot=gg,
               width=6,
               height=9)
    }
}
#hypothetical_plot(saveit=TRUE)

##################################################
# Visualize UCTs in hypothetical analysis
##################################################
hypothetical_uct_plot <- function(ext='png', saveit=FALSE){
    dset <- bind_rows(analysis1(specificity=0.97),
                      analysis1(specificity=0.98),
                      analysis1(specificity=0.99),
                      analysis1(specificity=1.0))
    dset <- dset %>% filter(prevalence == 0.01,
                            localization == 0.8)
    dset <- dset %>% select(-prevalence, -localization, -CD, -CDP)
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
#hypothetical_uct_plot(saveit=TRUE)

##################################################
# Visualize CDs in hypothetical analysis
##################################################
hypothetical_cd_plot <- function(ext='png', saveit=FALSE){
    dset <- analysis1(specificity=0.99)
    dset <- dset %>% ungroup()
    dset <- dset %>% filter(localization == 0.8)
    dset <- dset %>% select(-specificity, -localization, -UCT, -CDP)
    dset <- dset %>% mutate(sensitivity=factor(sensitivity,
                                               levels=unique(sensitivity),
                                               labels=sprintf('%2.0f%%', 100*unique(sensitivity))))
    gg_theme()
    gg <- ggplot(data=dset)
    gg <- gg+geom_line(aes(x=prevalence,
                           y=CD,
                           alpha=sensitivity),
                       size=0.6)
    gg <- gg+scale_x_continuous('\nPrevalence of cancers A and B',
                                labels=percent_format(accuracy=0.1),
                                breaks=seq(0, 0.01, by=0.002))
    gg <- gg+scale_y_continuous('Expected cancers detected per 1,000 women\n',
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
#hypothetical_cd_plot(saveit=TRUE)

##################################################
# Visualize outcomes in empirical analysis
##################################################
empirical_age_plot <- function(dset, ext='png', sensitivity=FALSE, saveit=FALSE){
    dset <- dset %>% mutate(UCT.CD=UCT/CD,
                            UCT.CDP=UCT/CDP)
    dset <- dset %>% select(-UCT, -CD, -CDP)
    dset <- dset %>% pivot_longer(cols=c(UCT.CD, UCT.CDP),
                                  names_to='outcome',
                                  values_to='value')
    dset <- dset %>% ungroup()
    dset <- dset %>% arrange(value)
    dset <- dset %>% mutate(age=sub('-[567]4', ' years', age),
                            outcome=factor(outcome,
                                           levels=c('UCT.CD', 'UCT.CDP'),
                                           labels=c('Unnecessary\nconfirmation\ntests per cancer\ndetected',
                                                    'Unnecessary\nconfirmation\ntests per cancer\ndeath prevented')),
                            site=factor(site, levels=unique(site)))
    if(sensitivity){
        dset <- dset %>% filter(outcome == 'Unnecessary\nconfirmation\ntests per cancer\ndeath prevented')
        height <- 4
    } else {
        height <- 6
    }
    gg_theme(axis.text.x=element_text(size=8, angle=90, vjust=0.5, hjust=1),
             axis.ticks.x=element_blank(),
             strip.text.y=element_text(size=10, angle=0))
    gg <- ggplot(data=dset)
    gg <- gg+geom_bar(aes(x=site, y=value),
                      stat='identity',
                      position='dodge')
    gg <- gg+geom_hline(aes(yintercept=0), colour='black')
    gg <- gg+facet_grid(outcome~age, scales='free_y')
    gg <- gg+scale_x_discrete(name='')
    gg <- gg+scale_y_continuous(name='',
                                expand=c(0, 0))
    print(gg)
    if(saveit){
        filename <- paste('empirical_age', datestamp, sep='_')
        filename <- paste(filename, ext, sep='.')
        ggsave(here('plots', filename),
               plot=gg,
               width=6,
               height=height)
    }
}
#empirical_age_plot(iset1, sensitivity=FALSE, saveit=TRUE)
#empirical_age_plot(iset1, sensitivity=TRUE, saveit=TRUE)

