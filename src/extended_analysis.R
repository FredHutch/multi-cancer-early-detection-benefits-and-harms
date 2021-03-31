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

#datestamp <- '2021-03-26'
#datestamp <- '2021-03-30'
datestamp <- '2021-03-31'

##################################################
# Project outcomes for k-cancer test
#
# Notation in paper         Variable name
# ---------------------------------------
# P_A                       prevalence (of A)
# P_A(T+)                   sensitivity (of A)
# L_A(T+)                   localization (of A)
# Sp                        overall specificity
# m                         mortality (of A)
# r                         mortality reduction (of A)
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
# Read merged SEER incidence and IBM data
##################################################
read_data <- function(filename){
    dset <- read_csv(here('data', filename), col_types='ccccddddddddd')
    dset <- dset %>% rename(sex='Sex', race='Race', age='Age', site='Site')
    dset <- dset %>% select(-Population,
                            -Diagnosis.Count,
                            -Death.Count,
                            -ends_with('Lower'),
                            -ends_with('Upper'))
    dset <- dset %>% rename(prevalence='Diagnosis.Rate',
                            mortality='Death.Rate')
    dset <- dset %>% group_by(site, sex)
    dset <- dset %>% mutate(absent=all(prevalence == 0))
    dset <- dset %>% filter(!absent)
    dset <- dset %>% select(-absent)
    dset <- dset %>% ungroup()
}

##################################################
# Quantify harms and benefits for initial site
##################################################
single_cancer_test <- function(dset, specificity=0.99){
    dset <- dset %>% pivot_longer(-site,
                                  names_to='feature',
                                  values_to='value')
    dset <- dset %>% group_by(site)
    dset <- dset %>% do(k_cancer_test(., specificity=specificity))
    dset <- dset %>% ungroup()
}

##################################################
# Control plot aesthetics
##################################################
gg_theme <- function(...){
    theme_set(theme_classic())
    theme_update(axis.text=element_text(size=10, colour='black'),
                 axis.line=element_line(colour='darkgray'),
                 axis.title=element_text(size=14),
                 axis.ticks=element_blank(),
                 axis.ticks.length=unit(0.15, 'cm'),
                 panel.border=element_rect(fill=NA, colour='darkgray'),
                 panel.background=element_rect(fill='lightgray'),
                 strip.background=element_rect(colour=NA, fill=NA),
                 strip.text=element_text(size=14))
    theme_update(...)
}

##################################################
# Visualize heatmap of SEER data
##################################################
seer_heatmap <- function(dset,
                         Race,
                         radix=1e5,
                         smidgen=0.1,
                         ext='png',
                         saveit=FALSE){
    dset <- dset %>% filter(race == Race)
    dset <- dset %>% mutate(age=factor(age),
                            #age=factor(sub('-[567]4', ' y', age)),
                            site=factor(site, levels=rev(sort(unique(site)))),
                            prevalence=ifelse(prevalence == 0,
                                              smidgen,
                                              prevalence*radix),
                            mortality=mortality*radix)
    gg_theme()
    gg <- ggplot(dset)
    gg <- gg+geom_tile(aes(x=age, y=site, fill=prevalence))
    gg <- gg+geom_point(aes(x=age, y=site, size=mortality))
    gg <- gg+facet_grid(.~sex)
    gg <- gg+scale_fill_viridis(name='Incidence rate\nper 100,000\n',
                                trans='log10',
                                option='cividis',
                                limits=c(0.1, 2000),
                                breaks=c(1, 10, 100, 1000),
                                guide=guide_colourbar(nbin=500,
                                                      title.hjust=0.5,
                                                      ticks=FALSE))
    gg <- gg+scale_size_continuous(name='15-year IBM rate\nper 100,000\n',
                                   limits=c(0.1, 1000),
                                   breaks=c(1, 10, 100, 1000),
                                   range=c(0, 4),
                                   guide=guide_legend(title.hjust=0.5,
                                                      reverse=TRUE))
    gg <- gg+geom_hline(yintercept=seq(1.5, nlevels(dset$site)-0.5),
                        colour='darkgray')
    gg <- gg+geom_vline(xintercept=seq(1.5, nlevels(dset$age)-0.5),
                        colour='darkgray')
    gg <- gg+scale_x_discrete(name='', expand=c(0, 0))
    gg <- gg+scale_y_discrete(name='', expand=c(0, 0))
    if(saveit){
        filename <- str_glue('heatmap_{tolower(Race)}_{datestamp}')
        filename <- paste(filename, ext, sep='.')
        ggsave(here('plots', filename),
               plot=gg,
               height=9,
               width=10)
    }
    return(gg)
}

##################################################
# Visualize single-cancer test outcomes
##################################################
single_cancer_plot <- function(dset, Race, ext='png', saveit=FALSE){
    dset <- dset %>% unnest(cols=data)
    dset <- dset %>% filter(race == Race)
    dset <- dset %>% select(-race, -UCT)
    dset <- dset %>% pivot_longer(-c(sex, age, site),
                                  names_to='outcome',
                                  values_to='value')
    dset <- dset %>% mutate(age=factor(sub('-[567]4', ' y', age)),
                            site=factor(site, levels=rev(sort(unique(site)))),
                            outcome=factor(outcome,
                                           levels=c('CD', 'LS'),
                                           labels=c('Cancers detected',
                                                    'Lives saved')))
    gg_theme(panel.grid.major.x=element_line(colour='darkgray'),
             axis.ticks.x=element_line(colour='black'))
    gg <- ggplot(data=dset)
    gg <- gg+geom_point(aes(x=value,
                            y=site,
                            alpha=outcome,
                            fill=outcome),
                        shape=21,
                        size=1.75,
                        show.legend=c(alpha=FALSE, fill=TRUE))
    gg <- gg+facet_grid(.~sex+age)
    gg <- gg+scale_x_continuous(name='',
                                limits=c(0, 10),
                                breaks=seq(0, 10, by=2))
    gg <- gg+scale_y_discrete(name='')
    gg <- gg+geom_hline(yintercept=seq(1.5, nlevels(dset$site)-0.5),
                        colour='darkgray')
    gg <- gg+scale_alpha_manual(values=c('Cancers detected'=1,
                                         'Lives saved'=0.65))
    gg <- gg+scale_fill_viridis(name='Expected outcomes\nper 1,000 persons',
                                discrete=TRUE,
                                begin=0.25,
                                end=0.75,
                                guide=guide_legend(title.hjust=0.5))
    print(gg)
    if(saveit){
        filename <- str_glue('single_{tolower(Race)}_{datestamp}')
        filename <- paste(filename, ext, sep='.')
        ggsave(here('plots', filename),
               plot=gg,
               width=12,
               height=10)
    }
}

##################################################
# Visualize MCED outcomes
##################################################
k_cancer_plot <- function(dset, Sex, Race, ext='png', saveit=FALSE){
    dset <- dset %>% filter(sex == Sex, race == Race)
    dset <- dset %>% arrange(desc(age), desc(CD))
    dset <- dset %>% mutate(age=factor(sub('-[567]4', ' y', age)),
                            site=factor(site, levels=rev(unique(site))))
    dset <- dset %>% group_by(age)
    dset <- dset %>% mutate(Minimum.UCT=log10(UCT),
                            Maximum.UCT=log10(cumsum(UCT)),
                            CD=cumsum(CD),
                            LS=cumsum(LS))
    dset <- dset %>% select(-sex, -race, -UCT)
    lset <- dset %>% pivot_longer(-c(age, site),
                                  names_to='outcome',
                                  values_to='value')
    lset <- lset %>% mutate(outcome=factor(outcome,
                                           levels=c('Minimum.UCT',
                                                    'Maximum.UCT',
                                                    'CD',
                                                    'LS'),
                                           labels=c('Minimum UCTs (log scale)',
                                                    'Maximum UCTs (log scale)',
                                                    'Cancers detected',
                                                    'Lives saved')))
    gg_theme(panel.grid.major.x=element_line(colour='darkgray'),
             axis.ticks.x=element_line(colour='black'))
    gg <- ggplot(data=lset)
    gg <- gg+geom_point(aes(x=value,
                            y=site,
                            fill=outcome),
                        size=1.75,
                        alpha=0.8,
                        shape=21)
    gg <- gg+geom_segment(data=dset,
                          aes(x=Minimum.UCT,
                              xend=Maximum.UCT,
                              y=site,
                              yend=site),
                          colour='#482576FF',
                          size=2,
                          alpha=0.4,
                          show.legend=FALSE)
    gg <- gg+facet_grid(.~age)
    gg <- gg+scale_x_continuous(name='',
                                limits=c(0, 20),
                                breaks=seq(0, 20, by=2))
    gg <- gg+scale_y_discrete(name='')
    gg <- gg+geom_hline(yintercept=seq(1.5, nlevels(dset$site)-0.5),
                        colour='darkgray')
    gg <- gg+scale_fill_manual(name='Expected outcomes\nper 1,000 persons',
                           values=c('Lives saved'='#BBDF27FF',
                                    'Cancers detected'='#21908CFF',
                                    'Minimum UCTs (log scale)'='#482576FF',
                                    'Maximum UCTs (log scale)'='#482576FF'),
                                guide=guide_legend(title.hjust=0.5))
    print(gg)
    if(saveit){
        filename <- str_glue('multi_{tolower(Sex)}_{tolower(Race)}_{datestamp}')
        filename <- paste(filename, ext, sep='.')
        ggsave(here('plots', filename),
               plot=gg,
               width=12,
               height=10)
    }
}

#################################################
# Read extended merged data
##################################################
#dset <- read_data('seer_merged_2000-2002_followup=15_extended_2021-03-26.csv')

##################################################
# Visualize heatmap of SEER data
##################################################
#seer_heatmap(dset, 'All', saveit=TRUE)
#seer_heatmap(dset, 'Black', saveit=TRUE)

##################################################
# Specify default sensitivity, localization probability
# and mortality reduction
##################################################
#dset <- dset %>% mutate(sensitivity=0.67,
#                        localization=0.90,
#                        effect=0.05)

##################################################
# Single-cancer test harms and benefits
##################################################
#sset <- dset %>% group_by(sex, race, age)
#sset <- sset %>% nest()
#sset <- sset %>% mutate(data=map(data, single_cancer_test))
#sset <- sset %>% ungroup()

##################################################
# Visualize single-cancer test harms and benefits
##################################################
single_cancer_plot(sset, 'All', saveit=TRUE)
single_cancer_plot(sset, 'Black', saveit=TRUE)

##################################################
# MCED test harms and benefits
##################################################
#mset <- sset %>% unnest(data)
#mset <- mset %>% ungroup()

##################################################
# Visualize MCED test harms and benefits
##################################################
k_cancer_plot(mset, 'Female', 'All', saveit=TRUE)
k_cancer_plot(mset, 'Male', 'All', saveit=TRUE)
#k_cancer_plot(mset, 'Female', 'Black', saveit=TRUE)
#k_cancer_plot(mset, 'Male', 'Black', saveit=TRUE)

