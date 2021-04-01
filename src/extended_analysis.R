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
k_cancer_test <- function(dset, specificity=0.99, size=100000){
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
    gg <- gg+scale_fill_manual(name='Expected outcomes\nper 1,000 persons',
                               values=c('Cancers detected'='#3B528BFF',
                                        'Lives saved'='#5DC863FF'),
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
# Visualize absolute MCED outcomes
##################################################
absolute_mced_plot <- function(dset, Sex, Outcome, ext='png', saveit=FALSE){
    dset <- dset %>% filter(sex == Sex)
    dset <- dset %>% arrange(race, desc(age), desc(!!sym(Outcome)))
    dset <- dset %>% mutate(age=factor(sub('-[567]4', ' y', age)),
                            site=factor(site, levels=rev(unique(site))))
    dset <- dset %>% group_by(race, age)
    dset <- dset %>% mutate(CD=cumsum(CD),
                            LS=cumsum(LS),
                            race=factor(race,
                                        levels=c('All', 'Black'),
                                        labels=c('All races', 'African Americans')))
    dset <- dset %>% select(-sex, -setdiff(c('UCT', 'CD', 'LS'), Outcome))
    if(Sex == 'Female'){
        ymin <- switch(Outcome,  UCT=1000, CD=0, LS=0)
        ymax <- switch(Outcome,  UCT=1050, CD=1000, LS=60)
        ystep <- switch(Outcome,  UCT=20, CD=200, LS=20)
    } else {
        ymin <- switch(Outcome,  UCT=1000, CD=0, LS=0)
        ymax <- switch(Outcome,  UCT=1100, CD=2200, LS=120)
        ystep <- switch(Outcome,  UCT=20, CD=500, LS=20)
    }
    gg_theme(panel.grid.major.x=element_line(colour='darkgray'),
             panel.spacing=unit(0.02, 'npc'),
             axis.ticks.x=element_line(colour='black'))
    gg <- ggplot(data=dset)
    gg <- gg+geom_point(aes_string(x=Outcome,
                                   y='site',
                                   fill='race'),
                        size=1.75,
                        alpha=0.8,
                        shape=21)
    gg <- gg+facet_grid(.~age)
    gg <- gg+scale_x_continuous(name='',
                                limits=c(ymin, ymax),
                                breaks=seq(ymin, ymax, by=ystep))
    gg <- gg+scale_y_discrete(name='')
    gg <- gg+geom_hline(yintercept=seq(1.5, nlevels(dset$site)-0.5),
                        colour='darkgray')
    gg <- gg+scale_fill_manual(name='Expected outcomes\nper 100,000 persons',
                               values=c('All races'='#482576FF',
                                        'African Americans'='#BBDF27FF'),
                               guide=guide_legend(title.hjust=0.5))
    print(gg)
    if(saveit){
        filename <- str_glue('mced_absolute_{tolower(Sex)}_{tolower(Outcome)}_{datestamp}')
        filename <- paste(filename, ext, sep='.')
        ggsave(here('plots', filename),
               plot=gg,
               width=12,
               height=10)
    }
}

##################################################
# Visualize relative MCED outcomes
##################################################
relative_mced_plot <- function(dset, Sex, Outcome, ext='png', saveit=FALSE){
    dset <- dset %>% filter(sex == Sex)
    dset <- dset %>% arrange(race, desc(age), desc(!!sym(Outcome)))
    dset <- dset %>% mutate(!!sym(Outcome):=cumsum(!!sym(Outcome)),
                            Outcome=UCT/!!sym(Outcome),
                            race=factor(race,
                                        levels=c('All', 'Black'),
                                        labels=c('All races', 'African Americans')))
    dset <- dset %>% mutate(age=factor(sub('-[567]4', ' y', age)),
                            site=factor(site, levels=rev(unique(site))))
    dset <- dset %>% select(-sex, -UCT, -CD, -LS)
    if(Sex == 'Female'){
        ymin <- switch(Outcome,  CD=0, LS=0)
        ymax <- switch(Outcome,  CD=4, LS=80)
        ystep <- switch(Outcome,  CD=1, LS=20)
    } else {
        ymin <- switch(Outcome,  CD=0, LS=0)
        ymax <- switch(Outcome,  CD=2, LS=40)
        ystep <- switch(Outcome,  CD=0.5, LS=10)
    }
    gg_theme(panel.grid.major.x=element_line(colour='darkgray'),
             panel.spacing=unit(0.02, 'npc'),
             axis.ticks.x=element_line(colour='black'))
    gg <- ggplot(data=dset)
    gg <- gg+geom_point(aes(x=Outcome,
                            y=site,
                            fill=race),
                        size=1.75,
                        alpha=0.8,
                        shape=21)
    gg <- gg+facet_grid(.~age)
    gg <- gg+scale_x_continuous(name='',
                                limits=c(ymin, ymax),
                                breaks=seq(ymin, ymax, by=ystep))
    gg <- gg+scale_y_discrete(name='')
    gg <- gg+geom_hline(yintercept=seq(1.5, nlevels(dset$site)-0.5),
                        colour='darkgray')
    gg <- gg+scale_fill_manual(name='Expected outcomes\nper 100,000 persons',
                               values=c('All races'='#482576FF',
                                        'African Americans'='#BBDF27FF'),
                               guide=guide_legend(title.hjust=0.5))
    print(gg)
    if(saveit){
        filename <- str_glue('mced_relative_{tolower(Sex)}_{tolower(Outcome)}_{datestamp}')
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
#single_cancer_plot(sset, 'All', saveit=TRUE)
#single_cancer_plot(sset, 'Black', saveit=TRUE)

##################################################
# MCED test harms and benefits
##################################################
#mset <- sset %>% unnest(data)
#mset <- mset %>% ungroup()

##################################################
# Visualize absolute MCED test harms and benefits
##################################################
#absolute_mced_plot(mset, 'Female', 'UCT', saveit=TRUE)
#absolute_mced_plot(mset, 'Female', 'CD', saveit=TRUE)
#absolute_mced_plot(mset, 'Female', 'LS', saveit=TRUE)

#absolute_mced_plot(mset, 'Male', 'UCT', saveit=TRUE)
#absolute_mced_plot(mset, 'Male', 'CD', saveit=TRUE)
#absolute_mced_plot(mset, 'Male', 'LS', saveit=TRUE)

##################################################
# Visualize relative MCED test harms and benefits
##################################################
#relative_mced_plot(mset, 'Female', 'CD', saveit=TRUE)
#relative_mced_plot(mset, 'Female', 'LS', saveit=TRUE)

#relative_mced_plot(mset, 'Male', 'CD', saveit=TRUE)
#relative_mced_plot(mset, 'Male', 'LS', saveit=TRUE)
