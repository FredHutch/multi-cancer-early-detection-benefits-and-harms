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

datestamp <- '2021-03-26'

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
# Age-specific outcomes from SEER data
##################################################
age_analysis <- function(bset, aset, specificity=0.99){
    dset <- bind_rows(bset, aset)
    dset <- dset %>% mutate(candidate=unique(bset$site))
    dset <- dset %>% pivot_longer(-c(age, site, candidate),
                                  names_to='feature',
                                  values_to='value')
    dset <- dset %>% group_by(age)
    dset <- dset %>% do(k_cancer_test(., specificity=specificity))
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
# Visualize outcomes in empirical analysis
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
    gg_theme(axis.ticks.x=element_line(colour='black'))
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
# Visualize outcomes in empirical analysis
##################################################
empirical_age_plot <- function(dset, figureno, ext='png', sensitivity=FALSE, saveit=FALSE){
    if(any(grepl('_low$', names(dset)))){
        dset <- dset %>% mutate(UCT.CD=UCT/CD,
                                UCT.LS_low=UCT/LS_low,
                                UCT.LS_high=UCT/LS_high)
        dset <- dset %>% select(-UCT, -CD, -LS_low, -LS_high)
        dset <- dset %>% pivot_longer(cols=c(UCT.CD, UCT.LS_low, UCT.LS_high),
                                      names_to='outcome',
                                      values_to='value')
        dset <- dset %>% separate(outcome, c('outcome', 'assumption'), sep='_', fill='right')
        dset <- dset %>% mutate(assumption=ifelse(is.na(assumption), 'point', assumption),
                                assumption=factor(assumption, levels=c('point', 'high', 'low')))
    } else {
        dset <- dset %>% mutate(UCT.CD=UCT/CD, UCT.LS=UCT/LS)
        dset <- dset %>% select(-UCT, -CD, -LS)
        dset <- dset %>% pivot_longer(cols=c(UCT.CD, UCT.LS),
                                      names_to='outcome',
                                      values_to='value')
        dset <- dset %>% mutate(assumption='high')
    }
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
    ymax <- switch(as.character(figureno), '2'=300, '3'=200, 'S1'=150, 'S2'=80)
    gg_theme(legend.position='none',
             axis.text.x=element_text(size=10, angle=90, vjust=0.5, hjust=1),
             axis.ticks.x=element_blank(),
             panel.spacing.y=unit(0.04, 'npc'),
             strip.text.y=element_text(size=10, angle=0))
    gg <- ggplot(data=dset)
    gg <- gg+geom_bar(aes(x=site, y=value, fill=assumption),
                      stat='identity',
                      position='dodge')
    gg <- gg+geom_hline(aes(yintercept=0), colour='black')
    if(!sensitivity)
        gg <- gg+geom_blank(data=dset %>% filter(outcome == 'Unnecessary\nconfirmation\ntests per cancer\ndetected'), aes(y=8))
    gg <- gg+geom_blank(data=dset %>% filter(outcome == 'Unnecessary\nconfirmation\ntests per life\nsaved'), aes(y=ymax))
    gg <- gg+facet_grid(outcome~age, scales='free_y')
    gg <- gg+scale_x_discrete(name='')
    if(!sensitivity)
        gg <- gg+labs(tag=paste('Figure', figureno))
    gg <- gg+scale_y_continuous(name='', expand=c(0, 0))
    gg <- gg+scale_fill_manual(name='', values=c(point='gray40',
                                                 high='gray40',
                                                 low='gray80'))
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
# Format table for realistic six-cancer test
##################################################
format_empirical <- function(dset, saveit=FALSE){
    dset <- dset %>% select(Test, age, UCT, CD, LS)
    dset <- dset %>% mutate(age=sub('-[567]4', '', age),
                            UCT.CD=UCT/CD,
                            UCT.LS=UCT/LS)
    dset <- dset %>% mutate(UCT=sprintf('%4.1f', UCT),
                            CD=sprintf('%3.1f', CD),
                            LS=sprintf('%3.1f', LS),
                            UCT.CD=sprintf('%3.1f', UCT.CD),
                            UCT.LS=sprintf('%3.1f', UCT.LS))
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
##################################################
#pset <- pset %>% mutate(effect=0.2)
#sset <- full_join(sset, pset, by='site')
#iset6 <- age_analysis_incremental(sset, setdiff(pset$site, 'Breast'))
#iset6 <- iset6 %>% mutate(Test='Pan-cancer')
#bset <- sset %>% filter(site == 'Breast')
#bset <- age_analysis(bset, bset[FALSE, ])
#bset <- bset %>% mutate(Test='Breast only')
#mset <- sset %>% filter(site == 'Breast')
#mset <- mset %>% mutate(sensitivity=0.869, localization=1)
#mset <- age_analysis(mset, mset[FALSE, ], specificity=0.889)
#mset <- mset %>% mutate(Test='Mammography')
#cset <- bind_rows(iset6, bset, mset)
#format_empirical(cset, saveit=TRUE)

##################################################
# Figure 2
##################################################
#lset <- sset %>% mutate(effect=0.05)
#iset1l <- age_analysis_incremental(lset, 'Breast')
#iset1h <- age_analysis_incremental(sset, 'Breast')
#iset1 <- full_join(iset1l, iset1h, by=c('site', 'age', 'UCT', 'CD'), suffix=c('_low', '_high'))
#empirical_age_plot(iset1, figureno=2, ext='pdf', saveit=TRUE)

##################################################
# Figure 3
##################################################
#iset2l <- age_analysis_incremental(lset, c('Breast', 'Lung'))
#iset2h <- age_analysis_incremental(sset, c('Breast', 'Lung'))
#iset2 <- full_join(iset2l, iset2h, by=c('site', 'age', 'UCT', 'CD'), suffix=c('_low', '_high'))
#empirical_age_plot(iset2, figureno=3, ext='pdf', saveit=TRUE)

##################################################
# Supplemental Figure 1
##################################################
#sset1s <- sset %>% mutate(effect=ifelse(site %in% c('Breast', 'Colorectal', 'Lung'), 0.1, 0.5))
#iset1s <- age_analysis_incremental(sset1s, 'Breast')
#empirical_age_plot(iset1s, figureno='S1', sensitivity=TRUE, saveit=TRUE)

