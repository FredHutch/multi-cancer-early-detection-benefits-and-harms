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

datestamp <- '2021-05-15'

##################################################
# Project outcomes for k-cancer test
#
# Notation in paper         Variable name
# ---------------------------------------
# P_A                       prevalence (of A)
# P_A(T+)                   sensitivity (of A)
# L_A(T+)                   localization (of A)
# Sp                        specificity
# m_A                       mortality (of A)
# R_A                       mortality reduction (of A)
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
                                     value[feature == 'sensitivity']*
                                     value[feature == 'localization']*
                                     value[feature == 'mortality'])
    lives <- with(lset, size*sum(lives))
    tibble(UCT=tests, CD=cancers, LS=lives)
}

##################################################
# Pan-cancer test for two hypothetical cancers
# across selected scenarios:
# 1. varying sensitivity of both cancers=0.5 to 0.9
# 2. fixed prevalence of cancer A=0.001
# 3. varying prevalence of cancer B=0 to 0.01
# 4. varying localization of both cancers=0.5 to 0.8
# 5. varying specificity=0.95 to 0.99
##################################################
hypothetical_test <- function(specificity){
    dset <- expand_grid(prevalence.a=0.005,
                        prevalence.b=c(0, 0.001, 0.005, 0.01, 0.02),
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
    dset <- read_csv(here('results', filename), col_types='ccddddddddd')
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
                                                 'Colorectal',
                                                 'Liver',
                                                 'Lung',
                                                 'Ovary',
                                                 'Pancreas')))
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
# Age-specific incremental impact of candidate cancer
##################################################
age_analysis_incremental <- function(dset, existing, specificity=0.99){
    aset <- dset %>% filter(site %in% existing)
    bset <- dset %>% filter(!site %in% existing)
    bset <- bset %>% group_by(site)
    rset <- bset %>% do(age_analysis(., aset, specificity=specificity))
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
hypothetical_uct_plot <- function(ext='svg', saveit=FALSE){
    dset <- bind_rows(hypothetical_test(specificity=0.97),
                      hypothetical_test(specificity=0.98),
                      hypothetical_test(specificity=0.99),
                      hypothetical_test(specificity=1.0))
    dset <- dset %>% filter(prevalence == 0.005,
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
                                limits=c(0, 35),
                                breaks=seq(0, 35, by=5),
                                expand=c(0, 0))
    gg <- gg+guides(alpha=guide_legend(title='Sensitivity for\ncancers A and B',
                                       keywidth=unit(1, 'cm'),
                                       title.theme=element_text(size=12),
                                       label.theme=element_text(size=12, angle=0)))
    print(gg)
    if(saveit){
        filename <- paste('hypothetical_uct', datestamp, sep='_')
        filename <- paste(filename, ext, sep='.')
        ggsave(here('results', filename),
               plot=gg,
               width=6,
               height=6)
    }
}

##################################################
# Visualize CDs in hypothetical analysis
##################################################
hypothetical_cd_plot <- function(ext='svg', saveit=FALSE){
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
                                breaks=seq(0, 0.02, by=0.005))
    gg <- gg+scale_y_continuous('Cancers detected per 1,000 women\n',
                                limits=c(0, 35),
                                breaks=seq(0, 35, by=5),
                                expand=c(0, 0))
    gg <- gg+guides(alpha=guide_legend(title='Sensitivity for\ncancers A and B',
                                       keywidth=unit(1, 'cm'),
                                       title.theme=element_text(size=12),
                                       label.theme=element_text(size=12, angle=0)))
    print(gg)
    if(saveit){
        filename <- paste('hypothetical_cd', datestamp, sep='_')
        filename <- paste(filename, ext, sep='.')
        ggsave(here('results', filename),
               plot=gg,
               width=6,
               height=6)
    }
}

##################################################
# Visualize outcomes in empirical analysis
##################################################
empirical_age_plot <- function(dset, figureno, ext='png', sensitivity=FALSE, saveit=FALSE){
    if(any(grepl('_low$', names(dset)))){
        dset <- dset %>% mutate(UCT.CD=UCT/CD,
                                UCT.LS_low=UCT/LS_low,
                                UCT.LS_mid=UCT/LS_mid,
                                UCT.LS_high=UCT/LS_high)
        dset <- dset %>% select(-UCT, -CD, -LS_low, -LS_mid, -LS_high)
        dset <- dset %>% pivot_longer(cols=c(UCT.CD, UCT.LS_low, UCT.LS_mid, UCT.LS_high),
                                      names_to='outcome',
                                      values_to='value')
        dset <- dset %>% separate(outcome, c('outcome', 'assumption'), sep='_', fill='right')
        dset <- dset %>% mutate(assumption=ifelse(is.na(assumption), 'point', assumption),
                                assumption=factor(assumption, levels=c('point', 'high', 'mid', 'low')))
        dset <- dset %>% arrange(desc(outcome), assumption, value)
    } else {
        dset <- dset %>% mutate(UCT.CD=UCT/CD, UCT.LS=UCT/LS)
        dset <- dset %>% select(-UCT, -CD, -LS)
        dset <- dset %>% pivot_longer(cols=c(UCT.CD, UCT.LS),
                                      names_to='outcome',
                                      values_to='value')
        dset <- dset %>% mutate(assumption='high')
        dset <- dset %>% arrange(age, desc(outcome), value)
    }
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
    ymax <- switch(as.character(figureno), '2'=80, '3'=50, 'S1'=100, 'S2'=30, 'S3'=40)
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
        gg <- gg+geom_blank(data=dset %>% filter(outcome == 'Unnecessary\nconfirmation\ntests per cancer\ndetected'), aes(y=2))
    gg <- gg+geom_blank(data=dset %>% filter(outcome == 'Unnecessary\nconfirmation\ntests per life\nsaved'), aes(y=ymax))
    gg <- gg+facet_grid(outcome~age, scales='free_y')
    gg <- gg+scale_x_discrete(name='')
    if(!sensitivity)
        gg <- gg+labs(tag=paste('Figure', figureno))
    gg <- gg+scale_y_continuous(name='', expand=c(0, 0))
    gg <- gg+scale_fill_manual(name='', values=c(point='gray40',
                                                 high='gray40',
                                                 mid='gray60',
                                                 low='gray80'))
    print(gg)
    if(saveit){
        filename <- str_glue('figure{figureno}_{datestamp}')
        filename <- paste(filename, ext, sep='.')
        ggsave(here('results', filename),
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
    dset <- dset %>% arrange(sensitivity,
                             localization,
                             specificity,
                             prevalence)
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
        write_csv(dset, here('results', filename))
    }
}

##################################################
# Format table of SEER incidence and mortality rates
##################################################
format_seer <- function(dset, fmt='%5.2f', saveit=FALSE){
    dset <- dset %>% select(site, age, prevalence, mortality)
    dset <- dset %>% mutate(site=factor(site, levels=c('Breast',
                                                       'Colorectal',
                                                       'Liver',
                                                       'Lung',
                                                       'Ovary',
                                                       'Pancreas')),
                            age=sub('-[567]4', '', age),
                            prevalence=sprintf(fmt, 100*prevalence),
                            mortality=sprintf(fmt, 100*mortality))
    dset <- dset %>% arrange(site, age)
    dset <- dset %>% rename('Tissue of origin'='site',
                            'Age, y'='age',
                            '5-year probability of diagnosis, %'=prevalence,
                            '15-year probability of death, %'=mortality)
    if(saveit){
        filename <- str_glue('table1_{datestamp}.csv')
        write_csv(dset, here('results', filename))
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
        write_csv(dset, here('results', filename))
    }
}

##################################################
# Format table for sensitivity analyses
##################################################
format_supplemental <- function(dset, tableno, saveit=FALSE){
    dset <- dset %>% select(age, site, UCT, CD, LS_low, LS_mid, LS_high)
    dset <- dset %>% mutate(age=sub('-[567]4', '', age),
                            site=factor(site, levels=c('Colorectal',
                                                       'Liver',
                                                       'Lung',
                                                       'Ovary',
                                                       'Pancreas')),
                            UCT=sprintf('%4.1f', UCT),
                            CD=sprintf('%3.1f', CD),
                            LS_low=sprintf('%3.1f', LS_low),
                            LS_mid=sprintf('%3.1f', LS_mid),
                            LS_high=sprintf('%3.1f', LS_high))
    dset <- dset %>% arrange(age, site)
    dset <- dset %>% rename('Tissue of origin'='site',
                            'Screening age, y'='age',
                            'Unnecessary confirmation tests, n'='UCT',
                            'Cancers detected, n'='CD',
                            'Lives saved (5%), n'='LS_low',
                            'Lives saved (10%), n'='LS_mid',
                            'Lives saved (20%), n'='LS_high')
    if(saveit){
        filename <- str_glue('supplemental_table{tableno}_{datestamp}.csv')
        write_csv(dset, here('results', filename))
    }
}

##################################################
# Read merged SEER incidence and mortality rates
##################################################
sset <- read_data('seer_merged_2000-2002_followup=15_2021-05-13.csv')

##################################################
# Table 1
##################################################
format_seer(sset, saveit=TRUE)

##################################################
# Table 2
##################################################
pset <- tribble(~site, ~sensitivity, ~localization,
                'Breast',     0.64, 0.96,
                'Colorectal', 0.74, 0.97,
                'Lung',       0.59, 0.92,
                'Ovary',      0.67, 0.96,
                'Pancreas',   0.78, 0.79,
                'Liver',      0.68, 0.72)
pset <- pset %>% mutate(marginal=sensitivity*localization,
                        effect=0.1/marginal)
sset <- full_join(sset, pset, by='site')

##################################################
# Table 3
##################################################
iset6 <- age_analysis_incremental(sset, setdiff(pset$site, 'Breast'))
iset6 <- iset6 %>% mutate(Test='Pan-cancer')
bset <- sset %>% filter(site == 'Breast')
bset <- age_analysis(bset, bset[FALSE, ])
bset <- bset %>% mutate(Test='Breast only')
mset <- sset %>% filter(site == 'Breast')
mset <- mset %>% mutate(sensitivity=0.869, localization=1)
mset <- age_analysis(mset, mset[FALSE, ], specificity=0.889)
mset <- mset %>% mutate(Test='Mammography')
cset <- bind_rows(iset6, bset, mset)
format_empirical(cset, saveit=TRUE)

##################################################
# Figure 1
##################################################
hypothetical_uct_plot(saveit=TRUE)
hypothetical_cd_plot(saveit=TRUE)

##################################################
# Figure 2
##################################################
lset <- sset %>% mutate(effect=0.05/marginal)
hset <- sset %>% mutate(effect=0.2/marginal)
fig2l <- age_analysis_incremental(lset, 'Breast')
fig2m <- age_analysis_incremental(sset, 'Breast')
fig2h <- age_analysis_incremental(hset, 'Breast')
fig2 <- full_join(fig2l, fig2m, by=c('site', 'age', 'UCT', 'CD'), suffix=c('_low', '_mid'))
fig2 <- full_join(fig2, fig2h, by=c('site', 'age', 'UCT', 'CD'))
fig2 <- fig2 %>% rename(LS_high=LS)
empirical_age_plot(fig2, figureno=2, ext='pdf', saveit=TRUE)

##################################################
# Figure 3
##################################################
fig3l <- age_analysis_incremental(lset, c('Breast', 'Lung'))
fig3m <- age_analysis_incremental(sset, c('Breast', 'Lung'))
fig3h <- age_analysis_incremental(hset, c('Breast', 'Lung'))
fig3 <- full_join(fig3l, fig3m, by=c('site', 'age', 'UCT', 'CD'), suffix=c('_low', '_mid'))
fig3 <- full_join(fig3, fig3h, by=c('site', 'age', 'UCT', 'CD'))
fig3 <- fig3 %>% rename(LS_high=LS)
empirical_age_plot(fig3, figureno=3, ext='pdf', saveit=TRUE)

##################################################
# Supplemental Table 1
##################################################
hset <- bind_rows(hypothetical_test(specificity=0.95),
                  hypothetical_test(specificity=0.99))
format_hypothetical(hset, saveit=TRUE)

##################################################
# Supplemental Table 2
##################################################
format_supplemental(fig2, tableno=2, saveit=TRUE)

##################################################
# Supplemental Table 3
##################################################
format_supplemental(fig3, tableno=3, saveit=TRUE)

##################################################
# Supplemental Figure 1
##################################################
figs1 <- age_analysis_incremental(sset, 'Breast', specificity=0.97)
empirical_age_plot(figs1, figureno='S1', sensitivity=TRUE, saveit=TRUE)

##################################################
# Supplemental Figure 2
##################################################
ssets2 <- sset %>% mutate(effect=ifelse(site %in% c('Breast', 'Colorectal', 'Lung'),
                                        0.1/marginal,
                                        0.5/marginal))
figs2 <- age_analysis_incremental(ssets2, 'Breast')
empirical_age_plot(figs2, figureno='S2', sensitivity=TRUE, saveit=TRUE)

##################################################
# Supplemental Figure 3
##################################################
ssets3 <- read_data('seer_merged_2000-2002_followup=10_2021-05-13.csv')
ssets3 <- full_join(ssets3, pset, by='site')
figs3 <- age_analysis_incremental(ssets3, 'Breast')
empirical_age_plot(figs3, figureno='S3', sensitivity=TRUE, saveit=TRUE)

