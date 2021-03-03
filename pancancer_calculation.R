library(dplyr)
library(readxl)

# For Table 1: Anlaysis 1: pan-cancer test for two hypothetical cancers a (fixed prevalence 0.001) and b (varying prevalence 0.0005-0.01),
# at varying levels of specifcity (0.95, 0.99)

# create a result table
result <- data.frame(prevalence = 1, specificity = 1, detect = 1, biopsy = 1, live = 1)

for (s in c(0.95, 0.99)) {
  for (p in c(0.0005, 0.001, 0.005, 0.01)) {
    
    #Sp
    specificity <- s

    #prevalence of a, b and neither
    #Pa
    p.a <- 0.001
    #Pb
    p.b <- p
    #1-Pa-Pb
    p.n <- 1 - p.a - p.b

    #population size: 1,000 women
    N <- 1000
    N.a <- N * p.a
    N.b <- N * p.b
    N.n <- N * p.n

    # assuming 50% of cancer patients would die without screening
    # death rate among cancer a patients
    death.a <- 0.5 * p.a
    # death rate among cancer b patients
    death.b <- 0.5 * p.b

    #overall sensitivity for A
    #PA (T+)
    sensitivity.a <- 0.5
    #overall sensitivity for B
    #PB (T+)
    sensitivity.b <- sensitivity.a

    # probability of correctly localizing b/o given b/o is present and test is positive
    # LA (A│T+)
    PAA <- 0.5
    # LB (B│T+)
    PBB <- PAA

    # fraction of false positives that are localized as cancer A
    # L(A│T-)
    PAN <- p.a / (p.a + p.b)

    # mortality reduction due to screening
    effect.a <- 0.2
    effect.b <- 0.2

    # LA (B│T+)
    PBA <- 1 - PAA
    
    # LB (A│T+)
    PAB <- 1 - PBB
    
    # LB (A│T+)
    PBN <- 1 - PAN

    # EC
    detect <- sensitivity.a * PAA * N.a + sensitivity.b * PBB * N.b

    # nUCT
    biopsy <- sensitivity.a * PBA * N.a + sensitivity.b * PAB * N.b + (1 - specificity) * N.n

    # LS (currently not presented in the paper)
    live <- N.a * effect.a * death.a + N.b * effect.b * death.b

    result <- rbind(
      result,
      data.frame(
        prevalence = p.b,
        specificity = specificity,
        detect = detect,
        biopsy = biopsy,
        live = live
      )
    )
  }
}

# final result
result <- result[-1, ] %>%
  mutate(ratio_bd = biopsy / detect, ratio_bl = biopsy / live) %>%
  dplyr::select(prevalence, specificity, biopsy, detect, ratio_bd)


# For Table 4. Anlaysis 2: pan-cancer test for breast cancer (indicated as b) + one of the other types of cancer (indicated as o)

setwd("~/Dropbox/Pan-cancer test/Figures")
seer <- read_excel("SEER CSR.xlsx")

result <- data.frame(prevalence = 1, specificity = 1, age = 1, detect = 1, biopsy = 1, live = 1, cancer = "b")

for (a in c(50, 60, 70)) {

  # create an input table
  seer_input <- seer %>%
    filter(age == a) %>%
    dplyr::select(cancer, prevalence, death)

  #overall sensitivity
  sensitivity <- c(0.74, 0.68, 0.78, 0.59, 0.64, 0.67)
  #probability of correct localization
  correct <- c(0.97, 0.72, 0.79, 0.92, 0.96, 0.96)
  
  input <- data.frame(seer_input, sensitivity = sensitivity, correct = correct)

  for (i in c(1:6)) {
    specificity <- 0.99

    #prevalence of breast cancer
    p.b <- input$prevalence[input$cancer == "Breast"]
    #prevalence of another cancer
    p.o <- input$prevalence[i]
    #prevalence of neither cancer
    p.n <- 1 - p.b - p.o

    N <- 1000
    N.b <- N * p.b
    N.o <- N * p.o
    N.n <- N * p.n

    death.b <- input$death[input$cancer == "Breast"]
    death.o <- input$death[i]

    #overall sensitivity
    sensitivity.b <- input$sensitivity[input$cancer == "Breast"]
    sensitivity.o <- input$sensitivity[i]

    #probability of correct localization
    PBB <- input$correct[input$cancer == "Breast"]
    POO <- input$correct[i]

    #fraction of false positives localized as breast cancer
    PBN <- p.b / (p.b + p.o)

    #mortality reduction
    effect.b <- 0.2
    effect.o <- 0.2

    #probability of incorrect localization
    POB <- 1 - PBB
    PBO <- 1 - POO
    
    #fraction of false positives localized as another cancer
    PON <- 1 - PBN
    
    #overall specificity 
    specificity <- 0.99

    #EC
    detect <- sensitivity.b * PBB * N.b + sensitivity.o * POO * N.o
    #nUCT
    biopsy <- sensitivity.b * POB * N.b + sensitivity.o * PBO * N.o + (1 - specificity) * N.n
    #LS (Note: here we use N, not N.a or N.b like in Analysis 1, 
    #because death.b refers to death rate in the overall population)
    live <- N * effect.b * death.b + N * effect.o * death.o

    result <- rbind(
      result,
      data.frame(
        prevalence = p.o,
        specificity = specificity,
        detect = detect,
        biopsy = biopsy,
        live = live,
        cancer = input$cancer[i],
        age = a
      )
    )
  }
}

result <- result[-1, ] %>%
  filter(cancer != "Breast") %>%
  mutate(ratio_bd = biopsy / detect, ratio_bl = biopsy / live) %>%
  dplyr::select(age, cancer, biopsy, detect, live, ratio_bd, ratio_bl)


# For Table 5: Analysis 3: pan-cancer test for breast cancer (indicated as b) + lung cancer (indicated as l) 
# + one of the other types of cancer (indicated as c)

result <- data.frame(prevalence = 1, specificity = 1, detect = 1, biopsy = 1, live = 1, cancer = "b", age = 1)

for (a in c(50, 60, 70)) {
  seer_input <- seer %>%
    filter(age == a) %>%
    dplyr::select(cancer, prevalence, death)

  sensitivity <- c(0.74, 0.68, 0.78, 0.59, 0.64, 0.67)
  correct <- c(0.97, 0.72, 0.79, 0.92, 0.96, 0.96)
  input <- data.frame(seer_input, sensitivity = sensitivity, correct = correct)

  for (i in 1:6) {
    #prevalence of breast cancer
    p.b <- input$prevalence[input$cancer == "Breast"]
    #prevalence of lung cancer
    p.l <- input$prevalence[input$cancer == "Lung and Bronchus"]
    #prevalence of another cancer
    p.c <- input$prevalence[i]
    #prevalence of no cancer
    p.n <- 1 - p.b - p.c - p.l

    N <- 1000
    N.b <- N * p.b
    N.l <- N * p.l
    N.c <- N * p.c
    N.n <- N * p.n

    #death rate of breast cancer
    death.b <- input$death[5]
    #death rate of lung cancer
    death.l <- input$death[4]
    #death rate of another cancer
    death.c <- input$death[i]
    
    #overall sensitivity of breast cancer
    sensitivity.b <- input$sensitivity[5]
    #overall sensitivity of lung cancer
    sensitivity.l <- input$sensitivity[4]
    #overall sensitivity of another cancer
    sensitivity.c <- input$sensitivity[i]

    #overall specificity
    specificity <- 0.99

    #overall sensitivity of breast cancer
    PBB <- input$correct[5]
    #overall sensitivity of lung cancer
    PLL <- input$correct[4]
    #overall sensitivity of another cancer
    PCC <- input$correct[i]

    #mortality reduction for breast cancer
    effect.b <- 0.2
    #mortality reduction for lung cancer
    effect.l <- 0.2
    #mortality reduction for another cancer
    effect.c <- 0.2

    #EC
    detect <-
      sensitivity.b * PBB * N.b +
      sensitivity.c * PCC * N.c +
      sensitivity.l * PLL * N.l

    #nUCT
    biopsy <-
      sensitivity.b * (1 - PBB) * N.b +
      sensitivity.c * (1 - PCC) * N.c +
      sensitivity.l * (1 - PLL) * N.l +
      (1 - specificity) * (N.n)

    #LS
    live <-
      N * effect.b * death.b +
      N * effect.c * death.c +
      N * effect.l * death.l

    result <- rbind(
      result,
      data.frame(
        prevalence = p.c,
        specificity = specificity,
        detect = detect,
        biopsy = biopsy,
        live = live,
        cancer = input$cancer[i],
        age = a
      )
    )
  }
}

result <- result[-1, ] %>%
  filter(cancer != "Lung and Bronchus" & cancer != "Breast") %>%
  mutate(ratio_bd = biopsy / detect, ratio_bl = biopsy / live) %>%
  dplyr::select(age, cancer, biopsy, detect, live, ratio_bd, ratio_bl)


# For Table 6: Analysis 4: pan-cancer test for breast cancer (indicated as b) + lung cancer (indicated as l) + ovary cancer (indicated as o) 
# + colon cancer (indicated as c) + liver cancer (indicated as li) + pancreatic cancer (indicated as p)

result <- data.frame(prevalence = 1, specificity = 1, detect = 1, biopsy = 1, live = 1, age = 1)

for (a in c(50, 60, 70)) {
  seer_input <- seer %>%
    filter(age == a) %>%
    dplyr::select(cancer, prevalence, death)

  sensitivity <- c(0.74, 0.68, 0.78, 0.59, 0.64, 0.67)
  correct <- c(0.97, 0.72, 0.79, 0.92, 0.96, 0.96)
  input <- data.frame(seer_input, sensitivity = sensitivity, correct = correct)

  #prevalence of breast cancer
  p.b <- input$prevalence[5]
  #prevalence of colon cancer
  p.c <- input$prevalence[1]
  #prevalence of liver cancer
  p.li <- input$prevalence[2]
  #prevalence of pancreatic cancer
  p.p <- input$prevalence[3]
  #prevalence of lung cancer
  p.l <- input$prevalence[4]
  #prevalence of ovarian cancer
  p.o <- input$prevalence[6]
  #prevalence of no cancer
  p.n <- 1 - p.b - p.c - p.l - p.li - p.p - p.o
  
  N <- 1000
  N.b <- N * p.b
  N.l <- N * p.l
  N.c <- N * p.c
  N.li <- N * p.li
  N.p <- N * p.p
  N.o <- N * p.o

  N.n <- N * p.n

  #death rate of each cancer in the overall population
  death.b <- input$death[5]
  death.c <- input$death[1]
  death.l <- input$death[4]
  death.li <- input$death[2]
  death.p <- input$death[3]
  death.o <- input$death[6]

  #overall sensitivity of each cancer in the overall population
  sensitivity.b <- input$sensitivity[5]
  sensitivity.c <- input$sensitivity[1]
  sensitivity.l <- input$sensitivity[4]
  sensitivity.li <- input$sensitivity[2]
  sensitivity.p <- input$sensitivity[3]
  sensitivity.o <- input$sensitivity[6]

  #overall specificity of each cancer
  specificity <- 0.99

  #probability of correct localiztion of each cancer 
  PBB <- input$correct[5]
  PCC <- input$correct[1]
  PLL <- input$correct[4]
  PLiLi <- input$correct[2]
  PPP <- input$correct[3]
  POO <- input$correct[6]

  #mortality reduction for each cancer
  effect.b <- 0.2
  effect.c <- 0.2
  effect.l <- 0.2
  effect.li <- 0.2
  effect.p <- 0.2
  effect.o <- 0.2

  #EC
  detect <-
    sensitivity.b * PBB * N.b +
    sensitivity.c * PCC * N.c +
    sensitivity.l * PLL * N.l +
    sensitivity.li * PLiLi * N.li +
    sensitivity.p * PPP * N.p +
    sensitivity.o * POO * N.o

  #uUCT
  biopsy <-
    sensitivity.b * (1 - PBB) * N.b +
    sensitivity.c * (1 - PCC) * N.c +
    sensitivity.l * (1 - PLL) * N.l +
    sensitivity.li * (1 - PLiLi) * N.li +
    sensitivity.p * (1 - PPP) * N.p +
    sensitivity.o * (1 - POO) * N.o +
    (1 - specificity) * (N.n)

  #LS
  live <-
    N * effect.b * death.b +
    N * effect.c * death.c +
    N * effect.l * death.l +
    N * effect.li * death.li +
    N * effect.p * death.p +
    N * effect.o * death.o

  result <- rbind(
    result,
    data.frame(
      prevalence = p.c,
      specificity = specificity,
      detect = detect,
      biopsy = biopsy,
      live = live,
      age = a
    )
  )
}

result <- result[-1, ] %>%
  mutate(ratio_bd = biopsy / detect, ratio_bl = biopsy / live) %>%
  dplyr::select(age, biopsy, detect, live, ratio_bd, ratio_bl)
