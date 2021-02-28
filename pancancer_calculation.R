library(dplyr)
library(readxl)

# Anlaysis 1: pan-cancer test for two hypothetical cancers b (fixed prevalence 0.001) and o (varying prevalence 0.0005-0.01),
# at varying levels of specifcity (0.95, 0.99)

# create a result table
result <- data.frame(prevalence = 1, specificity = 1, detect = 1, biopsy = 1, live = 1)

for (s in c(0.95, 0.99)) {
  for (p in c(0.0005, 0.001, 0.005, 0.01)) {
    specificity <- s

    #prevalence of b, o and neither
    p.b <- 0.001
    p.o <- p
    p.n <- 1 - p.b - p.o

    N <- 1000
    N.b <- N * p.b
    N.o <- N * p.o
    N.n <- N * p.n

    # assuming 50% of cancer patients would die without screening
    death.b <- 0.5 * p.b
    death.o <- 0.5 * p.o

    sensitivity.b <- 0.5
    sensitivity.o <- sensitivity.b

    # probability of correctly localizing b/o given b/o is present and test is positive
    PBB <- 0.5
    POO <- PBB

    # probability of misclassifying as b given neither cancer is present and test is positive
    PBN <- p.b / (p.b + p.o)

    # mortality reduction due to screening
    effect.b <- 0.2
    effect.o <- 0.2

    POB <- 1 - PBB
    PBO <- 1 - POO
    PON <- 1 - PBN

    # number of cancers detected
    detect <- sensitivity.b * PBB * N.b + sensitivity.o * POO * N.o

    # number of patients who have at least one unnecessary confirmation tests
    biopsy <- sensitivity.b * POB * N.b + sensitivity.o * PBO * N.o + (1 - specificity) * N.n

    # number of lives saved
    live <- N * effect.b * death.b + N * effect.o * death.o

    result <- rbind(
      result,
      data.frame(
        prevalence = p.o,
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


# Anlaysis 2: pan-cancer test for breast cancer (indicated as b) + one of the other types of cancer (indicated as o)

setwd("~/Dropbox/Pan-cancer test/Figures")
seer <- read_excel("SEER CSR.xlsx")

result <- data.frame(prevalence = 1, specificity = 1, age = 1, detect = 1, biopsy = 1, live = 1, cancer = "b")

for (a in c(50, 60, 70)) {

  # create an input table
  seer_input <- seer %>%
    filter(age == a) %>%
    dplyr::select(cancer, prevalence, death)

  sensitivity <- c(0.74, 0.68, 0.78, 0.59, 0.64, 0.67)
  correct <- c(0.97, 0.72, 0.79, 0.92, 0.96, 0.96)
  input <- data.frame(seer_input, sensitivity = sensitivity, correct = correct)

  for (i in c(1:6)) {
    specificity <- 0.99

    p.b <- input$prevalence[input$cancer == "Breast"]
    p.o <- input$prevalence[i]
    p.n <- 1 - p.b - p.o

    N <- 1000
    N.b <- N * p.b
    N.o <- N * p.o
    N.n <- N * p.n

    death.b <- input$death[input$cancer == "Breast"]
    death.o <- input$death[i]

    sensitivity.b <- input$sensitivity[input$cancer == "Breast"]
    sensitivity.o <- input$sensitivity[i]

    PBB <- input$correct[input$cancer == "Breast"]
    POO <- input$correct[i]

    PBN <- p.b / (p.b + p.o)

    effect.b <- 0.2
    effect.o <- 0.2

    POB <- 1 - PBB
    PBO <- 1 - POO
    PON <- 1 - PBN

    detect <- sensitivity.b * PBB * N.b + sensitivity.o * POO * N.o
    biopsy <- sensitivity.b * POB * N.b + sensitivity.o * PBO * N.o + (1 - specificity) * N.n
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


# Analysis 3: pan-cancer test for breast cancer (indicated as b) + lung cancer (indicated as l) 
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
    p.b <- input$prevalence[input$cancer == "Breast"]
    p.l <- input$prevalence[input$cancer == "Lung and Bronchus"]
    p.c <- input$prevalence[i]
    p.n <- 1 - p.b - p.c - p.l

    N <- 1000
    N.b <- N * p.b
    N.l <- N * p.l
    N.c <- N * p.c
    N.n <- N * p.n

    death.b <- input$death[5]
    death.l <- input$death[4]
    death.c <- input$death[i]

    sensitivity.b <- input$sensitivity[5]
    sensitivity.l <- input$sensitivity[4]
    sensitivity.c <- input$sensitivity[i]

    specificity <- 0.99

    PBB <- input$correct[5]
    PLL <- input$correct[4]
    PCC <- input$correct[i]

    effect.b <- 0.2
    effect.c <- 0.2
    effect.l <- 0.2

    detect <-
      sensitivity.b * PBB * N.b +
      sensitivity.c * PCC * N.c +
      sensitivity.l * PLL * N.l

    biopsy <-
      sensitivity.b * (1 - PBB) * N.b +
      sensitivity.c * (1 - PCC) * N.c +
      sensitivity.l * (1 - PLL) * N.l +
      (1 - specificity) * (N.n)

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
  filter(cancer != "Lung and Bronchus" | cancer != "Breast") %>%
  mutate(ratio_bd = biopsy / detect, ratio_bl = biopsy / live) %>%
  dplyr::select(age, cancer, biopsy, detect, live, ratio_bd, ratio_bl)


# Analysis 4: pan-cancer test for breast cancer (indicated as b) + lung cancer (indicated as l) + ovary cancer (indicated as o) 
# + colon cancer (indicated as c) + liver cancer (indicated as li) + pancreatic cancer (indicated as p)

result <- data.frame(prevalence = 1, specificity = 1, detect = 1, biopsy = 1, live = 1, age = 1)

for (a in c(50, 60, 70)) {
  seer_input <- seer %>%
    filter(age == a) %>%
    dplyr::select(cancer, prevalence, death)

  sensitivity <- c(0.74, 0.68, 0.78, 0.59, 0.64, 0.67)
  correct <- c(0.97, 0.72, 0.79, 0.92, 0.96, 0.96)
  input <- data.frame(seer_input, sensitivity = sensitivity, correct = correct)

  p.b <- input$prevalence[5]
  p.c <- input$prevalence[1]
  p.li <- input$prevalence[2]
  p.p <- input$prevalence[3]
  p.l <- input$prevalence[4]
  p.o <- input$prevalence[6]

  p.n <- 1 - p.b - p.c - p.l - p.li - p.p - p.o
  N <- 1000
  N.b <- N * p.b
  N.l <- N * p.l
  N.c <- N * p.c
  N.li <- N * p.li
  N.p <- N * p.p
  N.o <- N * p.o

  N.n <- N * p.n

  death.b <- input$death[5]
  death.c <- input$death[1]
  death.l <- input$death[4]
  death.li <- input$death[2]
  death.p <- input$death[3]
  death.o <- input$death[6]

  sensitivity.b <- input$sensitivity[5]
  sensitivity.c <- input$sensitivity[1]
  sensitivity.l <- input$sensitivity[4]
  sensitivity.li <- input$sensitivity[2]
  sensitivity.p <- input$sensitivity[3]
  sensitivity.o <- input$sensitivity[6]

  specificity <- 0.99

  PBB <- input$correct[5]
  PCC <- input$correct[1]
  PLL <- input$correct[4]
  PLiLi <- input$correct[2]
  PPP <- input$correct[3]
  POO <- input$correct[6]

  effect.b <- 0.2
  effect.c <- 0.2
  effect.l <- 0.2
  effect.li <- 0.2
  effect.p <- 0.2
  effect.o <- 0.2

  detect.singleforone <-
    sensitivity.b * PBB * N.b +
    sensitivity.c * PCC * N.c +
    sensitivity.l * PLL * N.l +
    sensitivity.li * PLiLi * N.li +
    sensitivity.p * PPP * N.p +
    sensitivity.o * POO * N.o

  biopsy <-
    sensitivity.b * (1 - PBB) * N.b +
    sensitivity.c * (1 - PCC) * N.c +
    sensitivity.l * (1 - PLL) * N.l +
    sensitivity.li * (1 - PLiLi) * N.li +
    sensitivity.p * (1 - PPP) * N.p +
    sensitivity.o * (1 - POO) * N.o +
    (1 - specificity) * (N.n)

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
