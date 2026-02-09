---
title: "BIOS 6624 Project 0:Cortisol and DHEA "
author: "Moj"
date: "2026-01-26"
header-includes: \usepackage{multirow}
output:
  html_document:
    df_print: paged
  pdf_document: default
urlcolor: blue
---
#Read data

```{r}
# Read data

library(dplyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(ggeffects)
library(binom)

cortisol <- read.csv("C:/Users/mojde/OneDrive/Desktop/Colorado SPH/PhD Courses/Spring 2026/BIOS 6624/Project 0/Project0_Clean_v2.csv", na.strings = c("", " ", "NA"))
```

## Explore Variables

# 
```{r}
colnames(cortisol)

```

```{r}
colSums(is.na(cortisol))

```
#Q1: Time difference in minutes
#Fix diary wake time to have values for each cell


```{r}
library(dplyr)
library(tidyr)

cortisol <- cortisol %>%
  arrange(SubjectID, `Collection.Date`, `Collection.Sample`) %>%  
  group_by(SubjectID, `Collection.Date`) %>%                    
  fill(`Sleep.Diary.reported.wake.time`, .direction = "down") %>% # fill down
  ungroup()
```


#prepare time
```{r}
head(cortisol$Booket..Clock.Time)
head(cortisol$MEMs..Clock.Time)
head(cortisol$Sleep.Diary.reported.wake.time)

```




```{r}
cortisol <- cortisol %>%
  mutate(
    wake_time_posix = as.POSIXct(
      Sleep.Diary.reported.wake.time,
      format = "%H:%M",
      tz = "UTC"
    ),
    booklet_time_posix = as.POSIXct(
      Booket..Clock.Time,
      format = "%H:%M",
      tz = "UTC"
    ),
    mems_time_posix = as.POSIXct(
      MEMs..Clock.Time,
      format = "%H:%M",
      tz = "UTC"
    )
  )
```


#difference between wake time and booklet time

```{r}
cortisol <- cortisol %>%
  mutate(
    booklet_minutes_since_wake =
      as.numeric(difftime(
        booklet_time_posix,
        wake_time_posix,
        units = "mins"
      )),

    mems_minutes_since_wake =
      as.numeric(difftime(
        mems_time_posix,
        wake_time_posix,
        units = "mins"
      ))
  )

```


```{r}
class(cortisol$booklet_minutes_since_wake)
# numeric

summary(cortisol$booklet_minutes_since_wake)
# values are in minutes
```


```{r}
library(ggplot2)

ggplot(cortisol, 
       aes(x = mems_minutes_since_wake, y = booklet_minutes_since_wake)) +
  geom_point(color = "blue", size = 2, alpha = 0.6) + 
  geom_smooth(method = "lm", color = "black", se = FALSE) +  # regression line
  labs(
    x = "MEMs minutes since wake",
    y = "Booklet minutes since wake",
    title = "Figure 1. Comparison of MEMs and Booklet Sample Timing",
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(face = "italic")
  )


```



```{r}
agreement_model <- lmer(
  booklet_minutes_since_wake ~ mems_minutes_since_wake + (1 | SubjectID),
  data = cortisol
)

summary(agreement_model)

```


#Q2: Adherence to +30 min and +10 hr protocol times
```{r}
library(dplyr)

# adherence variables using booklet times
cortisol <- cortisol %>%
  mutate(
    target_time = case_when(
      Collection.Sample == 2 ~ 30,    # 30 minutes
      Collection.Sample == 4 ~ 600,   # 600 minutes (10 hours)
      TRUE ~ NA_real_
    ),
    # deviation from target in minutes
    deviation = booklet_minutes_since_wake - target_time,
    # Mark adherent if within ±7.5 minutes
    adherent1 = !is.na(deviation) & abs(deviation) <= 7.5,
    # Mark adherent if within ±15 minutes
    adherent2 = !is.na(deviation) & abs(deviation) <= 15
  )

# Summarise adherence by sample
adherence_summary <- cortisol %>%
  filter(Collection.Sample %in% c(2, 4)) %>%
  group_by(Collection.Sample) %>%
  summarise(
    n = n(),
    adherent_count1 = sum(adherent1, na.rm = TRUE),
    proportion_adherent1 = adherent_count1 / n,
    adherent_count2 = sum(adherent2, na.rm = TRUE),
    proportion_adherent2 = adherent_count2 / n,
    .groups = "drop"
  )

adherence_summary
```


#Q3. 

#flag outliers based on PI information
#Cutoffs (nmol/L):
#Cortisol > 80
#DHEA > 5.205

```{r}
table(cortisol$SubjectID[cortisol$DHEA..nmol.L. >= 5.205])
```

```{r}
cortisol_clean <- cortisol %>%
  filter(Cortisol..nmol.L. <= 80,
         DHEA..nmol.L. <= 5.205)
```


#for participants who have multiple DHEA values at the detection limit (5.205)



#Basic hormone summaries

```{r}
summary(cortisol_clean$Cortisol..nmol.L.)
summary(cortisol_clean$DHEA..nmol.L.)

```

#cortisol boxplot

```{r}
ggplot(cortisol_clean,
       aes(x = factor(Collection.Sample),
           y = Cortisol..nmol.L.)) +
  geom_boxplot() +
  labs(x = "Sample number", y = "Cortisol (nmol/L)")
```


#DHEA boxplot

```{r}
ggplot(cortisol_clean,
       aes(x = factor(Collection.Sample),
           y = DHEA..nmol.L.)) +
  geom_boxplot() +
  labs(x = "Sample number", y = "DHEA (nmol/L)")

```


#log-transform hormones
```{r}
cortisol_clean <- cortisol_clean %>%
  mutate(
    log_cortisol = log(Cortisol..nmol.L.),
    log_dhea = log(DHEA..nmol.L.)
  )

```



#Log-scale plots

```{r}
ggplot(cortisol_clean,
       aes(x = factor(Collection.Sample),
           y = log_cortisol)) +
  geom_boxplot() +
  labs(x = "Sample number", y = "log Cortisol")

ggplot(cortisol_clean,
       aes(x = factor(Collection.Sample),
           y = log_dhea)) +
  geom_boxplot() +
  labs(x = "Sample number", y = "log DHEA")
```


```{r}
summary(cortisol_clean$booklet_minutes_since_wake)
```




```{r}
cortisol_clean <- cortisol_clean %>%
  mutate(
    time_min = booklet_minutes_since_wake, #renaming variable
    after_30 = ifelse(time_min > 30, time_min - 30, 0) #subtract 30 so that the second slope measures change after the breakpoint, not change from the start.
  )
```

#piecewise linear mixed-effects model for cortisol

```{r}
cortisol_model <- lmer(
  log_cortisol ~ time_min + after_30 + (1 | SubjectID),
  data = cortisol_clean
)

summary(cortisol_model)

```
#Piecewise regression for DHEA

```{r}
dhea_model <- lmer(
  log_dhea ~ time_min + after_30 + (1 | SubjectID),
  data = cortisol_clean
)

summary(dhea_model)
```

```{r}

library(dplyr)
library(tidyr)

spaghetti_data <- cortisol_clean %>%
  select(SubjectID, time_min, log_cortisol, log_dhea) %>%
  pivot_longer(
    cols = c(log_cortisol, log_dhea),
    names_to = "Hormone",
    values_to = "Value"
  ) %>%
  mutate(
    Hormone = recode(Hormone,
                     log_cortisol = "Cortisol",
                     log_dhea = "DHEA")
  )
```


```{r}
ggplot(spaghetti_data,
       aes(x = time_min,
           y = Value,
           group = SubjectID,
           color = Hormone)) +
  
  geom_line(alpha = 0.25, linewidth = 0.4) +
  
  facet_wrap(
    ~ Hormone,
    scales = "free_y",
    labeller = as_labeller(c(
      Cortisol = "log Cortisol",
      DHEA = "log DHEA"
    ))
  ) +
  
  labs(
    x = "Minutes since waking",
    y = "Log hormone level",
    title = "Figure 2. Spaghetti plot of individual cortisol and DHEA trajectories"
  ) +
  
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 12)
  )

```


R.version.string

This repository contains analysis code and information for BIOS6624.