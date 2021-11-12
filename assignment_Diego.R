############################################################################################################
# PHDS courses. Digital and Computational Demography 4
# C:/Statistiques/MPIDR/PHDS/Digital_Computational/4_Microsimulation_Alburez/Microsimulation.R

# Digital and Computational Demography. Day 4. Microsimulation. Diego Alburez

# Created by Liliana Calderon on 11-11-2021
# Last modified by Liliana Calderon on 12-11-2021

############################################################################################################

# You will need two different datasets, which you can find as .csv files in the folder “Assignment/Data”: 
# a. A subset of the Familinx database. 
# This particular version contains close to one million individual-level records representing the Swedish population. 
# b. The output of a SOCSIM microsimulation using historical Swedish mortality and fertility rates as input.  


# Cleaning the work space
rm(list=ls(all=T))

# Loading packages
library(tidyverse)  
library(viridis)

# Load the data
SWE_familinx <- read_csv("sweden_familinx.csv")
SWE_socsim <- read_csv("sweden_socsim.csv")

#----------------------------------------------------------------------------------------

#### Function to calculate the population alive at each calendar year ####
### Adapted from https://github.com/alburezg/socsim_example/blob/main/R/functions.R

# group_by_sex and group_by_age_sex are defined F by default. 
# If we want to get the counts only by sex or by sex and age, 
# we should set to true: group_by_sex = T or group_by_age_sex = T
# When using the function to calculate some indicators

## Pop Census

pop_census <- function(y, df, return_ids = F, group_by_sex = F, group_by_age_sex = F) {
  
  df$census <- y
  
  df$death_temp <- df$death
  
  year_df <- df %>% 
    dplyr::filter(birth <= y & death_temp >= y) %>% 
    dplyr::filter(!is.na(birth), !is.na(death))
  
  if(return_ids) {
    out <- year_df$profileid
  } else {
    if(group_by_sex) {
      out <- year_df %>% 
        dplyr::count(gender, census)
    } else if(group_by_age_sex) {
      out <- year_df %>% 
        dplyr::mutate(age_at_census = census - birth) %>% 
        dplyr::count(gender, age_at_census, census) %>% 
        tidyr::complete(gender, age_at_census, census, fill = list(n = 0))
    }
    else{
      out <- year_df
    }
    
  }
  
  return(out)
  
}


#-------------------------------------------------------------------------------------------------------------------------------
## Exercise 1: Annual sex ratios 

# In demography, sex ratios are traditionally defined as the proportion of men relative to women. 
# The annual sex ratio for in a given population in calendar year y is: 
# SR(y) number of men alive between time y and y+1 / number of women alive between time y and y+1


# a. Estimate period sex ratios for the 1850-1950 period using the Familinx and SOCSIM datasets.   

#### Function to calculate SR 

get_sr_period <- function(df, y_breaks, y_range){
  
  # Calculating the population who was ever alive during the year.
  # It includes all those who were alive at the start of a given year and did not die that year.
  # This is consistent with HMD's definition of "alive on January 1st". 
  
  opop2 <- df %>% 
    select(profileid, gender, birth = birth_year, death = death_year) %>% 
    filter(!is.na(gender), !is.na(birth), !is.na(death))
  
  yearly_pop_age_sex <- lapply(y_range, pop_census,
                               df = opop2, 
                               return_ids = F,
                               group_by_sex = T,
                               group_by_age_sex = F) 
  
  # 1. Numerator - number of men alive between time y and y+1 (by age)
  
  numerator <-
    data.frame(do.call(bind_rows, yearly_pop_age_sex)) %>%
    filter(gender == "male") %>%
    select(Year = census, everything()) %>%
    group_by(Year) %>%
    summarise(n_num = sum(n)) %>%
    arrange(Year)
  
  # 2. Denominator - number of women alive between time y and y+1 (by age)
  
  denominator <-
    data.frame(do.call(bind_rows, yearly_pop_age_sex)) %>%
    filter(gender == "female") %>%
    select(Year = census, everything()) %>%
    group_by(Year) %>%
    summarise(n_den = sum(n)) %>%
    arrange(Year)
  
  # 3. Rate (ASSR)
  
  sr <- numerator %>%
    full_join(denominator, by = c("Year")) %>%
    mutate(SR = n_num / n_den, 
           SR = ifelse(is.na(SR), 0, SR)) 
  
  return(sr)  
}


# Define some parameters

y_min <- 1850
y_max <- 1950
y_range <- y_min:y_max
y_breaks <- seq(y_min, y_max, 1)

## Calculate SR for each data base ####

Socsim_sr <- get_sr_period(df = SWE_socsim,
                           y_breaks = y_breaks, 
                           y_range = y_range) %>% 
              mutate(Source = "SOCSIM")

Familinx_sr <- get_sr_period(df = SWE_familinx,
                             y_breaks = y_breaks, 
                             y_range = y_range) %>% 
              mutate(Source = "Familinx")

# b. Plot the annual sex ratio over time for the two sources

SR <- Familinx_sr %>% bind_rows(Socsim_sr)

SR %>% 
  ggplot(aes(x = Year, y = SR, colour = Source)) +
  geom_line(size = 1)+
  theme_bw() +
  scale_color_viridis(option = "H", discrete=TRUE)

# ggsave(file="SR.jpeg", width=16, height=8, dpi=400)

#-------------------------------------------------------------------------------------------------------------------------------
## Exercise 2: Age-specific sex ratios
# We define age-specific sex ratios (ASSR) as the proportion of men aged ‘a’ relative to women aged ‘a’. T
# The ASSR for age ‘a’in a given population in calendar year y is
# ASSR(y) number of men in age range a to a+n alive between time y and y+1 
# by number of women in age range a to a+n alive between time y and y+1
# A value of 1 indicates an equivalent number of women and men in a given age group.a.

# a. Estimate age-specific sex ratios (ASSR) for the 1850-1950 period using the Familinx and SOCSIM datasets.  
# Use the following age groups: 0-1, 2-4, 5-29, 30-64, 65+

#### Function to calculate ASSR 

get_assr_period <- function(df, age_breaks, age_labels, y_breaks, y_range){
  
  # Calculating the population who was ever alive during the year.
  # It includes all those who were alive at the start of a given year and did not die that year.
  # This is consistent with HMD's definition of "alive on January 1st". 
  
  opop2 <- df %>% 
    select(profileid, gender, birth = birth_year, death = death_year) %>% 
    filter(!is.na(gender), !is.na(birth), !is.na(death))
  
  yearly_pop_age_sex <- lapply(y_range, pop_census,
                               df = opop2, 
                               return_ids = F,
                               group_by_sex = F,
                               group_by_age_sex = T) 
  
  # 1. Numerator - number of men alive between time y and y+1 (by age)
  
  numerator <-
    data.frame(do.call(bind_rows, yearly_pop_age_sex)) %>%
    filter(gender == "male") %>%
    select(Year = census, everything()) %>%
    mutate(Age = cut(age_at_census,breaks = age_breaks, labels = age_labels,include.lowest=TRUE,right=FALSE)) %>%
    group_by(Year, Age) %>%
    summarise(n_num = sum(n)) %>%
    arrange(Year, Age)
  
  # 2. Denominator - number of women alive between time y and y+1 (by age)
  
  denominator <-
    data.frame(do.call(bind_rows, yearly_pop_age_sex)) %>%
    filter(gender == "female") %>%
    select(Year = census, everything()) %>%
    mutate(Age = cut(age_at_census,breaks = age_breaks, labels = age_labels,include.lowest=TRUE,right=FALSE)) %>%
    group_by(Year, Age) %>%
    summarise(n_den = sum(n)) %>%
    arrange(Year, Age)
  
  # 3. Rate (ASSR)
  
  assr <- numerator %>%
    full_join(denominator, by = c("Year", "Age")) %>%
    mutate(ASSR = n_num / n_den, 
           ASSR = ifelse(is.na(ASSR), 0, ASSR)) 
  
  return(assr)  
}


# Define some parameters

y_min <- 1850
y_max <- 1950
y_range <- y_min:y_max
y_breaks <- seq(y_min, y_max, 1)

# Age groups for ASSR
age_breaks <- c(0, 2, 5, 30, 65, 100)
age_labels <- age_breaks[-length(age_breaks)]


## Calculate ASSR for each data base ####

Socsim_assr <- get_assr_period(df = SWE_socsim,
                               age_breaks = age_breaks,
                               age_labels = age_labels,
                               y_breaks = y_breaks, 
                               y_range = y_range) %>% 
                mutate(Source = "SOCSIM")

Familinx_assr <- get_assr_period(df = SWE_familinx,
                                 age_breaks = age_breaks,
                                 age_labels = age_labels,
                                 y_breaks = y_breaks, 
                                 y_range = y_range) %>% 
                mutate(Source = "Familinx")

## Check Missing values for Age in Familinx dataset

Familinx_assr %>% 
  filter(is.na(Age)) %>% 
  ungroup() %>% 
  summarise(men = sum(n_num), 
            women = sum(n_den))
## With 100 years as upper age bound, 
# 2681 men and 1984 women have no age at the census, 
# i.e. they are supposed to be older. 
# To allow comparability with the Socsim data (bounded at 100), 
# we will take out the data from this age group when plotting

# b.Plot the age-specific sex ratios over time for the two sources

ASSR <- Familinx_assr %>% 
  filter(!is.na(Age)) %>% 
  bind_rows(Socsim_assr)

ASSR %>% 
  ggplot(aes(x = Year, y = ASSR, colour = Age, group = Age)) +
  facet_wrap(~Source)+
  geom_line(size =1)+
  theme_bw()+
  scale_color_viridis(option = "H", discrete=TRUE)

# ggsave(file="ASSR.jpeg", width=16, height=8, dpi=400)
  
  
  
  