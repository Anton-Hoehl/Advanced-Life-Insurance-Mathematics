#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
# Advanced Life insurance mathematics: #
#            Assignment 2              #
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

#### 0. Settings ####
library(tidyverse)
library(demography)
library(forecast)
library(readr)
library(textreadr)
library(dplyr)

## Read the data
## -----------------------------------------------------------------------------

# read table seems to work better for me than read_delim
german_total_1x1 <- read_table(
  file = "./6_2_data/Germany_Life_Table_Total.txt",
  col_types = "nnnnnnnnnn",
  col_names = FALSE,
  skip = 1
)

names(german_total_1x1) <-
  c("Year", "Age", "mx", "qx", "ax", "lx", "dx", "Lx", "Tx", "ex")


## Descriptives
## -----------------------------------------------------------------------------
min(german_total_1x1$Year)
max(german_total_1x1$Year)

min(german_total_1x1$Age)
max(german_total_1x1$Age)

# The range of the available Years is small enough that we dont need to cut it
## -----------------------------------------------------------------------------

german_total <-
  german_total_1x1 %>% select(Year, Age, mx, qx, lx, dx, ex)

# Hello this is a tescomment