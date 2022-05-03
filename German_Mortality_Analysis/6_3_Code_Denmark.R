#### 0. Settings ####

packages <- c("tidyverse", 
              "demography", 
              "forecast", 
              "readr", 
              "textreadr", 
              "dplyr", 
              "ggplot2",
              "viridis")
suppressMessages(packages <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x)
    library(x, character.only = TRUE)
  }
}))

#---------------------------------------------------------------------------
#### 1. read & filter the data ####

DNK_bltper_1x1 <- read_table(
  file = "./6_2_data/bltper_1x1/DNK.bltper_1x1.txt",
  col_names = T,
  skip = 1
)

DNK_bltp <- as.tibble(DNK_bltper_1x1) %>%
  select(Year, Age, mx, qx, lx, dx, Lx, Tx, ex) %>%
  filter(Year >= 1950) %>%
  mutate(Age = as.numeric(Age)) %>%
  mutate(Age = replace_na(Age, 110))

str(DNK_bltp)
attach(DNK_bltp)

#----------------------------------------------------------------------------
#### 2. visualize metrics as functions of time / age ####

theme_custom <- function() {
 
  theme_bw() +
  
    theme(plot.title = element_text(
                       hjust = 0.5,
                       size = 18,
                       face = "bold",
                       vjust = 2),
          
          plot.subtitle = element_text(
                          hjust = 0.5,
                          size = 16,
                          face = "bold"),
          
          panel.grid.minor = element_line(
                             linetype = 2),
          
          axis.title.x = element_text(size = 12),
          
          axis.title.y = element_text(size = 12)
        )
}


logmx_b_viz <- DNK_bltp %>%
  ggplot(aes(Age, log(mx), group = Year)) +
  geom_line(aes(colour = Year), lwd = 0.5) +
  scale_colour_gradientn(colors = rainbow(20)) +
  scale_x_continuous(breaks = seq(min(Age), max(Age), 10)) +
  scale_y_continuous(breaks = seq(-12, 0, 1)) +
  labs( x = "Age (x)",
        y = expression("log" ~ m[x]),
        title = "Danish unisex death rates by age",
        subtitle = "Years 1950 - 2021",
        caption = "Source: HMD") +
 theme_custom()
