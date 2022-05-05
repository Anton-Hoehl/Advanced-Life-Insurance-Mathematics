#### 0. Settings ####

packages <- c("tidyverse", 
              "demography", 
              "forecast", 
              "readr", 
              "textreadr", 
              "dplyr", 
              "ggplot2")
suppressMessages(packages <- lapply(packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x)
    library(x, character.only = TRUE)
  }
}))

#---------------------------------------------------------------------------
#### 1. read & filter the data ####

DNK_bltper_1x1 <- 
  read_table(file = "./6_2_data/bltper_1x1/DNK.bltper_1x1.txt",
             col_names = T,
             skip = 1
            )

DNK_exp_1x1 <-
  read_table(file = "./6_2_data/Exposures_1x1/DNK.Exposures_1x1.txt",
             col_names = T,
             skip = 1
             )

DNK_bltp <- 
  as.tibble(DNK_bltper_1x1) %>%
  select(Year, Age, mx, qx, lx, dx, Lx, Tx, ex) %>%
  filter(Year >= 1950) %>%
  mutate(Age = as.numeric(Age)) %>%
  mutate(Age = replace_na(Age, 110))

DNK_texp <-
  as.tibble(DNK_exp_1x1) %>%
  select(Year, Age, Total) %>%
  filter(Year >= 1950) %>%
  mutate(Age = as.numeric(Age)) %>%
  mutate(Age = replace_na(Age, 110))

names(DNK_texp) <- c("Year_exp","Age_exp","Total_exp")

str(DNK_bltp)
str(DNK_texp)
attach(DNK_bltp)
attach(DNK_texp)
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
                  
DNK_metrics_summaries <- as.tibble( 
  data.frame(Year_int = 
             as.factor(c("1950-1961",
                          "1970-1981",
                          "1990-2001",
                          "2010-2021"))))

j = -1
repeat {    
  j = j + 1
    for (i in 0:3) {
        year_start <- 1950
        year_end <- 1961
        x <- vector(mode = "numeric", length = 4)
        x[i + 1] <- DNK_bltp %>%
        filter(Age == j, Year >= year_start + 20*i & Year <= year_end + 20*i) %>%
        summarize(mean_mx = mean(mx)) %>%
        pull() }
    next
      
  if (j == 110) {
      break
    } else {
      x <- as.matrix(x)
      names(x) <- paste("mean_m",j, sep = "")
      DNK_metrics_summaries <- DNK_metrics_summaries %>% bind_cols(x)
    }
    
  }



logmx_b_viz <- DNK_bltp %>%
            
               ggplot(aes(Age, log(mx), group = Year)) +
            
               geom_line(aes(colour = Year), lwd = 0.5) +
            
               scale_colour_gradientn(colors = rainbow(20)) +
               scale_x_continuous(breaks = seq(min(Age), max(Age), 10)) +
               scale_y_continuous(breaks = seq(-12, 0, 1)) +
            
               labs(x = "Age (x)",
                    y = expression("log" ~ m[x]),
                    title = "Danish unisex death rates",
                    subtitle = "Years 1950 - 2021",
                    caption = "Source: HMD") +
            
               theme_custom()

qx_b_viz <- DNK_bltp %>%
  
               ggplot(aes(Age, qx, group = Year)) +
  
               geom_line(aes(colour = Year), lwd = 0.5) +
               
               scale_colour_gradientn(colors = rainbow(20)) +
               scale_x_continuous(breaks = seq(min(Age), max(Age), 10)) +
                
               labs(x = "Age (x)",
                     y = expression(q[x]),
                     title = "Danish unisex death probabilities",
                     subtitle = "Years 1950 - 2021",
                     caption = "Source: HMD") +
                
               theme_custom()

ex0_b_viz <- DNK_bltp %>%
             filter(Age == 0) %>%
            
            ggplot(aes(Year, ex)) +
            
            geom_line() +
  
            scale_x_continuous(breaks = seq(min(Year), max(Year), 10)) +
            scale_y_continuous(breaks = seq(70, 85, 5)) +
            
            labs(x = "Time(t)",
                 y = expression( e[0]),
                 title = "Evolution of unisex life expectancies at birth in Denmark",
                 caption = "Source: HMD") +
            
            theme_custom()