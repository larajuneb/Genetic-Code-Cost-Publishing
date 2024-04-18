library(tidyverse)
library(hrbrthemes)
library(plotly)
library(viridis)


df <- read_csv(
  "matrices/SeqPredNN_cost_matrix_NORM.csv"
               )

df

df <- df %>%
  mutate(
    "Initial Codon"=`...1`,
    .keep = "unused",
    .before = 1
  ) %>%
  pivot_longer(
    cols = !`Initial Codon`,
    names_to = "Final Codon",
    values_to = "Cost"
  ) %>%
  filter(
    Cost != "-" &  Cost != "*" 
    ) %>%
  mutate(
    Cost = as.double(
      Cost
      )
    ) %>%
  separate_wider_position(
    cols = `Initial Codon`,
    widths = c("Initial Position 1"=1, 
               "Initial Position 2"= 1,
               "Initial Position 3" =1
               )
  ) %>%
  separate_wider_position(
    cols = `Final Codon`,
    widths = c("Final Position 1"=1, 
               "Final Position 2"= 1,
               "Final Position 3"=1
    )
  ) %>%
  mutate(
    Position = case_when(
      `Initial Position 1`!=`Final Position 1` ~ 1,
      `Initial Position 2`!=`Final Position 2` ~ 2,
      `Initial Position 3`!=`Final Position 3` ~ 3,
      ),
    .before = 1
    ) %>%
  mutate(
    "Initial Codon" = paste(
      `Initial Position 1`, 
      `Initial Position 2`, 
      `Initial Position 3`,
      sep = ""
      )
    ) %>%
  
  mutate(
    "Initial" = paste(
      "Initial Position", 
      Position
      ),
    "Final" = paste(
      "Final Position", 
      Position
      ),
    .before = 1
    ) %>%
  rowwise(
    ) %>%
  mutate(
    "Initial" = get(
      Initial),
    "Final" = get(
      Final
      )
    ) %>%
  ungroup(
    )

df

plot <- df %>%
  ggplot(
    aes(
      x = Final,
      y = `Initial Codon`,
      fill = Cost
      )
    ) +
  geom_raster(
    )+
  theme_ipsum()+
  scale_fill_viridis()+
  scale_x_discrete(
    limits=c(
      "U", "G", "C", "A"
      )
    )+
  facet_grid(
    factor(
      `Initial Position 1`,
      levels=c(
        "U", "G", "C", "A"
        )
      ) ~`Position`, 
    scales="free_y"
    )
plot
  #facet_grid(
  #  cols = vars(
  #    `Position`
  #  ), 
  #)+
  #coord_equal()
#plot
