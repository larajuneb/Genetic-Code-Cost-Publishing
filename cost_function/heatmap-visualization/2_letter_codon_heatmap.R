library(tidyverse)
library(hrbrthemes)
library(viridis)


df <- read_csv(
  "matrices/Higgs_primordial_cost_matrix_NORM.csv"
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
    )%>%
  separate_wider_position(
    cols = `Initial Codon`,
    widths = c("Initial Position 1"=1, 
               "Initial Position 2"= 1
               )
  ) %>%
  separate_wider_position(
    cols = `Final Codon`,
    widths = c("Final Position 1"=1, 
               "Final Position 2"= 1
    )
  ) %>%
  mutate(
    Position = case_when(
      `Initial Position 1`!=`Final Position 1` ~ 1,
      `Initial Position 2`!=`Final Position 2` ~ 2
      ),
    .before = 1
    ) %>%

  mutate(
    "Initial Codon" = paste(
      `Initial Position 1`, 
      `Initial Position 2`,
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
    ) %>%
  mutate(
    `Final` = factor(
      `Final`,
      levels=c(
        "U", "G", "C", "A"
        )
      ),
    ) 

  
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
  facet_grid(
    cols = vars(
      `Position`
      ), 
    )+
  coord_equal()


plot

