# --------------------------------------------------------------------------------------- #
# - FILE NAME:   Figures.R         
# - DATE:        20/07/2021
# - DESCRIPTION: Represent the results. 
# - AUTHORS:     Pol Capdevila Lanzaco (pcapdevila.pc@gmail.com)
# --------------------------------------------------------------------------------------- #

rm(list=ls(all=TRUE)) #remove everything

# Libraries

library(tidyverse)
library(dplyr)
library(tidybayes)
library(bayestestR)
library(ggthemes)
library(cowplot)
library(magrittr)
library(purrr)
library(forcats)
library(tidyr)
library(modelr)
library(ggdist)
library(rstan)
library(ggrepel)
library(ggsci)
library(wesanderson)
library(RColorBrewer)
library(posterior)
library(brms)
library(brmstools)
library(rphylopic)
library(RCurl)
library(png)
library(patchwork)

# Set default ggplot theme

theme_set(theme_minimal()+
            theme(axis.title.x = element_text(size=15, margin = margin(t = 10, r = 0, b = 0, l = 0)), 
                  axis.title.y = element_text(size=15, margin = margin(t = 0, r = 10, b = 0, l = 0)),
                  axis.line.x = element_line(color="black", size = 0.5),
                  axis.line.y = element_line(color="black", size = 0.5),
                  panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  axis.text.x = element_text(color="black", size = 12),
                  axis.text.y = element_text(color="black", size = 12),
                  strip.text.x = element_text(size = 12),
                  axis.ticks = element_line(color="black"),
                  plot.margin = margin(0.5,0.5,0.5,0.5, unit = "cm")))

# Set working directory

path <- gsub("Code", "", dirname(rstudioapi::getActiveDocumentContext()$path))

CodePath <- paste0(path,"Code")
DataPath <-  paste0(path,"Data")
ResultPath <-  paste0(path, "Results") 

# Load data 

load(paste0(ResultPath, "/Models.RData"))

# Figure 1: Global map #########################################################

# First we summarise the Number of threats by country #

lpi_map <- lpi %>% 
  group_by(Country) %>% 
  summarise(mean_threats=mean(n.threat, na.rm=T)) %>% 
  rename(region = Country)  # Rename columns

# Now by each 5 of latitude

lpi_lat <- lpi %>% 
  mutate(Lat=cut(Latitude, 
                 breaks = seq(-80,85, by=5), include.lowest = T, 
                 labels = as.character(seq(-80,80, by=5)))) %>% 
  group_by(Lat) %>% 
  summarise(n=n(), 
            mean_threats=mean(n.threat, na.rm=T),
            se_threats=plotrix::std.error(n.threat, na.rm = T))

# World map data

world_map <- map_data("world")

# Data on countries 

lpi_map <- lpi_map %>%
  mutate(region = ifelse(region=="United Kingdom", "UK",region),
         region = ifelse(region=="French Southern Territories", 
                         "French Southern and Antarctic Lands", region),
         region = ifelse(region=="Antigua and Barbuda", "Antigua", region),
         region = ifelse(region=="Bahamas, The", "Bahamas", region),
         region = ifelse(region=="Cote d'Ivoire", "Ivory Coast", region),
         region = ifelse(region=="Congo, Dem. Rep.", 
                         "Democratic Republic of the Congo", region),
         region = ifelse(region=="Congo, Rep.", "Republic of Congo", region),
         region = ifelse(region=="Cabo Verde", "Cape Verde", region),
         region = ifelse(region=="Egypt, Arab Rep.", "Egypt", region),
         region = ifelse(region=="Falkland Islands (Malvinas)", 
                         "Falkland Islands", region),
         region = ifelse(region=="Heard Island And McDonald Islands", 
                         "Heard Island", region),
         region = ifelse(region=="Iran, Islamic Rep.", "Iran", region),
         region = ifelse(region=="St. Kitts and Nevis", "Nevis", region),
         region = ifelse(region=="Korea, Rep.", "South Korea", region),
         region = ifelse(region=="Lao PDR", "Laos", region),
         region = ifelse(region=="St. Lucia", "Saint Lucia", region),
         region = ifelse(region=="North Macedoni", "Macedonia", region),
         region = ifelse(region=="Russian Federation", "Russia", region),
         region = ifelse(region=="South Georgia And The South Sandwich Islands", 
                         "South Georgia", region),
         region = ifelse(region=="Saint Helena, Ascension And Tristan Da Cunha", 
                         "Saint Helena", region),
         region = ifelse(region=="Slovak Republic", "Slovakia", region),
         region = ifelse(region=="Trinidad and Tobago", "Trinidad", region),
         region = ifelse(region=="Taiwan, Province Of China", "Taiwan", region),
         region = ifelse(region=="United States", "USA", region),
         region = ifelse(region=="Venezuela, RB", "Venezuela", region),
         region = ifelse(region=="Virgin Islands (U.S.)", "Virgin Islands", region)) 
#Change name of some countries to match the ones in world map
  

# Join with global data 

lpi_map <- world_map %>% left_join(lpi_map, by = "region")

# Plot it

(g1 <- ggplot(lpi_map, aes(long, lat, group = group))+
  geom_polygon(aes(fill = mean_threats), 
               color = "black")+
  scale_fill_viridis_c("Number of threats", 
                       option = "A",direction = -1) +
  labs(y="Latitude", x="Longitude") +
    scale_y_continuous(expand = expansion(mult = c(0, 0)))+
    scale_x_continuous(expand = expansion(mult = c(0, 0)))+
  theme_bw() +
  theme(panel.border = element_blank(),
        legend.position = "none"))

(g2 <- lpi_lat %>% ggplot(aes(y=mean_threats, x=Lat, 
                              fill=mean_threats))+
    geom_bar(stat = "Identity") + 
    geom_errorbar(aes(ymin= mean_threats-se_threats, 
                      ymax = mean_threats+se_threats), width = 0.5)+
    geom_text(aes(y = 2.7, 
                  label=paste("n=", n)),
              hjust=0.8)+
    scale_fill_viridis_c(option = "A",
                         direction = -1, breaks=c(0,1,2,3),
                         limits=c(0,3)) +
    scale_y_continuous(limits=c(0,3), 
                       expand = expansion(mult = c(0, 0)))+
    scale_x_discrete(breaks=c(80,50,0,-50,-80))+
    labs(x="", y="Number of threats")+
    coord_flip()+
    theme(legend.position = "none"))

legend <- get_legend(g1+theme(legend.position = "bottom",
                              legend.direction = "horizontal"))

row1 <- plot_grid(g1, g2,
          align = "hv",
          labels = c('a', 'b'),
          label_size = 12, rel_widths = c(1,0.5))

(figure1 <- plot_grid(row1, legend, 
                     ncol = 1, 
                     rel_heights = c(1,0.1)))

ggsave("Figure1.pdf", figure1,
       width = 16, height = 8,
       path = ResultPath)


# Figure 2: Missing values ##################################################### 

# Panel a: Proportion of missing values ---------------------------------------- 

(g2 <- lpi %>% ungroup() %>%
    distinct(ID, .keep_all=T) %>% 
    select(bm_g,pop_dens,
           Habitat_breadth_IUCN,
           Trophic_level, 
           Latitude) %>%  # replace to your needs
    rename(BodyMass=bm_g,
           HabitatBreadth=Habitat_breadth_IUCN,
           TrophicLevel=Trophic_level,
           PopulationDensity=pop_dens) %>% 
    gather(key = "key", value = "val") %>%
    mutate(val = na_if(val, NaN),
           isna = is.na(val),
           key=factor(key,levels=c("Latitude", 
                                   "BodyMass","PopulationDensity", 
                                   "HabitatBreadth", "TrophicLevel"))) %>%
    group_by(key) %>%
    mutate(total = n()) %>%
    group_by(key, total, isna) %>%
    summarise(num.isna = n()) %>%
    mutate(pct = num.isna / total * 100) %>% 
    ggplot() +
    geom_bar(aes(x = key, 
                 y = pct, fill=isna), 
             stat = 'identity', alpha=0.8) +
    scale_fill_manual(name = "", 
                      values = c('steelblue', 'tomato3'), 
                      labels = c("Present", "Missing")) +
    coord_flip() +
    theme(legend.position = "none")+
    scale_x_discrete(expand = c(0,0))+
    scale_y_continuous(expand = c(0,0))+
    labs(x ='', y = "% of missing values"))

# Panel b: Raster with combinations of missigness ------------------------------ 

(g3 <-  lpi %>%ungroup() %>% 
    distinct(ID, .keep_all=T) %>% 
    mutate(id=row_number()) %>% 
    select(id, bm_g,pop_dens,
           Habitat_breadth_IUCN,
           Trophic_level, Latitude) %>% 
    rename(BodyMass=bm_g,
           HabitatBreadth=Habitat_breadth_IUCN,
           TrophicLevel=Trophic_level,
           PopulationDensity=pop_dens) %>%
    gather(-id, key = "key", value = "val") %>%
    mutate(val = na_if(val, NaN),
           isna = is.na(val),
           key=factor(key,levels=c("Latitude",
                                   "BodyMass", "PopulationDensity", 
                                   "HabitatBreadth", "TrophicLevel"))) %>% 
    ggplot(aes(key, reorder(id, desc(isna)), fill = isna)) +
    geom_raster(alpha=0.8) +
    scale_fill_manual(name = "",
                      values = c('steelblue', 'tomato3'),
                      labels = c("Present", "Missing")) +
    labs(x = "",
         y = "Population") +
    coord_flip()+
    scale_x_discrete(expand = c(0,0))+
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()))

# Combine figures 

legend <- get_legend(g3+theme(legend.text = element_text(size = 14)))

(figure2 <- plot_grid(g2+theme(plot.margin = margin(5,0,0,0)), 
                      g3+theme(legend.position = "none", 
                               plot.margin = margin(5,0,0,0)), 
                      nrow = 1, align = "h", axis = "b", 
                      labels = "auto"))

(figure2 <- plot_grid(figure2, legend, rel_widths = c(1, 0.1)))

ggsave(figure2,filename = "Fig2.pdf", 
       width = 12, height = 4,
       path = ResultPath)

# Figure 3: Factors ##################################################
# Panel a: System -----------------------------------------------------

# Sample size

sample_size <- lpi %>% 
  distinct(ID, .keep_all=T) %>% 
  group_by(System) %>% 
  summarise(n=n())

med <- sys_1 %>%
  gather_draws(`b_.*`, regex = TRUE) %>%
  mutate(System = gsub("b_System", "", .variable), 
         System = gsub('[[:digit:]]+', '', System)) %>% 
  group_by(System) %>%
  summarise(med=quantile(.value,probs = 0.99))

sample_size<- sample_size %>% 
  left_join(med)

# Plot intervals 

dat <- sys_1 %>%
  gather_draws(`b_.*`, regex = TRUE)%>%
  mutate(System = gsub("b_System", "", .variable), 
         System = gsub('[[:digit:]]+', '', System)) 

(g3a <- sys_1 %>%
    gather_draws(`b_.*`, regex = TRUE)%>%
    median_qi(.width = c(.95, .8, .5)) %>%
    mutate(System = gsub("b_System", "", .variable), 
           System = gsub('[[:digit:]]+', '', System)) %>% 
    ggplot(aes(y = reorder(System,-.value), x = .value)) +
    stat_slab(data=dat, alpha=0.2, aes(fill=System)) +
    scale_fill_manual(values=c("#A1D6E2","#336B87", "#CC954E"))+
    geom_pointinterval(aes(xmin = .lower, xmax = .upper),
                       interval_size_range = c(0.5, 2), 
                       colour="grey40") +
    geom_text(data = sample_size, aes(x=med, y=System, 
                                      label = paste0("n=", n)),
              vjust   = -1.5, hjust=-.1)+
    scale_x_continuous(breaks = c(0,1,2),                        limits = c(0,2.5))+
    labs(x="Number of threats", y = "System") +
    theme(legend.position = "none"))

# Panel b: Class -----------------------------------------------------

# Sample size

sample_size <- lpi %>% 
  distinct(ID, .keep_all=T) %>% 
  group_by(Class) %>% 
  summarise(n=n())

med <- clas_1 %>%
  gather_draws(`b_.*`, regex = TRUE) %>%
  mutate(Class = gsub("b_Class", "", .variable), 
         Class = gsub('[[:digit:]]+', '', Class),
         Class = gsub("CartilaginousFish", "Cartilaginous Fish", Class),
         Class = gsub("BonyFish", "Bony Fish", Class)) %>% 
  group_by(Class) %>%
  summarise(med=quantile(.value,probs = 0.99))

sample_size<- sample_size %>% 
  left_join(med)

# Plot intervals 

dat <- clas_1 %>%
  gather_draws(`b_.*`, regex = TRUE)%>%
  mutate(Class = gsub("b_Class", "", .variable), 
         Class = gsub('[[:digit:]]+', '', Class),
         Class = gsub("CartilaginousFish", "Cartilaginous Fish", Class),
         Class = gsub("BonyFish", "Bony Fish", Class))   
(g3b <- clas_1 %>%
    gather_draws(`b_.*`, regex = TRUE)%>%
    median_qi(.width = c(.95, .8, .5)) %>%
    mutate(Class = gsub("b_Class", "", .variable), 
           Class = gsub('[[:digit:]]+', '', Class),
           Class = gsub("CartilaginousFish", "Cartilaginous Fish", Class),
           Class = gsub("BonyFish", "Bony Fish", Class)) %>% 
    ggplot(aes(y = reorder(Class,-.value), x = .value)) +
    stat_slab(data=dat, alpha=0.2, aes(fill=Class)) +
    scale_fill_jco(name="Class")+
    geom_pointinterval(aes(xmin = .lower, xmax = .upper),
                       interval_size_range = c(0.5, 2), 
                       colour="grey40") +
    geom_text(data = sample_size, aes(x=med, y=Class, 
                                      label = paste0("n=", n)),
              vjust   = -1.5, hjust=-.1)+
    scale_x_continuous(breaks = c(0,1,2),                        limits = c(0,2.5))+
    labs(x="Number of threats", y = "Class") +
    theme(legend.position = "none"))


# Panel c: Trophic level -------------------------------------------------------

sample_size <- lpi %>% 
  distinct(ID, .keep_all=T) %>% 
  group_by(Trophic_level) %>% 
  summarise(n=n())

med <- tl_1 %>%
  gather_draws(`b_.*`, regex = TRUE) %>%
  mutate(Trophic_level = gsub("b_Trophic_level", "", .variable)) %>%
  group_by(Trophic_level) %>%
  summarise(med=quantile(.value,probs = 0.99))

sample_size<- sample_size %>% 
  left_join(med) %>% 
  drop_na(med)

# Plot intervals 

dat <- tl_1 %>%
  gather_draws(`b_.*`, regex = TRUE) %>%
  mutate(Trophic_level = gsub("b_Trophic_level", "", .variable)) 

(g3c <- tl_1 %>%
    gather_draws(`b_.*`, regex = TRUE) %>%
    median_qi(.width = c(.95, .8, .5)) %>%
    mutate(Trophic_level = gsub("b_Trophic_level", "", .variable))  %>%
    ggplot(aes(y = reorder(Trophic_level, -.value), x = .value)) +
    stat_slab(data=dat, alpha=0.5, aes(fill=.variable)) +
    geom_pointinterval(aes(xmin = .lower, xmax = .upper),
                       interval_size_range = c(0.5, 2), 
                       colour="grey40") +
    geom_text(data = sample_size, aes(x=med, y=Trophic_level, 
                                      label = paste0("n=", n)),
              vjust   = -1.5, hjust=-.1)+
    scale_fill_manual(values = wes_palette("Royal1"))+
    labs(x="Number of threats", y = "Trophic level") +
    scale_x_continuous(breaks = c(0,1,2),                        limits = c(0,2.5))+
    theme(legend.position = "none"))

# Panel d: Body mass ------------------------------------------------------

bm_seq <- tibble(bm_g= seq(min(log(lpi$bm_g+1), na.rm = T),
                           max(log(lpi$bm_g+1), na.rm = T),
                           length.out = 100))

dat <- fitted(bm_1, 
              newdata = bm_seq,
              re_formula =NA,
              summary = F) %>% 
  data.frame() %>%
  pivot_longer(everything()) %>% 
  bind_cols(expand(bm_seq, iter=1:7200, 
                   nesting(bm_g))) 

mean_dat <- dat %>% 
  group_by(bm_g) %>% 
  summarise(Estimate=mean(value))

(g3d <- dat %>% 
    filter(iter%in%c(1:250)) %>% 
    ggplot(aes(x = log(bm_g+1), 
               y = value, group=iter)) +
    geom_line(alpha=0.05, color="#216894")+
    geom_line(data = mean_dat, aes(y=Estimate, 
                                   x=log(bm_g+1)),
              colour="#216894", size=1)+
        scale_y_continuous(breaks = c(0,1,2),                        limits = c(0,2.5))+
    labs(x = "Log(Body mass)",
         y = "Number of threats")) 

# Panel e: Latitude ------------------------------------------------------

lat_seq <- tibble(Latitude= seq(min(abs(lpi$Latitude), 
                                    na.rm = T),
                                max(abs(lpi$Latitude), 
                                    na.rm = T),
                                length.out = 100))

dat <- fitted(lat_1, 
              newdata = lat_seq,
              re_formula =NA,
              summary = F) %>% 
  data.frame() %>%
  pivot_longer(everything()) %>% 
  bind_cols(expand(lat_seq, iter=1:7200, 
                   nesting(Latitude))) 

mean_dat <- dat %>% 
  group_by(Latitude) %>% 
  summarise(Estimate=mean(value))

(g3e <- dat %>% 
    filter(iter%in%c(1:250)) %>% 
    ggplot(aes(x = Latitude, 
               y = value, group=iter)) +
    geom_line(alpha=0.05, color="#949228")+
    geom_line(data = mean_dat, aes(y=Estimate, 
                                   x=Latitude),
              colour="#949228", size=1)+
    scale_y_continuous(breaks = c(0,1,2),
                       limits = c(0,2.5))+
    labs(x = "Latitude",
         y = "Number of threats")) 

# Panel f: Population density --------------------------------------------------

pd_seq <- tibble(pop_dens= seq(min(lpi$pop_dens, na.rm = T),
                               max(lpi$pop_dens, na.rm = T),
                               length.out = 100))

dat <- fitted(pop_dens_1, 
              newdata = pd_seq,
              re_formula =NA,
              summary = F) %>% 
  data.frame() %>%
  pivot_longer(everything()) %>% 
  bind_cols(expand(pd_seq, iter=1:7200, nesting(pop_dens))) 

mean_dat <- dat %>% 
  group_by(pop_dens) %>% 
  summarise(Estimate=mean(value))

(g3f <- dat %>% 
    filter(iter%in%c(1:250)) %>% 
    ggplot(aes(x = pop_dens, y = value, group=iter)) +
    geom_line(alpha=0.05, color="#0396A6")+
    geom_line(data = mean_dat, aes(y=Estimate, x=pop_dens),
              colour="#0396A6", size=1)+
        scale_y_continuous(breaks = c(0,1,2),                        limits = c(0,2.5))+
    labs(x = "Human population density",
         y = "Number of threats")) 

# Panel g: Habitat breath ------------------------------------------------------

hab_seq <- tibble(Habitat_breadth_IUCN= seq(min(lpi$Habitat_breadth_IUCN, na.rm = T),
                                            max(lpi$Habitat_breadth_IUCN, na.rm = T),
                                            length.out = 100))

dat <- fitted(hab_1, 
              newdata = hab_seq,
              re_formula =NA,
              summary = F) %>% 
  data.frame() %>%
  pivot_longer(everything()) %>% 
  bind_cols(expand(hab_seq, iter=1:7200, nesting(Habitat_breadth_IUCN))) 

mean_dat <- dat %>% 
  group_by(Habitat_breadth_IUCN) %>% 
  summarise(Estimate=mean(value))

(g3g <- dat %>% 
    filter(iter%in%c(1:250)) %>% 
    ggplot(aes(x = Habitat_breadth_IUCN, 
               y = value, group=iter)) +
    geom_line(alpha=0.05, color="#A6323B")+
    geom_line(data = mean_dat, aes(y=Estimate, x=Habitat_breadth_IUCN),
              colour="#A6323B", size=1)+
        scale_y_continuous(breaks = c(0,1,2),                        limits = c(0,2.5))+
    labs(x = "Habitat breadth",
         y = "Number of threats")) 

# Combine figure 3 -------------------------------------------------------------

(row1 <- g3a+g3b+g3c+plot_layout(nrow = 1)+plot_annotation(tag_levels = "a")& 
   theme(plot.tag = element_text(face = "bold")))
(row2 <- g3d+g3e+ylab("")+g3f+ylab("")+g3g+ylab("")&
    plot_annotation(tag_levels = list(c("d", "e", "f", "g")))&
    plot_layout(nrow=1)& 
    theme(plot.tag = element_text(face = "bold")))

 (figure3 <- plot_grid(row1,
                        row2, 
                        nrow = 2))
 #                       align = "hv", axis = "bl",
 #                       labels = "auto")) 

ggsave(figure3, filename = "Figure3.pdf",
       height = 8, width = 12,path = ResultPath)

# Figure 4 #####################################################################

# Create a sequence for each of them 

bm_seq <- tibble(bm_g= seq(min(log(lpi$bm_g+1), na.rm = T),
                           max(log(lpi$bm_g+1), na.rm = T),
                           length.out = 100))

# Do summaries for each of them

mu <- fitted(bm_amp_fre,
             newdata = bm_seq,
             re_formula = NA) %>%
  as_tibble() %>%
  bind_cols(bm_seq) %>% 
  mutate(system="Freshwater",
         taxon="Amphibians") %>% 
  bind_rows(bm_amp_ter %>% fitted(newdata = bm_seq,
                                   re_formula = NA) %>%
              as_tibble() %>%
              bind_cols(bm_seq) %>% 
              mutate(system="Terrestrial",
                     taxon="Amphibians")) %>% 
  bind_rows(bm_bird_fre %>% fitted(newdata = bm_seq,
                                    re_formula = NA) %>%
              as_tibble() %>%
              bind_cols(bm_seq) %>% 
              mutate(system="Freshwater",
                     taxon="Birds")) %>% 
  bind_rows(bm_bird_mar %>% fitted(newdata = bm_seq,
                                    re_formula = NA) %>%
              as_tibble() %>%
              bind_cols(bm_seq) %>% 
              mutate(system="Marine",
                     taxon="Birds")) %>% 
  bind_rows(bm_bird_ter %>% fitted(newdata = bm_seq,
                                    re_formula = NA) %>%
              as_tibble() %>%
              bind_cols(bm_seq) %>% 
              mutate(system="Terrestrial",
                     taxon="Birds")) %>% 
  bind_rows(bm_car_mar %>% fitted(newdata = bm_seq,
                                   re_formula = NA) %>%
              as_tibble() %>%
              bind_cols(bm_seq) %>% 
              mutate(system="Marine",
                     taxon="Cartilaginous Fish")) %>% 
  bind_rows(bm_fish_fre %>% fitted(newdata = bm_seq,
                                    re_formula = NA) %>%
              as_tibble() %>%
              bind_cols(bm_seq) %>% 
              mutate(system="Freshwater",
                     taxon="Bony Fish")) %>% 
  bind_rows(bm_fish_mar %>% fitted(newdata = bm_seq,
                                    re_formula = NA) %>%
              as_tibble() %>%
              bind_cols(bm_seq) %>% 
              mutate(system="Marine",
                     taxon="Bony Fish")) %>% 
  bind_rows(bm_mam_fre %>% fitted(newdata = bm_seq,
                                   re_formula = NA) %>%
              as_tibble() %>%
              bind_cols(bm_seq) %>% 
              mutate(system="Freshwater",
                     taxon="Mammals")) %>% 
  bind_rows(bm_mam_mar %>% fitted(newdata = bm_seq,
                                   re_formula = NA) %>%
              as_tibble() %>%
              bind_cols(bm_seq) %>% 
              mutate(system="Marine",
                     taxon="Mammals")) %>% 
  bind_rows(bm_mam_ter %>% fitted(newdata = bm_seq,
                                   re_formula = NA) %>%
              as_tibble() %>%
              bind_cols(bm_seq) %>% 
              mutate(system="Terrestrial",
                     taxon="Mammals")) %>% 
  bind_rows(bm_rep_fre %>% fitted(newdata = bm_seq,
                                   re_formula = NA) %>%
              as_tibble() %>%
              bind_cols(bm_seq) %>% 
              mutate(system="Freshwater",
                     taxon="Reptiles")) %>% 
  bind_rows(bm_rep_mar %>% fitted(newdata = bm_seq,
                                   re_formula = NA) %>%
              as_tibble() %>%
              bind_cols(bm_seq) %>% 
              mutate(system="Marine",
                     taxon="Reptiles")) %>% 
  bind_rows(bm_rep_ter %>% fitted(newdata = bm_seq,
                                   re_formula = NA) %>%
              as_tibble() %>%
              bind_cols(bm_seq) %>% 
              mutate(system="Terrestrial",
                     taxon="Reptiles")) 

(figure4 <- mu %>%
    as.data.frame() %>% 
    ggplot(aes(x = bm_g, y = Estimate, 
               group=taxon)) +
    geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5, 
                    fill=taxon), 
                alpha=0.1) +
    geom_line(aes(y=Estimate, colour=taxon),
              size=1) +
    #ylim(-0.25,4)+
    labs(y="Number of threats", 
         x="log(Body mass)") +
    facet_wrap(.~ system)+   
    scale_color_jco(name="Class")+
    scale_fill_jco(name="Class")+
    theme(legend.position = "bottom", 
          legend.direction = "horizontal"))

ggsave(figure4, file="Figure4.pdf",
       width = 10, height = 6,
       path = ResultPath)

# Figure 5 #####################################################################

# Create a sequence for each of them 

lat_seq <- tibble(Latitude= seq(min(abs(lpi$Latitude), na.rm = T),
                           max(abs(lpi$Latitude), na.rm = T),
                           length.out = 100))

# Do summaries for each of them

mu <- fitted(lat_amp_fre,
             newdata = lat_seq,
             re_formula = NA) %>%
  as_tibble() %>%
  bind_cols(lat_seq) %>% 
  mutate(system="Freshwater",
         taxon="Amphibians") %>% 
  bind_rows(lat_amp_ter %>% fitted(newdata = lat_seq,
                                  re_formula = NA) %>%
              as_tibble() %>%
              bind_cols(lat_seq) %>% 
              mutate(system="Terrestrial",
                     taxon="Amphibians")) %>% 
  bind_rows(lat_bird_fre %>% fitted(newdata = lat_seq,
                                   re_formula = NA) %>%
              as_tibble() %>%
              bind_cols(lat_seq) %>% 
              mutate(system="Freshwater",
                     taxon="Birds")) %>% 
  bind_rows(lat_bird_mar %>% fitted(newdata = lat_seq,
                                   re_formula = NA) %>%
              as_tibble() %>%
              bind_cols(lat_seq) %>% 
              mutate(system="Marine",
                     taxon="Birds")) %>% 
  bind_rows(lat_bird_ter %>% fitted(newdata = lat_seq,
                                   re_formula = NA) %>%
              as_tibble() %>%
              bind_cols(lat_seq) %>% 
              mutate(system="Terrestrial",
                     taxon="Birds")) %>% 
  bind_rows(lat_car_mar %>% fitted(newdata = lat_seq,
                                  re_formula = NA) %>%
              as_tibble() %>%
              bind_cols(lat_seq) %>% 
              mutate(system="Marine",
                     taxon="Cartilaginous Fish")) %>% 
  bind_rows(lat_fish_fre %>% fitted(newdata = lat_seq,
                                   re_formula = NA) %>%
              as_tibble() %>%
              bind_cols(lat_seq) %>% 
              mutate(system="Freshwater",
                     taxon="Bony Fish")) %>% 
  bind_rows(lat_fish_mar %>% fitted(newdata = lat_seq,
                                   re_formula = NA) %>%
              as_tibble() %>%
              bind_cols(lat_seq) %>% 
              mutate(system="Marine",
                     taxon="Bony Fish")) %>% 
  bind_rows(lat_mam_fre %>% fitted(newdata = lat_seq,
                                  re_formula = NA) %>%
              as_tibble() %>%
              bind_cols(lat_seq) %>% 
              mutate(system="Freshwater",
                     taxon="Mammals")) %>% 
  bind_rows(lat_mam_mar %>% fitted(newdata = lat_seq,
                                  re_formula = NA) %>%
              as_tibble() %>%
              bind_cols(lat_seq) %>% 
              mutate(system="Marine",
                     taxon="Mammals")) %>% 
  bind_rows(lat_mam_ter %>% fitted(newdata = lat_seq,
                                  re_formula = NA) %>%
              as_tibble() %>%
              bind_cols(lat_seq) %>% 
              mutate(system="Terrestrial",
                     taxon="Mammals")) %>% 
  bind_rows(lat_rep_fre %>% fitted(newdata = lat_seq,
                                  re_formula = NA) %>%
              as_tibble() %>%
              bind_cols(lat_seq) %>% 
              mutate(system="Freshwater",
                     taxon="Reptiles")) %>% 
  bind_rows(lat_rep_mar %>% fitted(newdata = lat_seq,
                                  re_formula = NA) %>%
              as_tibble() %>%
              bind_cols(lat_seq) %>% 
              mutate(system="Marine",
                     taxon="Reptiles")) %>% 
  bind_rows(lat_rep_ter %>% fitted(newdata = lat_seq,
                                  re_formula = NA) %>%
              as_tibble() %>%
              bind_cols(lat_seq) %>% 
              mutate(system="Terrestrial",
                     taxon="Reptiles")) 

(figure5 <- mu %>%
  as.data.frame() %>% 
  ggplot(aes(x = Latitude, y = Estimate, 
             group=taxon)) +
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5, 
                  fill=taxon), 
              alpha=0.1) +
  geom_line(aes(y=Estimate, colour=taxon),
            size=1) +
  ylim(-0.25,4)+
  labs(y="Number of threats", 
       x="Absolut latitude") +
  facet_wrap(.~ system)+   
  scale_color_jco(name="Class")+
  scale_fill_jco(name="Class")+
    theme(legend.position = "bottom", 
          legend.direction = "horizontal"))

ggsave(figure5, file="Figure5.pdf",
       width = 12, height = 6,
       path = ResultPath)

