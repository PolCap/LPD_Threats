# --------------------------------------------------------------------------------------- #
# - FILE NAME:   SupFigures.R         
# - DATE:        03/08/2021
# - DESCRIPTION: Supplementary figures 
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
library(ggdist)
library(rstan)
library(ggsci)
library(wesanderson)
library(RColorBrewer)
library(brms)

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

# Phylogenetic signal ##############

load(paste0(ResultPath, "/PSModels.RData"))

# Create a common data frame

psign_res <- rbind(data.frame(taxon= "Amhibians", system="Terrestrial", 
                              dist=psign_amp_ter$samples$H1), 
                   # data.frame(taxon= "Amhibians", system="Freshwater", 
                   #            dist=psign_amp_fre$samples$H1), 
                   data.frame(taxon= "Birds", system="Terrestrial", 
                              dist=psign_bird_ter$samples$H1),
                   data.frame(taxon= "Birds", system="Marine", 
                              dist=psign_bird_mar$samples$H1),
                   data.frame(taxon= "Birds", system="Freshwater", 
                              dist=psign_bird_fre$samples$H1),
                   data.frame(taxon= "Mammals", system="Terrestrial", 
                              dist=psign_mam_ter$samples$H1),
                   data.frame(taxon= "Mammals", system="Marine", 
                              dist=psign_mam_mar$samples$H1),
                   data.frame(taxon= "Mammals", system="Freshwater", 
                              dist=psign_mam_fre$samples$H1),
                   data.frame(taxon= "Reptiles", system="Terrestrial", 
                              dist=psign_rep_ter$samples$H1),
                   data.frame(taxon= "Reptiles", system="Freshwater", 
                              dist=psign_rep_fre$samples$H1),
                   data.frame(taxon= "Bony Fishes", system="Freshwater", 
                              dist=psign_fish_fre$samples$H1),
                   data.frame(taxon= "Bony Fishes", system="Marine", 
                              dist=psign_fish_mar$samples$H1),
                   data.frame(taxon= "Cartilaginous Fishes", system="Marine", 
                              dist=psign_car_mar$samples$H1))

# Calculate the median 

res_median <- psign_res %>%
  group_by(taxon, system) %>% 
  summarise(median=median(dist))

# Plot 

(ga <- psign_res %>% 
    ggplot(aes(y=taxon, x=dist, color=taxon)) +
    stat_halfeye(aes(color = taxon,
                     fill=after_scale(colorspace::lighten(color, .3))),
                 shape = 18,
                 point_size = 3,
                 interval_size = 1.8,
                 adjust = .5) +
    geom_text(data=res_median,
              aes(x = median, label = format(round(median, 2), nsmall = 2)),
              stat = "unique",
              color = "black",
              fontface = "bold",
              size = 3.4,
              nudge_y = .15)+
    scale_color_manual(values = wesanderson::wes_palette("Cavalcanti1", 
                                                         type= "continuous",
                                                         n = 6))+
    facet_wrap(.~system)+
    labs(x="Phylogenetic signal", y="")+
    theme(legend.position = "none")+
    xlim(-0.1,0.3))

# Save it

ggsave("FigureS1.pdf", ga, 
       width = 9, height = 6,
       path = ResultPath)



# # Panel b: Red list --------------------------------------------------------------
# 
# # Sample size
# 
# sample_size <- lpi %>% 
#   distinct(ID, .keep_all=T) %>% 
#   group_by(Red_list_category) %>% 
#   summarise(n=n())
# 
# med <- red_1 %>%
#   gather_draws(`b_.*`, regex = TRUE) %>%
#   mutate(Red_list_category = gsub("b_Red_list_category", "", .variable),
#          Red_list_category = gsub("LRDlc", "LR/lc", Red_list_category),
#          Red_list_category = gsub("LRDnt", "LR/nt", Red_list_category)) %>%
#   group_by(Red_list_category) %>%
#   summarise(med=quantile(.value,probs = 0.99))
# 
# sample_size<- sample_size %>% 
#   left_join(med) 
# 
# # Plot intervals 
# 
# dat <-  red_1 %>%
#   gather_draws(`b_.*`, regex = TRUE) %>%
#   mutate(Red_list_category = gsub("b_Red_list_category", "", .variable),
#          Red_list_category = gsub("LRDlc", "LR/lc", Red_list_category),
#          Red_list_category = gsub("LRDnt", "LR/nt", Red_list_category)) 
# 
# (g6b <- red_1 %>%
#     gather_draws(`b_.*`, regex = TRUE) %>%
#     median_qi(.width = c(.95, .8, .5)) %>%
#     mutate(Red_list_category = gsub("b_Red_list_category", "", .variable),
#            Red_list_category = gsub("LRDlc", "LR/lc", Red_list_category),
#            Red_list_category = gsub("LRDnt", "LR/nt", Red_list_category)) %>%
#     ggplot(aes(y = reorder(Red_list_category, -.value), 
#                x = .value)) +
#     #geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
#     stat_slab(data=dat, alpha=0.5, aes(fill=.variable)) +
#     geom_pointinterval(aes(xmin = .lower, xmax = .upper),
#                        interval_size_range = c(0.5, 2), 
#                        colour="grey40") +
#     geom_text(data = sample_size, aes(x=med, y=Red_list_category,
#                                       label = paste0("n=", n)),
#               vjust   = 0, hjust=-.1)+
#     # scale_fill_manual(values = wes_palette("Chevalier1",14,
#     #                                        type="continuous"))+
#     labs(x="Number of threats", y = "Red list category") +
#     theme(legend.position = "none"))

# (g5c <- fitted(pop_dens_1, 
#        newdata = pd_seq,
#        re_formula =NA,
#        probs = c(.05, .95)) %>%
#   data.frame() %>%
#   bind_cols(pd_seq) %>% 
#   ggplot(aes(x = pop_dens, y = Estimate)) +
#     geom_smooth(aes(y = Estimate, ymin = Q5, ymax = Q95),
#                 stat = "identity",
#                 alpha = 1/4, size = 1/2) +
#     labs(x = "Population density (standardized)",
#        y = "Number of threats"))
# # Panel b: Realm --------------------------------------------------------------
# 
# # Sample size
# 
# sample_size <- lpi %>% 
#   distinct(ID, .keep_all=T) %>% 
#   group_by(Realm) %>% 
#   summarise(n=n())
# 
# med <- realm_1 %>%
#   gather_draws(`b_.*`, regex = TRUE) %>%
#   mutate(Realm = gsub("b_Realm", "", .variable),
#          Realm=ifelse(Realm=="Atlanticnorthtemperate", "Atlantic north temperate",
#                       ifelse(Realm=="Atlantictropicalandsubtropical", 
#                              "Atlantic tropical and subtropical",
#                              ifelse(Realm=="IndoMMalayan", 
#                                     "Indo-Malayan", 
#                                     ifelse(Realm=="Pacificnorthtemperate",
#                                            "Pacific north temperate", 
#                                            ifelse(Realm=="SouthtemperateandAntarctic",
#                                                   "South temperate and Antarctic", 
#                                                   ifelse(Realm=="TropicalandsubtropicalIndoMPacific",
#                                                          "Tropical and subtropical Indo-Pacific",
#                                                          Realm))))))) %>%
#   group_by(Realm) %>%
#   summarise(med=quantile(.value,probs = 0.99))
# 
# sample_size<- sample_size %>% 
#   left_join(med)
# 
# # Plot intervals 
# 
# dat <- realm_1 %>%
#   gather_draws(`b_.*`, regex = TRUE) %>%
#   mutate(Realm = gsub("b_Realm", "", .variable),
#          Realm=ifelse(Realm=="Atlanticnorthtemperate", "Atlantic north temperate",
#                       ifelse(Realm=="Atlantictropicalandsubtropical", 
#                              "Atlantic tropical and subtropical",
#                              ifelse(Realm=="IndoMMalayan", 
#                                     "Indo-Malayan", 
#                                     ifelse(Realm=="Pacificnorthtemperate",
#                                            "Pacific north temperate", 
#                                            ifelse(Realm=="SouthtemperateandAntarctic",
#                                                   "South temperate and Antarctic", 
#                                                   ifelse(Realm=="TropicalandsubtropicalIndoMPacific",
#                                                          "Tropical and subtropical Indo-Pacific",
#                                                          Realm)))))))
# (g5b <- realm_1 %>%
#     gather_draws(`b_.*`, regex = TRUE) %>%
#     median_qi(.width = c(.95, .8, .5)) %>%
#     mutate(Realm = gsub("b_Realm", "", .variable),
#            Realm=ifelse(Realm=="Atlanticnorthtemperate", "Atlantic north temperate",
#                         ifelse(Realm=="Atlantictropicalandsubtropical", 
#                                "Atlantic tropical and subtropical",
#                                ifelse(Realm=="IndoMMalayan", 
#                                       "Indo-Malayan", 
#                                       ifelse(Realm=="Pacificnorthtemperate",
#                                              "Pacific north temperate", 
#                                              ifelse(Realm=="SouthtemperateandAntarctic",
#                                                     "South temperate and Antarctic", 
#                                                     ifelse(Realm=="TropicalandsubtropicalIndoMPacific",
#                                                            "Tropical and subtropical Indo-Pacific",
#                                                            Realm))))))) %>%
#     ggplot(aes(y = reorder(Realm, -.value), x = .value)) +
#     geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
#     stat_slab(data=dat, alpha=0.5, aes(fill=.variable)) +
#     geom_pointinterval(aes(xmin = .lower, xmax = .upper),
#                        interval_size_range = c(0.5, 2), 
#                        colour="grey40") +
#     geom_text(data = sample_size, aes(x=med, y=Realm, 
#                                       label = paste0("n=", n)),
#               vjust   = 0, hjust=-.1)+
#     scale_fill_manual(values = wes_palette("Chevalier1",14,
#                                              type="continuous"))+
#     xlim(-1, 2.2)+
#     labs(x="Number of threats", y = "Realm") +
#     theme(legend.position = "none"))


# Panel d: Protection ----------------------------------------------------------

# # Sample size
# 
# sample_size <- lpi %>% 
#   distinct(ID, .keep_all=T) %>% 
#   group_by(Protected_status) %>% 
#   summarise(n=n())
# 
# med <- prot_1 %>%
#   gather_draws(`b_.*`, regex = TRUE) %>%
#   mutate(Protected_status = gsub("b_Protected_status", "", .variable)) %>%
#   group_by(Protected_status) %>%
#   summarise(med=quantile(.value,probs = 0.99))
# 
# sample_size<- sample_size %>% 
#   left_join(med)
# 
# # Plot intervals 
# 
# dat <- prot_1 %>%
#   gather_draws(`b_.*`, regex = TRUE) %>%
#   mutate(Protected_status = gsub("b_Protected_status", "", .variable)) 
# 
# (g5d <- prot_1 %>%
#     gather_draws(`b_.*`, regex = TRUE) %>%
#     median_qi(.width = c(.95, .8, .5)) %>%
#     mutate(Protected_status = gsub("b_Protected_status", "", .variable))  %>%
#     ggplot(aes(y = reorder(Protected_status, -.value), x = .value)) +
#     stat_slab(data=dat, alpha=0.5, aes(fill=.variable)) +
#     geom_pointinterval(aes(xmin = .lower, xmax = .upper),
#                        interval_size_range = c(0.5, 2), 
#                        colour="grey40") +
#     geom_text(data = sample_size, aes(x=med, y=Protected_status, 
#                                       label = paste0("n=", n)),
#               vjust   = -1.5, hjust=-.1)+
#      scale_fill_manual(values = wes_palette("Royal1"))+
#     labs(x="Number of threats", y = "Protection status") +
#     theme(legend.position = "none"))
# 
# Panel e: Management ----------------------------------------------------------

# Sample size

# sample_size <- lpi %>% 
#   distinct(ID, .keep_all=T) %>% 
#   group_by(Managed) %>% 
#   summarise(n=n()) %>% 
#   mutate(Managed=ifelse(Managed==0, "Unmanaged", "Managed"))
# 
# med <- man_1 %>%
#   gather_draws(`b_.*`, regex = TRUE) %>%
#   mutate(Managed = gsub("b_as.factorManaged", "", .variable),
#          Managed=ifelse(Managed==0, "Unmanaged", "Managed")) %>%
#   group_by(Managed) %>%
#   summarise(med=quantile(.value,probs = 0.99))
# 
# 
# sample_size<- sample_size %>% 
#   left_join(med)
# 
# # Plot intervals 
# 
# dat <- man_1 %>%
#   gather_draws(`b_.*`, regex = TRUE) %>%
#   mutate(Managed = gsub("b_as.factorManaged", "", .variable),
#          Managed=ifelse(Managed==0, "Unmanaged", "Managed"))
# 
# (g5f <- man_1 %>%
#     gather_draws(`b_.*`, regex = TRUE) %>%
#     median_qi(.width = c(.95, .8, .5)) %>%
#     mutate(Managed = gsub("b_as.factorManaged", "", .variable),
#            Managed=ifelse(Managed==0, "Unmanaged", "Managed")) %>%
#     ggplot(aes(y = reorder(Managed, -.value), x = .value)) +
#     stat_slab(data=dat, alpha=0.5, aes(fill=.variable)) +
#     geom_pointinterval(aes(xmin = .lower, xmax = .upper),
#                        interval_size_range = c(0.5, 2), 
#                        colour="grey40") +
#     geom_text(data = sample_size, aes(x=med, y=Managed, 
#                                       label = paste0("n=", n)),
#               vjust   = -1.5, hjust=-.1)+
#     scale_fill_manual(values = wes_palette("Chevalier1"))+
#     labs(x="Number of threats", y = "Managed") +
#     theme(legend.position = "none"))

# Figure x #####################################################################

# Create a sequence for each of them 

bm_seq <- expand.grid(bm_g= seq(min(log10(lpi$bm_g+1), na.rm = T),
                                max(log10(lpi$bm_g+1), na.rm = T),
                                length.out = 100),
                      Trophic_level=c("Carnivore",
                                      "Herbivore",
                                      "Omnivore")) %>% 
  as.tibble() %>% 
  mutate(Trophic_level=as.character(Trophic_level))

# Do summaries for each of them

mu <- fitted(bm_amp_fre,
             newdata = bm_seq,
             re_formula = NA) %>%
  as_tibble() %>%
  bind_cols(bm_seq) %>% 
  mutate(system="Freshwater",
         taxon="Amphibians") %>% 
  bind_rows(bm_amp_ter %>% fitted(newdata = bm_seq%>% filter(Trophic_level!="Herbivore"),
                                  re_formula = NA) %>%
              as_tibble() %>%
              bind_cols(bm_seq%>% filter(Trophic_level!="Herbivore")) %>% 
              mutate(system="Terrestrial",
                     taxon="Amphibians")) %>% 
  bind_rows(bm_bird_fre %>% fitted(newdata = bm_seq %>% 
                                     filter(Trophic_level!="Omnivore"),
                                   re_formula = NA) %>%
              as_tibble() %>%
              bind_cols(bm_seq %>% filter(Trophic_level!="Omnivore")) %>% 
              mutate(system="Freshwater",
                     taxon="Birds")) %>% 
  bind_rows(bm_bird_fre %>% 
              fitted(newdata = bm_seq%>% filter(Trophic_level!="Omnivore"),
                     re_formula = NA) %>%
              as_tibble() %>%
              bind_cols(bm_seq %>% filter(Trophic_level!="Omnivore")) %>% 
              mutate(system="Freshwater",
                     taxon="Birds")) %>% 
  # bind_rows(bm_bird_mar %>% fitted(newdata = bm_seq%>% filter(Trophic_level!="Omnivore"),
  #                                  re_formula = NA) %>%
  #             as_tibble() %>%
  #             bind_cols(bm_seq%>% filter(Trophic_level!="Omnivore")) %>% 
  #             mutate(system="Marine",
  #                    taxon="Birds")) %>% 
  bind_rows(bm_bird_ter %>% fitted(newdata = bm_seq,
                                   re_formula = NA) %>%
              as_tibble() %>%
              bind_cols(bm_seq) %>% 
              mutate(system="Terrestrial",
                     taxon="Birds")) %>% 
  bind_rows(bm_mam_fre %>% fitted(newdata = bm_seq%>% filter(Trophic_level!="Omnivore"),
                                  re_formula = NA) %>%
              as_tibble() %>%
              bind_cols(bm_seq%>% filter(Trophic_level!="Omnivore")) %>% 
              mutate(system="Freshwater",
                     taxon="Mammals")) %>% 
  
  bind_rows(bm_mam_ter %>% fitted(newdata = bm_seq,
                                  re_formula = NA) %>%
              as_tibble() %>%
              bind_cols(bm_seq) %>% 
              mutate(system="Terrestrial",
                     taxon="Mammals")) %>% 
  bind_rows(bm_rep_ter %>% fitted(newdata = bm_seq,
                                  re_formula = NA) %>%
              as_tibble() %>%
              bind_cols(bm_seq) %>% 
              mutate(system="Terrestrial",
                     taxon="Reptiles")) %>% 
  bind_rows(bm_rep_ter %>% fitted(newdata = bm_seq,
                                  re_formula = NA) %>%
              as_tibble() %>%
              bind_cols(bm_seq) %>% 
              mutate(system="Terrestrial",
                     taxon="Reptiles")) %>% 
  bind_rows(bm_fish_fre %>% 
              fitted(newdata = bm_seq,
                     re_formula = NA) %>%
              as_tibble() %>%
              bind_cols(bm_seq) %>% 
              mutate(system="Freshwater",
                     taxon="Bony Fish")) %>% 
  bind_rows(bm_fish_mar %>% 
              fitted(newdata = bm_seq,
                     re_formula = NA) %>%
              as_tibble() %>%
              bind_cols(bm_seq) %>% 
              mutate(system="Marine",
                     taxon="Bony Fish")) %>% 
  bind_rows(bm_car_mar %>% 
              fitted(newdata = bm_seq,
                     re_formula = NA) %>%
              as_tibble() %>%
              bind_cols(bm_seq) %>% 
              mutate(system="Marine",
                     taxon="Cartilaginous Fish"))

(figurex <- mu %>%
    as.data.frame() %>% 
    ggplot(aes(x = bm_g, y = Estimate, 
               group=Trophic_level)) +
    geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5, 
                    fill=Trophic_level), 
                alpha=0.1) +
    geom_line(aes(y=Estimate, colour=Trophic_level),
              size=1) +
    ylim(-0.25,4)+
    labs(y="Number of threats", 
         x="Log10(Body mass+1)") +
    facet_grid(taxon ~ system)+   
    scale_color_jco(name="Trophic level")+
    scale_fill_jco(name="Trophic level"))

ggsave(figurex, file="Figurex.pdf",
       width = 10, height = 6,
       path = ResultPath)

# Figure xx ####################################################################

# Create a sequence for each of them 

hab_seq <- expand.grid(Habitat_breadth_IUCN= seq(min(lpi$Habitat_breadth_IUCN, na.rm = T),
                                                 max(lpi$Habitat_breadth_IUCN, na.rm = T),
                                                 length.out = 100),
                       Trophic_level=c("Carnivore",
                                       "Herbivore",
                                       "Omnivore")) %>% 
  as.tibble() %>% 
  mutate(Trophic_level=as.character(Trophic_level))

# Do summaries for each of them

mu <- fitted(hab_amp_fre,
             newdata = hab_seq,
             re_formula = NA) %>%
  as_tibble() %>%
  bind_cols(hab_seq) %>% 
  mutate(system="Freshwater",
         taxon="Amphibians") %>% 
  bind_rows(hab_amp_ter %>% fitted(newdata = hab_seq%>% filter(Trophic_level!="Herbivore"),
                                   re_formula = NA) %>%
              as_tibble() %>%
              bind_cols(hab_seq%>% filter(Trophic_level!="Herbivore")) %>% 
              mutate(system="Terrestrial",
                     taxon="Amphibians")) %>% 
  bind_rows(hab_bird_fre %>% fitted(newdata = hab_seq %>% 
                                      filter(Trophic_level!="Omnivore"),
                                    re_formula = NA) %>%
              as_tibble() %>%
              bind_cols(hab_seq %>% filter(Trophic_level!="Omnivore")) %>% 
              mutate(system="Freshwater",
                     taxon="Birds")) %>% 
  bind_rows(hab_bird_fre %>% 
              fitted(newdata = hab_seq%>% filter(Trophic_level!="Omnivore"),
                     re_formula = NA) %>%
              as_tibble() %>%
              bind_cols(hab_seq %>% filter(Trophic_level!="Omnivore")) %>% 
              mutate(system="Freshwater",
                     taxon="Birds")) %>% 
  # bind_rows(hab_bird_mar %>% fitted(newdata = hab_seq%>% filter(Trophic_level!="Omnivore"),
  #                                  re_formula = NA) %>%
  #             as_tibble() %>%
  #             bind_cols(hab_seq%>% filter(Trophic_level!="Omnivore")) %>% 
  #             mutate(system="Marine",
  #                    taxon="Birds")) %>% 
  bind_rows(hab_bird_ter %>% fitted(newdata = hab_seq,
                                    re_formula = NA) %>%
              as_tibble() %>%
              bind_cols(hab_seq) %>% 
              mutate(system="Terrestrial",
                     taxon="Birds")) %>% 
  bind_rows(hab_mam_fre %>% fitted(newdata = hab_seq%>% filter(Trophic_level!="Omnivore"),
                                   re_formula = NA) %>%
              as_tibble() %>%
              bind_cols(hab_seq%>% filter(Trophic_level!="Omnivore")) %>% 
              mutate(system="Freshwater",
                     taxon="Mammals")) %>% 
  
  bind_rows(hab_mam_ter %>% fitted(newdata = hab_seq,
                                   re_formula = NA) %>%
              as_tibble() %>%
              bind_cols(hab_seq) %>% 
              mutate(system="Terrestrial",
                     taxon="Mammals")) %>% 
  bind_rows(hab_rep_ter %>% fitted(newdata = hab_seq,
                                   re_formula = NA) %>%
              as_tibble() %>%
              bind_cols(hab_seq) %>% 
              mutate(system="Terrestrial",
                     taxon="Reptiles")) %>% 
  bind_rows(hab_rep_ter %>% fitted(newdata = hab_seq,
                                   re_formula = NA) %>%
              as_tibble() %>%
              bind_cols(hab_seq) %>% 
              mutate(system="Terrestrial",
                     taxon="Reptiles")) %>% 
  bind_rows(hab_fish_fre %>% 
              fitted(newdata = hab_seq,
                     re_formula = NA) %>%
              as_tibble() %>%
              bind_cols(hab_seq) %>% 
              mutate(system="Freshwater",
                     taxon="Fishes")) %>% 
  bind_rows(hab_fish_mar %>% 
              fitted(newdata = hab_seq,
                     re_formula = NA) %>%
              as_tibble() %>%
              bind_cols(hab_seq) %>% 
              mutate(system="Marine",
                     taxon="Fishes")) 

(figurexx <- mu %>%
    as.data.frame() %>% 
    ggplot(aes(x = Habitat_breadth_IUCN, y = Estimate, 
               group=Trophic_level)) +
    geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5, 
                    fill=Trophic_level), 
                alpha=0.1) +
    geom_line(aes(y=Estimate, colour=Trophic_level),
              size=1) +
    ylim(-0.25,5)+
    labs(y="Number of threats", 
         x="Habitat breath") +
    facet_grid(taxon ~ system)+   
    scale_color_jco(name="Trophic level")+
    scale_fill_jco(name="Trophic level"))

ggsave(figurexx, file="Figurexx.pdf",
       width = 10, height = 6,
       path = ResultPath)

# Figure 6: Interactive effects of factors with body mass ######################
# Freshwater amphibians --------------------------------------------------------

dat <- in_amp_fre %>%
  gather_draws(`b_.*`, regex = TRUE) %>%
  mutate(Factor = gsub("b_", "", .variable), 
         Factor =gsub("scalelog10bm_gP1", "Body Mass", Factor),
         Factor = gsub("Trophic_level", "", Factor),
         Factor = gsub("Habitat_breadth_IUCN", "Habitat breadth", Factor),
         Factor = gsub("scaleabs", "", Factor),
         Factor = gsub("pop_dens", "Human population density", Factor)) 

(g1 <- in_amp_fre %>%
    gather_draws(`b_.*`, regex = TRUE) %>%
    median_qi(.value, .width = c(.95, .8, .5)) %>% 
    mutate(Factor = gsub("b_", "", .variable), 
           Factor =gsub("scalelog10bm_gP1", "Body Mass", Factor),
           Factor = gsub("Trophic_level", "", Factor),
           Factor = gsub("Habitat_breadth_IUCN", "Habitat breadth", Factor),
           Factor = gsub("scaleabs", "", Factor),
           Factor = gsub("pop_dens", "Human population density", Factor), 
           Factor = factor(Factor, levels = c("Intercept", 
                                              "Body Mass:Latitude",
                                              "Body Mass:Human population density",
                                              "Body Mass:Habitat breadth",
                                              "Body Mass:Carnivore", 
                                              "Body Mass:Herbivore",
                                              "Body Mass:Omnivore"))) %>% 
    ggplot(aes(y = Factor, x = .value)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
    #stat_slab(data=dat, aes(fill=Factor)) +
    geom_pointinterval(aes(xmin = .lower, xmax = .upper),
                       interval_size_range = c(0.5, 2),
                       colour="grey40") +
    scale_x_continuous(
      labels = scales::number_format(accuracy = 0.01))+
    labs(x="", y = "", title = "Freshwater amphibians") + 
    theme(legend.position = "none"))

# Terrestrial amphibians -------------------------------------------------------

dat <- in_amp_ter %>%
  gather_draws(`b_.*`, regex = TRUE) %>%
  mutate(Factor = gsub("b_", "", .variable), 
         Factor =gsub("scalelog10bm_gP1", "Body Mass", Factor),
         Factor = gsub("Trophic_level", "", Factor),
         Factor = gsub("Habitat_breadth_IUCN", "Habitat breadth", Factor),
         Factor = gsub("scaleabs", "", Factor),
         Factor = gsub("pop_dens", "Human population density", Factor)) 

(g2 <- in_amp_ter %>%
    gather_draws(`b_.*`, regex = TRUE) %>%
    median_qi(.value, .width = c(.95, .8, .5)) %>% 
    mutate(Factor = gsub("b_", "", .variable), 
           Factor =gsub("scalelog10bm_gP1", "Body Mass", Factor),
           Factor = gsub("Trophic_level", "", Factor),
           Factor = gsub("Habitat_breadth_IUCN", "Habitat breadth", Factor),
           Factor = gsub("scaleabs", "", Factor),
           Factor = gsub("pop_dens", "Human population density", Factor), 
           Factor = factor(Factor, levels = c("Intercept", 
                                              "Body Mass:Latitude",
                                              "Body Mass:Human population density",
                                              "Body Mass:Habitat breadth",
                                              "Body Mass:Carnivore", 
                                              "Body Mass:Herbivore",
                                              "Body Mass:Omnivore"))) %>% 
    ggplot(aes(y = Factor, x = .value)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
    #stat_slab(data=dat, aes(fill=Factor)) +
    geom_pointinterval(aes(xmin = .lower, xmax = .upper),
                       interval_size_range = c(0.5, 2),
                       colour="grey40") +
    scale_x_continuous(
      labels = scales::number_format(accuracy = 0.01))+
    labs(x="", y = "", title = "Terrestrial amphibians") + 
    theme(legend.position = "none"))

# Freshwater birds -------------------------------------------------------

dat <- in_bird_fre %>%
  gather_draws(`b_.*`, regex = TRUE) %>%
  mutate(Factor = gsub("b_", "", .variable), 
         Factor =gsub("scalelog10bm_gP1", "Body Mass", Factor),
         Factor = gsub("Trophic_level", "", Factor),
         Factor = gsub("Habitat_breadth_IUCN", "Habitat breadth", Factor),
         Factor = gsub("scaleabs", "", Factor),
         Factor = gsub("pop_dens", "Human population density", Factor), 
         Factor = factor(Factor, levels = c("Intercept", 
                                            "Body Mass:Latitude",
                                            "Body Mass:Human population density",
                                            "Body Mass:Habitat breadth",
                                            "Body Mass:Carnivore", 
                                            "Body Mass:Herbivore",
                                            "Body Mass:Omnivore")))

(g3 <- in_bird_fre %>%
    gather_draws(`b_.*`, regex = TRUE) %>%
    median_qi(.value, .width = c(.95, .8, .5)) %>% 
    mutate(Factor = gsub("b_", "", .variable), 
           Factor =gsub("scalelog10bm_gP1", "Body Mass", Factor),
           Factor = gsub("Trophic_level", "", Factor),
           Factor = gsub("Habitat_breadth_IUCN", "Habitat breadth", Factor),
           Factor = gsub("scaleabs", "", Factor),
           Factor = gsub("pop_dens", "Human population density", Factor), 
           Factor = factor(Factor, levels = c("Intercept", 
                                              "Body Mass:Latitude",
                                              "Body Mass:Human population density",
                                              "Body Mass:Habitat breadth",
                                              "Body Mass:Carnivore", 
                                              "Body Mass:Herbivore",
                                              "Body Mass:Omnivore"))) %>% 
    ggplot(aes(y = Factor, x = .value)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
    #stat_slab(data=dat, aes(fill=Factor)) +
    geom_pointinterval(aes(xmin = .lower, xmax = .upper),
                       interval_size_range = c(0.5, 2),
                       colour="grey40") +
    scale_y_discrete(drop=FALSE)+
    scale_x_continuous(
      labels = scales::number_format(accuracy = 0.01))+
    labs(x="", y = "", title = "Freshwater birds") + 
    theme(legend.position = "none"))

# Terrestrial birds -------------------------------------------------------

dat <- in_bird_ter %>%
  gather_draws(`b_.*`, regex = TRUE) %>%
  mutate(Factor = gsub("b_", "", .variable), 
         Factor =gsub("scalelog10bm_gP1", "Body Mass", Factor),
         Factor = gsub("Trophic_level", "", Factor),
         Factor = gsub("Habitat_breadth_IUCN", "Habitat breadth", Factor),
         Factor = gsub("scaleabs", "", Factor),
         Factor = gsub("pop_dens", "Human population density", Factor)) 

(g4 <- in_bird_ter %>%
    gather_draws(`b_.*`, regex = TRUE) %>%
    median_qi(.value, .width = c(.95, .8, .5)) %>% 
    mutate(Factor = gsub("b_", "", .variable), 
           Factor =gsub("scalelog10bm_gP1", "Body Mass", Factor),
           Factor = gsub("Trophic_level", "", Factor),
           Factor = gsub("Habitat_breadth_IUCN", "Habitat breadth", Factor),
           Factor = gsub("scaleabs", "", Factor),
           Factor = gsub("pop_dens", "Human population density", Factor), 
           Factor = factor(Factor, levels = c("Intercept", 
                                              "Body Mass:Latitude",
                                              "Body Mass:Human population density",
                                              "Body Mass:Habitat breadth",
                                              "Body Mass:Carnivore", 
                                              "Body Mass:Herbivore",
                                              "Body Mass:Omnivore"))) %>% 
    ggplot(aes(y = Factor, x = .value)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
    #stat_slab(data=dat, aes(fill=Factor)) +
    geom_pointinterval(aes(xmin = .lower, xmax = .upper),
                       interval_size_range = c(0.5, 2),
                       colour="grey40") +
    scale_x_continuous(
      labels = scales::number_format(accuracy = 0.01))+
    labs(x="", y = "", title = "Terrestrial birds") + 
    theme(legend.position = "none"))

# Marine bony fishes -------------------------------------------------------

dat <- in_fish_mar %>%
  gather_draws(`b_.*`, regex = TRUE) %>%
  mutate(Factor = gsub("b_", "", .variable), 
         Factor =gsub("scalelog10bm_gP1", "Body Mass", Factor),
         Factor = gsub("Trophic_level", "", Factor),
         Factor = gsub("Habitat_breadth_IUCN", "Habitat breadth", Factor),
         Factor = gsub("scaleabs", "", Factor),
         Factor = gsub("pop_dens", "Human population density", Factor)) 

(g5 <- in_fish_mar %>%
    gather_draws(`b_.*`, regex = TRUE) %>%
    median_qi(.value, .width = c(.95, .8, .5)) %>% 
    mutate(Factor = gsub("b_", "", .variable), 
           Factor =gsub("scalelog10bm_gP1", "Body Mass", Factor),
           Factor = gsub("Trophic_level", "", Factor),
           Factor = gsub("Habitat_breadth_IUCN", "Habitat breadth", Factor),
           Factor = gsub("scaleabs", "", Factor),
           Factor = gsub("pop_dens", "Human population density", Factor), 
           Factor = factor(Factor, levels = c("Intercept", 
                                              "Body Mass:Latitude",
                                              "Body Mass:Human population density",
                                              "Body Mass:Habitat breadth",
                                              "Body Mass:Carnivore", 
                                              "Body Mass:Herbivore",
                                              "Body Mass:Omnivore"))) %>%     
    ggplot(aes(y = Factor, x = .value)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
    #stat_slab(data=dat, aes(fill=Factor)) +
    geom_pointinterval(aes(xmin = .lower, xmax = .upper),
                       interval_size_range = c(0.5, 2),
                       colour="grey40") +
    #scale_fill_manual(values = wesanderson::wes_palette("Cavalcanti1", n = 5))+
    scale_x_continuous(
      labels = scales::number_format(accuracy = 0.01))+
    labs(x="", y = "", title = "Marine bony fishes") + 
    theme(legend.position = "none"))

# Freshwater bony fishes -------------------------------------------------------

dat <- in_fish_fre %>%
  gather_draws(`b_.*`, regex = TRUE) %>%
  mutate(Factor = gsub("b_", "", .variable), 
         Factor =gsub("scalelog10bm_gP1", "Body Mass", Factor),
         Factor = gsub("Trophic_level", "", Factor),
         Factor = gsub("Habitat_breadth_IUCN", "Habitat breadth", Factor),
         Factor = gsub("scaleabs", "", Factor),
         Factor = gsub("pop_dens", "Human population density", Factor)) 

(g6 <- in_fish_fre %>%
    gather_draws(`b_.*`, regex = TRUE) %>%
    median_qi(.value, .width = c(.95, .8, .5)) %>% 
    mutate(Factor = gsub("b_", "", .variable), 
           Factor =gsub("scalelog10bm_gP1", "Body Mass", Factor),
           Factor = gsub("Trophic_level", "", Factor),
           Factor = gsub("Habitat_breadth_IUCN", "Habitat breadth", Factor),
           Factor = gsub("scaleabs", "", Factor),
           Factor = gsub("pop_dens", "Human population density", Factor), 
           Factor = factor(Factor, levels = c("Intercept", 
                                              "Body Mass:Latitude",
                                              "Body Mass:Human population density",
                                              "Body Mass:Habitat breadth",
                                              "Body Mass:Carnivore", 
                                              "Body Mass:Herbivore",
                                              "Body Mass:Omnivore"))) %>%   
    ggplot(aes(y = Factor, x = .value)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
    #stat_slab(data=dat, aes(fill=Factor)) +
    geom_pointinterval(aes(xmin = .lower, xmax = .upper),
                       interval_size_range = c(0.5, 2),
                       colour="grey40") +
    #scale_fill_manual(values = wesanderson::wes_palette("Cavalcanti1", n = 5))+
    scale_x_continuous(
      labels = scales::number_format(accuracy = 0.01))+
    labs(x="", y = "", title = "Freshwater bony fishes") + 
    theme(legend.position = "none"))

# Freshwater mammals -------------------------------------------------------

dat <- in_mam_fre %>%
  gather_draws(`b_.*`, regex = TRUE) %>%
  mutate(Factor = gsub("b_", "", .variable), 
         Factor =gsub("scalelog10bm_gP1", "Body Mass", Factor),
         Factor = gsub("Trophic_level", "", Factor),
         Factor = gsub("Habitat_breadth_IUCN", "Habitat breadth", Factor),
         Factor = gsub("scaleabs", "", Factor),
         Factor = gsub("pop_dens", "Human population density", Factor)) 

(g7 <- in_mam_fre %>%
    gather_draws(`b_.*`, regex = TRUE) %>%
    median_qi(.value, .width = c(.95, .8, .5)) %>% 
    mutate(Factor = gsub("b_", "", .variable), 
           Factor =gsub("scalelog10bm_gP1", "Body Mass", Factor),
           Factor = gsub("Trophic_level", "", Factor),
           Factor = gsub("Habitat_breadth_IUCN", "Habitat breadth", Factor),
           Factor = gsub("scaleabs", "", Factor),
           Factor = gsub("pop_dens", "Human population density", Factor), 
           Factor = factor(Factor, levels = c("Intercept", 
                                              "Body Mass:Latitude",
                                              "Body Mass:Human population density",
                                              "Body Mass:Habitat breadth",
                                              "Body Mass:Carnivore", 
                                              "Body Mass:Herbivore",
                                              "Body Mass:Omnivore"))) %>% 
    ggplot(aes(y = Factor, x = .value)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
    #stat_slab(data=dat, aes(fill=Factor)) +
    geom_pointinterval(aes(xmin = .lower, xmax = .upper),
                       interval_size_range = c(0.5, 2),
                       colour="grey40") +
    #scale_fill_manual(values = wesanderson::wes_palette("Cavalcanti1", n = 5))+
    scale_x_continuous(
      labels = scales::number_format(accuracy = 0.01))+
    labs(x="", y = "", title = "Freshwater mammals") + 
    theme(legend.position = "none"))

# Terrestrial mammals -------------------------------------------------------

dat <- in_mam_ter %>%
  gather_draws(`b_.*`, regex = TRUE) %>%
  mutate(Factor = gsub("b_", "", .variable), 
         Factor =gsub("scalelog10bm_gP1", "Body Mass", Factor),
         Factor = gsub("Trophic_level", "", Factor),
         Factor = gsub("Habitat_breadth_IUCN", "Habitat breadth", Factor),
         Factor = gsub("scaleabs", "", Factor),
         Factor = gsub("pop_dens", "Human population density", Factor)) 

(g8 <- in_mam_ter %>%
    gather_draws(`b_.*`, regex = TRUE) %>%
    median_qi(.value, .width = c(.95, .8, .5)) %>% 
    mutate(Factor = gsub("b_", "", .variable), 
           Factor =gsub("scalelog10bm_gP1", "Body Mass", Factor),
           Factor = gsub("Trophic_level", "", Factor),
           Factor = gsub("Habitat_breadth_IUCN", "Habitat breadth", Factor),
           Factor = gsub("scaleabs", "", Factor),
           Factor = gsub("pop_dens", "Human population density", Factor), 
           Factor = factor(Factor, levels = c("Intercept", 
                                              "Body Mass:Latitude",
                                              "Body Mass:Human population density",
                                              "Body Mass:Habitat breadth",
                                              "Body Mass:Carnivore", 
                                              "Body Mass:Herbivore",
                                              "Body Mass:Omnivore"))) %>% 
    ggplot(aes(y = Factor, x = .value)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
    #stat_slab(data=dat, aes(fill=Factor)) +
    geom_pointinterval(aes(xmin = .lower, xmax = .upper),
                       interval_size_range = c(0.5, 2),
                       colour="grey40") +
    #scale_fill_manual(values = wesanderson::wes_palette("Cavalcanti1", n = 5))+
    scale_x_continuous(
      labels = scales::number_format(accuracy = 0.01))+
    labs(x="", y = "", title = "Terrestrial mammals") + 
    theme(legend.position = "none"))

# Terrestrial reptiles -------------------------------------------------------

dat <- in_rep_ter %>%
  gather_draws(`b_.*`, regex = TRUE) %>%
  mutate(Factor = gsub("b_", "", .variable), 
         Factor =gsub("scalelog10bm_gP1", "Body Mass", Factor),
         Factor = gsub("Trophic_level", "", Factor),
         Factor = gsub("Habitat_breadth_IUCN", "Habitat breadth", Factor),
         Factor = gsub("scaleabs", "", Factor),
         Factor = gsub("pop_dens", "Human population density", Factor)) 

(g9 <- in_rep_ter %>%
    gather_draws(`b_.*`, regex = TRUE) %>%
    median_qi(.value, .width = c(.95, .8, .5)) %>% 
    mutate(Factor = gsub("b_", "", .variable), 
           Factor =gsub("scalelog10bm_gP1", "Body Mass", Factor),
           Factor = gsub("Trophic_level", "", Factor),
           Factor = gsub("Habitat_breadth_IUCN", "Habitat breadth", Factor),
           Factor = gsub("scaleabs", "", Factor),
           Factor = gsub("pop_dens", "Human population density", Factor), 
           Factor = factor(Factor, levels = c("Intercept", 
                                              "Body Mass:Latitude",
                                              "Body Mass:Human population density",
                                              "Body Mass:Habitat breadth",
                                              "Body Mass:Carnivore", 
                                              "Body Mass:Herbivore",
                                              "Body Mass:Omnivore"))) %>% 
    ggplot(aes(y = Factor, x = .value)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour="grey50") +
    #stat_slab(data=dat, aes(fill=Factor)) +
    geom_pointinterval(aes(xmin = .lower, xmax = .upper),
                       interval_size_range = c(0.5, 2),
                       colour="grey40") +
    #scale_fill_manual(values = wesanderson::wes_palette("Cavalcanti1", n = 5))+
    scale_x_continuous(
      labels = scales::number_format(accuracy = 0.01))+
    labs(x="", y = "", title = "Terrestrial reptiles") + 
    theme(legend.position = "none"))

# Combined figure --------------------------------------------------------------

# (figure6 <- plot_grid(g1,g2+ theme(axis.text.y = element_blank()),
#                       g3+ theme(axis.text.y = element_blank()),
#                      g4,g5+ theme(axis.text.y = element_blank()),
#                      g6+ theme(axis.text.y = element_blank()),
#                      g7,g8+ theme(axis.text.y = element_blank()),
#                      g9+ theme(axis.text.y = element_blank()),
#                      ncol=3, align = "hv",axis = "tbl"))

(figure6 <- (g1|g2+theme(axis.text.y = element_blank())|g3+theme(axis.text.y = element_blank()))/
   (g4|g5+theme(axis.text.y = element_blank())|g6+theme(axis.text.y = element_blank()))/
   (g7|g8+theme(axis.text.y = element_blank())|g9+theme(axis.text.y = element_blank())))

ggsave(figure6, filename = "Figure6.pdf",
       height = 10, width = 12,path = ResultPath)

