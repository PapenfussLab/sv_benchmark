out_of_four <- 
    ensemble_df %>% 
    filter(CallSet == "All calls", 
          ncallers == 4, minhits %in% c(2,3),
          str_detect(callers, "gridss")) %>%
    rename(callers_4 = callers) %>%
    mutate(
        callers = str_replace_all(callers_4, "(,gridss)|(^gridss,)", ""))



out_of_four <- out_of_four %>% left_join(ensemble_df, by = c("callers", "CallSet")) %>% filter(minhits.y == 2, ncallers.y == 3)

qplot(tp.y, precision.y, data = out_of_four) + theme_cowplot() + theme(aspect.ratio = 1) + geom_segment(aes(xend = tp.x, yend = precision.x, color = ensemble.x), arrow = arrow(type = "closed", length = unit(5, "pt"))) + scale_color_brewer(palette = "Set1")