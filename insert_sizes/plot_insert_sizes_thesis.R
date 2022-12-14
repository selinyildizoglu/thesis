.libPaths('/data/Blizard-Rakyan/Fran/R/packages_4.1.1/')

library(tidyverse)

insert_sizes <- bind_rows(
  lapply(
    X=list.files('/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/insert_sizes/'),
    FUN=function(x) {
      read_table(paste0('/data/Blizard-Rakyan/Selin/projects/MEF_SD_DEC/RRBS/Data_Novogene/X204SC22010869-Z01-F001/insert_sizes/', x),
                 col_names='insert_size',
                 col_types=cols()) %>%
        group_by(insert_size) %>%
        summarise(n=n()) %>%
        ungroup() %>%
        mutate(f=n/sum(n)) %>%
        mutate(sample=gsub('.txt', '', x)) %>%
        relocate(sample)}))

insert_sizes %>% 
  write_csv('insert_sizes_f.csv')

(insert_sizes %>%
  ggplot(aes(x=insert_size, y=f)) +
  geom_line(aes(group=sample), alpha=0.5) +
  theme_bw() +
  labs(x='Insert Size', y='Frequency') +
  theme(
    aspect.ratio=1,
    panel.grid=element_blank())) %>% 
ggsave(filename='insert_sizes.png',
       width=5, height=5)
