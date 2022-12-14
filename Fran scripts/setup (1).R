library(tidyverse)
library(ggpubr)
library(patchwork)

# Variables

allele_colours <- c(
  'A'='orangered3',
  'C'='gold1',
  'G'='forestgreen',
  'T'='steelblue2'
)

allele_colours_cb <- c(
  'A'='#D55E00',
  'C'='#F0E442',
  'G'='#009E73',
  'T'='#56B4E9'
)

haplotype_colours <- c(
  'C_T_T_A'='chartreuse4',
  'C_T_A'='grey65',
  'CTA'='grey65',
  'C_T_C_A'='coral3',
  'C_C_A'='grey35',
  'CCA'='grey35',
  'A_G_T_A'='royalblue3',
  'A_T_A'='royalblue3',
  'ATA'='royalblue3',
  'A_G_T_G'='goldenrod3',
  'A_T_G'='goldenrod3',
  'ATG'='goldenrod3',
  'A_G_T_G_-178'='goldenrod4',
  'A_G_T_G_1495'='goldenrod3',
  'A_G_T_G_6832'='goldenrod2',
  'other'='grey50',
  'Other'='grey50',
  'Non-ATA'='grey50',
  'Ambiguous'='grey50',
  'A'='green4',
  'C'='grey50',
  'Sim'='red'
)

strains_colours <- c(
  '129S'='orange',
  'C3H'='dodgerblue', 
  'BL6N'='red',
  'BL6J'='springgreen4', 
  'CAST'='yellow',  
  'AJ'='purple'
)

strains_colours_cb <- c(
  '129S'='#E69F00',
  'C3H'='#56B4E9', 
  'BL6N'='#009E73',
  'BL6J'='#D55E00', 
  'CAST'='#F0E442',  
  'AJ'='#CC79A7'
)

strains_shapes <- c(
  '129S'=15,
  'C3H'=13, 
  'BL6N'=28,
  'BL6J'=11, 
  'CAST'=23,  
  'AJ'=5
)

coding_regions_m <- tribble(
  ~region,~start,~end,
  '18S',4008,5877,
  '5.8S',6878,7034,
  '28S',8123,12852)

coding_regions_h <- tribble(
  ~region,~start,~end,
  '18S',3658,5526,
  '5.8S',6597,6753,
  '28S',7921,12971)

data_fd <- '../../../data/'

# Auxiliary functions

adjust_pos <- function(pos, loop_size=3008) {
  if_else(pos<=loop_size, pos-loop_size-1, pos-loop_size)
}

adjust_haplo <- function(haplotype)  {
  case_when(
      haplotype=='A_G_T_A' ~ 'A_T_A',
      haplotype=='A_G_T_G' ~ 'A_T_G',
      haplotype=='C_T_C_A' ~ 'C_C_A',
      haplotype=='C_T_T_A' ~ 'C_T_A',
      TRUE ~ haplotype)
}

read_data <- function(fn, fd='', the_data_fd='') {
  the_data_fd <- ifelse(the_data_fd=='', data_fd, the_data_fd)
  read_csv(paste0(the_data_fd, fd, fn), col_types=cols())
}

read_blink <- function(fd='', loop_size=3008, the_data_fd='') {
  the_data_fd <- ifelse(the_data_fd=='', data_fd, the_data_fd)
  blink_fd <- paste0(the_data_fd, fd, 'blink/')
  bind_rows(
    lapply(X=list.files(blink_fd),
           FUN=function(x) read_csv(
             paste0(blink_fd, x, '/', x, '_link.csv'),
             col_types=cols()) %>%
             mutate(snp_pos=adjust_pos(snp_pos, loop_size)) %>%
             mutate(meth_pos=adjust_pos(meth_pos, loop_size)) %>%
             mutate(sample=x)))
}

read_alleles <- function(fd='', loop_size=3008, the_data_fd='') {
  the_data_fd <- ifelse(the_data_fd=='', data_fd, the_data_fd)
  alleles_fd <- paste0(the_data_fd, fd, 'snp_extract/')
  bind_rows(
    lapply(X=list.files(alleles_fd),
           FUN=function(x) read_csv(
             paste0(alleles_fd, x, '/', x, '_alleles.csv'),
             col_types=cols()) %>%
             mutate(snp_pos=adjust_pos(snp_pos, loop_size)) %>%
             mutate(sample=x)))
}

read_snps <- function(fd='', loop_size=3008, the_data_fd='') {
  the_data_fd <- ifelse(the_data_fd=='', data_fd, the_data_fd)
  blink_fd <- paste0(the_data_fd, fd, 'blink/')
  bind_rows(
    lapply(X=list.files(blink_fd),
           FUN=function(x) read_csv(
             paste0(blink_fd, x, '/', x, '_snps.csv'),
             col_types=cols()) %>%
             mutate(snp_pos=adjust_pos(snp_pos, loop_size)) %>%
             mutate(sample=x)))
}

read_bcf_all <- function(fn, loop_size=3008) {
  read_tsv(fn, comment = "#", col_names=F, col_types=cols()) %>%
    select(pos=X2, ref=X4, alt=X5, aux=X10) %>%
    mutate(alt=if_else(alt=="<*>", '', alt)) %>%
    filter(nchar(ref)==1) %>%
    separate(alt, into=c('alt'), sep=',', extra='drop') %>%
    filter(nchar(alt)<=1) %>%
    separate(aux, into=c('PL', 'AD'), sep=':') %>%
    select(-PL) %>%
    mutate(AD=paste(AD, 0, sep=',')) %>%
    separate(AD, into=c('num_ref', 'num_alt'), sep=',',
             extra='drop') %>%
    mutate(num_ref=as.numeric(num_ref),
           num_alt=as.numeric(num_alt)) %>%
    mutate(pos_ref=adjust_pos(pos, loop_size)) %>%
    select(pos_ref, everything(), -pos)
}


read_bcf <- function(fn, loop_size=3008) {
  read_tsv(fn, comment = "#", col_names=F, col_types=cols()) %>%
    select(pos=X2, ref=X4, alt=X5, aux=X10) %>%
    filter(alt!="<*>", nchar(ref)==1) %>%
    separate(alt, into=c('alt'), sep=',', extra='drop') %>%
    filter(nchar(alt)==1) %>%
    separate(aux, into=c('PL', 'AD'), sep=':') %>%
    select(-PL) %>%
    separate(AD, into=c('num_ref', 'num_alt'), sep=',', extra='drop') %>%
    mutate(num_ref=as.numeric(num_ref),
          num_alt=as.numeric(num_alt)) %>%
    mutate(pos_ref=adjust_pos(pos, loop_size)) %>%
    select(pos_ref, everything(), -pos)
}

read_mpileup <- function(fd='', loop_size=3008, the_data_fd='', all_pos=F) {
  the_data_fd <- ifelse(the_data_fd=='', data_fd, the_data_fd)
  mpileup_fd <- paste0(the_data_fd, fd, 'mpileup/')
  if (!all_pos) {
    bind_rows(
      lapply(X=list.files(mpileup_fd),
             FUN=function(x) read_bcf(
               paste0(mpileup_fd, x, '/', x, '.bcf'), loop_size) %>%
               mutate(sample=x)))
  } else {
    bind_rows(
      lapply(X=list.files(mpileup_fd),
             FUN=function(x) read_bcf_all(
               paste0(mpileup_fd, x, '/', x, '.bcf'), loop_size) %>%
               mutate(sample=x)))
  }
}

read_lofreq <- function(fd, fn='snps', the_data_fd='', loop_size=3008) {
  the_data_fd <- ifelse(the_data_fd=='', data_fd, the_data_fd)
  fd <- paste0(the_data_fd, '/', fd, '/')
  bind_rows(
    lapply(
      X=list.files(fd),
      FUN=function(x) { 
        read_tsv(
          paste0(fd, x, "/", x, ".", fn, ".vcf"), 
          comment = "#",
          col_names=F,
          col_types=cols()) %>%
          mutate(sample=x)})) %>%
    select(sample, pos=X2, ref=X4, alt=X5, info=X8) %>%
    separate(info, into=c('dp', 'af'), sep=';', extra='drop') %>%
    mutate(af=as.numeric(gsub('AF=', '', af)),
           dp=as.numeric(gsub('DP=', '', dp))) %>%
    mutate(pos_ref=adjust_pos(pos, loop_size=loop_size))
}

read_snps_lofreq <- function(fd, fn='snps', data_fd='../data') {
  fd <- paste0(data_fd, '/', fd, '/')
  bind_rows(
    lapply(
      X=list.files(fd),
      FUN=function(x) { 
        read_tsv(
          paste0(fd, x, "/", x, ".", fn, ".vcf"), 
          comment = "#",
          col_names=F,
          col_types=cols()) %>%
          mutate(sample=x)})) %>%
    select(sample, pos=X2, ref=X4, alt=X5, info=X8) %>%
    separate(info, into=c('dp', 'af'), sep=';', extra='drop') %>%
    mutate(af=as.numeric(gsub('AF=', '', af)),
           dp=as.numeric(gsub('DP=', '', dp)))
}

# Plotting functions

plot_dp <- function(df, condition, comparisons, method='wilcox.test',
                    values='meth', ylabel='Methylation Level',
                    binwidth=0.05, label_y=0.95, 
                    title='', subtitle='') {
  
  p <- df %>% ggplot(aes(x=get(condition), y=get(values), color=haplotype)) +
    geom_dotplot(aes(fill=haplotype),
      binaxis = "y", stackdir = "center", alpha=0.7, stroke=2, binwidth=binwidth)
  
  if (length(comparisons)==1) {
    p <- p + stat_compare_means(
      aes(label=paste0("p = ", ..p.format..)), comparisons=comparisons,
      method=method, label.y=label_y, label.x.npc=0.25, vjust=-0.25,
      size=3)  
  } else {
    p <- p + stat_compare_means(
      aes(label=paste0("p = ", ..p.format..)), comparisons=comparisons,
      method=method, size=3) 
  }
  
  p <- p + facet_grid(cols=vars(haplotype)) +
    theme_bw() +
    scale_color_manual(name='',
                       values=haplotype_colours) +
    scale_fill_manual(name='',
                      values=haplotype_colours) +
    scale_y_continuous(limits=c(0, 1), breaks=c(0, 0.5, 1),
                       labels=c('', 0.5, 1)) +
    labs(x='', y=ylabel) +
    theme(axis.text.x=element_text(angle=45, hjust=1),
          strip.text.y=element_text(angle=0),
          legend.position='none',
          panel.grid=element_blank())
  
  if(title != '' | subtitle != '') {
    p + ggtitle(title, subtitle)
  } else {
    p
  }
  
}

plot_cnt_covs <- function(df, the_strain, the_mouse, max_y=8000) {
  df %>%
    filter(strain==the_strain, antibody!='IGG', mouse==the_mouse) %>%
    ggplot() +
    geom_rect(data=coding_regions_m,
              aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf),
              fill='grey10', alpha=0.1) +
    geom_vline(xintercept=0, alpha=0.3) +
    geom_ribbon(aes(x=pos_ref, ymin=0, ymax=cov, fill=antibody)) +
    geom_text(aes(label=antibody), x=35000, y=(7/8)*max_y, size=4,
              data=data.frame(antibody=c('H3K27me3', 'H3K9me3'))) + 
    facet_grid(rows=vars(antibody)) +
    theme_bw() +
    labs(x='Position in rDNA reference', y='Read Depth') +
    scale_fill_manual(name='Antibody',
                      values=c('IGG'='gray70',
                               'H3K27me3'='darkorange3',
                               'H3K9me3'='darkseagreen4')) +
    scale_y_continuous(limits=c(0, max_y),
                       breaks=c(0, max_y/2, max_y),
                       labels=c('', max_y/2, max_y)) +
    scale_x_continuous(labels=c('TSS', seq(from=10000, to=40000, by=10000))) +
    theme(strip.background=element_blank(),
          strip.text=element_blank(),
          legend.position='none',
          panel.grid = element_blank(),
          panel.border=element_rect(size=1.1),
          panel.spacing=unit(0.2, 'lines'),
          axis.text=element_text(size=12),
          axis.title=element_text(size=14)) +
    ggtitle(the_mouse)
}
