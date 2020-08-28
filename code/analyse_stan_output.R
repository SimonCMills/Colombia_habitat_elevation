## Analyse output of stan model. 


## functions ----
f <- function(b0, b1, b2, b3, r_scale, r_mid=0, habitat, ele, response = TRUE) {
    psi <- exp(r_scale - ifelse(ele<r_mid, b2*habitat, b3*habitat)) * -(ele-r_mid)^2 + 
        b1*habitat + b0   
    if(response) boot::inv.logit(psi) else psi
}

dtc <- function(psi, d0, d1, time) {
    psi * boot::inv.logit(d0 + d1*time)
}

x <- seq(-2, 2, len=100)
get_range_margin <- function(b0, habitat = c(0, .5), b1, b2, b3, r_scale, r_mid, pct = c(.1, .9)) {
    y_1 <- f(b0, b1, b2, b3, r_scale, r_mid, habitat[1], ele = x)
    y_2 <- f(b0, b1, b2, b3, r_scale, r_mid, habitat[2], ele = x)
    matrix(c(HQ_lwr = x[max(which(cumsum(y_1)/sum(y_1) < pct[1]))], 
             HQ_upr = x[min(which(cumsum(y_1)/sum(y_1) > pct[2]))],
             LQ_lwr = x[max(which(cumsum(y_2)/sum(y_2) < pct[1]))],
             LQ_upr = x[min(which(cumsum(y_2)/sum(y_2) > pct[2]))]), ncol=4)
}

range_unscaling <- function(ele_scaled, Upper, Lower) {
    scale <- (Upper-Lower)/2
    mid <- (Upper+Lower)/2
    ele_scaled *scale + mid
}

range_scaling <- function(ele, Upper, Lower) {
    scale <- (Upper-Lower)/2
    mid <- (Upper+Lower)/2
    (ele - mid)/(scale)
}

range_basefun <- function(ele, a0, a21, a22, b1, b21, b22, h) {
    d <- abs(ele)
    a0 + ifelse(ele <0, a21, a22)*d^2 + b1*h + ifelse(ele<0, b21, b22)*h*d
}


## packages ----
library(dplyr); library(ggplot2); library(posterior)

## read draws ----
variable <- "pct_tc"
# draws <- readRDS("Y:/edwards_lab1/User/bo1scm/Colombia_rangeLims/output/draws_habit")
stan_data <- readRDS("Y:/edwards_lab1/User/bo1scm/Colombia_rangeLims/output/stanData_habitatModel_reparam2_pct_tc_D2020-08-04_T08-52.rds")
draws1 <- readRDS("Y:/edwards_lab1/User/bo1scm/Colombia_rangeLims/output/draws_habitatModel_reparam2_pct_tc_D2020-08-04_T08-52.rds")
draws2 <- readRDS("Y:/edwards_lab1/User/bo1scm/Colombia_rangeLims/output/draws_habitatModel_reparam2_pct_tc_D2020-08-07_T00-17.rds")
draws <- bind_draws(draws1, draws2, along="chain")
# draws <- readRDS("Y:/edwards_lab1/User/bo1scm/Colombia_rangeLims/output/draws_habitatModel_reparam2_pct_tc_D2020-08-04_T08-52.rds")
subset_draws(draws, "mu", regex=T) %>% summarise_draws()
# stan_data <- readRDS("output/stanData_habitatModel_asymm_2020-07-13 03_45_31.rds")

a0 <- subset_draws(draws, "a0\\[", regex = T) %>%
    as_draws_df() %>%
    select(iter=.draw, everything(), -.iteration, -.chain) %>%
    reshape2::melt(id.vars="iter") %>%
    as_tibble %>%
    mutate(id_sp = gsub("^.*\\[(.*)\\]$", "\\1", variable)) %>%
    select(id_sp, iter, a0 = value)

a21 <- subset_draws(draws, "a21\\[", regex = T) %>%
    as_draws_df() %>%
    select(iter=.draw, everything(), -.iteration, -.chain) %>%
    reshape2::melt(id.vars="iter") %>%
    as_tibble %>%
    mutate(id_sp = gsub("^.*\\[(.*)\\]$", "\\1", variable)) %>%
    select(id_sp, iter, a21 = value)


a22 <- subset_draws(draws, "a22\\[", regex = T) %>%
    as_draws_df() %>%
    select(iter=.draw, everything(), -.iteration, -.chain) %>%
    reshape2::melt(id.vars="iter") %>%
    as_tibble %>%
    mutate(id_sp = gsub("^.*\\[(.*)\\]$", "\\1", variable)) %>%
    select(id_sp, iter, a22 = value)



d0 <- subset_draws(draws, "d0\\[", regex = T) %>%
    as_draws_df() %>%
    select(iter=.draw, everything(), -.iteration, -.chain) %>%
    reshape2::melt(id.vars="iter") %>%
    as_tibble %>%
    mutate(id_sp = gsub("^.*\\[(.*)\\]$", "\\1", variable)) %>%
    select(id_sp, iter, d0 = value)

d1 <- subset_draws(draws, "d1\\[", regex = T) %>%
    as_draws_df() %>%
    select(iter=.draw, everything(), -.iteration, -.chain) %>%
    reshape2::melt(id.vars="iter") %>%
    as_tibble %>%
    mutate(id_sp = gsub("^.*\\[(.*)\\]$", "\\1", variable)) %>%
    select(id_sp, iter, d1 = value)

b1 <- subset_draws(draws, "b1\\[", regex = T) %>%
    as_draws_df() %>%
    select(iter=.draw, everything(), -.iteration, -.chain) %>%
    reshape2::melt(id.vars="iter") %>%
    as_tibble %>%
    mutate(id_sp = gsub("^.*\\[(.*)\\]$", "\\1", variable)) %>%
    select(id_sp, iter, b1 = value)

b21 <- subset_draws(draws, "b21\\[", regex = T) %>%
    as_draws_df() %>%
    select(iter=.draw, everything(), -.iteration, -.chain) %>%
    reshape2::melt(id.vars="iter") %>%
    as_tibble %>%
    mutate(id_sp = gsub("^.*\\[(.*)\\]$", "\\1", variable)) %>%
    select(id_sp, iter, b21 = value)

b22 <- subset_draws(draws, "b22\\[", regex = T) %>%
    as_draws_df() %>%
    select(iter=.draw, everything(), -.iteration, -.chain) %>%
    reshape2::melt(id.vars="iter") %>%
    as_tibble %>%
    mutate(id_sp = gsub("^.*\\[(.*)\\]$", "\\1", variable)) %>%
    select(id_sp, iter, b22 = value)

all <- full_join(a21, a22) %>%
    full_join(., a0) %>%
    full_join(., b1) %>%
    full_join(., b21) %>%
    full_join(., b22) %>%
    full_join(., d0) %>%
    full_join(., d1) %>%
    mutate(id_sp = as.numeric(id_sp)) %>%
    group_by(id_sp, iter)# %>%
    # mutate(range_margin = get_range_margin(b0 = b0, b1 = b1, b2 = b2, b3 = b3, 
    #                                         r_scale = r_scale, r_mid = r_mid, pct = c(.05, .95))) %>%
    # ungroup() 

# all <- all %>%
#     mutate(lwr_100 = range_margin[,1], 
#            lwr_50 = range_margin[,3], 
#            upr_100 = range_margin[,2], 
#            upr_50 = range_margin[,4],
#            lwr_diff = abs(range_margin[,1]) - abs(range_margin[,3]),
#            upr_diff = abs(range_margin[,2]) - abs(range_margin[,4]), 
#            rWidth_100 = upr_100 - lwr_100, 
#            rWidth_50 = upr_50 - lwr_50)

## parameter estimates by species ----
var_lookup <- c(b21 = "b21: hab-ele lower", 
                b22 = "b22: hab-ele upper", 
                a0 = "a0: 'commonness'",
                b1 = "b1: habitat intercept effect", 
                a21 = "a21: lower-range scaling",
                a22 = "a22:upper-range scaling",
                d0 = "d0: detection intercept", 
                d1 = "d1: detection time-of-day", 
                upr_diff="Upper difference", 
                lwr_diff = "Lower difference")
var_subset <- c("a0", "a21", "a22", "b1", "b21", "b22")
var_order <- c("a0", "a21", "a22", "b21", "b22", "b1")


# get order of 'commonness'
all_summary <- all %>% 
    select(-iter) %>%
    reshape2::melt(., id.vars="id_sp") %>%
    group_by(id_sp, variable) %>%
    summarise_all(list(mid=mean, lwr=~quantile(., probs=.05), upr=~quantile(., probs=.95))) %>%
    left_join(., stan_data$species_lookup) %>%
    filter(!variable == "iter") %>%
    mutate(variable = factor(variable, levels=var_order)) %>% 
    group_by(species, variable) %>%
    mutate(xmin=min(lwr), xmax=max(upr))

levels_Sp <- all_summary %>% filter(variable=="b21") %>%
    arrange(desc(mid)) %>% pull(species) %>% as.character

# get the panel limits
x_lims <- all_summary %>%
    group_by(variable) %>%
    summarise(xmin=min(lwr), xmax=max(upr)) %>%
    filter(variable %in% var_subset)

all_summary <- all_summary %>% mutate(species = factor(species, levels=levels_Sp))


# get
mu_summ <- subset_draws(draws, "mu", regex = T) %>%
    summarise_draws() %>%
    as_tibble %>%
    mutate(variable = gsub("mu_", "", variable),
           variable = factor(variable, levels=var_order)) %>%
    filter(!is.na(variable))

breaks_fun <- function(x) {
    if(max(x)-min(x) < 2) {
        seq(-100, 100, .5)
        } else if(max(x)-min(x) < 4) {
        seq(-100, 100, 1)
        } else if(max(x)-min(x) < 8) {
            seq(-100, 100, 2)
        } else if(max(x)-min(x) < 16) {
            seq(-100, 100, 4)}
}

all_summary %>%
    filter(variable %in% c("a21", "b21", "a22", "b22")) %>%
    dcast(., id_sp~variable, value.var = "mid") %>%
    select(-1) %>%
    pairs()


ggplot(all_summary, aes()) + geom_point()

p1 <- all_summary %>%
    ungroup %>%
    filter(variable %in% var_subset, 
           species %in% levels_Sp[1:20]) %>%
    ggplot(aes(mid, species)) +
    geom_point() +
    geom_errorbarh(aes(xmin=lwr, xmax=upr), height=0) +
    scale_x_continuous(breaks=breaks_fun) +
    facet_wrap(~variable, ncol=8, scale="free_x",
               labeller = labeller(
                   variable = var_lookup
               )) +
    geom_vline(xintercept=0, lty=2) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(colour="black"),
          axis.title.x = element_blank(),
          strip.text=element_text(hjust=0, face="bold"),
          strip.background = element_blank()) +
    geom_blank(data=x_lims, aes(x = xmin), inherit.aes=F) +
    geom_blank(data=x_lims, aes(x = xmax), inherit.aes=F)

p2 <- all_summary %>%
    ungroup %>%
    filter(variable %in% var_subset, 
           species %in% levels_Sp[101:120]) %>%
    ggplot(aes(mid, species)) +
    geom_point() +
    geom_errorbarh(aes(xmin=lwr, xmax=upr), height=0) +
    scale_x_continuous(breaks=breaks_fun) +
    facet_wrap(~variable, ncol=8, scale="free_x",
               labeller = labeller(
                   variable = var_lookup
               )) +
    geom_vline(xintercept=0, lty=2) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(colour="black"),
          axis.title.x = element_blank(),
          strip.text=element_text(hjust=0, face="bold"),
          strip.background = element_blank()) +
    geom_blank(data=x_lims, aes(x = xmin), inherit.aes=F) +
    geom_blank(data=x_lims, aes(x = xmax), inherit.aes=F)

p3 <- all_summary %>%
    ungroup %>%
    filter(variable %in% var_subset, 
           species %in% levels_Sp[172:191]) %>%
    ggplot(aes(mid, species)) +
    geom_point() +
    geom_errorbarh(aes(xmin=lwr, xmax=upr), height=0) +
    scale_x_continuous(breaks=breaks_fun) +
    facet_wrap(~variable, ncol=8, scale="free_x",
               labeller = labeller(
                   variable = var_lookup
               )) +
    geom_vline(xintercept=0, lty=2) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(colour="black"),
          axis.title.x = element_blank(),
          strip.text=element_text(hjust=0, face="bold"),
          strip.background = element_blank()) +
    geom_blank(data=x_lims, aes(x = xmin), inherit.aes=F) +
    geom_blank(data=x_lims, aes(x = xmax), inherit.aes=F)


p4 <- mu_summ %>%
    filter(variable %in% var_subset) %>%
    ggplot(aes(y=variable, x=mean, xmin=q5, xmax=q95)) +
    geom_point() + geom_errorbarh(height=0) +
    facet_wrap(~variable, scales="free", ncol=8) +
    scale_x_continuous(breaks=breaks_fun) +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          strip.text=element_blank(),
          axis.text=element_text(colour="black")) +
    geom_blank(data=x_lims, aes(x = xmin), inherit.aes=F) +
    geom_blank(data=x_lims, aes(x = xmax), inherit.aes=F) +
    geom_vline(xintercept=0, lty=2) +
    labs(x = "Parameter estimate (\u00B1 95% CI)")

p_both <- egg::ggarrange(p3, p2, p1, p4, ncol=1, heights=c(1,1,1,.05))

ggsave("figures/plot_range_limits.png",plot=p_both, width=300, height=200, units="mm")

## Elevation range by species ----
h90 <- quantile(stan_data$habitat, .9)
h10 <- quantile(stan_data$habitat, .1)

## rescale elevations 
species_lookup <- read.csv("data/initial_species_list.csv", as.is=T) 
colnames(species_lookup) <- tolower(colnames(species_lookup))
species_lookup <- species_lookup %>%
    mutate(species_clements = gsub(" ", "_", ebird), 
           species_hbw = gsub(" ", "_", hbw)) %>%
    select(-x)

ele_range <- read.csv("data/elevational_ranges_Quinones.csv", as.is=T) %>%
    mutate(species_clements = paste0(genus, "_", species)) %>%
    mutate(lower = ifelse(is.na(lower), 0, lower)) %>%
    left_join(species_lookup, .) %>%
    mutate(species = species_hbw) %>% 
    as_tibble %>%
    left_join(., stan_data$species_lookup) %>%
    filter(!is.na(id_sp))

# common_sp <- ele_dat %>% filter(Q==1) %>% group_by(id_sp) %>% summarise(N=n()) %>%
#     filter(N > 3)

ele_pred <- seq(-1.5, 1.5, len=40)
all_sub  <- all %>% ungroup %>% filter(id_sp %in% 1:191)
expit_pars <- c(paste0(c("p100"), c("_mid", "_lwr1", "_upr1")), 
               paste0(c("p50"), c("_mid", "_lwr1", "_upr1")))

ele_preds <- replicate(40, all_sub, simplify=F) %>%
    bind_rows(., .id = "id") %>%
    left_join(., ele_range) %>%
    mutate(id = as.numeric(id), 
           ele= ele_pred[id],
           midpoint = (upper+lower)/2,
           ele_raw = range_unscaling(ele, mid = (upper+lower)/2, 
                                    scale = (upper-lower)/2), 
           p100 = range_basefun(a0 = a0, a21=a21, a22 = a22, b1=b1, b21=b21, b22=b22, 
                                h=h90, ele=ele),
           p50 = range_basefun(a0 = a0, a21=a21, a22 = a22, b1=b1, b21=b21, b22=b22, 
                               h=h10, ele=ele),
           p_diff = p100 - p50) %>%
    group_by(id_sp, ele, ele_raw, midpoint) %>%
    summarise_at(.vars=vars("p100", "p50", "p_diff"), 
                 list(mid=mean, 
                      # lwr2=~quantile(., probs=.05), 
                      lwr1=~quantile(., probs=.25),
                      # upr2=~quantile(., probs=.95),
                      upr1=~quantile(., probs=.75))) %>%
    mutate_at(all_of(expit_pars), boot::inv.logit) %>%
    mutate_at(contains("p_diff", vars = names(.)), exp)

plot_ele_sp1 <- ele_preds %>%
    ggplot(aes(ele, p100_mid, group=id_sp)) + 
    scale_colour_viridis_c() +
    geom_line(alpha=.5) +
    # xlim(200, 3900) +
    geom_hline(yintercept=1, lty=2, col="red") 

plot_ele_sp2 <- ele_preds %>%
    ggplot(aes(ele, p_diff_mid, group=id_sp)) + 
    scale_colour_viridis_c() +
    geom_line(alpha=.5) +
    # xlim(200, 3900) +
    geom_hline(yintercept=1, lty=2, col="red") +
    scale_x_continuous(expand = c(0,0), breaks=seq(-100, 5, .5)) +
    scale_y_continuous(breaks=c(1, 10, 20, 30, 40)) +
    geom_hline(yintercept=1, lty=2) +
    labs(x = "Scaled elevation",
         y = "Odds ratio", 
         title = bquote("(c) Odds ratio ("~psi[90]~", "~psi[10]~")")) +
    theme_bw() +
    theme(axis.text=element_text(colour="black"), 
          plot.title = element_text(size=12),
          panel.background = element_blank(),
          panel.grid =element_blank(), 
          axis.ticks.length = unit(1.5, "mm"), axis.line = element_line())

## work out the elevation histogram
# ele_dat <- tibble(dist2 = stan_data$dist2, 
#                   half = stan_data$half, 
#                   id_sp = stan_data$id_sp,
#                   Q = stan_data$Q, 
#                   ele = ifelse(half == 1, -sqrt(dist2), sqrt(dist2))) %>%
    

# p_ele_sp <- ggplot(ele_dat, aes(ele, fill=factor(Q), group=factor(Q))) + 
#     geom_histogram(binwidth=.05, col="black") +
#     labs(x="Scaled elevation") +
#     scale_x_continuous(limits=c(-1,1), expand = c(0,0)) +
#     scale_y_continuous(trans="reverse", expand=expansion(mult=0.0),breaks=seq(0,500,250)) +
#     scale_fill_manual(values=c("grey90", "grey50")) +
#     theme(panel.background = element_blank(), 
#           # axis.line.y = element_line(),
#           panel.grid = element_blank(),
#           axis.text=element_text(colour="black"), 
#           axis.title.x = element_blank(),
#           axis.text.x=element_blank(),
#           axis.ticks.x=element_blank()) +
#     guides(fill=FALSE)


samp_sp <- sample(unique(ele_preds$id_sp), 8)
ele_preds_samp <- ele_preds %>% filter(id_sp %in% samp_sp)
ele_preds %>%
    ggplot(aes(ele_raw, p100_mid, group=id_sp)) + 
    # scale_colour_viridis_c() +
    geom_line(alpha=.3) +
    geom_line(data=ele_preds_samp, aes(col=factor(id_sp)), size=1) +
    guides(col=F) + 
    xlim(200, 3900) +
    geom_hline(yintercept=1, lty=2, col="red") +
    coord_fixed(1000)


p1 <- ggplot(ele_preds_samp, aes(ele_raw, group=id_sp, p100_mid, col=factor(id_sp))) +
    geom_line(size=.7) +
    geom_line(data=ele_preds, #%>% rename(id2 = id_sp), 
               aes(ele_raw, p100_mid, group=id_sp), alpha=.1, col="black") +
    # theme(aspect.ratio = .7) + 
    xlim(0, 3900) +
    guides(col=F) +ylim(0, 1)

p2 <- ggplot(ele_preds_samp, aes(ele_raw, group=id_sp, p_diff_mid, col=factor(id_sp))) +
    geom_line(size=.7) +
    geom_line(data=ele_preds, #%>% rename(id2 = id_sp), 
              aes(ele_raw, p_diff_mid, group=id_sp), alpha=.1, col="black") +
    # theme(aspect.ratio = .7) + 
    xlim(0, 3900) +
    guides(col=F) #+ylim(0, 1)


# p2 <- ggplot(ele_preds_samp, aes(ele_raw, group=id_sp, p50_mid, col=factor(id_sp))) +
#     geom_line(size=.7) +
#     geom_line(data=ele_preds, #%>% rename(id2 = id_sp), 
#               aes(ele_raw, p50_mid, group=id_sp), alpha=.1, col="black") +
#     theme(aspect.ratio = .7) +ylim(0, 1)

p3 <- ggplot(ele_preds_samp, aes(ele_raw, group=id_sp, p100_mid, col=factor(id_sp))) +
    geom_line(size=.7) +
    geom_line(aes(y=p50_mid), lty=2) +
    geom_ribbon(aes(ymin=p100_lwr1, ymax=p100_upr1), alpha=.1, col=NA) +
    geom_ribbon(aes(ymin=p50_lwr1, ymax=p50_upr1), alpha=.1, col=NA) +
    # geom_line(data=ele_preds, #%>% rename(id2 = id_sp),
    #           aes(ele_raw, p50_mid, group=id_sp), alpha=.1, col="black") +
    theme(#aspect.ratio = .7, 
          strip.background = element_blank(), 
          strip.text = element_text(face="bold", hjust=0)) + 
    # ylim(0, 1) + 
    facet_wrap(~id_sp, nrow=2, scales="free") +
    guides(col=F)

p4 <- ggplot(ele_preds_samp, aes(ele_raw, group=id_sp, p_diff_mid, col=factor(id_sp))) +
    geom_line(size=.7) +
    geom_ribbon(aes(ymin=p_diff_lwr1, ymax=p_diff_upr1), alpha=.1, col=NA) +
    # geom_line(aes(y=p50_mid), lty=2) +
    # geom_line(data=ele_preds, #%>% rename(id2 = id_sp),
    #           aes(ele_raw, p50_mid, group=id_sp), alpha=.1, col="black") +
    theme(#aspect.ratio = .7, 
          strip.background = element_blank(), 
          strip.text = element_text(face="bold", hjust=0)) + 
    facet_wrap(~id_sp, nrow=2, scales="free") +
    guides(col=F)


egg::ggarrange(p1, p2, p3, p4, ncol=2, widths=c(1, 1), heights=c(1.2,1))    
#+
    # geom_line(aes(y=p50_mid, col=id_sp), lty=2) #+
    # facet_grid(~id_sp) 
p1
# p1 <- ggplot(ele_preds_samp, aes(ele_raw, group=id_sp, p100_mid)) +
#     geom_line(col="red") +
#     geom_line(aes(y=p50_mid), col="red", lty=2) +
#     facet_grid(~id_sp)

p2 <- ggplot(ele_preds_samp, aes(ele_raw, group=id_sp, p_diff_mid)) +
    geom_line() +
    # geom_line(aes(y=p50_mid), lty=2) +
    facet_grid(~id_sp)

egg::ggarrange(p1, p2)

p_ele_all_50 <- ele_preds %>%
    left_join(all_sub_sum) %>%
    group_by(id_sp) %>%
    mutate(max_ele = ele_raw[which.max(p100_mid)]) %>%
    mutate(p50_mid = boot::inv.logit(p50_mid)) %>%
    ggplot(aes(ele, p50_mid, group=id_sp, col=max_ele)) + 
    scale_colour_viridis_c() +
    # facet_grid(~cut(max_ele, seq(0, 4000, 500))) +
    geom_line(alpha=.8) + #alpha=p100_upr1-p100_lwr1)) +
    # scale_y_continuous(trans="log") +
    # ylim(0, 100) +
    # xlim(200, 3900) +
    # geom_line(aes(y=boot::logit(p100_mid)), lty=2) +    # facet_grid(cut(max_ele, 5)~.) +
    geom_line(data=ele_mu_preds %>%  mutate(p50_mid = boot::inv.logit(p50_mid)), 
              aes(ele, p50_mid), size=2, col="red", inherit.aes=F) +
    geom_hline(yintercept=1, lty=2, col="red") +
    xlim(c(-1.1, 1.1))

p_ele_all_diff <- ele_preds %>%
    left_join(all_sub_sum) %>%
    group_by(id_sp) %>%
    mutate(max_ele = ele_raw[which.max(p100_mid)]) %>%
    ggplot(aes(ele, exp(p100_mid-p50_mid), group=id_sp, col=max_ele)) + 
    scale_colour_viridis_c() +
    # facet_grid(~cut(max_ele, seq(0, 4000, 500))) +
    geom_line(alpha=.8) + #alpha=p100_upr1-p100_lwr1)) +
    # scale_y_continuous(trans="log") +
    # ylim(0, 100) +
    # xlim(200, 3900) +
    # geom_line(aes(y=boot::logit(p100_mid)), lty=2) +    # facet_grid(cut(max_ele, 5)~.) +
    geom_line(data=ele_mu_preds, aes(ele, exp(p100_mid-p50_mid)), size=2, col="red", inherit.aes=F) +
    geom_hline(yintercept=1, lty=2, col="red") +
    xlim(c(-1.1, 1.1)) +
    ylim(c(0, 30))


egg::ggarrange(p_ele_all_100, p_ele_all_50, p_ele_all_diff)


ele_preds %>%
    left_join(all_sub_sum) %>%
    group_by(id_sp) %>%
    mutate(max_ele = ele_raw[which.max(p100_mid)]) %>%
    ggplot(aes(ele_raw, exp(p100_mid-p50_mid), group=id_sp, col=max_ele)) + 
    scale_colour_viridis_c() +
    # facet_grid(~cut(max_ele, seq(0, 4000, 500))) +
    geom_line(alpha=.8) + #alpha=p100_upr1-p100_lwr1)) +
    # scale_y_continuous(trans="log") +
    # ylim(0, 100) +
    # xlim(200, 3900) +
    # geom_line(aes(y=boot::logit(p100_mid)), lty=2) +    # facet_grid(cut(max_ele, 5)~.) +
    # geom_line(data=ele_mu_preds, aes(ele, exp(p100_mid-p50_mid)), size=2, col="red", inherit.aes=F) +
    geom_hline(yintercept=1, lty=2, col="red") +
    # xlim(c(-1.1, 1.1)) +
    ylim(c(0, 30))

   
hist(pars_species$a0_mid + pars_species$b1_mid*h90)
abline(v=pars_mu$a0_mid + pars_mu$b1_mid*h90, col="red")

hist(pars_species$a0_mid + pars_species$b1_mid*h10)
abline(v=pars_mu$a0_mid + pars_mu$b1_mid*h10, col="red")

range_basefun
    
ggplotly(p_ele_all_100)
   stan_data$species_lookup %>% filter(id_sp == 160)
p_ele_all_50 <- ele_preds %>%
        left_join(all_sub_sum) %>%
        group_by(id_sp) %>%
        mutate(max_ele = ele_raw[which.max(p100_mid)]) %>%
        ggplot(aes(ele_raw, p50_mid, group=id_sp, col=max_ele)) + 
        scale_colour_viridis_c() +
        # geom_ribbon(aes(ymin = p_relDiff_lwr1, ymax = p_relDiff_upr1), alpha=.5) +
        geom_line(alpha=.8) + #alpha=p100_upr1-p100_lwr1)) +
        scale_y_continuous() +
        xlim(200, 3900) +
        # facet_grid(cut(max_ele, 5)~.) +
        # geom_line(data=ele_mu_preds, col="red") +
        geom_hline(yintercept=1, lty=2, col="red")# +
    facet_wrap(~id_sp)
   
    ncols <- length(levels(cut(ele_preds$ele_raw, seq(0, 3600, 100), right=F)))
    cols <- colorRampPalette(c("#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", "#FFFFBF", 
                               "#FEE090", "#FDAE61", "#F46D43", "#D73027"))(ncols) 
    set.seed(101)
    ele_preds %>%
        left_join(all_sub_sum) %>%
        ungroup %>%
        filter(id_sp %in% sample(id_sp, 20)) %>%
        group_by(id_sp) %>%
        mutate(max_ele = ele_raw[which.max(p100_mid)], 
               max_ele_f = cut(max_ele, seq(0, 3600, 100), right=F)) %>%
        ggplot(aes(ele_raw, p50_mid, group=id_sp, col=max_ele_f)) + 
        scale_colour_manual(values=cols) +
        scale_fill_manual(values=cols) +
        # scale_fill_viridis_c() +
        geom_ribbon(aes(ymin = p50_lwr1, ymax = p50_upr1, fill=max_ele_f), alpha=.2, col=NA) +
        geom_ribbon(aes(ymin = p100_lwr1, ymax = p100_upr1, fill=max_ele_f), alpha=.2, col=NA) +
        geom_line(alpha=.8, lty=2) + 
            geom_line(aes(y=p100_mid), lty=1) + #alpha=p100_upr1-p100_lwr1)) +
        scale_y_continuous() +
        xlim(200, 3900) + facet_wrap(~id_sp, scales="free_y")
    

egg::ggarrange(p_ele_all_100, p_ele_all_50, ncol=2)

ele_preds %>%
    left_join(all_sub_sum) %>%
    ggplot(aes(ele_raw, p100_mid)) + 
    geom_line(aes(group=id_sp, col=b21_mid/a21_mid), alpha=.2) + #alpha=p100_upr1-p100_lwr1)) +
    # scale_y_continuous(trans="log") +
    geom_line(data=ele_mu_preds, col="red") +
    geom_hline(yintercept=1, lty=2, col="red") 





ggplot(ele_preds, aes(ele, p100_mid)) + geom_line(aes(group=id_sp)) +
    geom_ribbon(aes(ymin=p100_lwr1, ymax=p100_upr1), alpha=.1) +
    # scale_y_continuous(trans="log") +
    facet_wrap(~id_sp, scales="free") +
    geom_line(aes(y=p50_mid), col="red") +
    geom_ribbon(aes(ymin=p50_lwr1, ymax=p50_upr1), fill="red", alpha=.1)
    # geom_hline(yintercept=1, lty=2, col="red")

left_join(ele_preds, all_sub_sum) %>%
    ggplot(aes(ele_raw, p_relDiff_mid, group=id_sp, col=a0_mid)) + 
    geom_line(alpha=.7) + 
    xlim(c(200, 3900)) +
    # ylim(c(0, 1)) +
    labs("Occupancy 80% treecover")  +
    scale_colour_gradient2()

p1
ggplot(ele_preds, aes(ele_raw, p_relDiff_mid, group=id_sp)) + 
    geom_line(alpha=.4) + 
    xlim(c(200, 3900)) +
    # ylim(c(0, 1)) +
    labs("Occupancy 80% treecover") 

p1
p2 <- ggplot(ele_preds, aes(ele_raw, p50_mid)) + 
    geom_line(alpha=.4, col="red") + 
    xlim(c(200, 3900)) +
    ylim(c(0, 1)) + 
    labs("Occupancy 20% treecover")

egg::ggarrange(p1, p2, ncol=2)
library(plotly)
ggplotly(p1)

plot(stan_data$Q[stan_data$id_sp == 41] ~ stan_data$dist2[stan_data$id_sp == 41])
# facet_wrap(~cut(midpoint, seq(0, 4000, 500)))

# geom_ribbon(aes(ymin=p100_lwr1, ymax=p100_upr1), alpha=.1) +
# scale_y_continuous(trans="log") +
# facet_wrap(~id_sp, scales="free") +
geom_line(aes(y=p50_mid, group=id_sp), col="red") 
# geom_ribbon(aes(ymin=p50_lwr1, ymax=p50_upr1), fill="red", alpha=.1)
# geom_hline(yintercept=1, lty=2, col="red")


ele_pred <- seq(-1, 1, len=3)
all_sub  <- all %>% ungroup %>% filter(id_sp %in% 1:100)
ele_preds <- all_sub %>%
    left_join(., ele_range) %>%
    mutate(id = as.numeric(id), 
           ele= ele_pred[id],
           midpoint = (upper+lower)/2,
           ele_raw = range_unscaling(ele, mid = (upper+lower)/2, scale = (upper-lower)/2), 
           p100_rlwr = range_basefun(a0 = a0, a21=a21, a22 = a22, b1=b1, b21=b21, b22=b22, 
                                h=h90, ele=ele),
           p50_rlwr = range_basefun(a0 = a0, a21=a21, a22 = a22, b1=b1, 
                               b21=b21, b22=b22, h=h10, 
                               ele=ele), 
           pdiff = p100/p50) %>%
    group_by(id_sp, ele,midpoint) %>%
    summarise_at(.vars=vars("p100", "p50", "pdiff"), 
                 list(mid=mean, 
                      # lwr2=~quantile(., probs=.05), 
                      lwr1=~quantile(., probs=.25),
                      # upr2=~quantile(., probs=.95),
                      upr1=~quantile(., probs=.75))) #%>%
    # group_by(id_sp, ele) %>%
    # summarise(pdiff_diff1 = pdiff_mid[1] - pdiff_mid[2],
    #           pdiff_diff2 = pdiff_mid[3] - pdiff_mid[2])

ggplot(ele_preds, aes(pdiff_diff1)) + geom_histogram()

geom_point(position=position_dodge(width=.5)) 


    geom_errorbar(aes(ymin=pdiff_lwr1, ymax=pdiff_upr1), 
                  position=position_dodge(width=.5), 
                  width=0) +
    scale_y_continuous(trans="log") +
    geom_line(data=ele_mu_preds, aes(ele, p_relDiff_mid), col="red", inherit.aes=F)

splits <- ele_preds %>% split(., .$ele)
names(splits[[1]])[-c(28, 30:32)]
full_join(splits[[1]], splits[[2]], by=names(splits[[1]])[-c(1, 28, 30:33)], 
          suffix= c("_rmin", "_rcent")) %>%
    full_join(., splits[[3]], by=names(splits[[1]])[-c(1, 28, 30:33)], 
              suffix = c("", "_rcent"))
    
           
p_lwrmar = (p100/(1-p100))/(p50/(1-p50))) %>%
    group_by(id_sp, ele, ele_raw, midpoint) %>%
    summarise_at(.vars=vars("p100", "p50", "p_relDiff"), 
                 list(mid=mean, 
                      # lwr2=~quantile(., probs=.05), 
                      lwr1=~quantile(., probs=.05),
                          # upr2=~quantile(., probs=.95),
                          upr1=~quantile(., probs=.95)))


# ,
#            sim100 = rbinom(n(), 1, p100),
#            sim50 = rbinom(n(), 1, p50)) #%>%
#     
# 
# ePred_summ <- ele_preds %>% group_by(ele_raw, iter) %>% 
#     summarise(SR100 = sum(sim100), 
#               SR50 = sum(sim50)) %>%
#     group_by(ele_raw) %>%
#     summarise(SR100_lwr = quantile(SR100, .1), 
#               SR100_mid = mean(SR100),
#               SR100_upr = quantile(SR100, .9), 
#               SR50_lwr = quantile(SR50, .1), 
#               SR50_mid = mean(SR50),
#               SR50_upr = quantile(SR50, .9))
# 
# ggplot(ePred_summ, aes(ele_raw, SR100_mid, ymin=SR100_lwr, ymax=SR100_upr)) + 
#     geom_point() + 
#     geom_ribbon(alpha=.1) +
#     geom_point(aes(y=SR50_mid), col="red") +
#     geom_ribbon(aes(ymin=SR50_lwr, ymax=SR50_upr), alpha=.1, fill="red")
# 
# 
# 
ele_speciesPred <- ele_preds_sp %>%
    group_by(id_sp, Species, ele_raw, Lower, Upper) %>%
    summarise(p100_lwr = quantile(p100, 0.1, na.rm=T),
              p100_mid = quantile(p100, 0.5, na.rm=T),
              p100_upr = quantile(p100, 0.9, na.rm=T),
              p50_lwr = quantile(p50, 0.1, na.rm=T),
              p50_mid = mean(p50, na.rm=T),
              p50_upr = quantile(p50, 0.9, na.rm=T)) 
# ele_speciesPred$p50_lwr

ele_speciesPred %>% 
    filter(Lower < 500) %>%
    # group_by(Species, ele_raw) %>%
    ggplot(aes(ele_raw, p100_mid)) +
    geom_line() +
    geom_ribbon(aes(ymax=p100_upr, ymin=p100_lwr), alpha=.1) +
    geom_line(aes(y=p50_mid), col="red") +
    geom_ribbon(aes(ymin=p50_lwr, ymax=p50_upr), fill="red", alpha=.1) +
    facet_wrap(~Species) +
    theme(strip.background=element_blank(), 
          strip.text = element_text(hjust=0, face="bold")) +
    labs(x="Elevation (m)", y="P(Occupancy)")

ggsave("figures/range_exampes.png")    
# 
# ggplot(ele_preds_sub, aes(ele, boot::inv.logit(p100_mid) - boot::inv.logit(p50_mid))) + 
#     geom_point() +
#     # geom_line(aes(y=boot::inv.logit(p100_upr))) + 
#     # geom_line(aes(y=boot::inv.logit(p100_lwr))) +
#     # geom_line(aes(y=boot::inv.logit(p50_upr)), col="red") + 
#     # geom_line(aes(y=boot::inv.logit(p50_lwr)), col="red") + 
#     facet_wrap(~Species)
# 
# 
# ele_preds_hyper <- replicate(20, mu, simplify=F) %>%
#     bind_rows(., .id = "id") %>%
#     as_tibble %>%
#     mutate(id = as.numeric(id), 
#            ele = ele_pred[id], 
#            p100 = -(ele - mu_r_mid)^2 * (1) + mu_b0,
#            p50 = -(ele - mu_r_mid)^2 * (1 + b2$b2*.5) + mu_b0 + mu_b1*.5) %>%
#     group_by(ele) %>%
#     summarise(p100_lwr = quantile(p100, 0.1), 
#               p100_mid = quantile(p100, 0.5), 
#               p100_upr = quantile(p100, 0.9), 
#               p50_lwr = quantile(p50, 0.1), 
#               p50_mid = quantile(p50, 0.5), 
#               p50_upr = quantile(p50, 0.9)) 
# 
# ggplot(ele_preds_hyper, aes(ele, boot::inv.logit(p100_mid))) +
#     geom_line(aes(y=boot::inv.logit(p100_mid))) + 
#     geom_ribbon(aes(ymax = boot::inv.logit(p100_upr), 
#                     ymin = boot::inv.logit(p100_lwr)), alpha=.2) +
#     geom_line(aes(y=boot::inv.logit(p50_mid)), col="red") + 
#     geom_ribbon(aes(ymax = boot::inv.logit(p50_upr), 
#                     ymin = boot::inv.logit(p50_lwr)), fill="red", alpha=.2) 
# 
# 
# 
# 
## Average elevation relationship ----
ele_dat <- tibble(dist2 = stan_data$dist2, 
                  half = stan_data$half, 
                  id_sp = stan_data$id_sp,
                  Q = stan_data$Q, 
                  ele = ifelse(half == 1, -sqrt(dist2), sqrt(dist2)))

mu <- subset_draws(draws, "mu", regex = T) %>%
    as_draws_df() %>%
    as_tibble %>%
    select(iter=.draw, everything(), -.iteration, -.chain) 
expit_pars <- c(paste0(c("p100"), c("_mid", "_lwr1", "_upr1", "_lwr2", "_upr2")), 
               paste0(c("p50"), c("_mid", "_lwr1", "_upr1", "_lwr2", "_upr2")))
ele_pred <- seq(-1, 1, len=100)
ele_mu_preds <- replicate(100, mu, simplify=F) %>%
    bind_rows(., .id = "id") %>%
    mutate(id = as.numeric(id), 
           ele = ele_pred[id], 
           p100 = range_basefun(a0 = mu_a0, a21=mu_a21, a22 = mu_a22, b1=mu_b1, 
                                b21=mu_b21, b22=mu_b22, h=h90, ele=ele),
           p50 = range_basefun(a0 = mu_a0, a21=mu_a21, a22 = mu_a22, b1=mu_b1, 
                               b21=mu_b21, b22=mu_b22, h=h10, ele=ele),
           p_diff = p100 - p50,  
           p_relDiff = p100/p50) %>%
    group_by(ele) %>%
    summarise_at(.vars=vars("p100", "p50", "p_diff", "p_relDiff"), 
                 list(mid=mean, 
                      lwr2=~quantile(., probs=.1), 
                      lwr1=~quantile(., probs=.25),
                      upr2=~quantile(., probs=.90),
                      upr1=~quantile(., probs=.75))) %>%
    mutate_at(all_of(transPars), boot::inv.logit) %>%
    mutate_at(contains("p_diff", vars = names(.)), exp)

y_range <- c(0, max(ele_mu_preds$p100_upr2))

p_ele_1 <- ggplot(ele_mu_preds, aes((ele), ymin=p100_lwr2, y=p100_mid, ymax=p100_upr2)) +
    geom_line(col="royalblue1") +
    geom_ribbon(alpha=.1, fill="royalblue1") +
    geom_ribbon(aes(ymin=p100_lwr1, ymax=p100_upr1), alpha=.1, fill="royalblue1") +
    geom_line(aes(y=p50_mid), col="indianred") +
    geom_ribbon(aes(ymin=p50_lwr2, ymax=p50_upr2), alpha=.2, fill="indianred") +
    geom_ribbon(aes(ymin=p50_lwr1, ymax=p50_upr1), alpha=.2, fill="indianred") +
    labs(x = "Scaled elevation",
         y = "P(Occupancy)", 
         title = "(a) Occupancy probability") +
    theme_bw() +
    theme(axis.text=element_text(colour="black"), 
          plot.title = element_text(size=12),
          panel.background = element_blank(),
          panel.grid =element_blank(), 
          axis.ticks.length = unit(1.5, "mm"), axis.line = element_line())

p_3 <- ele_mu_preds %>%
    mutate(abs_cat = ifelse(p100_mid-p50_mid > .01, 1, 2), 
           abs_cat_grp = abs_cat * sign(ele), 
           abs_cat = as.factor(abs_cat)) %>%
    ggplot(aes(ele, ymin=p_diff_lwr2, y=p_diff_mid, ymax=p_diff_upr2)) +
    geom_line() +
    geom_ribbon(alpha=.1) +
    geom_ribbon(aes(ymin=p_diff_lwr1, ymax=p_diff_upr1), alpha=.2) +
    scale_x_continuous(expand = c(0,0), breaks=seq(-100, 5, .5)) +
    scale_y_continuous(breaks=c(1, 10, 20, 30, 40)) +
    geom_hline(yintercept=1, lty=2) +
    labs(x = "Scaled elevation",
         y = "Odds ratio", 
         title = bquote("(b) Odds ratio ("~psi[90]~", "~psi[10]~")")) +
    theme_bw() +
    theme(axis.text=element_text(colour="black"), 
          plot.title = element_text(size=12),
          panel.background = element_blank(),
          panel.grid =element_blank(), 
          axis.ticks.length = unit(1.5, "mm"), axis.line = element_line())

p_ele_20 <- ggplot(ele_dat, aes(ele, fill=factor(Q), group=factor(Q))) + 
    geom_histogram(binwidth=.05, col="black") +
    labs(x="Scaled elevation") +
    scale_x_continuous(limits=c(-1,1), expand = c(0,0)) +
    scale_y_continuous(trans="reverse", expand=expansion(mult=0.0),breaks=seq(0,500,250)) +
    scale_fill_manual(values=c("grey90", "grey50")) +
    theme(panel.background = element_blank(), 
          # axis.line.y = element_line(),
          panel.grid = element_blank(),
          axis.text=element_text(colour="black"), 
          axis.title.x = element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    guides(fill=FALSE)

library(egg)
mu_plots <- ggarrange(p_ele_1, plot_ele_sp1, p_3, plot_ele_sp2, p_ele_20, ncol=2, heights=c(1, 1, .3))

ggsave("figures/average_occupancy&detection_elevation_pct_tc_reparam.png", mu_plots, 
       width=150, height=250, units="mm")



## JUNK ----
## JUNK ----
# 
# 
# 
# 
# 
# mu <- subset_draws(draws, "mu", regex = T) %>%
#     as_draws_df() %>%
#     as_tibble %>%
#     select(iter=.draw, everything(), -.iteration, -.chain) 
# summarise_draws(mu)
# observed_sites <- stan_data$ele[stan_data$Q == 1]
# ele_range <- c(-1.5, 1.5)#range(observed_sites)
# ele_pred <- seq(-1.5, 1.5, len=40)
# ele_preds <- replicate(40, mu, simplify=F) %>%
#     bind_rows(., .id = "id") %>%
#     mutate(id = as.numeric(id), 
#            ele = ele_pred[id], 
#            p100 = f(b0 = mu_b0, b1 = mu_b1, b2 = mu_b2, b3 = mu_b3, r_scale = mu_r_scale, 
#                     r_mid=mu_r_mid, habitat=quantile(stan_data$habitat, .9), ele=ele),
#            p50 = f(b0 = mu_b0, b1 = mu_b1, b2 = mu_b2, b3 = mu_b3, r_scale = mu_r_scale, 
#                    r_mid = mu_r_mid,
#                    habitat=quantile(stan_data$habitat, .1), ele=ele), 
#            p100_lin = f(b0 = mu_b0, b1 = mu_b1, b2 = mu_b2, b3 = mu_b3, r_scale = mu_r_scale, 
#                     r_mid=mu_r_mid,habitat=quantile(stan_data$habitat, .9), ele=ele, response = F),
#            p50_lin = f(b0 = mu_b0, b1 = mu_b1, b2 = mu_b2, b3 = mu_b3, r_scale = mu_r_scale, 
#                    r_mid = mu_r_mid,
#                    habitat=quantile(stan_data$habitat, .1), ele=ele, response = F), 
#            p_diff_abs = p100 - p50,
#            p50_dtc = dtc(p50, mu_d0, mu_d1, -1.53),
#            p100_dtc = dtc(p100, mu_d0, mu_d1, -1.53)) %>%
#     group_by(ele) %>%
#     summarise_at(.vars=vars("p100", "p50", "p100_dtc", "p50_dtc", "p100_lin", 
#                             "p50_lin", "p_diff_abs"), 
#                  list(mid=mean, lwr=~quantile(., probs=.05), upr=~quantile(., probs=.95))) 
# 
# y_range <- c(0, max(ele_preds$p100_upr))
# 
# p_ele_1 <- ggplot(ele_preds, aes(ele, ymin=p100_lwr, y=p100_mid, ymax=p100_upr)) +
#     geom_line() +
#     geom_ribbon(alpha=.1) +
#     geom_line(aes(y=p50_mid), col="indianred") +
#     geom_ribbon(aes(ymin=p50_lwr, ymax=p50_upr), alpha=.2, fill="indianred") +
#     scale_x_continuous(limits=ele_range, expand = c(0,0), breaks=seq(-5, 5, .5)) +
#     scale_y_continuous(limits=y_range, expand=expansion(.04), breaks=seq(0, 1, .05)) +
#     # xlim(ele_range) +
#     # coord_flip()+
#     labs(x = "Scaled elevation",
#          y = "P(Occupancy)", 
#          title = "(a) Occupancy probability") +
#     theme_bw() +
#     theme(axis.text=element_text(colour="black"), 
#           plot.title = element_text(size=12),
#           panel.background = element_blank(),
#           panel.grid =element_blank(), 
#           axis.ticks.length = unit(1.5, "mm"), axis.line = element_line())
# 
# library(ggthemes)
# p_ele_2 <- ggplot() + 
#     geom_histogram(aes(x=stan_data$ele[stan_data$Q == 1]), binwidth=.05, col="black", fill="grey95") +
#     labs(x="Scaled elevation") +
#     scale_x_continuous(limits=ele_range, expand = c(0,0)) +
#     scale_y_continuous(trans="reverse", expand=expansion(mult=0.0),breaks=seq(0,80,20)) +
#     # geom_rangeframe(data = data.frame(y=c(0,60))) +
#     # ggthemes::theme_tufte() +
#     theme(panel.background = element_blank(), 
#           # axis.line.y = element_line(),
#           panel.grid = element_blank(),
#           axis.text=element_text(colour="black"), 
#           axis.title.x = element_blank(),
#           axis.text.x=element_blank(),
#           axis.ticks.x=element_blank())
# p_ele_2
# # p_dtc_1 <- ggplot(ele_preds, aes(ele, ymin=p100_dtc_lwr, y=p100_dtc_mid, ymax=p100_dtc_upr)) +
# #     geom_line(lty=2) +
# #     geom_ribbon(alpha=.1) +
# #     geom_line(aes(y=p50_dtc_mid), col="indianred", lty=2) +
# #     geom_ribbon(aes(ymin=p50_dtc_lwr, ymax=p50_dtc_upr), alpha=.2, fill="indianred") +
# #     scale_x_continuous(limits=ele_range, expand = c(0,0), breaks=seq(-5, 5, .5)) +
# #     scale_y_continuous(limits=y_range, expand=expansion(.04), breaks=seq(0, 1, .05)) +
# #     # xlim(ele_range) +
# #     # coord_flip()+
# #     labs(x = "Scaled elevation",
# #          y = "P(Occupancy)", 
# #          title = "(b) Detection probability") +
# #     theme_bw() +
# #     theme(axis.text=element_text(colour="black"), 
# #           plot.title = element_text(size=12),
# #           panel.background = element_blank(),
# #           panel.grid =element_blank(), 
# #           axis.ticks.length = unit(1.5, "mm"), axis.line = element_line())
# 
# p_ele_lin <- ggplot(ele_preds, aes(ele, ymin=p100_lin_lwr, y=p100_lin_mid, ymax=p100_lin_upr)) +
#     geom_line() +
#     geom_ribbon(alpha=.1) +
#     geom_line(aes(y=p50_lin_mid), col="indianred") +
#     geom_ribbon(aes(ymin=p50_lin_lwr, ymax=p50_lin_upr), alpha=.2, fill="indianred") +
#     scale_x_continuous(limits=ele_range, expand = c(0,0), breaks=seq(-100, 5, .5)) +
#     #scale_y_continuous(expand=expansion(.04), breaks=seq(0, 1, .05)) +
#     # xlim(ele_range) +
#     # coord_flip()+
#     labs(x = "Scaled elevation",
#          y = "P(Occupancy)", 
#          title = "(a) Occupancy probability") +
#     theme_bw() +
#     theme(axis.text=element_text(colour="black"), 
#           plot.title = element_text(size=12),
#           panel.background = element_blank(),
#           panel.grid =element_blank(), 
#           axis.ticks.length = unit(1.5, "mm"), axis.line = element_line()) +
#     theme_bw() +
#     theme(text = element_blank(), 
#           panel.grid=element_blank(), 
#           axis.ticks = element_blank(), plot.background = element_blank()) +
#     guides(lty=F)
# 
# p_ele_1_both <- p_ele_1 + 
#     annotation_custom(
#         ggplotGrob(p_ele_lin), 
#         xmin = .45, xmax = 1.5, ymax = y_range[2], ymin= y_range[2] *2/3
#     )
# 
# plot_diff <- ggplot(ele_preds, aes(ele, ymin=p_diff_abs_lwr, y=p_diff_abs_mid, ymax=p_diff_abs_upr)) +
#     geom_line() +
#     geom_ribbon(alpha=.1) +
#     # geom_line(aes(y=p50_lin_mid), col="indianred") +
#     # geom_ribbon(aes(ymin=p50_lin_lwr, ymax=p50_lin_upr), alpha=.2, fill="indianred") +
#     scale_x_continuous(limits=ele_range, expand = c(0,0), breaks=seq(-100, 5, .5)) +
#     #scale_y_continuous(expand=expansion(.04), breaks=seq(0, 1, .05)) +
#     # xlim(ele_range) +
#     # coord_flip()+
#     geom_hline(yintercept=0, lty=2) +
#     labs(x = "Scaled elevation",
#          y = bquote(psi[90] - psi[10]), 
#          title = "(b) Difference in occupancy") +
#     theme_bw() +
#     theme(axis.text=element_text(colour="black"), 
#           plot.title = element_text(size=12),
#           panel.background = element_blank(),
#           panel.grid =element_blank(), 
#           axis.ticks.length = unit(1.5, "mm"), axis.line = element_line())
# 
# mu_plots <- ggarrange(p_ele_1_both, plot_diff, p_ele_2, ncol=1, heights=c(1,1,.3))
# ggsave("figures/average_occupancy&detection_elevation_dist_from_edge.png", mu_plots, 
#        width=150, height=200, units="mm")

ggplot(all_summary, aes(x=mid, xmin=lwr, xmax=upr, y=id_sp, col=variable)) +
    geom_point() +
    # geom_boxplot(data=ele_forPlot, aes(Species_Clements, ele), inherit.aes=F) +
    geom_errorbarh(width=0)  +
    facet_wrap(~variable, )

geom_point(aes(y=midpoint), col="black") +
    geom_point(aes(y=Lower), col="black") +
    geom_point(aes(y=Upper), col="black") +
    # theme(legend.position = "bottom") +
    labs(y="Elevation (m)") +
    coord_flip()
sp_r_mid_order <- all_summary %>% ungroup %>% 
    mutate(midpoint = (upper+lower)/2) %>% select(Species, midpoint) %>%
    unique %>%
    arrange(midpoint) %>% pull(Species)

# rLims_summary <- all_summary %>%
#     filter(variable %in% c("lwr_100", "upr_100", "r_mid", "rWidth_100", "rWidth_50")) %>%
#     group_by(Species, variable) %>%
#     summarise_at(vars("mid", "lwr", "upr"), list(~range_unscaling(., upper, lower))) %>%
#     ungroup %>%
#     left_join(., lookup)%>%
#     mutate(midpoint = (Upper+Lower)/2) %>%
#     mutate(Species = factor(Species, levels=sp_r_mid_order))
# 
# 
# rLims_summary %>%
#     filter(!variable %in% c("rWidth_100", "rWidth_50")) %>%
#     ggplot(aes(x=mid, xmin=lwr, xmax=upr, y=Species, col=variable)) + 
#     geom_point() + 
#     geom_errorbarh(height=0) +
#     geom_point(aes(x=midpoint), col="black") +
#     geom_point(aes(x=Lower), col="black") +
#     geom_point(aes(x=Upper), col="black") +
#     # theme(legend.position = "bottom") +
#     labs(x="Elevation (m)") 
# 
# rLims_summary %>%
#     filter(!variable %in% c("rWidth_100", "rWidth_50")) %>%
#     ggplot(aes(y=mid, ymin=lwr, ymax=upr, x=Species, col=variable)) + 
#     geom_point() + 
#     geom_boxplot(data=ele_forPlot, aes(Species_Clements, ele), inherit.aes=F) +
#     geom_errorbar(width=0) +
#     geom_point(aes(y=midpoint), col="black") +
#     geom_point(aes(y=Lower), col="black") +
#     geom_point(aes(y=Upper), col="black") +
#     # theme(legend.position = "bottom") +
#     labs(y="Elevation (m)") +
#     coord_flip()
#     #geom_
ggsave("figures/elevational_limits_bySpecies.png", height=250,width=200, units="mm")

ele_forPlot <- eles %>% 
    filter(Species_Clements %in% lookup$Species) %>%
    group_by(Point, Species) %>%
    mutate(Q = max(`1`, `2`, `3`, `4`)) %>%
    filter(Q == 1) %>%
    mutate(Species_Clements = factor(Species_Clements, levels=levels(rLims_summary$Species))) #%>%
ggplot(aes(Species_Clements, ele)) + geom_boxplot() + coord_flip()
# rLims_summary %>%
#     filter(grepl("rWidth_100", variable)) %>%
#     mutate(lwr_at_0 = ifelse(Lower == 0, T, F)) %>%
#     ggplot(aes(x=midpoint, ymin=lwr , ymax=upr, y = mid, group=Species, col=lwr_at_0)) + 
#     geom_point(position=position_dodge(width=50)) + 
#     geom_errorbar(width=0, position = position_dodge(width=50)) +
#     coord_equal() +
#     geom_abline(lty=2) +
#     theme(legend.position = "bottom") +
#     labs(x="Quinones' midpoint (m)", 
#          y = "Estimated range width (m)")

p1_rLims <- rLims_summary %>%
    filter(grepl("upr_100", variable)) %>%
    mutate(lwr_at_0 = ifelse(Lower == 0, T, F)) %>%
    ggplot(aes(x=Upper, ymin=lwr , ymax=upr, y = mid, group=Species, col=lwr_at_0)) + 
    geom_point(position=position_dodge(width=50)) + 
    geom_errorbar(width=0, position = position_dodge(width=50)) +
    coord_equal() +
    geom_abline(lty=2) +
    theme(legend.position = "bottom") +
    labs(x="Quinones' upper elevation margin (m)", 
         y = "Estimated 95% margin (m)")

p2_rLims <- rLims_summary %>%
    filter(grepl("lwr_100", variable)) %>%
    mutate(lwr_at_0 = ifelse(Lower == 0, T, F)) %>%
    ggplot(aes(x=Lower, ymin=lwr , ymax=upr, y = mid, group=Species, col=lwr_at_0)) + 
    geom_point(position=position_dodge(width=50)) + 
    geom_errorbar(width=0, position = position_dodge(width=50)) +
    coord_equal() +
    geom_abline(lty=2) +
    theme(legend.position = "bottom") +
    labs(x="Quinones' upper elevation margin (m)", 
         y = "Estimated 5% margin (m)")


p_both <- egg::ggarrange(p1_rLims, p2_rLims, ncol=2)
ggsave("figures/estimated_range_margins.png",p_both,  width=210, units="mm")


## parameter covariance ----
pars_species <- all_sub_sum %>% dplyr::select(contains("mid")) %>%
    dplyr::select(-iter_mid)
pars_mu <- mu %>% 
    as_draws_df() %>%
    summarise_at(vars(mu_a0:mu_d0), list(mid=mean)) %>%
    replicate(10, simplify=F, .) %>% bind_rows
names(pars_mu) <- gsub("mu_", "", names(pars_mu))

pars_both <- bind_rows(pars_species, pars_mu, .id="id")
library(GGally)

ggpairs(pars_both, aes(col=id))
ggsave("figures/pars_covariance.png")

ggplot()
