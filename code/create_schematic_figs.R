library(dplyr); library(ggplot2)

f <- function(ele, hab, mid, scale, b0, b1, b2, b3) {
    b0 + b1*hab + (scale - ifelse(ele < mid, b2*hab, b3*hab)) * -(ele-mid)^2
}

df <- expand.grid(ele = seq(-1, 1, len=100), 
            scale = 4, 
            mid = c(-.9, 0, .9), 
            b0 = 0,
            b1 = -.2, #c(-.2, 0, .2), 
            b2 = c(-1, 0, 1), 
            #b3 = c(-1, 0, 1),
            hab =c(-1,1)) %>%
    as_tibble() %>%
    mutate(b3 = b2, 
           scale = ifelse(mid == 0, 4, 2.5),
           #h = ifelse(h_sc==0, 1, .4),
           p = boot::inv.logit(f(ele, hab, mid, scale, b0, b1, b2, b3)), 
           mid_lab = paste0("list(mid == ", mid, ", b[2] ==", scale, ")"), 
           b11_lab = paste0("b[3]~','~b[4] == ", b2), 
           size = ifelse(hab==1, 1, 1.01), 
           x_nudge = ifelse(hab==1, x, x+.05))

p1 <- ggplot(df, aes(ele, p, lty=factor(hab))) +
    #annotate(geom="rect", xmin=-1, xmax=1, ymin=-Inf, ymax=Inf, alpha=.2) +
    geom_line(size = .5) + facet_grid(b11_lab ~ mid_lab, labeller = label_parsed) +
    # theme(strip.text = element_text(face="bold", hjust=0), strip.background = element_blank()) +
    scale_size_continuous(limits = c(1, 1.5)) +
    theme(strip.text = element_text(face="bold", hjust=0), strip.text.y = element_text(angle=0),
          strip.background = element_blank(), 
          axis.text = element_text(colour="black"), legend.position = "bottom", 
          aspect.ratio=1) +
    ylim(0, .7) +
    labs(y = "p(Z=1)", x="Elevational range position", 
         lty="% treecover in 200m") + guides(lty="none")
p1

pg <- ggplotGrob(p1)
for(i in which(grepl("strip-r", pg$layout$name))){
    pg$grobs[[i]]$layout$clip <- "off"
}

grid::grid.draw(pg)
ggsave("figures/example_panel.png")


# df_3eg <- tibble(x = rep(seq(-1.5, 1.5, len=200), 3), 
#        scale = rep(c(.75, 3, 1), each=200), 
#        mid = rep(c(1, 0, -1.3), each=200), 
#        b1 = rep(c(-.1, -.5, -.25), each=200),  
#        b11 = rep(c(2, 0, -.5), each=200), 
#        b0 = rep(c(-1, -1, 0), each=200)) %>%
#     bind_rows(., .) %>%
#     mutate(h_sc = rep(c(0, 1), each = 600), 
#            h = ifelse(h_sc == 0, 1, .4),
#            p = boot::inv.logit(f(x, mid, scale, b1, b11, b0, h=h_sc)), 
#            id = paste0("mid=", mid,"; scale=",scale, "; a0=", b0,"; a1=", b1, "; a2=", b11))
# 
# ggplot(df_3eg, aes(x, p, lty=factor(h, levels=c(1, 0.4)))) +
#     annotate(geom="rect", xmin=-1, xmax=1, ymin=-Inf, ymax=Inf, alpha=.2) +
#     geom_line(size = .5) + 
#     facet_wrap(~id) +
#     theme(strip.text = element_text(face="bold", hjust=0), strip.background = element_blank()) +
#     scale_size_continuous(limits = c(1, 1.5)) +
#     theme(strip.text = element_text(face="bold", hjust=0), strip.text.y = element_text(angle=0),
#           strip.background = element_blank(), 
#           axis.text = element_text(colour="black"), legend.position = "bottom", 
#           aspect.ratio=1) +
#     ylim(0, .7) +
#     labs(y = "p(Z=1)", x="Elevational range position", 
#          lty="% treecover in 200m") 

# ggsave("figures/3examples.png")



# main figure ----
df <- expand.grid(ele = seq(-1, 1, len=100), id=1:3, hab=c(-1,1)) %>% 
    left_join(., tibble(id = 1:3, 
                        mid = rep(0, 3), 
                        scale = rep(4.5, 3), 
                        b0 = rep(-.2, 3), 
                        b1 = rep(.2, 3), 
                        b2 = c(0, 1, 0), 
                        b3 = c(0, 1, 1))) %>%
    mutate(psi_logit = f(ele, hab, mid, scale, b0, b1, b2, b3), 
           psi = boot::inv.logit(psi_logit), 
           hab_f = factor(hab, levels=c(1,-1)))

df_wide <- bind_cols(df[df$hab==1,], psi2 = df$psi[df$hab==-1])


# plots ----
plot_lims <- c(0, max(df_wide$psi/df_wide$psi2))
# plot 1: constant
plot_1 <- df %>%
    filter(id == 1) %>%
    ggplot(aes(ele, psi, lty=hab_f)) + 
    geom_line() +
    scale_x_continuous(expand = expansion(add = .02), breaks = c(-1,0,1), 
                       labels=c("Lower range edge", "Range centre", "Upper range limit")) +
    scale_y_continuous(expand = expansion(add = .01), breaks = c(0.025, .5), 
                       labels=c("Low", "High")) +
    # scale_colour_manual(values=c("black", "indianred")) +
    scale_linetype_manual(values=c(1, 2)) +
    # theme_bw() +
    theme(axis.ticks.length = unit(2, "mm"),
          #panel.grid = element_blank(),
          panel.border = element_rect(colour="black", fill=NA),
          axis.text = element_text(colour="black"), 
          axis.ticks = element_line(size=.5)) +
    guides(lty="none") +
    labs(x = "Elevational range position", 
         y = "Occupancy probability")

plot_1_sub <- df %>%
    filter(id == 1) %>%
    ggplot(aes(ele, psi, lty=hab_f)) + 
    geom_line() +
    scale_x_continuous(expand = expansion(add = .02), breaks=c(-1,0,1)) +
    scale_y_continuous(trans="logit", expand = expansion(mult=.1), breaks=c(0.025, .5)) +
    # scale_colour_manual(values=c("black", "indianred")) +
    scale_linetype_manual(values=c(1, 2)) +
    # theme_bw() +
    theme(panel.grid.minor = element_blank(), 
          panel.border = element_rect(colour="black", fill=NA),
          plot.background = element_blank(),
          axis.ticks =element_blank(), 
          axis.title=element_blank(), 
          #axis.text = element_text(colour="black"), 
          axis.text = element_blank()) +
    guides(lty="none") +
    labs(x = "Elevational range position", 
         y = "Occupancy probability")

plot_1_both <- plot_1 + 
    annotation_custom(
        ggplotGrob(plot_1_sub), 
        xmin = .45, xmax = 1.02, ymax = .51, ymin= .35
    )

# symmetric
plot_2 <- df %>%
    filter(id == 2) %>%
    ggplot(aes(ele, psi, lty=hab_f)) + 
    geom_line() +
    scale_x_continuous(expand = expansion(add = .02), breaks = c(-1,0,1), 
                       labels=c("Lower range edge", "Range centre", "Upper range limit")) +
    scale_y_continuous(expand = expansion(add = .01), breaks = c(0.025, .5), 
                       labels=c("Low", "High")) +
    # scale_colour_manual(values=c("black", "indianred")) +
    scale_linetype_manual(values=c(1, 2)) +
    # theme_bw() +
    theme(axis.ticks.length = unit(2, "mm"),
          #panel.grid = element_blank(),
          panel.border = element_rect(colour="black", fill=NA),
          axis.text = element_text(colour="black"), 
          axis.ticks = element_line(size=.5)) +
    guides(lty="none") +
    labs(x = "Elevational range position", 
         y = "Occupancy probability")

plot_2_sub <- df %>%
    filter(id == 2) %>%
    ggplot(aes(ele, psi, lty=hab_f)) + 
    geom_line() +
    scale_x_continuous(expand = expansion(add = .02), breaks=c(-1,0,1)) +
    scale_y_continuous(trans="logit", expand = expansion(mult=.1), breaks=c(0.025, .5)) +
    # scale_colour_manual(values=c("black", "indianred")) +
    scale_linetype_manual(values=c(1, 2)) +
    # theme_bw() +
    theme(panel.grid.minor = element_blank(), 
          panel.border = element_rect(colour="black", fill=NA),
          plot.background = element_blank(),
          axis.ticks =element_blank(), 
          axis.title=element_blank(), 
          #axis.text = element_text(colour="black"), 
          axis.text = element_blank()) +
    guides(lty="none") +
    labs(x = "Elevational range position", 
         y = "Occupancy probability")

plot_2_both <- plot_2 + 
    annotation_custom(
        ggplotGrob(plot_3_sub), 
        xmin = .45, xmax = 1.02, ymax = .51, ymin= .35
    )

# plot: asymmetric
plot_3 <- df %>%
    filter(id == 3) %>%
    ggplot(aes(ele, psi, lty=hab_f)) + 
    geom_line() +
    scale_x_continuous(expand = expansion(add = .02), breaks = c(-1,0,1), 
                       labels=c("Lower range limit", "Range centre", "Upper range limit")) +
    scale_y_continuous(expand = expansion(add = .01), breaks = c(0.025, .5), 
                       labels=c("Low", "High")) +
    # scale_colour_manual(values=c("black", "indianred")) +
    scale_linetype_manual(values=c(1, 2)) +
    # theme_bw() +
    theme(axis.ticks.length = unit(2, "mm"),
          #panel.grid = element_blank(),
          panel.border = element_rect(colour="black", fill=NA),
          axis.text = element_text(colour="black"), 
          axis.ticks = element_line(size=.5)) +
    guides(lty="none") +
    labs(x = "Elevational range position", 
         y = "Occupancy probability")

plot_3_sub <- df %>%
    filter(id == 3) %>%
    ggplot(aes(ele, psi, lty=hab_f)) + 
    geom_line() +
    scale_x_continuous(expand = expansion(add = .02), breaks=c(-1,0,1)) +
    scale_y_continuous(trans="logit", expand = expansion(mult=.1), breaks=c(0.025, .5)) +
    # scale_colour_manual(values=c("black", "indianred")) +
    scale_linetype_manual(values=c(1, 2)) +
    # theme_bw() +
    theme(panel.grid.minor = element_blank(), 
        panel.border = element_rect(colour="black", fill=NA),
        plot.background = element_blank(),
        axis.ticks =element_blank(), 
        axis.title=element_blank(), 
        #axis.text = element_text(colour="black"), 
        axis.text = element_blank()) +
    guides(lty="none") +
    labs(x = "Elevational range position", 
         y = "Occupancy probability")

plot_3_both <- plot_3 + 
    annotation_custom(
        ggplotGrob(plot_3_sub), 
        xmin = .45, xmax = 1.02, ymax = .51, ymin= .35
    )

## combine into single figure
p_all <- egg::ggarrange(plot_1_both +
                            labs(title = "(a) No differential effect of habitat") +
                            theme(plot.title = element_text(face = "plain", size=10), 
                                  axis.text.x = element_blank(), 
                                  axis.title.x = element_blank(), 
                                  axis.title.y=element_blank(), 
                                  panel.grid.minor=element_blank(), 
                                  plot.margin = unit(c(5,30,5,5), "pt")), 
                        plot_2_both +
                            labs(title = "(b) Differential response: symmetrical")+
                            theme(plot.title = element_text(face = "plain", size=10), 
                                  axis.text.x = element_blank(), 
                                  axis.title.x = element_blank(), 
                                  panel.grid.minor=element_blank()), 
                        plot_3_both +
                            labs(title = "(c) Differential response: asymmetrical") +
                            theme(plot.title = element_text(face = "plain", size=10), 
                                  axis.title.y=element_blank(), 
                                  panel.grid.minor=element_blank()), ncol=1) 


png("figures/schematic_3_examples.png", height=200, width=100, units="mm", res=200)
p_all
dev.off()
