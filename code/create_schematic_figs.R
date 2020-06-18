library(dplyr); library(ggplot2)

f <- function(x, mid, scale, b1, b11, b0, h) {
    (scale + b11*h) * -(x-mid)^2 + .5 + b1*h + b0
}

df <- expand.grid(x = seq(-1.2, 1.2, len=200), 
            scale = .5, 
            mid = c(-1.2, 0, 1.2), 
            b1 = 0, 
            b11 = c(-1, 0, 2), 
            h_sc=c(0,1)) %>%
    as_tibble() %>%
    mutate(scale = ifelse(mid == 0, 3, 1.5),
           h = ifelse(h_sc==0, 1, .4),
           p = boot::inv.logit(f(x, mid, scale, b1, b11, 0, h_sc)), 
           mid_lab = paste0("list(mid == ", mid, ", scale ==", scale, ")"), 
           b11_lab = paste0("alpha == ", b11), 
           size = ifelse(h==1, 1, 1.01), 
           x_nudge = ifelse(h==1, x, x+.05))

p1 <- ggplot(df, aes(x, p, lty=factor(h, levels=c(1, .4)))) +
    annotate(geom="rect", xmin=-1, xmax=1, ymin=-Inf, ymax=Inf, alpha=.2) +
    geom_line(size = .5) + facet_grid(b11_lab ~ mid_lab, labeller = label_parsed) +
    theme(strip.text = element_text(face="bold", hjust=0), strip.background = element_blank()) +
    scale_size_continuous(limits = c(1, 1.5)) +
    theme(strip.text = element_text(face="bold", hjust=0), strip.text.y = element_text(angle=0),
          strip.background = element_blank(), 
          axis.text = element_text(colour="black"), legend.position = "bottom", 
          aspect.ratio=1) +
    ylim(0, .7) +
    labs(y = "p(Z=1)", x="Elevational range position", 
         lty="% treecover in 200m") #+ guides(fill=FALSE)

p1
pg <- ggplotGrob(p1)
for(i in which(grepl("strip-r", pg$layout$name))){
    pg$grobs[[i]]$layout$clip <- "off"
}

grid::grid.draw(pg)
ggsave("figures/example_panel.png")


df_3eg <- tibble(x = rep(seq(-1.5, 1.5, len=200), 3), 
       scale = rep(c(.75, 3, 1), each=200), 
       mid = rep(c(1, 0, -1.3), each=200), 
       b1 = rep(c(-.1, -.5, -.25), each=200),  
       b11 = rep(c(2, 0, -.5), each=200), 
       b0 = rep(c(-1, -1, 0), each=200)) %>%
    bind_rows(., .) %>%
    mutate(h_sc = rep(c(0, 1), each = 600), 
           h = ifelse(h_sc == 0, 1, .4),
           p = boot::inv.logit(f(x, mid, scale, b1, b11, b0, h=h_sc)), 
           id = paste0("mid=", mid,"; scale=",scale, "; a0=", b0,"; a1=", b1, "; a2=", b11))

ggplot(df_3eg, aes(x, p, lty=factor(h, levels=c(1, 0.4)))) +
    annotate(geom="rect", xmin=-1, xmax=1, ymin=-Inf, ymax=Inf, alpha=.2) +
    geom_line(size = .5) + 
    facet_wrap(~id) +
    theme(strip.text = element_text(face="bold", hjust=0), strip.background = element_blank()) +
    scale_size_continuous(limits = c(1, 1.5)) +
    theme(strip.text = element_text(face="bold", hjust=0), strip.text.y = element_text(angle=0),
          strip.background = element_blank(), 
          axis.text = element_text(colour="black"), legend.position = "bottom", 
          aspect.ratio=1) +
    ylim(0, .7) +
    labs(y = "p(Z=1)", x="Elevational range position", 
         lty="% treecover in 200m") 

ggsave("figures/3examples.png")
