suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))

df = fread("chicken_specific_egene_nonegene_GO_BP.txt")
p=ggplot(df %>% filter(Type == "eGene"), aes(x=reorder(term_name, negative_log10_of_adjusted_p_value), y=negative_log10_of_adjusted_p_value))+
    geom_bar(fill="#F4A582",stat = "identity")+
    xlab("")+
    ylab(expression("-log"[10]*italic(P)))+
    coord_flip()+
    #scale_x_discrete(position = "top")+
    theme_classic(base_line_size = 0.4)+
    theme(legend.title = element_blank(), legend.position="none",
          axis.text.x = element_text(color="black", size=8),
          axis.title.x = element_text(color="black", size=10),
          axis.text.y = element_text(color="black", size=8),
          axis.title.y = element_text(color="black", size=10),
          #strip.background = element_rect(color="black", fill="#C0C0C0", size=0.2, linetype="solid"),
          axis.ticks= element_line(color="black"),
          axis.line= element_line(color="black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
ggsave("chicken_specific_egene_go.pdf", p, width=4, height=2)
