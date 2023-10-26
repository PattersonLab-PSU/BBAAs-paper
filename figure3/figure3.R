library(tidyverse)
library(readxl)
library(ggsci)
library(viridis)
library(RColorBrewer)

theme_clean <- function(){
  theme_bw() %+replace%    #replace elements we want to change
    theme(
      panel.grid = element_blank(),
      strip.background = element_blank()
    )
}

### Data Table
path <- "dataTables/DataTables.xlsx"
fig3_data <- path %>% 
  excel_sheets() %>% 
  set_names() %>%
  map(read_excel, path=path)


### Figure3-a: Shannon Diversity
s_div <- fig3_data$`FIgure3-a`

s_div %>%
  ggplot(aes(x=visit, y=shannon, color=visit, fill=visit)) +
  geom_violin() +
  geom_jitter(width = 0.15, size=0.9) +
  scale_color_viridis(end=0.8, discrete = T) +
  scale_fill_viridis(end=0.8, discrete = T, alpha=0.4) +
  scale_y_continuous(limits=c(0,7)) +
  guides(color="none") +
  theme_clean() +
  theme(
    aspect.ratio=1
  ) +
  labs(
    x="",
    y="Shannon diversity"
  )



## Figure3b: PCOA Plot
theme_clean_pcoa <- function(){
  theme_bw() %+replace%    #replace elements we want to change
    theme(
      axis.text=element_blank(),
      axis.ticks=element_blank(),
      panel.grid = element_blank(),
      aspect.ratio = 1,
      strip.background = element_blank()
    )
}

pcoA_data <- fig3_data$`Figure3-b`   %>%
  mutate(visit=as.factor(visit))

pcoA_data %>%
  ggplot(aes(x=Axis.1,y=Axis.2,color=visit)) +
  geom_point(size=2) +
  theme_clean_pcoa() +
  scale_color_manual(values=viridis(3, end=0.8)[c(2,3)]) +
  theme(
    legend.position = "bottom"
  ) +
  labs(color="")



### Figure 3c: Proportions plot
props_toPlot <- fig3_data$`Figure3-c`
props_toPlot %>%
  #write.table(file="igram_patterson_barplot_1M_12M_data.txt", quote=F, row.names=F, sep='\t')
  mutate(Taxa = factor(Taxa)) %>%
  mutate(Taxa = fct_relevel(Taxa, "Other", after=Inf)) %>%
  mutate(Taxa = fct_rev(Taxa)) %>%
  
  ggplot(aes(x=visit, y=props, fill=Taxa)) +
  geom_histogram(stat="identity") +
  #scale_fill_d3(palette = "category10") +
  scale_fill_manual(values=c("#CFCFCF", pal_nejm()(9))) +
  scale_y_continuous(label=scales:::percent) +
  theme_clean() +
  labs(
    x="", fill="",
    y="Relative abundance"
  )



## Figure3d: BBAAs plot
IGRAM_data <- fig3_data$`Figure3-d`
### Extract meatadata
IGRAM_sample_metadata<-IGRAM_data %>% select(c("Label","Sample","visit","feeding_type")) 
### Bile acids metadata
conjugated_CAs <- c("Ala.CA","Arg.CA","His.CA","Ile.Leu.CA","Phe.CA","Thr.CA","Trp.CA","Tyr.CA","Ser.CA")
conjugated_DCAs <- c("His.DCA","Glu.DCA","Tyr.DCA","Met.DCA","Trp.DCA","Ile.Leu.DCA","Phe.DCA")
conjugated_CDCAs <- c("Glu.CDCA","His.CDCA","Ile.Leu.CDCA","Met.CDCA","Phe.CDCA","Trp.CDCA","Tyr.CDCA")
Tauro_BAs <- c("TCA","THCA","TUDCA","THDCA","TCDCA","TDCA")
Gly_BAs <- c("GCA","GCDCA","GDCA","GHDCA","GLCA")


unconjugated <- c("CDCA","DCA","HDCA","UDCA")
microbial_conjugated_BAs <- c(conjugated_CAs,conjugated_DCAs,conjugated_CDCAs)
host_conjugated <- c(Tauro_BAs,Gly_BAs)
all_BAs <- c(microbial_conjugated_BAs, host_conjugated, unconjugated)


IGRAM_long <- IGRAM_data %>% #mutate(Visit=as.factor(gsub("^.*-","",Label))) %>% 
  gather(key="Bile Acid",value="Concentration",
         -c(Sample,Label,feeding_type, Included, Proxy, Weight,Visits, visit)) %>%
  mutate(`Bile Acid`=gsub("-",".",`Bile Acid`)) %>%
  mutate(
    AAs = case_when(
      `Bile Acid` %in% microbial_conjugated_BAs ~ "BBAAs",
      `Bile Acid` %in% host_conjugated ~ "Host Conjugated",
      `Bile Acid` %in% unconjugated ~ "Unconjugated"
    )
  ) %>%
  mutate(Concentration=as.numeric(Concentration),
         visit=factor(visit,levels=c("0-4D","1M","12M")),
         AAs=factor(AAs,levels=c("Host Conjugated","Unconjugated","BBAAs"))
  )

IGRAM_long %>% select(`Bile Acid`,AAs) %>% unique()



viridis_colors <- c("#440154FF","#2C788F","#5DC863FF")

IGRAM_long %>% filter(Sample!=209) %>% group_by(Label,AAs,visit) %>% 
  summarise(Concentration = sum(Concentration)) %>%
  filter(!is.na(AAs)) %>% 
  ggplot(aes(x=AAs,y=Concentration,fill=visit,color=visit)) +
  #geom_boxplot(alpha=0,outlier.shape=NA,lwd=0.2,width=0.2,position=position_dodge(0.9)) + 
  geom_point(aes(color=visit),size=1,shape=19,
             position=position_jitterdodge(dodge.width = 0.9,jitter.width = 0.25)) +
  geom_violin(alpha=0.5,lwd=0.48,scale = "width") + 
  scale_fill_manual(values=viridis_colors) +
  scale_color_manual(values=viridis_colors) +
  ylab("Bile Acids(nM/mg)") +
  #scale_y_continuous(breaks=c(1e-1,1e+1,1e+3,1e+5),limits=c(1e-1,1e+5)) +
  scale_y_log10(breaks=c(1e-1,1e+1,1e+3,1e+5),limits=c(1e-1,1e+5)) +
  #facet_wrap(~BileType) +
  theme_classic() +
  theme(
    panel.border = element_rect(color="black",fill=NA),
    strip.background = element_blank(),
    axis.line = element_line(linewidth = 0),
    axis.ticks = element_line(linewidth=0.2),
    text = element_text(size=14,colour = "black"),
    axis.title.x = element_blank()
  ) 


### Figure3e: Correlation plot
corr_data <- fig3_data$`Figure3-e`
corr_data %>%
  ggplot(aes(x=Conjugation, y=visit, fill=estimate)) +
  geom_tile(fill = "white", color="black") +
  geom_point(shape = 21, aes(size = abs(estimate)), color="white") +
  geom_text(aes(label=p_label, color=label_color), size=3.5, vjust=0.76) +
  facet_grid(Taxa~., space="free", scales="free") +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  scale_fill_gradient2(low = brewer.pal(11, 'RdBu')[11], high = brewer.pal(11, 'RdBu')[1], mid=brewer.pal(11, 'RdBu')[6], midpoint=0) +
  scale_color_manual(values=c("#D3D3D3", "#000000")) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    strip.text.y = element_text(angle = 0),
    strip.text.x = element_text(angle = 90),
    strip.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA)
  ) +
  #coord_equal() +
  guides(size="none", color="none") +
  labs(
    x="", fill="Correlation\ncoefficient",
    y=""
  )






