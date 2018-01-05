#listgene2.hpa.r

library(readr)

pro_class<- read_tsv("./protein_class_Predicted.5480.tsv.gz")

#subcellular localization
pro_class%>%
  filter(grepl("Plasma membrane",`Subcellular location`)) %>%
  rbind(filter(pro_class,is.na(`Subcellular location`)))->pro_sub.plasm.in

#reliable
#pro_class %>%
 # filter(`Reliability (IH)`=="Uncertain"|`Reliability (Mouse Brain)`=="Uncertain"|`Reliability (IF)`=="Uncertain")->pro_uncertain

pro_sub.plasm.in %>%
  filter(Evidence=="Evidence at protein level",!is.na(Antibody))->pro_reliable


#locate by cell organs
cell_plasma <- read_tsv("./subcell_location_Plasma.tsv.gz")
cell_junc <- read_tsv("./subcell_location_Cell.tsv.gz")
pro_subcell <- cell_junc%>%
  rbind(cell_plasma)%>%
  unique()

membrane.protein.hpa <- pro_reliable%>%
  rbind(pro_subcell)%>%
  unique()



memb.1 <- pro_reliable%>%
  anti_join(pro_subcell)%>%
  count(`Subcellular location`)
  #filter(grepl("Predicted membrane proteins",`Protein class`))

memb.2 <- pro_subcell %>%
  anti_join(pro_reliable)%>%
  filter(grepl("Predicted membrane proteins",`Protein class`))
  



#test dataset
pro_subcell%>%
  filter(Gene=="ADGRE2")
pro_reliable%>%
  filter(Gene=="ADGRE2")

membrane.protein.hpa%>%
  filter(Gene=="CD44")

pro_class%>%
  filter(Gene=="ADGRE2")->test

genelist.membrane%>%
  filter(Gene=="ADGRE2")

#car-T target&camcidate
#Car_t <- c("CD19","CLEC12A","CD44","CLEC12A","FOLR2","FUT3","CD33","CCR1","CD70","LILRB4","ADGRE2","LILRA2")

