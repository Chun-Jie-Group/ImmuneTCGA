#listgene2.hpa.r

library(readr)

pro_class<- read_tsv("./protein_class_Predicted.5480.tsv.gz")

#subcellular localization
pro_class%>%
  filter(grepl("Plasma membrane",`Subcellular location`)) %>%
  rbind(filter(pro_class,is.na(`Subcellular location`)))->pro_sub.plasm.in

#reliable
pro_class %>%
  filter(`Reliability (IH)`=="Uncertain"|`Reliability (Mouse Brain)`=="Uncertain"|`Reliability (IF)`=="Uncertain")->pro_uncertain

pro_sub.plasm.in %>%
  filter(Evidence=="Evidence at protein level",!is.na(Antibody))%>%
  anti_join(pro_uncertain)->pro_reliable

#empirical

pro_reliable %>%
  filter(grepl("ATP",`Gene description`))->pro_atp.related
pro_reliable %>%
  filter(grepl("Carbonic",`Gene description`))->pro_c.related
pro_reliable %>%
  filter(grepl("Calcium",`Gene description`))->pro_ca.related
pro_reliable%>%
  filter(grepl("channel",`Gene description`))->pro_channel.rela
pro_reliable%>%
  filter(grepl("kinase",`Gene description`))->pro_kinase.rela




#car-T target&camcidate
Car_t <- c("CD19","CLEC12A","CD44","CLEC12A","FOLR2","FUT3","CD33","CCR1","CD70","LILRB4","ADGRE2","LILRA2")

