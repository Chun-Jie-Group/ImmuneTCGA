#listgene.hpa_protein.classes.r
#the whole protein.classes consists of 3489 proteins collected from HPA.The human tissue.membrane proteome

library(dplyr)


protein.classes <- read_tsv("./protein_class_Predicted.tsv.gz")

#Subcellular
protein.classes%>%
  filter(`Subcellular location`=="Plasma membrane")->pro_sub.plasm
protein.classes%>%
  filter(is.na(`Subcellular location`))->pro_sub.na
protein.classes%>%
  filter(`Subcellular location`=="Cell Junctions")->pro_sub.celljunc
protein.classes%>%
  filter(grepl("Plasma membrane",`Subcellular location`))->pro_sub.plasm.in

pro_sub.na%>%
  rbind(pro_sub.celljunc)%>%
  rbind(pro_sub.plasm)%>%
  rbind(pro_sub.plasm.in)->pro_sub


#At first,divided by evidence,then Antibody
protein.classes%>%
  filter(Evidence=="No evidence",!is.na(Antibody))%>%
  count( `Reliability (IH)`,`Reliability (Mouse Brain)`,`Reliability (IF)`)



protein.classes%>%
filter(Evidence=="Evidence at transcript level",!is.na(Antibody))%>%
  count( `Reliability (IH)`,`Reliability (Mouse Brain)`,`Reliability (IF)`)


protein.classes %>%
 filter(Evidence=="Evidence at protein level",!is.na(Antibody))%>%
  mutate(Antibody.validation=ifelse(is.na(Antibody),"No.antibody","Yes.antibody")) %>%
  select(Gene,Evidence,Antibody.validation,
         `Reliability (IH)`,`Reliability (Mouse Brain)`,`Reliability (IF)`)->protein.reliability

#Draw a mosaic plot

library(grid)
library(vcd)

protein.reliability$`Reliability (IH)`[is.na(protein.reliability$`Reliability (IH)`)] <- "Non"
protein.reliability$`Reliability (Mouse Brain)`[is.na(protein.reliability$`Reliability (Mouse Brain)`)] <- "Non"
protein.reliability$`Reliability (IF)`[is.na(protein.reliability$`Reliability (IF)`)] <- "Non"

colnames(protein.reliability) <- c("Gene","Evi","An","Ih","Mb","If")

structable(Hi~If+Mb,
           protein.reliability)%>%
  mosaic(shade=T,
       legend=T)

#get genlist from protein.reliability

protein.classes%>%
  filter(`Reliability (IH)`=="Enhanced"|`Reliability (Mouse Brain)`=="Enhanced"|`Reliability (IF)`=="Enhanced")->protein.enhanced

protein.reliability%>%
  filter(Ih=="Uncertain"|Mb=="Uncertain"|If=="Uncertain")->protein.uncertain


protein.classes %>%
  filter(Evidence=="Evidence at protein level",!is.na(Antibody))%>%
anti_join(protein.uncertain)->pro_except.unc

intersect(pro_except.unc,pro_sub)->overlap

