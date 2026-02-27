library(dplyr)
library(ggplot2)
library(stringr)

npeaks = read.table("npeaks.v2.tab")%>%
  mutate(condition = str_extract(.$V2, "AR|antago|ctrl"))%>%
  mutate(condition = str_replace(condition, "antago", "antago+AR"))%>%
  mutate(experiment = str_extract(.$V2, "GSM\\d*(?=\\d{2})"))

t.test(npeaks[npeaks$condition=="ctrl","V1"],
       npeaks[npeaks$condition=="AR","V1"]) # p-value = 0.2664
t.test(npeaks[npeaks$condition=="antago+AR","V1"],
       npeaks[npeaks$condition=="AR","V1"]) # p-value = 0.09845


gsm="GSM13289"
t.test(npeaks[npeaks$condition=="ctrl" & npeaks$experiment==gsm,"V1"],
       npeaks[npeaks$condition=="AR" & npeaks$experiment==gsm,"V1"])
## No 0.3109
t.test(npeaks[npeaks$condition=="antago+AR" & npeaks$experiment==gsm,"V1"],
       npeaks[npeaks$condition=="AR" & npeaks$experiment==gsm,"V1"])
## No 0.3109

gsm="GSM44626"
t.test(npeaks[npeaks$condition=="ctrl" & npeaks$experiment==gsm,"V1"],
       npeaks[npeaks$condition=="AR" & npeaks$experiment==gsm,"V1"])
## Yes 0.01173
# t.test(npeaks[npeaks$condition=="antago+AR" & npeaks$experiment==gsm,"V1"],
#        npeaks[npeaks$condition=="AR" & npeaks$experiment==gsm,"V1"])
# ## not enough obs

gsm="GSM66154"
t.test(npeaks[npeaks$condition=="antago+AR" & npeaks$experiment==gsm,"V1"],
       npeaks[npeaks$condition=="AR" & npeaks$experiment==gsm,"V1"])
## No 0.07919

# gsm="GSM15270"
# t.test(npeaks[npeaks$condition=="ctrl" & npeaks$experiment==gsm,"V1"],
#        npeaks[npeaks$condition=="AR" & npeaks$experiment==gsm,"V1"])
# ## not enough obs

for(gsm in unique(npeaks$experiment)){
  npeaks%>%
    filter(experiment==gsm)%>%
    group_by(condition)%>%
    summarise(sd = sd(V1), V1 = sum(V1), mean = mean(V1))%>%
    ggplot(., aes(y=V1, x=condition))+
      geom_bar(stat = 'identity')+
      geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd))+
      ylab("Number of peaks")+
      ggtitle(paste0("ChIP-seq binding of AR-linked TF: ",gsm,"--"))
  ggsave(paste0("npeaks_",gsm,".png"), width = 6, height = 6)
}

npeaks%>%
  group_by(condition)%>%
  summarise(sd = sd(V1), V1 = sum(V1), mean = mean(V1))%>%
  ggplot(., aes(y=V1, x=condition))+
  geom_bar(stat = 'identity')+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd))+
  ylab("Number of peaks")+
  ggtitle("ChIP-seq binding of AR-linked TF")
ggsave("npeaks.png", width = 6, height = 6)


