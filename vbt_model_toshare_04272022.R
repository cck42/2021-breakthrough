library(tidyverse)
library(lubridate)
library(broom)
library(splines)
library(plyr)
library(dplyr)
library(ggplot2)
library(tibble)
library(lme4)
library(janitor)
library(flextable)
library(cowplot)

data <- read.csv("allsamples_clean_010922.csv") %>%
  as_tibble()

#Collapse sublineage names and add each as a binary variable
data2 <- data %>% 
  mutate(collection_date = yday(as.Date(Collection.date))) %>%
  mutate(vax_to_collection=as.Date(latest_vax.x)-as.Date(Collection.date)) %>%
  mutate(delta=fct_collapse(data$lineagename,'TRUE'='Delta',other_level ='FALSE')) %>%
  mutate(alpha=fct_collapse(data$lineagename,'TRUE'='Alpha',other_level ='FALSE')) %>%
  mutate(bet=fct_collapse(data$lineagename,'TRUE'='Beta',other_level ='FALSE')) %>%
  mutate(gam=fct_collapse(data$lineagename,'TRUE'='Gamma',other_level ='FALSE')) %>%
  mutate(iota=fct_collapse(data$lineagename,'TRUE'='Iota',other_level ='FALSE')) %>%
  mutate(mu=fct_collapse(data$lineagename,'TRUE'='Mu',other_level ='FALSE')) %>%
  mutate(nonvoc=fct_collapse(data$lineagename,'TRUE'=c('Unnamed',"Lambda","Kappa","Eta","Epsilon","Beta"),other_level ='FALSE'))%>%
  mutate(delta=fct_relevel(delta,"FALSE",after=0)) %>%
  mutate(alpha=fct_relevel(alpha,"FALSE",after=0)) %>%
  mutate(bet=fct_relevel(bet,"FALSE",after=0)) %>%
  mutate(gam=fct_relevel(gam,"FALSE",after=0)) %>%
  mutate(iota=fct_relevel(iota,"FALSE",after=0)) %>%
  mutate(mu=fct_relevel(mu,"FALSE",after=0)) %>%
  mutate(nonvoc=fct_relevel(nonvoc,"FALSE",after=0))

#add a column for time since vaccination, grouped
vax_group <- case_when(
  data2$vax_status == "None" ~ "No vax",
  data2$vax_status == "Interim" ~ "Interim",
  is.na(data2$latest_vax.x) ~ "No date info",
  abs(data2$vax_to_collection)<=60 ~ "<60 days",
  60<abs(data2$vax_to_collection) & abs(data2$vax_to_collection)<=120 ~ "61-120 days",
  120<abs(data2$vax_to_collection) & abs(data2$vax_to_collection)<=180 ~ "121-180 days",
  abs(data2$vax_to_collection)>180 ~ "> 180 days",
  TRUE ~ "other"
)
vax_group <- as_tibble(vax_group)
data3 <- tibble(data2,vax_group) %>%
rename(vaxgroup=value)

agegroup <- case_when(
  data$birthyear <1956 ~ ">=65",
  data$birthyear >2005 ~ "<16",
  TRUE ~ ">=16 and <65"
)

analyze_data <- data3 %>%
  bind_cols(age_group=as_factor(agegroup)) %>%
  filter(age_group!="<16") %>%
  filter(vax_status!="Interim")
write.csv(analyze_data,"all_2022updated.csv")

 ##----------------------------------
##Delta numbers
deltacounts <- analyze_data %>%
  ddply(.(delta,vax_status),nrow) %>%
  bind_rows(
    analyze_data %>%
      ddply(.(delta,county),nrow)# %>%
  )%>%
  bind_rows(
    analyze_data %>%
      ddply(.(delta),summarize,
            V1=mean(birthyear,na.rm=T),
            cat="birthyear")
  )%>%
  bind_rows(
    analyze_data %>%
      ddply(.(delta),summarize,
            V1=mean(collection_date,na.rm=T),
            cat="collection date")
  )%>%
  bind_rows(
    analyze_data %>%
      ddply(.(delta,age_group),nrow)
  )%>%
  pivot_wider(values_from = V1,names_from = delta) %>%
  pivot_longer(cols=c(vax_status,county,cat,age_group),names_to = "category",values_drop_na = TRUE)


##----------------------------------
##Delta model
## df to use = 3
deltamod <- analyze_data %>%
  glm(formula=delta ~ vax_status + age_group + county + bs(collection_date,df=3),family="binomial"(link="logit"),model=TRUE)
deltamod_unadjusted <- analyze_data %>%
  glm(formula=delta ~ vax_status,family="binomial"(link="logit"),model=TRUE)%>% 
  tidy(exponentiate=TRUE,conf.int=TRUE) %>%
  mutate(across(where(is.numeric),round,digits=2))
deltanums <- deltamod %>% tidy(exponentiate=TRUE,conf.int=TRUE) %>%
  mutate(across(where(is.numeric),round,digits=2))
deltamodlinear <- analyze_data %>%
  glm(formula=delta ~ vax_status + birthyear + county + collection_date,family="binomial"(link="logit"),model=TRUE)
deltanumslinear <- deltamodlinear %>% tidy(exponentiate=TRUE,conf.int=TRUE) %>%
  mutate(across(where(is.numeric),round,digits=2))
deltamodplot <- analyze_data %>%
  filter(!(X %in% deltamodlinear$na.action)) %>%
  ggplot() +
  stat_count(aes(x=collection_date,fill=delta),position="fill")+
  geom_line(aes(x=collection_date,y=deltamodlinear$fitted.values))
plot(deltamodplot)

##----------------------------------
##alpha numbers
alphanums <- analyze_data %>%
  ddply(.(alpha,vax_status),nrow) %>%
  bind_rows(
    analyze_data %>%
      ddply(.(alpha,county),nrow)
  )%>%
  bind_rows(
    analyze_data %>%
      ddply(.(alpha),summarize,
            V1=mean(birthyear,na.rm=T),
            cat="birthyear")
  )%>%
  bind_rows(
    analyze_data %>%
      ddply(.(alpha),summarize,
            V1=mean(collection_date,na.rm=T),
            cat="collection date")
  )%>%
  bind_rows(
    analyze_data %>%
      ddply(.(alpha,age_group),nrow)
  ) %>%
  pivot_wider(values_from = V1,names_from = alpha) %>%
  pivot_longer(cols=c(vax_status,county,cat,age_group),names_to = "category",values_drop_na = TRUE)

## ----------------------------------------------------------
## Alpha model
##DF to use = 4
alphamod <- analyze_data %>%
  glm(formula=alpha ~ vax_status + age_group + county + bs(collection_date,df=4),family="binomial"(link="logit"))
alphamod_unadjusted <- analyze_data %>%
  glm(formula=alpha ~ vax_status,family="binomial"(link="logit"),model=TRUE)%>% 
  tidy(exponentiate=TRUE,conf.int=TRUE) %>%
  mutate(across(where(is.numeric),round,digits=2))
alphamodnums<- alphamod %>% tidy(exponentiate=TRUE,conf.int=TRUE)
alphaplot <- analyze_data %>%
  filter(!(X %in% alphamod$na.action)) %>%
  ggplot() +
  stat_count(aes(x=collection_date,fill=alpha),position="fill")+
  geom_line(aes(x=collection_date, y=alphamod$fitted.values))
plot(alphaplot)

## ----------------------------------------------------------
## Gamma numbers
gammacounts <- analyze_data %>%
  ddply(.(gam,vax_status),nrow) %>%
  bind_rows(
    analyze_data %>%
      ddply(.(gam,county),nrow)
  )%>%
  bind_rows(
    analyze_data %>%
      ddply(.(gam),summarize,
            V1=mean(birthyear,na.rm=T),
            cat="birthyear")
  )%>%
  bind_rows(
    analyze_data %>%
      ddply(.(gam),summarize,
            V1=mean(collection_date,na.rm=T),
            cat="collection date")
  )%>%
  bind_rows(
    analyze_data %>%
      ddply(.(gam,age_group),nrow)
  ) %>%
  pivot_wider(values_from = V1,names_from = gam) %>%
  pivot_longer(cols=c(vax_status,county,cat,age_group),names_to = "category",values_drop_na = TRUE)

## ----------------------------------------------------------
## gamma model
#df to use = 3
gammamod <- analyze_data %>%
  glm(formula=gam ~ vax_status + age_group + county + bs(collection_date,df=3),family="binomial"(link="logit"))
gammamod_unadjusted <- analyze_data %>%
  glm(formula=gam ~ vax_status,family="binomial"(link="logit"),model=TRUE)%>% 
  tidy(exponentiate=TRUE,conf.int=TRUE) %>%
  mutate(across(where(is.numeric),round,digits=2))
gammamodnums <- gammamod %>% tidy(exponentiate=TRUE,conf.int=TRUE)
gammaplot <- analyze_data %>%
  filter(!(X %in% gammamod$na.action)) %>%
  ggplot() +
  stat_count(aes(x=collection_date,fill=gam),position="fill")+
  geom_line(aes(x=collection_date,y=gammamod$fitted.values))
plot(gammaplot)

## ----------------------------------------------------------
## Iota numbers
iotacounts <- analyze_data %>%
  ddply(.(iota,vax_status),nrow) %>%
  bind_rows(
    analyze_data %>%
      ddply(.(iota,county),nrow)
  )%>%
  bind_rows(
    analyze_data %>%
      ddply(.(iota),summarize,
            V1=mean(birthyear,na.rm=T),
            cat="birthyear")
  )%>%
  bind_rows(
    analyze_data %>%
      ddply(.(iota),summarize,
            V1=mean(collection_date,na.rm=T),
            cat="collection date")
  )%>%
  bind_rows(
    analyze_data %>%
      ddply(.(iota,age_group),nrow)
  ) %>%
  pivot_wider(values_from = V1,names_from = iota) %>%
  pivot_longer(cols=c(vax_status,county,cat,age_group),names_to = "category",values_drop_na = TRUE)

## ----------------------------------------------------------
## Iota model
## df to use = 4
iotamod <- analyze_data %>%
  glm(formula=iota ~ vax_status + age_group + county + bs(collection_date,df=4),family="binomial"(link="logit"))
iotamod_unadjusted <- analyze_data %>%
  glm(formula=iota ~ vax_status,family="binomial"(link="logit"),model=TRUE)%>% 
  tidy(exponentiate=TRUE,conf.int=TRUE) %>%
  mutate(across(where(is.numeric),round,digits=2))
iotamodnums <- iotamod %>% tidy(exponentiate=TRUE,conf.int=TRUE)
iotaplot <- analyze_data %>%
  filter(!(X %in% iotamod$na.action)) %>%
  ggplot() +
  stat_count(aes(x=collection_date,fill=iota),position="fill")+
  geom_line(aes(x=collection_date,y=iotamod$fitted.values))
plot(iotaplot)

## ----------------------------------------------------------
## Mu numbers
mucounts <- analyze_data %>%
  ddply(.(mu,vax_status),nrow) %>%
  bind_rows(
    analyze_data %>%
      ddply(.(mu,county),nrow)
  )%>%
  bind_rows(
    analyze_data %>%
      ddply(.(mu),summarize,
            V1=mean(birthyear,na.rm=T),
            cat="birthyear")
  )%>%
  bind_rows(
    analyze_data %>%
      ddply(.(mu),summarize,
            V1=mean(collection_date,na.rm=T),
            cat="collection date")
  )%>%
  bind_rows(
    analyze_data %>%
      ddply(.(mu,age_group),nrow)
  ) %>%
  pivot_wider(values_from = V1,names_from = mu) %>%
  pivot_longer(cols=c(vax_status,county,cat,age_group),names_to = "category",values_drop_na = TRUE)

## ----------------------------------------------------------
## Mu model
## df to use = 3
mumod <- analyze_data %>%
  glm(formula=mu ~ vax_status + age_group + county + bs(collection_date,df=3),family="binomial"(link="logit"))
mumod_unadjusted <- analyze_data %>%
  glm(formula=mu ~ vax_status,family="binomial"(link="logit"),model=TRUE)%>% 
  tidy(exponentiate=TRUE,conf.int=TRUE) %>%
  mutate(across(where(is.numeric),round,digits=2))
mumodnums <- mumod %>% tidy(exponentiate=TRUE,conf.int=TRUE)
muplot <- analyze_data %>%
  filter(!(X %in% mumod$na.action)) %>%
  ggplot() +
  stat_count(aes(x=collection_date,fill=mu),position="fill")+
  geom_line(aes(x=collection_date,y=mumod$fitted.values))
plot(muplot)

## ----------------------------------------------------------
## VOC vs. non-VOC numbers
nonvoccounts <- analyze_data %>%
  ddply(.(nonvoc,vax_status),nrow) %>%
  bind_rows(
    analyze_data %>%
      ddply(.(nonvoc,county),nrow)
  )%>%
  bind_rows(
    analyze_data %>%
      ddply(.(nonvoc),summarize,
            V1=mean(birthyear,na.rm=T),
            cat="birthyear")
  )%>%
  bind_rows(
    analyze_data %>%
      ddply(.(nonvoc),summarize,
            V1=mean(collection_date,na.rm=T),
            cat="collection date")
  )%>%
  bind_rows(
    analyze_data %>%
      ddply(.(nonvoc,age_group),nrow)
  ) %>%
  pivot_wider(values_from = V1,names_from = nonvoc) %>%
  pivot_longer(cols=c(vax_status,county,cat,age_group),names_to = "category",values_drop_na = TRUE)

## ----------------------------------------------------------
## All VOC vs. Non-VOC model
## df to use = 4
allvocmodel <- analyze_data %>%
  glm(formula=nonvoc ~ vax_status + age_group + county + bs(collection_date,df=4),family="binomial"(link="logit"))
allvocmod_unadjusted <- analyze_data %>%
  glm(formula=nonvoc ~ vax_status,family="binomial"(link="logit"),model=TRUE)%>% 
  tidy(exponentiate=TRUE,conf.int=TRUE) %>%
  mutate(across(where(is.numeric),round,digits=2))
allvocmodelnums <- allvocmodel %>% tidy(exponentiate=TRUE,conf.int=TRUE)
allvocplot <- analyze_data %>%
  filter(!(X %in% allvocmodel$na.action)) %>%
  ggplot() +
  stat_count(aes(x=collection_date,fill=allvoc),position="fill")+
  geom_line(aes(x=collection_date,y=allvocmodel$fitted.values))
plot(allvocplot)

## ----------------------------------------------------------
## Bind model results together

modelnums <- as_tibble(alphamodnums) %>%
  bind_rows(list(deltanums,gammamodnums,iotamodnums,mumodnums,allvocmodelnums),.id="lineage") %>%
  mutate(lineage=fct_recode(lineage,"Alpha"="1","Delta" = "2","Gamma"="3","Iota"="4","Mu"="5","Non-VOC" = "6"))
modelnums_unadjusted<-as_tibble(alphamod_unadjusted) %>%
  bind_rows(list(deltamod_unadjusted,gammamod_unadjusted,iotamod_unadjusted,mumod_unadjusted,allvocmod_unadjusted),.id="lineage") %>%
  mutate(lineage=fct_recode(lineage,"Alpha"="1","Delta" = "2","Gamma"="3","Iota"="4","Mu"="5","Non-VOC" = "6"))

## ----------------------------------------------------------
## Bind counts together

allcounts <- as_tibble(alphanums) %>%
  bind_rows(list(deltacounts,gammacounts, iotacounts,mucounts,nonvoccounts),.id="lineage")%>%
  mutate(lineage=fct_recode(lineage,"Alpha"="1","Delta" = "2","Gamma"="3","Iota"="4","Mu"="5","Non-VOC" = "6")) %>%
  mutate(name=str_c(category,value)) %>%
  mutate(total=`FALSE` + `TRUE`) %>%
  left_join(modelnums,by=c("lineage"="lineage",name="term")) %>%
  left_join(modelnums_unadjusted,by=c("lineage"="lineage",name="term")) %>%
  filter(category=="vax_status")%>%
  #select(!(FALSE|category))
  select(lineage,"value","TRUE","estimate.y",p.value.y,conf.low.y,conf.high.y,"estimate.x",p.value.x,conf.low.x,conf.high.x)

#Just vaccination status numbers
vaxnums <- as_tibble(allcounts) %>%
  #filter(category=="vax_status") %>%
  mutate(estimate.y=round(estimate.y,digits = 2))%>%
  mutate(p.value.y=round(p.value.y,digits = 2))%>%
  mutate(conf.low.y=round(conf.low.y,digits = 2))%>%
  mutate(conf.high.y=round(conf.high.y,digits = 2))%>%
  mutate(estimate.x=round(estimate.x,digits = 2))%>%
  mutate(p.value.x=round(p.value.x,digits = 2))%>%
  mutate(conf.low.x=round(conf.low.x,digits = 2))%>%
  mutate(conf.high.x=round(conf.high.x,digits = 2))%>%
  #select(!category) %>%
  qflextable() %>%
  add_header_row(
    top=TRUE,
    values = c("Lineage",
               "Vaccination status",
               "N",
               "Unadjusted odds ratio of infection with this variant vs. all others",
               "",
               "95% CI",
               "",
               "Adjusted odds ratio of infection with this variant vs. all others",
               "",
               "95% CI",
               ""
    )) %>%
  set_header_labels(
    lineage = "",
    "value" = "",
    "TRUE" = "",
    "estimate.y"="Estimate",
    p.value.y="P-value",
    conf.low.y="Lower",
    conf.high.y = "Upper",
    "estimate.x"="Estimate",
    p.value.x="P-value",
    conf.low.x="Lower",
    conf.high.x = "Upper"
  ) %>%
  merge_at(i=1,j=4:5,part="header")%>%
  merge_at(i=1,j=6:7,part="header")%>%
  merge_at(i=1,j=8:9,part="header")%>%
  merge_at(i=1,j=10:11,part="header")%>%
  merge_at(i=1:2,j=1,part="header")%>%
  merge_at(i=1:2,j=2,part="header")%>%
  merge_at(i=1:2,j=3,part="header")%>%
  merge_at(i=1:2,j=1,part="body") %>%
  merge_at(i=3:4,j=1,part="body") %>%
  merge_at(i=5:6,j=1,part="body") %>%
  merge_at(i=7:8,j=1,part="body") %>%
  merge_at(i=9:10,j=1,part="body") %>%
  merge_at(i=11:12,j=1,part="body") %>%
  border_remove() %>%
  theme_booktabs()%>%
  align(align="center",part="all") %>%
  bold(i=1,bold=TRUE,part="header") %>%
  vline(part="all",j=3)%>%
  vline(part="all",j=7)%>%
  font(part="all",fontname="times new roman") %>%
  fontsize(part="all",size=10) %>%
  hline(i=2,part="body")%>%
  hline(i=4,part="body")%>%
  hline(i=6,part="body")%>%
  hline(i=8,part="body")%>%
  hline(i=10,part="body")%>%
  fix_border_issues()
  
## Save this table in a doc
save_as_pptx(vaxnums,path='lineagemodel_010922.pptx')
  

## ----------------------------------------------------------
## Make some figures

# Vaccine status x lineage
colorvals=c(Alpha="#ed4d18", "Delta" = "#09769A", "Gamma"="#BF0D40","Iota"= "#B18195","Mu"="#958ECD","Non-VOC" ="#7C7C7C")
categoricals=c("vax_statusInterim","vax_statusPostseries")
catplots <- modelnums %>%
  filter(term %in% categoricals) %>%
  ggplot() +
  geom_hline(aes(yintercept=1, color="grey"),linetype="dotted")+
  geom_pointrange(aes(x=lineage,y=estimate,ymax=conf.high,ymin=conf.low,color=lineage)) +
  scale_color_manual(values=colorvals)+
  labs(x="Lineage",y="Adjusted Odds Ratio",title="Odds ratio of infection with each lineage\n in cases among fully vaccinated \nvs. never vaccinated individuals")+
  theme_cowplot(12)+
  theme(legend.position="none")
plot(catplots)

#Breakdown of variant frequencies by week
Variant_freq_by_week <- analyze_data %>%
  mutate(lineagename=fct_collapse(lineagename,
    "Alpha" = "Alpha","Delta"="Delta","Gamma"="Gamma",
    "Iota"="Iota","Mu"="Mu",other_level = "Non-VOC"
  ))%>%
  mutate(col_date=as_date(Collection.date)) %>%
  ggplot() +
  aes(x=col_date,fill=lineagename) +
  geom_area(stat="bin",position="fill",bins=33)+
  facet_wrap(vars(vax_status),nrow=3,scales="free_y")+
  labs(y="Frequency",x="Date positive specimen collected",title="Frequency of variants over time, by vaccination status") +
  scale_fill_manual(values=colorvals,name="Lineage")+
  scale_x_date(date_breaks="1 months",date_labels="%b %d")+
  theme_cowplot(12)+
  theme(legend.position="none")+
  scale_y_continuous(labels=scales::percent)
plot(Variant_freq_by_week)

variants_by_vax <- analyze_data %>%
  mutate(lineagename=fct_collapse(lineagename,
    "Alpha" = "Alpha","Delta"="Delta","Gamma"="Gamma",
    "Iota"="Iota","Mu"="Mu",other_level = "Non-VOC"
  )) %>%
  ggplot() +
  aes(x=vax_status,fill=lineagename) +
  geom_bar(stat="count",position="fill") +
  labs(y="Frequency",x="Vaccination status",title="Variant frequencies by \n vaccination status") +
  scale_fill_manual(values=colorvals,name="Lineage")+
  theme_cowplot(12)+
  theme(legend.position="left")
plot(variants_by_vax)

Variants_by_week <- analyze_data %>%
  mutate(lineagename=fct_collapse(lineagename,
    "Alpha" = "Alpha","Delta"="Delta","Gamma"="Gamma",
    "Iota"="Iota","Mu"="Mu",other_level = "Non-VOC"
  )) %>%
  mutate(col_date=as_date(Collection.date)) %>%
  ggplot() +
  aes(x=col_date,fill=lineagename) +
  geom_histogram(binwidth=7,position="stack") +
  labs(y="Number of cases",x="Date positive specimen collected",title="Number of cases over time, by variant") +
  scale_fill_manual(values=colorvals,name="Lineage")+
  scale_x_date(date_breaks="1 months",date_labels="%b %d")+
  theme_cowplot(12)+
  theme(legend.position="none")
plot(Variants_by_week)

#Combine these plots

fig1a <- align_plots(Variants_by_week,Variant_freq_by_week,align = "h",axis="lr")
ggdraw(fig1a[[1]]) + draw_plot(fig1a[[2]])
fig1bc <- plot_grid(Variants_by_week,Variant_freq_by_week,labels=c('C','D'),nrow=2,rel_heights = c(.4,.6))
fig1bc
fig1b <- plot_grid(variants_by_vax,catplots,labels=c('A','B'),rel_widths = c(.4,.6))
fig1b
fig1 <- plot_grid(fig1b,fig1bc,labels=c('',''),nrow=2,rel_heights = c(.3,.6))
fig1
