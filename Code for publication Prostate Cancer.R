#Load Packages ####
require("pacman")
require("dplyr")
require("ggplot2")
require("data.table")
require("bit64")
require("lubridate")
require("survival")
require("survminer")
require("Rmisc")
require("RODBC")
require("bit64")

#Import Data ####
#Import Preprocessed Cancer Data
MasterTab<- fread("All.csv")


#Identification of Diabetic PAtients Using Each method in Prostate Cancer

Comorb_Prostate<- MasterTab %>% dplyr::filter(grepl("C61", ICD10CDS)==T)
#Preprocess All Cancer Data ####

#Prediagnosis

Diabetic_Prostate_Tab<- Comorb_Prostate %>% 
  dplyr::mutate( HESDiabetic = ifelse(T1Diabetes == 1|
                                        T2Diabetes ==1|
                                        OtherDiabetes == 1, 1, 0),
                 BloodDiabetic = Diagnostic)%>% # identification of patients with pre cancer clinical coding for diabetes
  dplyr::mutate(Diagnosis_Year = as.numeric(strtrim(as.character(DiagnosisDate), width = 4))) %>%
  dplyr::filter(Diagnosis_Year>2004) %>% # remove patients diagnosed before 2005 as bloods available from 2004 onwards
  dplyr::arrange(StudyId, DiagnosisDate)%>%
  dplyr::filter(!duplicated(`StudyId`)) # select only the first cancer event for the patients in the cohort

BloodNumber <- sum(Diabetic_Prostate_Tab$BloodDiabetic)#Number of prediagnosis Diagnostic HBA1c
HesNumber<- sum(Diabetic_Prostate_Tab$HESDiabetic)# Number of coded diabetics
AllNumber<- sum(Diabetic_Prostate_Tab$AnyDiabetes) # Number of diabetics including HBA1c and Coding
BothNumber<- nrow(dplyr::filter(Diabetic_Prostate_Tab, HESDiabetic ==1 & BloodDiabetic ==1)) # Number with both blood and coding
AllCancer<- nrow(Diabetic_Prostate_Tab) # Number of PAtients in total

#create summary table for patient numbers pre cancer diagnosis
DiabeticSumTab<- as.data.frame(cbind(c("Non-Diabetic", "Coded as Diabetic & Diagnostic HBA1c","Clinical Coding Only", "Abnormal HBA1c Only",
                                       "Hybrid Definition", "Total Patients"),
                                     c(AllCancer-AllNumber, BothNumber, HesNumber, BloodNumber,  AllNumber, AllCancer )))

colnames(DiabeticSumTab)<- c("Group", "NumberOfPatients")

DiabeticSumTab<- DiabeticSumTab%>% dplyr::mutate(Percentage = round((as.numeric(as.character(NumberOfPatients))/AllCancer)*100, digits = 2), # include percentages
                                                 TimePoint = "Before Cancer Diagnosis") # add time label


# Post Diangosis
#Identification of Diabetic PAtients Using Each method in Prostate Cancer

Diabetic_Prostate_Tab<- Diabetic_Prostate_Tab %>% dplyr::mutate( Late.HESDiabetic = ifelse( Late.T1Diabetes == 1|
                                                                                              Late.T2Diabetes ==1|
                                                                                              Late.OtherDiabetes == 1, 1, 0),
                                                                 Late.BloodDiabetic = Late.Diagnostic) # add label for late effect clinical coding and bloods

Late.BloodNumber <- sum(Diabetic_Prostate_Tab$Late.BloodDiabetic, na.rm = T) #number of late blood diabetics
Late.HesNumber<- sum(Diabetic_Prostate_Tab$Late.HESDiabetic, na.rm = T) # number of late coded diabetics
Late.AllNumber<- sum(Diabetic_Prostate_Tab$Late.AnyDiabetes) # total number of late diabetics
Late.BothNumber<- nrow(dplyr::filter(Diabetic_Prostate_Tab, Late.HESDiabetic ==1 & Late.BloodDiabetic ==1)) # number of diabetics with both
AllCancer<- nrow(Diabetic_Prostate_Tab) # total number of patients

#summary Table for late diabetes numbers
Late.DiabeticSumTab<- as.data.frame(cbind(c("Non-Diabetic", "Coded as Diabetic & Diagnostic HBA1c","Clinical Coding Only", "Abnormal HBA1c Only",
                                            "Hybrid Definition", "Total Patients"),
                                          c(AllCancer-Late.AllNumber, Late.BothNumber, Late.HesNumber, Late.BloodNumber,  Late.AllNumber, AllCancer )))

colnames(Late.DiabeticSumTab)<- c("Group", "NumberOfPatients")

Late.DiabeticSumTab<- Late.DiabeticSumTab%>% dplyr::mutate(Percentage = round((as.numeric(as.character(NumberOfPatients))/AllCancer)*100, digits = 2), # add percentage
                                                           TimePoint = "After Cancer Diagnosis") # add time label




#Any Time Point repeat of above but for any time before or after cancer diagnosis
Diabetic_Prostate_Tab<- Diabetic_Prostate_Tab %>% dplyr::mutate( TotalHESDiabetic = ifelse( T1Diabetes == 1|
                                                                                              T2Diabetes ==1|
                                                                                              OtherDiabetes == 1 |
                                                                                              Late.T1Diabetes ==1| 
                                                                                              Late.T2Diabetes ==1|
                                                                                              Late.OtherDiabetes ==1, 1, 0),
                                                                 TotalBloodDiabetic = ifelse(Diagnostic ==1 |
                                                                                               Late.Diagnostic ==1, 1,0),
                                                                 TotalAnyDiabetes = ifelse(AnyDiabetes ==1|
                                                                                             Late.AnyDiabetes==1, 1, 0),
                                                                 Late.Diagnostic =  ifelse(is.na(Late.Diagnostic) == T, 0, Late.Diagnostic))

TotalBloodNumber <- sum(Diabetic_Prostate_Tab$TotalBloodDiabetic, na.rm = T)
TotalHesNumber<- sum(Diabetic_Prostate_Tab$TotalHESDiabetic)
TotalAllNumber<- sum(Diabetic_Prostate_Tab$TotalAnyDiabetes)
TotalBothNumber<- nrow(dplyr::filter(Diabetic_Prostate_Tab, TotalHESDiabetic ==1 & TotalBloodDiabetic ==1))


TotalDiabeticSumTab<- as.data.frame(cbind(c("Non-Diabetic", "Coded as Diabetic & Diagnostic HBA1c","Clinical Coding Only", "Abnormal HBA1c Only",
                                            "Hybrid Definition", "Total Patients"),
                                          c(AllCancer-TotalAllNumber, TotalBothNumber,TotalHesNumber, TotalBloodNumber,  TotalAllNumber, AllCancer )))

colnames(TotalDiabeticSumTab)<- c("Group", "NumberOfPatients")

TotalDiabeticSumTab<- TotalDiabeticSumTab%>% dplyr::mutate(Percentage = round((as.numeric(as.character(NumberOfPatients))/AllCancer)*100, digits = 2),
                                                           TimePoint = "Before and After Cancer Diagnosis")



#Merge Summary Tables to creat all numbers before and after in one table

OverallDiabSumTab<- rbind(DiabeticSumTab[3:5,], Late.DiabeticSumTab[3:5,],TotalDiabeticSumTab[3:5,]) %>% 
  dplyr::mutate(NumberOfPatients = as.numeric(as.character(NumberOfPatients)) ) %>%
  dplyr::mutate(TimePoint = factor(TimePoint, levels = c("Before Cancer Diagnosis", "After Cancer Diagnosis", 
                                                         "Before and After Cancer Diagnosis")))
# Plot PAtient numbers and percentages by time point
Numbers_Bar_plot<- ggplot(OverallDiabSumTab) + 
  geom_bar(aes(y = NumberOfPatients, x = reorder(Group, NumberOfPatients), fill = TimePoint),  
           stat = "identity") +
  theme(text = element_text(size=30),
        axis.text.x = element_text(angle = 90), legend.position = "none") +
  ggtitle("Comparison of Diabetic Prostate Cancer Patient Numbers by Method of Identification and Time Period") +
  xlab("Method of Identification") + ylab("Number of Patients") +
  facet_wrap( ~TimePoint)

Perc_Bar_plot<- ggplot(OverallDiabSumTab) + 
  geom_bar(aes(y = Percentage, x = reorder(Group, Percentage), fill = TimePoint),  
           stat = "identity") +
  theme(text = element_text(size=30),
        axis.text.x = element_text(angle = 90), legend.position = "none") +
  ggtitle("Comparison of Diabetic Prostate Cancer Patient Percentage by Method of Identification and Time Period") +
  xlab("Method of Identification") + ylab("Percentage of Patients") +
  facet_wrap( ~TimePoint)



# Create Tables for patient numbers by method of identification 

AllBloodOnly<- nrow(dplyr::filter(Diabetic_Prostate_Tab, (Diagnostic == 1 | Late.Diagnostic == 1) &
                                    T1Diabetes == 0 & T2Diabetes ==0 & OtherDiabetes ==0 &
                                    Late.T1Diabetes ==0 & Late.T2Diabetes ==0 & Late.OtherDiabetes ==0))

AllCodedOnly<- nrow(dplyr::filter(Diabetic_Prostate_Tab,(T1Diabetes == 1 | T2Diabetes ==1 | OtherDiabetes ==1 |
                                                           Late.T1Diabetes ==1 | Late.T2Diabetes ==1 | Late.OtherDiabetes ==1) &
                                    Diagnostic == 0 & Late.Diagnostic == 0 ))

AllBoth <- nrow(dplyr::filter(Diabetic_Prostate_Tab,(T1Diabetes == 1 | T2Diabetes ==1 | OtherDiabetes ==1 |
                                                       Late.T1Diabetes ==1 | Late.T2Diabetes ==1 | Late.OtherDiabetes ==1) &
                                (Diagnostic == 1 | Late.Diagnostic == 1 )))


PreBloodOnly<- nrow(dplyr::filter(Diabetic_Prostate_Tab, (Diagnostic == 1 ) &
                                    T1Diabetes == 0 & T2Diabetes ==0 & OtherDiabetes ==0))

PreCodedOnly<- nrow(dplyr::filter(Diabetic_Prostate_Tab,(T1Diabetes == 1 | T2Diabetes ==1 | OtherDiabetes ==1 ) &
                                    Diagnostic == 0 ))

PreBoth <- nrow(dplyr::filter(Diabetic_Prostate_Tab,(T1Diabetes == 1 | T2Diabetes ==1 | OtherDiabetes ==1 ) &
                                (Diagnostic == 1 )))


PostBloodOnly<- nrow(dplyr::filter(Diabetic_Prostate_Tab, (Late.Diagnostic == 1) &
                                     Late.T1Diabetes ==0 & Late.T2Diabetes ==0 & Late.OtherDiabetes ==0))

PostCodedOnly<- nrow(dplyr::filter(Diabetic_Prostate_Tab,( Late.T1Diabetes ==1 | Late.T2Diabetes ==1 | Late.OtherDiabetes ==1) &
                                     Late.Diagnostic == 0 ))

PostBoth <- nrow(dplyr::filter(Diabetic_Prostate_Tab, (Late.T1Diabetes ==1 | Late.T2Diabetes ==1 | Late.OtherDiabetes ==1) &
                                 Late.Diagnostic == 1 ))

TotalPre<- sum(PreBloodOnly, PreCodedOnly, PreBoth)
TotalPost<- sum(PostBloodOnly, PostCodedOnly, PostBoth)
TotalAll<- sum(AllBloodOnly, AllCodedOnly, AllBoth)

Diab_Pie_Tab<- as.data.frame(cbind(rep(c("Uniquely Identified by HBA1c", "Uniquely Identified by Clinical Coding", "Universally Identified"), 3),
                                   c(PreBloodOnly, PreCodedOnly, PreBoth, PostBloodOnly, PostCodedOnly, PostBoth, AllBloodOnly, AllCodedOnly, AllBoth ),
                                   rep(c("Before Cancer Diagnosis", "After Cancer Diagnosis", "Before or After Cancer Diagnosis"), each =3),
                                   rep(c(TotalPre, TotalPost, TotalAll), each = 3)))

colnames(Diab_Pie_Tab)<- c("Identification", "Number", "TimePoint", "Total")

Diab_Pie_Tab<- Diab_Pie_Tab%>% 
  dplyr::mutate(Percentage = round((as.numeric(as.character(Number))/as.numeric(as.character(Total)))*100, digits =1))%>%
  dplyr::mutate(Identification = factor(Identification, levels = c("Uniquely Identified by HBA1c", "Uniquely Identified by Clinical Coding", "Universally Identified")))%>%
  dplyr::mutate(TimePoint = factor(TimePoint, levels = c("Before Cancer Diagnosis", "After Cancer Diagnosis", 
                                                         "Before or After Cancer Diagnosis"))) %>%
  group_by(TimePoint) %>%
  dplyr::mutate(pos = cumsum(Percentage)- Percentage/2) 

Identification_Pie<- ggplot(Diab_Pie_Tab, aes(y=Percentage, x=factor(1))) + 
  geom_bar(aes(fill= Identification), stat="identity") +
  facet_wrap(~TimePoint)+ 
  coord_polar("y", start=0) +
  geom_text(aes(label = paste0(Percentage, "%"), y = 100-pos), size = 10) +
  xlab("")+
  theme(text = element_text(size = 30),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank())


#Population Visualisation BAsed on Identification Method




Diabetic_Prostate_Tab<- Diabetic_Prostate_Tab %>%
  dplyr::mutate(Sex_Label = ifelse(Sex == 781, "Female", "Male"))

HBA1C_Pre<- Diabetic_Prostate_Tab%>% dplyr::filter(Diagnostic == 1)
Coded_Pre<- Diabetic_Prostate_Tab%>% dplyr::filter(T1Diabetes == 1 | T2Diabetes ==1 | OtherDiabetes == 1,
                                                   is.na(Sex_Label)==F)
Any_Pre<- Diabetic_Prostate_Tab%>% dplyr::filter(AnyDiabetes == 1,
                                                 is.na(Sex_Label)==F)

Demographics_Comp_tab<- data.frame(cbind(c("Median Age", "Standard Deviation Age", "Median IMD" ),
                                         c(median(HBA1C_Pre$Age), sd(HBA1C_Pre$Age), median(HBA1C_Pre$IMDQuintile, na.rm = T)),
                                         c(median(Coded_Pre$Age), sd(Coded_Pre$Age),median(Coded_Pre$IMDQuintile, na.rm = T)),
                                         c(median(Any_Pre$Age), sd(Any_Pre$Age), median(Any_Pre$IMDQuintile, na.rm = T))))



Unique_HBA1c<- Diabetic_Prostate_Tab %>%
  dplyr::filter(Diagnostic==1 &
                  (T1Diabetes == 0&
                     T2Diabetes ==0 &
                     OtherDiabetes == 0),is.na(Sex_Label)==F) %>%
  dplyr::mutate(Diabetes_Route = "HBA1c")

Unique_Coding<- Diabetic_Prostate_Tab %>%
  dplyr::filter(Diagnostic==0 &
                  (T1Diabetes == 1 |
                     T2Diabetes == 1 |
                     OtherDiabetes == 1),is.na(Sex_Label)==F)%>%
  dplyr::mutate(Diabetes_Route = "Coding")

Both_Identified<- Diabetic_Prostate_Tab %>%
  dplyr::filter(Diagnostic==1 &
                  (T1Diabetes == 1 |
                     T2Diabetes == 1 |
                     OtherDiabetes == 1),is.na(Sex_Label)==F)%>%
  dplyr::mutate(Diabetes_Route = "Both")

Demographics_Unique_tab<- data.frame(cbind(c("Median Age", "Standard Deviation Age", "Median IMD"),
                                           c(median(Unique_HBA1c$Age), sd(Unique_HBA1c$Age), median(Unique_HBA1c$IMDQuintile, na.rm = T)),
                                           c(median(Unique_Coding$Age), sd(Unique_Coding$Age),median(Unique_Coding$IMDQuintile, na.rm = T)),
                                           c(median(Both_Identified$Age), sd(Both_Identified$Age), median(Both_Identified$IMDQuintile, na.rm = T))))

colnames(Demographics_Unique_tab)<- c("Variable", "HBA1c Only", "Coding Only", "Both")
#HBA1c VS Both

HB<- as.data.frame(rbind(Unique_HBA1c, Both_Identified)) %>%
  dplyr:: mutate(Label = ifelse(Diabetes_Route =="Both", 1, 0),
                 Sex_Label_Numeric = ifelse(Sex_Label =="Male", 1, 0))

HB_Age<- wilcox.test(Age5 ~ Label, data=HB) 
HB_IMD<- wilcox.test(IMDQuintile ~ Label, data=HB)

#Coding VS Both
CB <-   as.data.frame(rbind(Unique_Coding, Both_Identified))%>%
  dplyr:: mutate(Label = ifelse(Diabetes_Route =="Both", 1, 0),
                 Sex_Label_Numeric = ifelse(Sex_Label =="Male", 1, 0))

CB_Age<- wilcox.test(Age5 ~ Label, data=CB) 
CB_IMD<- wilcox.test(IMDQuintile ~ Label, data=CB)
#HBA1c vs Coding

HC<- as.data.frame(rbind(Unique_HBA1c, Unique_Coding))%>%
  dplyr:: mutate(Label = ifelse(Diabetes_Route =="Coding", 1, 0),
                 Sex_Label_Numeric = ifelse(Sex_Label =="Male", 1, 0))

HC_Age<- wilcox.test(Age5 ~ Label, data=HC) 
HC_IMD<- wilcox.test(IMDQuintile ~ Label, data=HC)

Stats_Tests<- cbind(c("HBA1c Only VS Coding Only", "HBA1c Only VS HBA1c and Coding", "Coding Only VS HBA1c and Coding"),
                    c(HC_Age$p.value, HB_Age$p.value, CB_Age$p.value),
                    c(HC_IMD$p.value, HB_IMD$p.value, CB_IMD$p.value))

colnames(Stats_Tests)<- c("Comparitor_Groups", "Age", "IMD Quintile")


AgeLabels<- c("", "0 - 4", "5 - 9", 
              "10 - 14", "15 - 19",
              "20 - 24", "25 - 29",
              "30 - 34", "35 - 39",
              "40 - 44", "45 - 49",
              "50 - 54", "55 - 59",
              "60 - 64", "65 - 69",
              "70 - 74", "75 - 79",
              "80 - 84", "85 - 89",
              "90 - 94", "95 - 99",
              "100 - 104", "105 - 109",
              "110 - 114", "115 - 119")

n1 <- ggplot() + 
  geom_bar(data = dplyr::filter(HBA1C_Pre, Sex_Label == "Female"), 
           aes(x = Age5, fill = Sex_Label), color = "black",  
           stat = "count")  + 
  geom_bar(data = dplyr::filter(HBA1C_Pre, Sex_Label == "Male"),  
           aes(x = Age5, fill = Sex_Label, y=..count..*(-1)), color = "black",
           stat = "count") + 
  scale_y_continuous(breaks = seq(-30000, 30000, 1000),
                     labels = c(seq(30000, 0, -1000), seq(1000, 30000, 1000))) + 
  scale_x_continuous(breaks = seq(-2, 122, 5),
                     labels = AgeLabels) +
  coord_flip() + 
  scale_fill_brewer(palette = "Set1") + 
  theme_bw() +
  xlab("Number of Patients") +
  ylab("Age Band") +
  ggtitle("Diagnostic HBA1c")

n2 <- ggplot() + 
  geom_bar(data = dplyr::filter(Coded_Pre, Sex_Label == "Female"), 
           aes(x = Age5, fill = Sex_Label), color = "black",  
           stat = "count")  + 
  geom_bar(data = dplyr::filter(Coded_Pre, Sex_Label == "Male"),  
           aes(x = Age5, fill = Sex_Label, y=..count..*(-1)), color = "black",
           stat = "count") + 
  scale_y_continuous(breaks = seq(-30000, 30000, 1000),
                     labels = c(seq(30000, 0, -1000), seq(1000, 30000, 1000))) + 
  scale_x_continuous(breaks = seq(-2, 122, 5),
                     labels = AgeLabels) +
  coord_flip() + 
  scale_fill_brewer(palette = "Set1") + 
  theme_bw()+
  xlab("Number of Patients") +
  ylab("Age Band") +
  ggtitle("Diabetic Coding")

n3 <- ggplot() + 
  geom_bar(data = dplyr::filter(Any_Pre, Sex_Label == "Female"), 
           aes(x = Age5, fill = Sex_Label), color = "black",  
           stat = "count")  + 
  geom_bar(data = dplyr::filter(Any_Pre, Sex_Label == "Male"),  
           aes(x = Age5, fill = Sex_Label, y=..count..*(-1)), color = "black",
           stat = "count") + 
  scale_y_continuous(breaks = seq(-30000, 30000, 1000),
                     labels = c(seq(30000, 0, -1000), seq(1000, 30000, 1000))) + 
  scale_x_continuous(breaks = seq(-2, 122, 5),
                     labels = AgeLabels) +
  coord_flip() + 
  scale_fill_brewer(palette = "Set1") + 
  theme_bw() +
  xlab("Number of Patients") +
  ylab("Age Band") +
  ggtitle("Diagnostic HBA1c or Diabetic Coded")

multiplot(n1, n2, n3)

#KM Plots


All_Surv <- Surv(time = Diabetic_Prostate_Tab$SurvivalTime/365, event = Diabetic_Prostate_Tab$SurvivalStatus)

HBA1c_KMFit<- survfit(All_Surv ~ Diagnostic  , data = Diabetic_Prostate_Tab)
HBA1c_KM<- ggsurvplot(HBA1c_KMFit, data = Diabetic_Prostate_Tab, risk.table = T,
                      legend.title = "",
                      legend.labs = c( "Non Diabetic", "Diabetic HBA1c"),
                      conf.int = T,
                      censor.shape= NA,
                      xlim= c(0,15),
                      ggtheme = theme_bw(),
                      break.time.by=1, xlab = "Time (Years)",
                      title = "Prostate Cancer Survival by Diabetic HBA1c Status" ) 

Coded_KMFit<- survfit(All_Surv ~ HESDiabetic  , data = Diabetic_Prostate_Tab)
Coded_KM<- ggsurvplot(Coded_KMFit, data = Diabetic_Prostate_Tab, risk.table = T,
                      legend.title = "",
                      legend.labs = c( "Non Diabetic", "Diabetic Coded"),
                      conf.int = T,
                      censor.shape= NA,
                      xlim= c(0,15),
                      ggtheme = theme_bw(),
                      break.time.by=1, xlab = "Time (Years)",
                      title = "Prostate Cancer Survival by Diabetic CodingStatus" ) 

Any_KMFit<- survfit(All_Surv ~ AnyDiabetes  , data = Diabetic_Prostate_Tab)
Any_KM<- ggsurvplot(Any_KMFit, data = Diabetic_Prostate_Tab, risk.table = T,
                    legend.title = "",
                    legend.labs = c( "Non Diabetic", "Diabetic by HBA1c or Coding"),
                    conf.int = T,
                    censor.shape= NA,
                    xlim= c(0,15),
                    ggtheme = theme_bw(),
                    break.time.by=1, xlab = "Time (Years)",
                    title = "Prostate Cancer Survival by Any Diabetic Status" ) 

surv_median(Coded_KMFit)
surv_median(HBA1c_KMFit)
surv_median(Any_KMFit)

HESGroup<- Diabetic_Prostate_Tab %>% dplyr::filter (HESDiabetic ==1) %>%
  dplyr::mutate(SurvivalComp = "HES Only")
BloodsGroup<- Diabetic_Prostate_Tab %>% dplyr::filter (Diagnostic ==1) %>%
  dplyr::mutate(SurvivalComp = "Bloods Only")
BothGroup<- Diabetic_Prostate_Tab %>% dplyr::filter (Diagnostic ==1 | HESDiabetic ==1) %>%
  dplyr::mutate(SurvivalComp = "Both")
AllCancers<- Diabetic_Prostate_Tab%>% dplyr::mutate(SurvivalComp = "All Cancer Patients")

Surv_Comp_Tab<- rbind(HESGroup, BloodsGroup, BothGroup, AllCancers)

Comp_Surv <- Surv(time = Surv_Comp_Tab$SurvivalTime/365, event = Surv_Comp_Tab$SurvivalStatus)
Comp_KMFit<- survfit(Comp_Surv ~ SurvivalComp  , data = Surv_Comp_Tab)
Comp_KM<- ggsurvplot(Comp_KMFit, data = Surv_Comp_Tab, risk.table = T,
                     legend.title = "",
                     conf.int = T,
                     censor.shape= NA,
                     xlim= c(0,12),
                     ggtheme = theme_bw(),
                     break.time.by=1, xlab = "Time (Years)",
                     title = "All Cancer Survival by Any Diabetic Status",
                     legend.labs = c("Total Prostate Cancer Population", "HBA1c Only Diabetics", "Hybrid Definition Diabetics ", "Clinical Coding Only Diabetics"),
                     surv.median.line = "hv"  ,
                     font.x = c(30),
                     font.y = c(30),
                     font.tickslab = c(30, "plain", "black"),
                     fontsize = c(10),
                     font.main = c(30),
                     font.legend = c(20),
                     base_size = c(30), tables.theme = theme_cleantable(),risk.table.ylab = c(20)) 

Comp_KM$table<- Comp_KM$table + theme(axis.text.y = element_text(size = 20)) 

survdiff(Comp_Surv ~ SurvivalComp  , data = Surv_Comp_Tab)

#Repeat Adjusting for Age Sex_Label and Deprivation

Cox_HES<- survfit(coxph(Surv(time = SurvivalTime/365, event = SurvivalStatus)~Age5 + IMDQuintile, data = HESGroup))
Cox_Bloods<- survfit(coxph(Surv(time = SurvivalTime/365, event = SurvivalStatus)~Age5 +IMDQuintile, data = BloodsGroup))
Cox_Both<- survfit(coxph(Surv(time = SurvivalTime/365, event = SurvivalStatus)~Age5 + IMDQuintile, data = BothGroup))
Cox_All<- survfit(coxph(Surv(time = SurvivalTime/365, event = SurvivalStatus)~Age5 + IMDQuintile, data = AllCancers))

Median_Table<- cbind(c("HES", "Blood", "Both", "All"),
                     c(surv_median(Cox_HES)$median, surv_median(Cox_Bloods)$median, surv_median(Cox_Both)$median, surv_median(Cox_All)$median),
                     c(surv_median(Cox_HES)$lower, surv_median(Cox_Bloods)$lower, surv_median(Cox_Both)$lower, surv_median(Cox_All)$lower),
                     c(surv_median(Cox_HES)$upper, surv_median(Cox_Bloods)$upper, surv_median(Cox_Both)$upper, surv_median(Cox_All)$upper))
colnames(Median_Table)<- c("Identification", "Median", "Lower", "Upper")

Cox_Comp<- ggsurvplot_combine(fit = list(Cox_HES, Cox_Bloods, Cox_Both, Cox_All), 
                              data = Diabetic_Prostate_Tab,
                              risk.table = T,
                              legend.title = "",
                              legend.labs = c("Clinical Coding Only Diabetics", "HBA1c Only Diabetics", "Hybrid Definition Diabetics ", "Total Prostate Cancer Population" ),
                              conf.int = T,
                              censor.shape= NA,
                              xlim= c(0,10),
                              ggtheme = theme_bw(),
                              break.time.by=1, xlab = "Time (Years)",
                              surv.median.line = c("hv"),
                              title = "Cancer Survival by Diabetic Definition", sep =" ",
                              font.x = c(30),
                              font.y = c(30),
                              font.tickslab = c(30, "plain", "black"),
                              fontsize = c(10),
                              font.main = c(30),
                              font.legend = c(20),
                              base_size = c(30), tables.theme = theme_cleantable(),risk.table.ylab = c(20)) 

Cox_Comp$table<- Cox_Comp$table + theme(axis.text.y = element_text(size = 20)) 

#Change in diagnosis date

BothDiab<- Diabetic_Prostate_Tab%>% dplyr::filter((Diagnostic == 1 | Late.Diagnostic ==1) & 
                                                    (T1Diabetes == 1| T2Diabetes ==1|OtherDiabetes == 1 |
                                                       Late.T1Diabetes == 1|Late.T2Diabetes ==1|Late.OtherDiabetes == 1)) %>%
  dplyr::mutate(PrePostHESDate = pmin( 
    as_date(as.character(T1DiabetesDate)),
    as_date(as.character(T2DiabetesDate)),
    as_date(as.character(OtherDiabetesDate)), na.rm = T)) %>%
  dplyr::mutate( DiabeticDateDiff = as.numeric(as_date(as.character(DiagnosticDate)) - as_date(as.character(PrePostHESDate)))) %>%
  dplyr:: mutate(ComorbShift = ifelse(HESDiabetic == 0 & Diagnostic ==1, 1, 0))

Perc_Shift<- sum(BothDiab$ComorbShift)/nrow(BothDiab)

TimeShiftPlot<- ggplot(BothDiab) + 
  geom_histogram(aes(x=DiabeticDateDiff/365, fill = as.factor(ComorbShift)), binwidth = 0.5, alpha = 0.6, color = "black") +
  xlab("Difference Between Diabetic Diagnosis by Clinical Coding and HBA1c Identification in Years") +
  ylab("Count") +
  geom_vline(xintercept = 0, size = 1.5 ) +
  scale_fill_discrete(name = "",labels = c("Pre/Post Cancer Diabetes Diagnosis Label 
Unchanged With Inclusion of HBA1c", paste("Post Cancer Diabetes Diagnosis Label
Changes With Inclusion of HBA1c (", round(Perc_Shift*100,digits = 2) , "%)", sep ="" ))) +
  theme(text = element_text(size = 30))

# Identification of high volume HBA1c areas for cancer patients


Volume_Locations<- fread("Volumelocations.csv")
Volume_Patients<- fread("VolumePatients.csv")

Volume_Prostate<- inner_join(Diabetic_Prostate_Tab , Volume_Patients, by = c("Study_ID"))

# Create Tables for patient numbers by method of identification 

AllBloodOnly2<- nrow(dplyr::filter(Volume_Prostate, (Diagnostic == 1 | Late.Diagnostic == 1) &
                                     T1Diabetes == 0 & T2Diabetes ==0 & OtherDiabetes ==0 &
                                     Late.T1Diabetes ==0 & Late.T2Diabetes ==0 & Late.OtherDiabetes ==0))

AllCodedOnly2<- nrow(dplyr::filter(Volume_Prostate,(T1Diabetes == 1 | T2Diabetes ==1 | OtherDiabetes ==1 |
                                                      Late.T1Diabetes ==1 | Late.T2Diabetes ==1 | Late.OtherDiabetes ==1) &
                                     Diagnostic == 0 & Late.Diagnostic == 0 ))

AllBoth2 <- nrow(dplyr::filter(Volume_Prostate,(T1Diabetes == 1 | T2Diabetes ==1 | OtherDiabetes ==1 |
                                                  Late.T1Diabetes ==1 | Late.T2Diabetes ==1 | Late.OtherDiabetes ==1) &
                                 (Diagnostic == 1 | Late.Diagnostic == 1 )))


PreBloodOnly2<- nrow(dplyr::filter(Volume_Prostate, (Diagnostic == 1 ) &
                                     T1Diabetes == 0 & T2Diabetes ==0 & OtherDiabetes ==0))

PreCodedOnly2<- nrow(dplyr::filter(Volume_Prostate,(T1Diabetes == 1 | T2Diabetes ==1 | OtherDiabetes ==1 ) &
                                     Diagnostic == 0 ))

PreBoth2 <- nrow(dplyr::filter(Volume_Prostate,(T1Diabetes == 1 | T2Diabetes ==1 | OtherDiabetes ==1 ) &
                                 (Diagnostic == 1 )))


PostBloodOnly2<- nrow(dplyr::filter(Volume_Prostate, (Late.Diagnostic == 1) &
                                      Late.T1Diabetes ==0 & Late.T2Diabetes ==0 & Late.OtherDiabetes ==0))

PostCodedOnly2<- nrow(dplyr::filter(Volume_Prostate,( Late.T1Diabetes ==1 | Late.T2Diabetes ==1 | Late.OtherDiabetes ==1) &
                                      Late.Diagnostic == 0 ))

PostBoth2 <- nrow(dplyr::filter(Volume_Prostate, (Late.T1Diabetes ==1 | Late.T2Diabetes ==1 | Late.OtherDiabetes ==1) &
                                  Late.Diagnostic == 1 ))

TotalPre2<- sum(PreBloodOnly2, PreCodedOnly2, PreBoth2)
TotalPost2<- sum(PostBloodOnly2, PostCodedOnly2, PostBoth2)
TotalAll2<- sum(AllBloodOnly2, AllCodedOnly2, AllBoth2)

Filtered_Pie_Tab<- as.data.frame(cbind(rep(c("Uniquely Identified by HBA1c", "Uniquely Identified by Clinical Coding", "Universally Identified"), 3),
                                       c(PreBloodOnly2, PreCodedOnly2, PreBoth2, PostBloodOnly2, PostCodedOnly2, PostBoth2, AllBloodOnly2, AllCodedOnly2, AllBoth2 ),
                                       rep(c("Before Cancer Diagnosis", "After Cancer Diagnosis", "Before or After Cancer Diagnosis"), each =3),
                                       rep(c(TotalPre2, TotalPost2, TotalAll2), each = 3)))

colnames(Filtered_Pie_Tab)<- c("Identification", "Number", "TimePoint", "Total")

Filtered_Pie_Tab<- Filtered_Pie_Tab%>% 
  dplyr::mutate(Percentage = round((as.numeric(as.character(Number))/as.numeric(as.character(Total)))*100, digits =1))%>%
  dplyr::mutate(Identification = factor(Identification, levels = c("Uniquely Identified by HBA1c", "Uniquely Identified by Clinical Coding", "Universally Identified")))%>%
  dplyr::mutate(TimePoint = factor(TimePoint, levels = c("Before Cancer Diagnosis", "After Cancer Diagnosis", 
                                                         "Before or After Cancer Diagnosis"))) %>%
  group_by(TimePoint) %>%
  dplyr::mutate(pos = cumsum(Percentage)- Percentage/2) 

FilteredIdentification_Pie<- ggplot(Filtered_Pie_Tab, aes(y=Percentage, x=factor(1))) + 
  geom_bar(aes(fill= Identification), stat="identity") +
  facet_wrap(~TimePoint)+ 
  coord_polar("y", start=0) +
  geom_text(aes(label = paste0(Percentage, "%"), y = 100-pos), size  = 10) +
  xlab("")+
  theme(text = element_text(size = 30),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank())


#identify Patients with coding only and no HBA1c tests ever

Filtered_Coding_Only<- dplyr::filter(Volume_Prostate,(T1Diabetes == 1 | T2Diabetes ==1 | OtherDiabetes ==1 |
                                                        Late.T1Diabetes ==1 | Late.T2Diabetes ==1 | Late.OtherDiabetes ==1) &
                                       (Diagnostic == 0 | Late.Diagnostic == 0 )) %>%
  dplyr::mutate(Study_ID = as.integer64(as.character(Study_ID)))

HBA1C_Import<- HBA1C_Import %>% dplyr::mutate(Study_ID = as.integer64(as.character(StudyId1)))
Filtered_Coding_Only_Joined<- left_join(Filtered_Coding_Only, HBA1C_Import, by = c("Study_ID"))

Test_Volume<- Filtered_Coding_Only_Joined %>%
  dplyr::filter(is.na(ReceivedDate)==F) %>%
  group_by(Study_ID)%>%
  summarise(dplyr::n())

Test_Volume_Tab<- as.data.frame(table(unlist(Test_Volume$`dplyr::n()`))) %>%
  dplyr::mutate(Var1= as.numeric(as.character(Var1)))
Test_Volume_Tab<- as.data.frame(rbind(c(0, nrow(dplyr::filter(Filtered_Coding_Only_Joined, is.na(ReceivedDate)==T))),
                                      Test_Volume_Tab))
Testing_Volume<- ggplot(Test_Volume_Tab) + 
  geom_bar(aes(x= Var1, y = Freq), 
           stat = "identity", color = "black", 
           fill = "blue", alpha = 0.4 ) + 
  xlab("Number of HBA1c blood Tests") + 
  ylab("Number of Patients") +
  geom_vline(xintercept = median(Test_Volume$`dplyr::n()`), color = "red" ) =
  theme(text = element_text(size = 20))



# Generate Data Required for Publications

Total_Diagnoses<- nrow( MasterTab %>% 
                          dplyr::mutate(Diagnosis_Year = as.numeric(strtrim(as.character(DiagnosisDate), width = 4))) %>%
                          dplyr::filter(Diagnosis_Year>2004)) # remove patients diagnosed before 2005 as bloods available from 2004 onwards

Total_Patients<- nrow(Diabetic_All_Tab)  


