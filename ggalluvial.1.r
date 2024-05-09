library(ggalluvial)

data(mtcars)

ggplot(mtcars,
  aes(axis1 = factor(gear, levels = c(5,4,3)),  # 调整 gear 的顺序为 3, 4, 5
      axis2 = factor(cyl, levels = c(8,6,4)),   # 调整 cyl 的顺序为 4, 6, 8
      y = mpg)) +
  geom_alluvium(aes(fill = factor(vs)), alpha = 0.7) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_minimal()


library(ggplot2)
library(ggalluvial)

ggplot(as.data.frame(UCBAdmissions), aes(y = Freq, axis1 = Gender, axis2 = Dept)) +
  geom_alluvium(aes(fill = Admit), width = 1/14) +
  geom_stratum(width = 1/8, fill = "green",color = "black") + #设置节点及边框颜色
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) + #设置节点标签
  scale_x_discrete(limits = c("Gender", "Dept"), expand = c(0, 0)) +
  # scale_fill_brewer(type = "qual",palette = 6) + #颜色填充
  ggtitle("UC Berkeley admissions and rejections, by sex and department")



ggplot(as.data.frame(HairEyeColor), aes(y = Freq,axis1 = Hair, axis2 = Eye, fill = Eye,axis3 = Sex,alpha=Freq))+
  geom_flow()+
  geom_alluvium(aes(fill = Eye), curve_type = "sine") +
  scale_fill_manual(values = c(Brown = "#70493D", Hazel = "#E2AC76",Green = "#3F752B", Blue  = "#81B0E4")) +
  guides(fill = "none") +
  geom_stratum(alpha = .2) +
  geom_text(stat = "stratum", size=3, aes(label = after_stat(stratum)), reverse = T) +
  scale_x_continuous(breaks = 1:3, expand = c(0,0), labels = c("Hair", "Eye", "Sex")) +
  ggtitle("Eye colors of 592 subjects, by sex and hair color")



titanic_wide <- data.frame(Titanic)
head(titanic_wide)
#>   Class    Sex   Age Survived Freq
#> 1   1st   Male Child       No    0
#> 2   2nd   Male Child       No    0
#> 3   3rd   Male Child       No   35
#> 4  Crew   Male Child       No    0
#> 5   1st Female Child       No    0
#> 6   2nd Female Child       No    0
ggplot(data = titanic_wide,
       aes(axis1 = Class, axis2 = Sex, axis3 = Age,
           y = Freq)) +
  scale_x_discrete(limits = c("Class", "Sex", "Age"), expand = c(.2, .05)) +
  xlab("Demographic") +
  geom_alluvium(aes(fill = Survived)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_minimal() +
  ggtitle("passengers on the maiden voyage of the Titanic",
          "stratified by demographics and survival")

