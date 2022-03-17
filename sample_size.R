
library(tidyverse)
library(moonBook)
library(pwr)

#fx for effect size
effect.size <- function (m1, m2, s1, s2) {
    (abs(m1-m2))/(((s1^2+s2^2)/2)^0.5)
}



# effect size calculation (mean & sd for ep from pilot study) 
ep <- read_csv("Tables/Early progression/ep.csv")
mytable(recur~.-X1, data= ep, digits=2)
d_1228 <- effect.size(0.9, -1.81, 0.88, 1.9)
d_33a <- effect.size(0.72, -2.22, 2.17, 2.84)
d_3200 <- effect.size(0.87, -2.11, 2.39, 2.22)
d_146a <- effect.size(-1.05, 1.14, 1.93, 1.77)
d_6815 <- effect.size(0.12, -2.31, 1.73, 3.22)
d_fcgr1a <- effect.size(1.94, -3.04, 2.27, 3.09)
d_pck1 <- effect.size(-1.39, 2.79, 3.15, 1.74)
d_group <- effect.size(3.66, -9.6, 3.24, 2.26)

# p < 0.05 power =0.8
pwr.t.test(d=d_1228, sig.level=0.05, power=0.8, type="two.sample") # 6 vs 6
pwr.t.test(d=d_33a, sig.level=0.05, power=0.8, type="two.sample") # 13 vs 13
pwr.t.test(d=d_3200, sig.level=0.05, power=0.8, type="two.sample") # 11 vs 11
pwr.t.test(d=d_146a, sig.level=0.05, power=0.8, type="two.sample") # 13 vs 13
pwr.t.test(d=d_6815, sig.level=0.05, power=0.8, type="two.sample") # 19 vs 19
pwr.t.test(d=d_fcgr1a, sig.level=0.05, power=0.8, type="two.sample") # 6 vs 6
pwr.t.test(d=d_pck1, sig.level=0.05, power=0.8, type="two.sample") # 7 vs 7
pwr.t.test(d=d_group, sig.level=0.05, power=0.8, type="two.sample") # 3 vs 3


# unequal sample size (if 7 of 45 underwent ep)
pwr.t2n.test(n1=7, n2=38, sig.level=0.01,  d=d_1228) #0.01 95
pwr.t2n.test(n1=7, n2=38,  d=d_33a) #0.05 79
pwr.t2n.test(n1=7, n2=38,  d=d_3200) #0.05 87
pwr.t2n.test(n1=7, n2=38,  d=d_146a) #0.05 80
pwr.t2n.test(n1=7, n2=38, d=d_6815) #0.05 61
pwr.t2n.test(n1=7, n2=38,  sig.level=0.01, d=d_fcgr1a) #0.01 96
pwr.t2n.test(n1=7, n2=38,  sig.level=0.01, d=d_pck1) #0.01 90
pwr.t2n.test(n1=7, n2=38,  sig.level=0.01, d=d_group) #0.01 1


# mets effect size from pilot study
dm <- read_csv("Tables/Distant metastasis/dm.csv")
mytable(dm~.-X1-stage, data= dm, digits=2)
d_new <- effect.size(1.17, -2.4, 2.92, 2.32)
d_605 <- effect.size(-0.51, 2.8, 2.98, 2.25)
d_6791 <- effect.size(-0.85, 2.18, 2.85, 2.46)
d_6826 <- effect.size(-1.58, 2.34, 2.08, 2.77)
d_6780a <- effect.size(-1.53, 0.72, 2.97, 2.47)
d_c1qtnf <- effect.size(-0.16, -3.09, 3.45, 3.43)
d_rbp3 <- effect.size(-1.35, 0.83, 2.23, 1.62)
d_group <- effect.size(-5.64, 10.44, 4.94, 4.7)

# p < 0.05 power =0.8
pwr.t.test(d=d_new, sig.level=0.05, power=0.8, type="two.sample") # 10 vs 10
pwr.t.test(d=d_605, sig.level=0.05, power=0.8, type="two.sample") # 12 vs 12
pwr.t.test(d=d_6791, sig.level=0.05, power=0.8, type="two.sample") # 14 vs 14
pwr.t.test(d=d_6826, sig.level=0.05, power=0.8, type="two.sample") # 8 vs 8
pwr.t.test(d=d_6780a, sig.level=0.05, power=0.8, type="two.sample") # 25 vs 25
pwr.t.test(d=d_c1qtnf, sig.level=0.05, power=0.8, type="two.sample") # 23 vs 23
pwr.t.test(d=d_rbp3, sig.level=0.05, power=0.8, type="two.sample") # 14 vs 14
pwr.t.test(d=d_group, sig.level=0.05, power=0.8, type="two.sample") # 3 vs 3

## extrapelvic mets - 11, others 34
pwr.t2n.test(n1=11, n2=34, sig.level = 0.01,  d=d_new) #p=0.01, power=0.88
pwr.t2n.test(n1=11, n2=34, sig.level = 0.01,  d=d_605) #p=0.01, power=0.82
pwr.t2n.test(n1=11, n2=34,  d=d_6791) #p=0.05, power=0.89
pwr.t2n.test(n1=11, n2=34, sig.level = 0.01,  d=d_6826) #p=0.01, power=0.97
pwr.t2n.test(n1=11, n2=34,  d=d_6780a) #p=0.05, power=0.64
pwr.t2n.test(n1=11, n2=34,  d=d_c1qtnf) #p=0.05, power=0.67
pwr.t2n.test(n1=11, n2=34,  d=d_rbp3) #p=0.05, power=0.88
pwr.t2n.test(n1=11, n2=34,  sig.level=0.01, d=d_group) #p=0.01, power=1



# ar  effect size from pilot study
ar <- read_csv("Tables/Acute tumor response/ar.csv")
mytable(ar_m~.-X1-ar, data= ar, digits=2)
d_3928 <- effect.size(-2.29, 1.54, 2.99, 2.98)
d_3960 <- effect.size(-4.07, -0.88, 2.17, 2.95)
d_92a <- effect.size(-2.57, 0.26, 3.32, 3.13)
d_574 <- effect.size(-0.45, 0.22, 0.54, 0.53)
d_DIABLO <- effect.size(-2.21, 0.1, 1.73, 0.99)
d_group <- effect.size(-10.81, 1.35, 4.2, 4.33)


# p < 0.05 power =0.8
pwr.t.test(d=d_3928, sig.level=0.05, power=0.8, type="two.sample") # 11 vs. 11
pwr.t.test(d=d_3960, sig.level=0.05, power=0.8, type="two.sample") # 12 vs. 12
pwr.t.test(d=d_92a, sig.level=0.05, power=0.8, type="two.sample") # 22 vs. 22
pwr.t.test(d=d_DIABLO, sig.level=0.05, power=0.8, type="two.sample") # 7 vs. 7
pwr.t.test(d=d_group, sig.level=0.05, power=0.8, type="two.sample") # 4 vs 4

## if 9 of 45 presents poor response
pwr.t2n.test(n1=9, n2=36, d=d_3928) #p=0.05, power=0.91
pwr.t2n.test(n1=9, n2=36, d=d_3960) #p=0.05, power=0.90
pwr.t2n.test(n1=9, n2=36,  d=d_92a) #p=0.05, power=0.63
pwr.t2n.test(n1=9, n2=36, sig.level = 0.01,  d=d_DIABLO) #p=0.01, power=0.95
pwr.t2n.test(n1=9, n2=36,  sig.level=0.01, d=d_group) #p=0.01, power=1

## linear correlation
summary(lm(ar~hsa.miR.3928.3p , data=ar))
pwr.f2.test(u=1, f2=0.1573^0.5, sig.level = 0.05, power=0.8) #22
pwr.f2.test(u=1, f2=0.1573^0.5, sig.level = 0.01, power=0.9) #41

summary(lm(ar~hsa.miR.3960 , data=ar))
pwr.f2.test(u=1, f2=0.1471^0.5, sig.level = 0.05, power=0.8) #23
pwr.f2.test(u=1, f2=0.1471^0.5, sig.level = 0.01, power=0.9) #43

summary(lm(ar~hsa.miR.202.5p  , data=ar))
pwr.f2.test(u=1, f2=0.1128^0.5, sig.level = 0.05, power=0.8) #26
pwr.f2.test(u=1, f2=0.1128^0.5, sig.level = 0.01, power=0.9) #48

summary(lm(ar~DIABLO , data=ar))
pwr.f2.test(u=1, f2=0.4506^0.5, sig.level = 0.05, power=0.8) #14
pwr.f2.test(u=1, f2=0.4506^0.5, sig.level = 0.01, power=0.9) #26

summary(lm(ar~hsa.miR.92a.1.5p  , data=ar))
pwr.f2.test(u=1, f2=0.2424^0.5, sig.level = 0.05, power=0.8) #19
pwr.f2.test(u=1, f2=0.2424^0.5, sig.level = 0.01, power=0.9) #34

summary(lm(ar~hsa.miR.574.3p  , data=ar))
pwr.f2.test(u=1, f2=0.1788^0.5, sig.level = 0.05, power=0.8) #21
pwr.f2.test(u=1, f2=0.1788^0.5, sig.level = 0.01, power=0.9) #39

# multiple linear regression
summary(lm(ar~group  , data=ar))
pwr.f2.test(u=5, f2=0.7299^0.5, sig.level = 0.05, power=0.8) #18
pwr.f2.test(u=5, f2=0.7299^0.5, sig.level = 0.01, power=0.9) #30



# Dx
dx <- read_csv("Tables/Diagnosis/dx.csv")
mytable(cancer~.-X1, data=dx, digits=2)
d_RGS18 <- effect.size(4.1, 5.87, 0.72, 0.82)
d_SNORA12 <- effect.size(6.28, 4.1, 2.73, 2.22)
d_SNORD97 <- effect.size(6.8, 1.34, 4.09, 2.14)


## SNORA12
pwr.t2n.test(n2=44, sig.level = 0.01, power=0.9,  d=d_RGS18) # n1 =4 

## SNORA12
pwr.t2n.test(n2=44, sig.level = 0.05, power=0.8,  d=d_SNORA12) # n1 =14
pwr.t2n.test(n2=44, sig.level = 0.05, power=0.9,  d=d_SNORA12) # n1 =21 
pwr.t2n.test(n2=44, sig.level = 0.01, power=0.9,  d=d_SNORA12) # n1 =38 


## SNORD97
pwr.t2n.test(n2=44, sig.level = 0.05, power=0.8,  d=d_SNORD97) # n1 =4
pwr.t2n.test(n2=44, sig.level = 0.01, power=0.9,  d=d_SNORD97) # n1 =7

