setwd('/home/luederm/Desktop/Brown-Assignment/pipeline')

baseCov = read.csv('base_anova.csv')

boxplot(Coverage ~ Base, data = baseCov)

a = aov(Coverage ~ Base, data = baseCov)
summary(a)
