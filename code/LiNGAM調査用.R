setwd("C:/Users/Rstudio/Desktop/graduation_thesis")

n <- 10000
eps1 <- rnorm(n,mean = 0) # ノイズ1
eps2 <- rnorm(n,mean = 0) # ノイズ2
z <- rnorm(n,mean = 0,sd = 1) # 未観測共通要因

# train1
x <- 0.3*z + eps1
y <- 0.8*x + 0.3*z + eps2
train1 <- data.frame(X = x,Y = y,label = x)

#train2
x <- 0.8*y + 0.3*z + eps1
y <- 0.3*z + eps2
train2 <- data.frame(X = x,Y = y,label = x)

#train3
x <- 0.9*z + eps1
y <- 0.9*z + eps2
train3 <- data.frame(X = x,Y = y,label = x)

# 作図
library(ggplot2)
g <- ggplot(train1,aes(X,Y)) + geom_point(colour = "blue") + geom_smooth(method = "lm")
g <- g + annotate("text",label="x=0.3z+eps1",x=2.5,y=-3,fontface="italic")  
g <- g + annotate("text",label="y=0.8x+0.3z+eps2",x=2.5,y=-3.4,fontface="italic")
ggsave("matomebox/data/gap_sample1.png", g)

g <- ggplot(train2,aes(X,Y)) + geom_point(colour = "blue") + geom_smooth(method = "lm")
g <- g + annotate("text",label="x=0.8y+0.3z+eps1",x=2.5,y=-3,fontface="italic")  
g <- g + annotate("text",label="y=0.3z+eps2",x=2.5,y=-3.4,fontface="italic")
ggsave("matomebox/data/gap_sample2.png", g)

g <- ggplot(train3,aes(X,Y)) + geom_point(colour = "blue") + geom_smooth(method = "lm")
g <- g + annotate("text",label="x=0.9z+eps1",x=2.5,y=-3,fontface="italic")  
g <- g + annotate("text",label="y=0.9z+eps2",x=2.5,y=-3.4,fontface="italic")
ggsave("matomebox/data/gap_sample3.png", g)