require(ggplot2)
require(reshape)
require(gridExtra)
library(reshape2)
library(plotly)
require(webshot)
library(htmlwidgets)


# clear everything
rm(list = ls())

set.seed(100)

setwd(
  "C:/wk/odrive/Amazon Cloud Drive/ajays_stuff/georgia_tech_ms/math_cse6644_interative_methods_for_system_of_equations/project"
)

sink("project_1_script_output.txt",
     append = FALSE,
     split = TRUE)
sink()


# read the matlab results into R dataframe
xDf <-
  read.csv(paste("project_1_x", ".csv", sep = ""),
           header = FALSE)

# read the matlab results into R dataframe
yDf <-
  read.csv(paste("project_1_y", ".csv", sep = ""),
           header = FALSE)


# read the matlab results into R dataframe
pDf <-
  read.csv(paste("project_1_perf", ".csv", sep = ""),
           header = FALSE)


colnames(pDf) <- c('T','Iterations','Time(Secs)','Residual Norm','h')

rownames(xDf) <- 1:length(xDf[['V1']])
rownames(yDf) <- 1:length(xDf[['V1']])  

xM <- data.matrix(xDf)
yM <- data.matrix(yDf)


# rows in each x,y
M <- as.vector((pDf$T/.1)+1)

# melt the performance data for gg plot of the data frame with pivot as id
xDf$id <- 1:length(xDf[['V1']])
yDf$id <- 1:length(xDf[['V1']])

# calculate V


# GG Plot of X,Y for each T
dat <- setNames(lapply(1:length(M), function(x) cbind(x=xDf[1:M[x],x], y=yDf[1:M[x],x])), as.vector(pDf$T))

list.names <- names(dat)
lns <- sapply(dat, nrow)
dat <- as.data.frame(do.call("rbind", dat))
dat$group <- rep(list.names, lns)

# Plot x(t), y(t)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

par(
  mfrow = c(1, 1),
  oma = c(0, 0, 2, 0),
  lab = c(2, 5, 3),
  lwd = 1,
  pch = 19
)

g1 <- ggplot(dat, aes(x = x, y = y)) +
  labs(x = "X(t)",
       y = "Y(t))",
       colour = "")  +
  geom_line(linetype="solid",aes(colour = group)) + 
  ggtitle(paste("X(t) Vs Y(t) by T",sep="")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(name="T",breaks=pDf$T,values=gg_color_hue(length(M)))

grid.arrange(g1)

dev.copy(filename = paste("proj_1_g_xy", ".png", sep = ""),
         device=png, width=1000,height=700)

dev.off()


# Plot individual plots for each T
for( i in 1:length(pDf$T)) {
  
  par(
    mfrow = c(1, 1),
    oma = c(0, 0, 2, 0),
    lab = c(2, 5, 3),
    lwd = 1,
    pch = 19
  )
  
  gi <- ggplot(dat[which(dat$group == pDf$T[i] ),], 
               aes(x = x, y = y),group=pDf$T[i]) +
    labs(x = "X(t)",
         y = "Y(t))",
         colour = "")  +
    geom_line(linetype="solid",aes(colour = group)) + 
    ggtitle(paste("X(t) Vs Y(t) for T=",pDf$T[i],sep="")) +
    theme(plot.title = element_text(hjust = 0.5),legend.position="none")
  
  grid.arrange(gi)
  
  dev.copy(filename = paste("proj_1_g_xy_t_",pDf$T[i], ".png", sep = ""),
           device=png, width=1000,height=700)
  
  dev.off()
  
}

# Plot x(t), y(t) as facet

par(
  mfrow = c(1, 1),
  oma = c(0, 0, 2, 0),
  lab = c(2, 5, 3),
  lwd = 1,
  pch = 19
)

# arrange facets in order
dat$group_f <- factor(dat$group,levels=pDf$T)

g2 <- ggplot(dat, aes(x=x, y=y, group = 1)) +
  geom_line(aes(colour = group)) +
  facet_grid(group_f~., scales = "free_y") +
  labs(x = "X(t)",
       y = "Y(t))",
       colour = "") +
  ggtitle(paste("X(t) Vs Y(t) by T",sep="")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(name="T",breaks=pDf$T,values=gg_color_hue(length(M)))


grid.arrange(g2)

dev.copy(filename = paste("proj_1_g_xyf", ".png", sep = ""),
         device=png, width=1000,height=700)

dev.off()




# Plot of the performance for the different T's

#melt with method name as pivot
molten <-
  melt(
    pDf[, c('T', 'Iterations', 'Time(Secs)', 'Residual Norm')],
    id.vars = c('T'),
    measure.vars = c('Iterations', 'Time(Secs)', 'Residual Norm'),
    variable_name = 'series'
  )

par(
  mfrow = c(1, 1),
  oma = c(0, 0, 2, 0),
  lab = c(2, 5, 3),
  lwd = 1,
  pch = 19
)

klabs <- c(pDf$T)


g3 <- ggplot(molten, aes(x=T, y=value, group = 1)) +
  geom_line(aes(colour = T)) +
  facet_grid(series~., scales = "free_y") +
  labs(x = "T",
       y = "",
       colour = "") +
  ggtitle(paste("Newtons Method Performance Vs T",sep="")) +
  theme(plot.title = element_text(hjust = 0.5),legend.position="none") 


grid.arrange(g3)

dev.copy(png,
         filename = "proj_1_g_perf.png",
         width = 1000,
         height = 700)

dev.off()



# 3D Plot of V(x,y)

# Use the one with the max x,y
#vI <- which(M==max(M))

# v(x,y) function
V_x_y <- function(x,y) {
  ((1 - x^2)^2)/4 + ((y + x^2 - 1)^2)/2 
}


for (vI in 1:length(M)) {
  
  
  mX <- min(xDf[vI])
  iX <- max(xDf[vI])
  
  mX <- -1
  iX <- 1
  
  mY <- min(yDf[vI])
  iY <- max(yDf[vI])
  
  xx <- seq(mX, iX, length.out = 100)  
  yy <- seq(mY, iY, length.out = 100)   
  
  
  zz <- data.frame(outer(xx, yy, V_x_y))  
  
  colnames(zz) <- yy
  rownames(zz) <- xx
  
  zz <- as.matrix(zz)
  # Create lists for axis properties
  f1 <- list(
    family = "Arial, sans-serif",
    size = 18,
    color = "red")
  
  f2 <- list(
    family = "Old Standard TT, serif",
    size = 14,
    color = "#ff9999")
  
  axis <- list(
    titlefont = f1,
    tickfont = f2,
    showgrid = T
  )
  
  xaxis <- c(list(title='x(t)'),axis)
  yaxis <- c(list(title='y(t)'),axis)
  zaxis <- c(list(title='V(x,y)'),axis)
  
  scene = list(
    xaxis = yaxis,
    yaxis = xaxis,
    zaxis = zaxis,
    camera = list(eye = list(x = -1.65, y = 1.25, z = 1.55)))
  
  
  
  ply <- plot_ly(y=xx,x=yy,z=zz, type = "surface", size = I(3),
                 width=1000,
                 height=800) %>% 
    layout(title = paste("V(x,y)"), scene = scene)  
  
  
  #zz <- as.matrix(100* V_x_y(xDf[vI],yDf[vI]))
  
  #rm(ply)
  
  #ply <- plot_ly( width=1000,height=800) %>% 
  #  add_trace(z=zz, x = xDf[vI], y = yDf[vI], mode = "lines+markers", type="scatter3d",opacity = 1, line = list(width = 6,  reverscale = FALSE))
  #%>% add_surface(z=z, x=x, y=y, size = I(3))
  #%>% layout(title = paste("V(x,y)"), scene = scene) 
  
  
  saveWidget(as.widget(ply), paste("project_1_g_5_",vI ,".html", sep = ""))
  #webshot("temp.html", file = paste("hw2_q_5_g_sax_", n, ".png", sep = ""),
  #        cliprect = "viewport")
  #export(ply, file = "test.png")
  
  #rm(ply)
}




