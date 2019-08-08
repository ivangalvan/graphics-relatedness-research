#
#
#    March 15, 2017  
#
#    Graphics for relatedness research     ######
#
#    Galvan-Femenia I., Graffelman J., Barcelo-i-Vidal C.
#
#
#
#    CPU: Linux 3.16.0-4-amd64 #1 SMP Debian 3.16.36-1+deb8u2 (2016-10-19) x86_64 GNU/Linux
#    Model name: Intel(R) Xeon(R) CPU E5-2609 v2 @ 2.50GHz
#    Memory: 32.00 GB
#    Cores: 8
#
#
#    Data from Rosenberg lab at Standford University: https://rosenberglab.stanford.edu/diversity.html
#
#    multiplot function from Peter Haschke: http://www.peterhaschke.com/r/2013/04/24/MultiPlot.html
#
#    stat_chull function from A. Kassambara : https://github.com/kassambara/ggpubr/blob/master/R/stat_chull.R
#



# install.packages("curl") # version 2.3
# install.packages("data.table") # version 1.10.4
# install.packages("ggplot2") # version 2.2.1
# install.packages("ggtern")  # version 2.2.0
# install.packages("dplyr")   # version 0.5.0
# install.packages("iterpc")  # version 0.3.0
# install.packages("compiler") # version 3.3.1
# install.packages("genetics") # version 1.3.8.1
# install.packages("Rsolnp")   # version 1.16
# install.packages("doParallel") # version 1.0.10
# install.packages("foreach")   # version 1.4.3
# install.packages("doSNOW")  # version 1.0.14


library(curl)
library(data.table) 
library(ggplot2) 
library(ggtern) 
library(dplyr)  
library(iterpc)  
library(compiler) 
library(genetics)
library(Rsolnp)   
library(doParallel) 
library(foreach)   
library(doSNOW)  


source("functions.R")


n_cores <- 4 # set the number of available cores in your PC for the IBD calculations

set.seed(080317)


X_data = fread("https://rosenberglab.stanford.edu/data/rosenbergEtAl2002/diversitydata.stru",na.strings = "-9")



#       The first five columns are: #####
#
#  https://rosenberglab.stanford.edu/data/rosenbergEtAl2002/diversityreadme.txt
#
#
# (V1) Individual code number assigned by CEPH.
# (V2) Population code number assigned by us.
# (V3) Population name.
# (V4) Geographic information about the population.
# (V5) Pre-defined region, as was used in the article.
#



### Subset of unrelated individuals #####


maya = X_data %>% filter(V3=="Maya")

maya[1:5,1:10]

related_individuals = c(867,874,866,865,869,878) # Documented by Rosenberg NA, 2006

Maya <- maya[!(maya$V1 %in% related_individuals),]

Maya <- onerowperAll(x = Maya, id = Maya$V1 , colmarkers = 6:382)
Maya[1:8,1:8]

dim(Maya) # 19 unrelated individuals

ibs_matrix <- allelesharing(X = Maya, colmarkers = 1:377)
ibs_matrix[1:4,1:8]

per <- percentages(ibs_matrix)

per[1:3,]

plot(per$P0,per$P2,asp=1,xlim = c(0,1),ylim = c(0,1),xlab="p0",ylab="p2",pch=16) # All individuals are unrelated



### Simulate family relationships #####

Maya_sim <- Maya

for(i in 1:4){
  
  Maya_sim <- children(Maya_sim,4*(i-1) + 1,4*(i-1) + 2) 
  Maya_sim <- children(Maya_sim,4*(i-1) + 1,4*(i-1) + 2)
  Maya_sim <- children(Maya_sim,4*(i-1) + 1,4*(i-1) + 2) 
  Maya_sim <- children(Maya_sim,4*(i-1) + 1,4*(i-1) + 2) 
  Maya_sim <- children(Maya_sim,4*(i-1) + 3,19 + 10*(i-1) + 1) 
  Maya_sim <- children(Maya_sim,4*(i-1) + 3,19 + 10*(i-1) + 1) 
  Maya_sim <- children(Maya_sim,4*(i-1) + 3,19 + 10*(i-1) + 1) 
  Maya_sim <- children(Maya_sim,4*(i-1) + 4,19 + 10*(i-1) + 2) 
  Maya_sim <- children(Maya_sim,4*(i-1) + 4,19 + 10*(i-1) + 2) 
  Maya_sim <- children(Maya_sim,4*(i-1) + 4,19 + 10*(i-1) + 2) 
  
  
}


dim(Maya_sim) # 40 simulated individuals


for(i in 20:59){
  Maya_sim[i,][grep("NA",Maya_sim[i,])] = NA} # convert to NA heterozygous individuals with X/NA or NA/X genotypes, where X is not missing

for(i in 20:59){
  Maya_sim[i,sample(1:377,sample(1:3,1))] = "0/0"} # convert from one to three genotypes to 0/0 as genotyping error



ibs_matrix_sim <- allelesharing(X = Maya_sim, colmarkers = 1:377)
ibs_matrix_sim[1:4,1:8]

per_sim <- percentages(ibs_matrix_sim)

plot(per_sim$P0,per_sim$P2,asp=1,xlim = c(0,1),ylim = c(0,1),xlab="p0",ylab="p2",pch=16)



### Label each simulated family relationship ####


### First cousins (FC) #####


a_ID1 = 24:26
b_ID2 = 27:29

c_ID1 = 34:36
d_ID2 = 37:39

e_ID1 = 44:46
f_ID2 = 47:49

g_ID1 = 54:56
h_ID2 = 57:59

names1 = rep(0,9)
names2 = rep(0,9)
names3 = rep(0,9)
names4 = rep(0,9)

count = 1

for(i in 1:3){
  for(j in 1:3){
    
    names1[count] = paste0("ID",a_ID1[i],"-ID",b_ID2[j])
    names2[count] = paste0("ID",c_ID1[i],"-ID",d_ID2[j])
    names3[count] = paste0("ID",e_ID1[i],"-ID",f_ID2[j])
    names4[count] = paste0("ID",g_ID1[i],"-ID",h_ID2[j])
    
    
    count <- count + 1
  }
  
}


per_fc = per_sim[which(rownames(per_sim) %in% c(names1,names2,names3,names4)),]


plot(per_fc$P0,per_fc$P2,asp=1,xlim = c(0,1),ylim = c(0,1),xlab="p0",ylab="p2",pch=16)




### Second degree #####


names_pairs = NULL

for(k in 0:3){
  
  
  a_ID1 = 21:23 + k*10
  b_ID2 = 24:26 + k*10
  
  c_ID1 = c(20,22,23) + k*10
  d_ID2 = 27:29 + k*10
  
  e_ID1 = 1 + (4*k)
  
  g_ID1 = 2 + (4*k)
  
  names1 = rep(0,9)
  names2 = rep(0,9)
  names3 = rep(0,3)
  names4 = rep(0,3)
  names5 = rep(0,3)
  names6 = rep(0,3)
  
  
  count = 1
  
  for(i in 1:3){
    for(j in 1:3){
      
      names1[count] = paste0("ID",a_ID1[i],"-ID",b_ID2[j])
      names2[count] = paste0("ID",c_ID1[i],"-ID",d_ID2[j])
      names3[count] = paste0("ID",e_ID1,"-ID",b_ID2[j])
      names4[count] = paste0("ID",e_ID1,"-ID",d_ID2[j])
      names5[count] = paste0("ID",g_ID1,"-ID",b_ID2[j])
      names6[count] = paste0("ID",g_ID1,"-ID",d_ID2[j])
      
      count <- count + 1
    }
  }
  
  names_pairs = c(names_pairs,names1,names2,names3,names4,names5,names6)
  
}



per_second = per_sim[which(rownames(per_sim) %in% names_pairs),]


plot(per_second$P0,per_second$P2,asp=1,xlim = c(0,1),ylim = c(0,1),xlab="p0",ylab="p2",pch=16)



### Full sibs (FS) #######



names_pairs = NULL

for(k in 0:3){
  
  
  a_ID1 = 20:23 + k*10
  b_ID1 = 24:26 + k*10
  c_ID1 = 27:29 + k*10
  
  names1 = rep(0,16)
  names2 = rep(0,9)
  names3 = rep(0,9)
  
  
  count = 1
  
  for(i in 1:4){
    for(j in 1:4){
      
      names1[count] = paste0("ID",a_ID1[i],"-ID",a_ID1[j])
      
      count <- count + 1
    }
  }
  
  
  count = 1
  
  for(i in 1:3){
    for(j in 1:3){
      
      names2[count] = paste0("ID",b_ID1[i],"-ID",b_ID1[j])
      names3[count] = paste0("ID",c_ID1[i],"-ID",c_ID1[j])
      
      count <- count + 1
    }
  }
  
  names_pairs = c(names_pairs,names1,names2,names3)
  
}



per_fs = per_sim[which(rownames(per_sim) %in% names_pairs),]


plot(per_fs$P0,per_fs$P2,asp=1,xlim = c(0,1),ylim = c(0,1),xlab="p0",ylab="p2",pch=16)




### Parent-offspring (PO) #####


names_pairs = NULL

for(k in 0:3){
  
  
  a_ID = 20:23 + k*10
  b_ID = 24:26 + k*10
  c_ID = 27:29 + k*10
  
  
  d_ID = 1 + (4*k)
  e_ID = 2 + (4*k)
  f_ID = 3 + (4*k)
  g_ID = 4 + (4*k)
  h_ID = 20 + k*10
  i_ID = 21 + k*10
  
  names1 = rep(0,4)
  names2 = rep(0,4)
  names3 = rep(0,3)
  names4 = rep(0,3)
  names5 = rep(0,3)
  names6 = rep(0,3)
  
  count = 1
  
  for(i in 1:4){
    
    names1[count] = paste0("ID",d_ID,"-ID",a_ID[i])
    names2[count] = paste0("ID",e_ID,"-ID",a_ID[i])
    
    count <- count + 1
  }
  
  
  
  count = 1
  
  for(i in 1:3){
    
    names3[count] = paste0("ID",f_ID,"-ID",b_ID[i])
    names4[count] = paste0("ID",g_ID,"-ID",c_ID[i])
    names5[count] = paste0("ID",h_ID,"-ID",b_ID[i])
    names6[count] = paste0("ID",i_ID,"-ID",c_ID[i])
    
    count <- count + 1
  }
  
  
  names_pairs = c(names_pairs,names1,names2,names3,names4,names5,names6)
  
}



per_po = per_sim[which(rownames(per_sim) %in% names_pairs),]


plot(per_po$P0,per_po$P2,asp=1,xlim = c(0,1),ylim = c(0,1),xlab="p0",ylab="p2",pch=16)



### Plot simulated family relationships with colors #####


per_sim$groups = "UN" # unrelated

per_sim[which(rownames(per_sim) %in% rownames(per_fc)),]$groups = "FC"

per_sim[which(rownames(per_sim) %in% rownames(per_second)),]$groups = "second"

per_sim[which(rownames(per_sim) %in% rownames(per_fs)),]$groups = "FS"

per_sim[which(rownames(per_sim) %in% rownames(per_po)),]$groups = "PO"

table(per_sim$groups)


plot(per_sim$P0,per_sim$P2,asp=1,xlim = c(0,1),ylim = c(0,1),xlab="p0",ylab="p2",pch=16,col=as.numeric(as.factor(per_sim$groups))+3)
title("(p0-p2)-plot for simulated data")




### 
### 
### 
### IBS Graphics #####
### 
### 
### 
### 



#
#
### Mean versus standard deviation plot #####
#
#



Maya <- onerowperAll(x = maya, id = maya$V1 , colmarkers = 6:382)
Maya[1:8,1:8]

dim(Maya) # 25 individuals

ibs_matrix <- allelesharing(X = Maya, colmarkers = 1:377)
ibs_matrix[1:4,1:8]


# related pairs IDs reported by Rosenberg, 2002 ####

fc <- c("858-867", "862-866", "865-873" ,"865-874" ,"873-874", "868-869", "859-865")
po <- c("858-866", "862-867")
fs <- c("876-878")
second <- c("854-874","866-867")

fc <- which(rownames(ibs_matrix) %in% fc)
po <- which(rownames(ibs_matrix) %in% po)
fs <- which(rownames(ibs_matrix) %in% fs)
second <- which(rownames(ibs_matrix) %in% second)


# Compute mean and sd for the Maya population and simulated data 

m_s = mean_sd(ibs_matrix)

rownames(m_s) = rownames(ibs_matrix)


m_s$groups = "UN"

m_s[po,]$groups = "PO"
m_s[fs,]$groups = "FS"
m_s[second,]$groups = "second"
m_s[fc,]$groups = "FC"

table(m_s$groups)

plot(m_s$mean,m_s$sd,xlim=c(0,2),ylim=c(0,1),pch=16,xlab="mean",ylab="sd",col=as.numeric(as.factor(m_s$groups))+3)
title("Mean vs sd plot for the Maya population")


m_s_sim = mean_sd(ibs_matrix_sim)

rownames(m_s_sim) = rownames(ibs_matrix_sim)

m_s_sim$groups = per_sim$groups

table(m_s_sim$groups)

plot(m_s_sim$mean,m_s_sim$sd,xlim=c(0,2),ylim=c(0,1),pch=16,xlab="mean",ylab="sd",col=as.numeric(as.factor(m_s_sim$groups))+3)
title("Mean vs sd plot for simulated data")


# joint Maya and simulated individuals

m_s_maya_sim = rbind(m_s_sim,m_s)

m_s_maya_sim$groups = as.factor(m_s_maya_sim$groups)

table(m_s_maya_sim$groups)


## compute the domain of the mean versus sd plot


I = iterpc(c(0,1,2), 100, replace=TRUE)
a = as.matrix(getall(I))

mmean2 <- apply(a, MARGIN = 1, function(x) mean(x,na.rm = T))
msd2 <- apply(a, MARGIN = 1, function(x) sd(x,na.rm = T))

m_s2 = data.frame(mean=mmean2-1,sd=msd2)

curve1 = 1
sum = 1
for(i in 1:100){ sum = sum + i; curve1 = c(curve1,sum)}

curve2 = NULL
sum = 1
for(i in 2:100){ sum = sum + i; curve2 = c(curve2,sum)}

curve3 = 5051:5151

aux1 = m_s2[curve1,]; aux1$groups = "curve"
aux2 = m_s2[c(1,curve2,5151),]; aux2$groups = "curve"
aux3 = m_s2[curve3,]; aux3$groups = "curve"


m_s_maya_sim = rbind(m_s_maya_sim,aux1)
m_s_maya_sim = rbind(m_s_maya_sim,aux2)
m_s_maya_sim = rbind(m_s_maya_sim,aux3)



pa = ggplot(m_s_maya_sim[1712:2011,], aes(mean, sd, col = groups,shape=groups)) + 
  geom_line(data = m_s_maya_sim[2012:2112,],aes(x=mean,y=sd),size=1,col="gray")+
  geom_line(data = m_s_maya_sim[2113:2213,],aes(x=mean,y=sd),size=1,col="gray")+
  geom_line(data = m_s_maya_sim[2214:2314,],aes(x=mean,y=sd),size=1,col="gray")+
  stat_chull(data=m_s_maya_sim[1:1711,],fill = NA,show_guide = FALSE) +
  scale_colour_manual(values =  c("darkgoldenrod2","dodgerblue2","firebrick2","darkorchid2","forestgreen","gray")) +
  scale_shape_manual(values = c(3,17, 3, 7, 15,18))+
  coord_fixed(ratio=1,xlim=c(0.09,2),ylim=c(0.05,1.03))+
  theme(panel.border=element_blank(), axis.line=element_line()) + 
  theme(axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(), 
        panel.background=element_blank(),
        axis.ticks = element_blank(),
        axis.title.x=element_text(size=12,vjust=-0.5),
        axis.title.y=element_text(size=12,vjust=1.5),
        plot.title = element_text(lineheight=1, face="bold",hjust=0.5),
        axis.line.x=element_line(),
        axis.line.y=element_line()) +
  guides(shape = guide_legend(reverse = TRUE), color = guide_legend(reverse = TRUE)) +
  ggtitle("A")+
  xlab("Mean")+ylab("Standard deviation")+
  geom_point(size=2,show_guide = FALSE)


pa





#
#
### (p0-p2)-plot ####
#
#




per_maya <- percentages(ibs_matrix)

per_maya$groups = m_s$groups

per_maya_sim = rbind(per_sim,per_maya)

per_maya_sim$groups = as.factor(per_maya_sim$groups)

table(per_maya_sim$groups)


x_lab =expression(paste(p["0"]))
y_lab= expression(paste(p["2"]))


pb = ggplot(per_maya_sim[1712:2011,], aes(P0, P2, col = groups,shape=groups))+
  stat_chull(data=per_maya_sim[1:1711,],fill = NA,show_guide = FALSE) +
  scale_colour_manual(values =  c("darkgoldenrod2","dodgerblue2","firebrick2","darkorchid2","forestgreen")) +
  scale_shape_manual(values = c(17, 3, 7, 15,18))+
  coord_fixed(ratio=1,xlim=c(0.03,1),ylim=c(0.03,1))+
  theme(panel.border=element_blank(), axis.line=element_line()) + 
  theme(axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(), 
        panel.background=element_blank(),
        axis.ticks = element_blank(),
        axis.title.x=element_text(size=12,vjust=-0.5),
        axis.title.y=element_text(size=12,vjust=1.5),
        plot.title = element_text(lineheight=1, face="bold",hjust=0.5),
        axis.line.x=element_line(),
        axis.line.y=element_line()) +
  guides(shape = guide_legend(reverse = TRUE), color = guide_legend(reverse = TRUE)) +
  ggtitle("B")+
  xlab(x_lab)+ylab(y_lab)+
  geom_point(size=2,show_guide = FALSE)+
  geom_abline(intercept = 1,slope=-1,col="gray") 



pb




#
#
### Ternary diagram ####
#
#



x_lab=expression(p[0])
y_lab=expression(p[1])
z_lab=expression(p[2])



pc = ggtern(data = per_maya_sim[1712:2011,],
            aes(P0, P1, P2,colour=groups,shape=groups)) + geom_mask()+
  
  theme_bw() +
  scale_colour_manual(values =  c("darkgoldenrod2","dodgerblue2","firebrick2","darkorchid2","forestgreen")) +
  scale_shape_manual(values = c(17, 3, 7, 15,18))+
  
  theme(plot.title = element_text( size=14,face="bold",hjust = 0.5)) +
  xlab(x_lab) + ylab(y_lab) + zlab(z_lab) +
  guides(color = "none", fill = "none", alpha = "none") + 
  labs(title="C") +
  
  geom_path(data=per_maya_sim[per_maya_sim$groups=="second",][c(chull(per_maya_sim[per_maya_sim$groups=="second",]),
                                                                chull(per_maya_sim[per_maya_sim$groups=="second",])[1]),],
            color="darkorchid2") +
  
  geom_path(data=per_maya_sim[per_maya_sim$groups=="FC",][c(chull(per_maya_sim[per_maya_sim$groups=="FC",]),
                                                            chull(per_maya_sim[per_maya_sim$groups=="FC",])[1]),],
            color="darkgoldenrod2") +
  
  geom_path(data=per_maya_sim[per_maya_sim$groups=="UN",][c(chull(per_maya_sim[per_maya_sim$groups=="UN",]),
                                                            chull(per_maya_sim[per_maya_sim$groups=="UN",])[1]),],
            color="forestgreen") +
  
  geom_path(data=per_maya_sim[per_maya_sim$groups=="PO",][c(chull(per_maya_sim[per_maya_sim$groups=="PO",]),
                                                            chull(per_maya_sim[per_maya_sim$groups=="PO",])[1]),],
            color="firebrick2") +  
  
  geom_path(data=per_maya_sim[per_maya_sim$groups=="FS",][c(chull(per_maya_sim[per_maya_sim$groups=="FS",]),
                                                            chull(per_maya_sim[per_maya_sim$groups=="FS",])[1]),],
            color="dodgerblue2") +
  
  geom_point(data=per_maya_sim[1712:2011,],aes(P0, P1, P2,
                                               colour=groups),show_guide = FALSE,size=2) 


pc




#
#
### ilr plot #####
#
#



z1 = 1/sqrt(2)*log(per_maya_sim$P2/per_maya_sim$P0)
z2 = 1/sqrt(6)*log((per_maya_sim$P0*per_maya_sim$P2)/(per_maya_sim$P1^2))

ilr_maya_sim = data.frame(z1=z1,z2=z2)

ilr_maya_sim$groups = per_maya_sim$groups




x_lab = expression(paste(z["11"]))

y_lab = expression(paste(z["12"]))


pd = ggplot(ilr_maya_sim[1712:2011,], aes(z1, z2, col = groups,shape=groups)) + 
  geom_point(size=2,show_guide = FALSE)+
  stat_chull(data=ilr_maya_sim[1:1711,],fill = NA,show_guide = FALSE) +
  scale_colour_manual(values =  c("darkgoldenrod2","dodgerblue2","firebrick2","darkorchid2","forestgreen")) +
  scale_shape_manual(values = c(17, 3, 7, 15,18))+
  theme(panel.border=element_blank(), axis.line=element_line()) + 
  theme(axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(), 
        panel.background=element_blank(),
        axis.ticks = element_blank(),
        axis.title.x=element_text(size=12,vjust=-0.5),
        axis.title.y=element_text(size=12,vjust=1.5),
        plot.title = element_text(lineheight=1, face="bold",hjust=0.5),
        axis.line.x=element_line(),
        axis.line.y=element_line()) +
  guides(shape = guide_legend(reverse = TRUE), color = guide_legend(reverse = TRUE)) +
  ggtitle("D")+
  xlab(x_lab)+ylab(y_lab)+coord_fixed(ratio=1)


pd



#
#
## legend ####
#
#



legend_data = data.frame(k0 = c(0,0.25,0.5,0.75,1),k1=c(1,0.5,0.5,0.25,0),k2=c(0,0.25,0,0,0))

cols = rev(c("forestgreen","darkgoldenrod2","darkorchid2","dodgerblue2","firebrick2"))

legend_data$family = cols

leg_lab = c("Unrelated (UN)","First cousins (FC)","Half sibs (HS), Avuncular (AV) \n or Grandparent-grandchild (GG)",
            "Full sibs (FS)","Parent-offspring (PO)") 

plegend = ggplot(legend_data,aes(x=legend_data[,1]+200,y=legend_data[,3]+200,xend=1,yend=1,
                                 shape=legend_data[,4], 
                                 colour=legend_data[,4], 
                                 group=legend_data[,4])) +
  geom_point(size=0)+ geom_blank()+
  xlab("") + 
  ylab("") +
  theme_bw() + theme(axis.ticks = element_blank()
                     , legend.position = c(.3,.5), legend.text=element_text(size=10),
                     legend.title=element_text(size=12,face = "bold"),
                     legend.key = element_rect(fill = 'white', colour = 'gray'),
                     panel.background = element_rect(fill = 'white', colour = 'white')
                     ,panel.grid.major = element_blank()
                     ,panel.grid.minor = element_blank()
                     ,panel.border = element_blank()
                     ,axis.line=element_blank(), axis.text = element_blank()) +
  scale_colour_manual(name="Family relationships",values=cols,labels=leg_lab,breaks=cols)+
  scale_colour_manual(name = "Family relationships",
                      labels = leg_lab,
                      values = rev(cols)) +   
  scale_shape_manual(name = "Family relationships",
                     labels = leg_lab,
                     values = c(18, 17, 15, 3,7))+
  guides(shape = guide_legend(reverse = TRUE), color = guide_legend(override.aes = list(size=4),
                                                                    reverse = TRUE))+
  scale_x_continuous(limits = c(0, 1))+
  scale_y_continuous(limits = c(0, 1))



plegend




#
#
#
### IBS panel: Figure 7 of the paper ####
#
#
#



plots = list()
plots[[1]] = pa
plots[[2]] = pb
plots[[3]] = pc
plots[[4]] = pd
plots[[5]] = plegend



lay_vec = 
  c(1,1,2,2,0,
    1,1,2,2,5,
    3,3,4,4,5,
    3,3,4,4,0)


layout <- matrix(lay_vec, nrow = 4,ncol=5, byrow = TRUE)


pdf("IBS_panel.pdf",width = 14,height = 10)
multiplot(plotlist = plots, layout = layout)
dev.off()



### 
### 
### 
### IBD Graphics #####
### 
### 
### 
### 




#
#
#
## (k0,k1)-plot ####
#
#
#


cluster = makeCluster(n_cores, type = "SOCK") # open n_cores in a cluster for parallel process


time2 <- Sys.time()

registerDoSNOW(cluster)

ibd_maya = foreach(i=1:300, .combine='rbind') %dopar% {       # around 3 minutes
  
  cotterman(X = Maya, pairs = i, colmarkers = 1:377, ibs = ibs_matrix )
  
  
}

stopCluster(cluster)

time1 <- Sys.time()


time2-time1





cluster = makeCluster(n_cores, type = "SOCK") # open n_cores in a cluster for parallel process

time4 <- Sys.time()

registerDoSNOW(cluster)

ibd_sim = foreach(i=1:1711, .combine='rbind') %dopar% {  # around 25 minutes
  
  cotterman(X = Maya_sim, pairs = i, colmarkers = 1:377, ibs = ibs_matrix_sim )
  
  
}

stopCluster(cluster)


time3 <- Sys.time()

time4-time3







ibd_maya_sim = rbind(ibd_sim,ibd_maya)


ibd_maya_sim$groups = per_maya_sim$groups

ibd_maya_sim$groups = as.factor(ibd_maya_sim$groups)

table(ibd_maya_sim$groups)



p <- seq(0, 1, by = 0.005) # inbreeding curve

inbreeding <- cbind(p^2, 2 * p * (1 - p), (1 - p)^2)

inbreeding = data.frame(k0=inbreeding[,1],k1=inbreeding[,2],k2=inbreeding[,3]) 

inbreeding$groups = "curve"

ibd_maya_sim = rbind(ibd_maya_sim,inbreeding)


x_lab =expression(paste(k["0"]))
y_lab= expression(paste(k["1"]))


pa_2 = ggplot(ibd_maya_sim[1712:2011,], aes(k0, k1, col = groups,shape=groups))+
  stat_chull(data=ibd_maya_sim[1:1711,],fill = NA,show_guide = FALSE) +
  scale_colour_manual(values =  c("darkgoldenrod2","dodgerblue2","firebrick2","darkorchid2","forestgreen")) +
  scale_shape_manual(values = c(3,17, 3, 7, 15,18))+
  coord_fixed(ratio=1,xlim=c(0.03,1),ylim=c(0.03,1))+
  theme(panel.border=element_blank(), axis.line=element_line()) + 
  theme(axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(), 
        panel.background=element_blank(),
        axis.ticks = element_blank(),
        axis.title.x=element_text(size=12,vjust=-0.5),
        axis.title.y=element_text(size=12,vjust=1.5),
        plot.title = element_text(lineheight=1, face="bold",hjust=0.5),
        axis.line.x=element_line(),
        axis.line.y=element_line()) +
  guides(shape = guide_legend(reverse = TRUE), color = guide_legend(reverse = TRUE)) +
  ggtitle("A")+
  xlab(x_lab)+ylab(y_lab)+
  geom_point(size=2,show_guide = FALSE)+
  geom_abline(intercept = 1,slope=-1,col="gray") +
  geom_line(data=ibd_maya_sim[2012:2212,],size=1,col="gray") 




pa_2



#
#
#
#### ternary diagram of k0,k1,k2 #####
#
#
#


x_lab=expression(k[0])
y_lab=expression(k[1])
z_lab=expression(k[2])



pb_2 = ggtern(data = ibd_maya_sim[1712:2011,],
              aes(k0, k1, k2,colour=groups,shape=groups)) + geom_mask()+
  
  theme_bw() +
  scale_colour_manual(values =  c("darkgoldenrod2","dodgerblue2","firebrick2","darkorchid2","forestgreen")) +
  scale_shape_manual(values = c(3,17, 3, 7, 15,18))+
  
  theme(plot.title = element_text( size=14,face="bold",hjust = 0.5)) +
  xlab(x_lab) + ylab(y_lab) + zlab(z_lab) +
  guides(color = "none", fill = "none", alpha = "none") + 
  labs(title="B") +
  
  geom_path(data=ibd_maya_sim[ibd_maya_sim$groups=="second",][c(chull(ibd_maya_sim[ibd_maya_sim$groups=="second",]),
                                                                chull(ibd_maya_sim[ibd_maya_sim$groups=="second",])[1]),],
            color="darkorchid2") +
  
  geom_path(data=ibd_maya_sim[ibd_maya_sim$groups=="FC",][c(chull(ibd_maya_sim[ibd_maya_sim$groups=="FC",]),
                                                            chull(ibd_maya_sim[ibd_maya_sim$groups=="FC",])[1]),],
            color="darkgoldenrod2") +
  
  geom_path(data=ibd_maya_sim[ibd_maya_sim$groups=="UN",][c(chull(ibd_maya_sim[ibd_maya_sim$groups=="UN",]),
                                                            chull(ibd_maya_sim[ibd_maya_sim$groups=="UN",])[1]),],
            color="forestgreen") +
  
  geom_path(data=ibd_maya_sim[ibd_maya_sim$groups=="PO",][c(chull(ibd_maya_sim[ibd_maya_sim$groups=="PO",]),
                                                            chull(ibd_maya_sim[ibd_maya_sim$groups=="PO",])[1]),],
            color="firebrick2") +  
  
  geom_path(data=ibd_maya_sim[ibd_maya_sim$groups=="FS",][c(chull(ibd_maya_sim[ibd_maya_sim$groups=="FS",]),
                                                            chull(ibd_maya_sim[ibd_maya_sim$groups=="FS",])[1]),],
            color="dodgerblue2") +
  
  geom_point(data=ibd_maya_sim[1712:2011,],aes(k0, k1, k2,
                                               colour=groups),show_guide = FALSE,size=2) + 
  geom_line(data=ibd_maya_sim[2012:2212,],size=1,color="gray") 



pb_2





#
#
#
### ilr plot of k0,k1,k2 ####
#
#
#



cluster = makeCluster(n_cores, type = "SOCK") # open n_cores in a cluster for parallel process


time2 <- Sys.time()

registerDoSNOW(cluster)

ilr_maya_max = foreach(i=1:300, .combine='rbind') %dopar% {       # around 3 minutes
  
  prob.rel = as.data.frame(probIBD(Maya,i,1:377,ibs_matrix))
  
  ibd.out <- optim(c(0,-0.6),ilr_maximum, method ="L-BFGS-B",lower=c(-Inf,-Inf),upper=c(Inf,-(sqrt(2)*log(2))/sqrt(3)))
  
  ibd.out$par  
  
}

stopCluster(cluster)

time1 <- Sys.time()


time2-time1





cluster = makeCluster(n_cores, type = "SOCK") # open n_cores in a cluster for parallel process

time4 <- Sys.time()

registerDoSNOW(cluster)

ilr_sim_max = foreach(i=1:1711, .combine='rbind') %dopar% {       # around 25 minutes
  
  prob.rel = as.data.frame(probIBD(Maya_sim,i,1:377,ibs_matrix_sim))
  
  ibd.out <- optim(c(0,-0.6),ilr_maximum, method ="L-BFGS-B",lower=c(-Inf,-Inf),upper=c(Inf,-(sqrt(2)*log(2))/sqrt(3)))
  
  ibd.out$par  
  
}

stopCluster(cluster)

time3 <- Sys.time()


time4-time3





ilr_maya_sim_max = rbind(ilr_sim_max,ilr_maya_max)


ilr_maya_sim_max = data.frame(ilr_maya_sim_max)

ilr_maya_sim_max$groups = per_maya_sim$groups

ilr_maya_sim_max$groups = as.factor(ilr_maya_sim_max$groups)

table(ilr_maya_sim_max$groups)




x_lab = expression(paste(z["11"]))

y_lab = expression(paste(z["12"]))



pc_2 = ggplot(ilr_maya_sim_max[1712:2011,], aes(-X1, X2, col = groups,shape=groups)) + 
  geom_point(size=2,show_guide = FALSE)+
  stat_chull(data=ilr_maya_sim_max[1:1711,],fill = NA,show_guide = FALSE) +
  scale_colour_manual(values =  c("darkgoldenrod2","dodgerblue2","firebrick2","darkorchid2","forestgreen")) +
  scale_shape_manual(values = c(17, 3, 7, 15,18))+
  coord_fixed(ratio=1)+
  theme(panel.border=element_blank(), axis.line=element_line()) + 
  theme(axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black"),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(), 
        panel.background=element_blank(),
        axis.ticks = element_blank(),
        axis.title.x=element_text(size=12,vjust=-0.5),
        axis.title.y=element_text(size=12,vjust=1.5),
        plot.title = element_text(lineheight=1, face="bold",hjust=0.5),
        axis.line.x=element_line(),
        axis.line.y=element_line()) +
  guides(shape = guide_legend(reverse = TRUE), color = guide_legend(reverse = TRUE)) +
  ggtitle("C")+
  xlab(x_lab)+ylab(y_lab)+ xlim(c(-15,1)) +
  geom_abline(intercept = -2/3*log(2),slope=0,col="gray")  


pc_2







#
#
#
## IBD panel: Figure 8 of the paper ####
#
#
#





plots = list()
plots[[1]] = pa_2
plots[[2]] = pb_2
plots[[3]] = pc_2
plots[[4]] = plegend




lay_vec = 
  c(1,1,1,1,2,2,2,2,
    3,3,3,3,3,3,4,4)


layout <- matrix(lay_vec, nrow = 2,ncol=8, byrow = TRUE)


pdf("IBD_panel.pdf",width = 14,height = 10)
multiplot(plotlist = plots, layout = layout)
dev.off()

