# =============================================================================
# EXAMPLE SCRIPT — not called by the pipeline
# Legacy standalone alpha diversity script with hardcoded paths.
# Retained here as a reference only. Use the Nextflow pipeline instead.
# =============================================================================

####################
# Data Import Done #
####################

library(phyloseq)
library(stringr)
library(data.table)
library(vegan)
library(ggplot2)
library(grid) #We need grid to draw the arrows

abund_table<-read.csv("../../data/Test1/feature-table.tsv",sep="\t",header=T,row.names=1,check.names=FALSE)
meta_table<-read.csv("../../data/Test1/metadata.csv",header=T,row.names=1)
OTU_taxonomy<-read.csv("../../data/Test1/taxonomy.csv",header=T,row.names=1,check.names=FALSE)

exists("abund_table")
exists("meta_table")
exists("OTU_taxonomy")

#PARAMETERS ###########################
which_level="Otus" #Otus Genus Family Order Class Phylum
text_size=16
point_size=5
point_opacity=0.8
axis_text_size=12
strip_text_size=16
increment_divider=2
p_adjust_method="none" #P-value adjustment: BH | bonferroni | holm | fdr | none
exclude_pvalues_text_from_drawing=FALSE
legends_position_bottom=TRUE
exclude_legends=TRUE
pairwise_text_size=7
number_of_rows=1
legend_text_size=12
axis_title_size=16
height_image=10
width_image=25
use_provided_colors=FALSE
colours <- c("#ffa172", "#81fc76", "#68aeff","#c320d8","#2BCE48","#FFCC99","#808080","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#FFFF00",grey.colors(1000));
#/PARAMETERS ###########################

# Hypothesis 1
label="xRun"
meta_table$Groups<-as.factor(as.character(meta_table$run))
meta_table$Type<-as.factor(as.character(meta_table$sweetening.flow))
meta_table$Type2<-as.factor(meta_table$BED)
meta_table$Connections<-NULL

abund_table = t(abund_table)

#PARAMETERS CHANGE THE GROUPING COLUMN AS YOU DESIRE############################

#Align samples — keep only IDs present in both tables
common_samps<-intersect(rownames(abund_table),rownames(meta_table))
if(length(common_samps)==0) stop("No shared sample IDs between feature table and metadata.")
abund_table<-abund_table[common_samps,,drop=FALSE]
meta_table<-meta_table[common_samps,,drop=FALSE]
#After adjustment, get rid of OTUs that are all empty
abund_table<-abund_table[,colSums(abund_table)>0]
#Adjust OTU taxonomy
OTU_taxonomy<-OTU_taxonomy[colnames(abund_table),]

#COLLATE OTUS AT A PARTICULAR LEVEL#######################################
new_abund_table<-NULL
if(which_level=="Otus"){
  new_abund_table<-abund_table
} else {
  list<-unique(OTU_taxonomy[,which_level])
  new_abund_table<-NULL
  for(i in list){
    tmp<-data.frame(rowSums(abund_table[,rownames(OTU_taxonomy)[OTU_taxonomy[,which_level]==i]]))
    if(i==""){colnames(tmp)<-c("__Unknowns__")} else {colnames(tmp)<-paste("",i,sep="")}
    if(is.null(new_abund_table)){new_abund_table<-tmp} else {new_abund_table<-cbind(tmp,new_abund_table)}
  }
}
new_abund_table<-as.data.frame(as(new_abund_table,"matrix"))
abund_table<-new_abund_table
#/COLLATE OTUS AT A PARTICULAR LEVEL#######################################

grouping_column<-"Groups"

#Calculate Richness
R<-vegan::rarefy(abund_table,min(rowSums(abund_table)))
df_R<-data.frame(sample=names(R),value=R,measure=rep("Richness",length(R)))

#Calculate Shannon entropy
H<-vegan::diversity(abund_table)
df_H<-data.frame(sample=names(H),value=H,measure=rep("Shannon",length(H)))

#Calculate Simpson diversity index
simp <- vegan::diversity(abund_table, "simpson")
df_simp<-data.frame(sample=names(simp),value=simp,measure=rep("Simpson",length(simp)))

#Calculate Fisher alpha
alpha <- vegan::fisher.alpha(abund_table)
df_alpha<-data.frame(sample=names(alpha),value=alpha,measure=rep("Fisher alpha",length(alpha)))

#Calculate Pielou's evenness
S <- vegan::specnumber(abund_table)
J <- H/log(S)
df_J<-data.frame(sample=names(J),value=J,measure=rep("Pielou's evenness",length(J)))

#Uncomment to retain everything
df<-rbind(df_R,df_H,df_simp,df_alpha,df_J)

#Write all the metrics in a file for further analyses elsewhere
data_to_write<-data.frame(df_R[,"value",drop=F],df_H[,"value",drop=F],df_simp[,"value",drop=F],df_alpha[,"value",drop=F],df_J[,"value",drop=F])
colnames(data_to_write)<-c("Richness","Shannon","Simpson","FisherAlpha","PielouEvenness")
write.csv(data_to_write,paste("Diversity_",which_level,"_",label,".csv",sep=""))
#/Write all the metrics in a file for further analyses else where


rownames(df)<-NULL

#Incorporate categorical data in df
df<-data.frame(df,meta_table[as.character(df$sample),])

#To do anova, we will convert our data.frame to data.table

#Since we can't pass a formula to data.table, I am creating
#a dummy column .group. so that I don't change names in the formula
dt<-data.table(data.frame(df,.group.=df[,grouping_column]))

#I am also specifying a p-value cutoff for the ggplot2 strips
pValueCutoff<-0.05
pval<-dt[, list(pvalue = sprintf("%.2g",
                                 tryCatch(summary(aov(value ~ .group.))[[1]][["Pr(>F)"]][1],error=function(e) NULL))),
         by=list(measure)]

#Filter out pvals that we don't want
pval<-pval[!pval$pvalue=="",]
pval<-pval[as.numeric(pval$pvalue)<=pValueCutoff,]

#I am using sapply to generate significances for pval$pvalue using the cut function.
pval$pvalue<-sapply(as.numeric(pval$pvalue),function(x){as.character(cut(x,breaks=c(-Inf, 0.001, 0.01, 0.05, Inf),label=c("***", "**", "*", "")))})

#Update df$measure to change the measure names if the grouping_column has more than three classes
if(length(unique(as.character(meta_table[,grouping_column])))>2){
  df$measure<-as.character(df$measure)
  if(dim(pval)[1]>0){
    for(i in seq(1:dim(pval)[1])){
      df[df$measure==as.character(pval[i,measure]),"measure"]=paste(as.character(pval[i,measure]),as.character(pval[i,pvalue]))
    }
  }
  df$measure<-as.factor(df$measure)
}

#Get all possible combination of values in the grouping_column
s<-combn(unique(as.character(df[,grouping_column])),2)

#df_pw will store the pair-wise p-values
df_pw<-NULL
for(k in unique(as.character(df$measure))){
  bas<-max(df[(df$measure==k),"value"])
  inc<-0.05*(bas-min(df[(df$measure==k),"value"]))
  bas<-bas+inc

  # First pass: collect raw p-values for every pair in this metric
  pw_rows<-lapply(seq_len(ncol(s)),function(l){
    sub<-df[(df$measure==k) & (df[,grouping_column]==s[1,l] | df[,grouping_column]==s[2,l]),]
    pv<-tryCatch(
      if(is.null(meta_table$Connections)){
        summary(aov(as.formula(paste("value ~",grouping_column)),data=sub))[[1]][["Pr(>F)"]][1]
      } else {
        summary(aov(as.formula(paste("value ~",grouping_column,"+",paste("Error(","Connections/",grouping_column,")",sep=""))),data=sub))[[1]][[1]][["Pr(>F)"]][1]
      },
      error=function(e) NA_real_
    )
    list(from=s[1,l],to=s[2,l],pvalue=pv)
  })

  # Apply multiple-testing correction across all pairs within this metric
  raw_pv<-sapply(pw_rows,`[[`,"pvalue")
  adj_pv<-if(p_adjust_method=="none") raw_pv else p.adjust(raw_pv,method=p_adjust_method)

  # Second pass: retain significant pairs and assign bracket y-positions
  for(l in seq_len(ncol(s))){
    if(!is.na(adj_pv[l]) && adj_pv[l]<pValueCutoff){
      tmp<-c(k,pw_rows[[l]]$from,pw_rows[[l]]$to,bas,
             paste(sprintf("%.2g",adj_pv[l]),"",sep=""))
      if(is.null(df_pw)){df_pw<-tmp}else{df_pw<-rbind(df_pw,tmp)}
      bas<-bas+inc
    }
  }
}

if(!is.null(df_pw)){
  if(sum(class(df_pw) %in% c("character"))>0){
    df_pw<-t(as.matrix(df_pw))
  }
  df_pw<-data.frame(row.names=NULL,df_pw)
  names(df_pw)<-c("measure","from","to","y","p")
}


#Draw the boxplots
p<-NULL
if(!"Type" %in% colnames(meta_table)){
    if(!"Type2" %in% colnames(meta_table)){
      p<-ggplot(aes_string(x=grouping_column,y="value",color=grouping_column,group=grouping_column),data=df)
    } else {
      p<-ggplot(aes_string(x=grouping_column,y="value",color=grouping_column,group=grouping_column,shape="Type2"),data=df)
    }
  } else {
  p<-ggplot(aes_string(x=grouping_column,y="value",color=grouping_column,group=grouping_column,shape="interaction(Type,Type2)"),data=df)
}

p<-p+geom_boxplot(outlier.size=0,show.legend=FALSE)#+geom_jitter(position = position_jitter(height = 0, width=0),show.legend=FALSE)
p<-p+theme_bw()
p<-p+geom_point(size=point_size,alpha=point_opacity)
p<-p+facet_wrap(~measure,scales="free_y",nrow=number_of_rows)+ylab("Observed Values")+xlab("Samples")

if(!is.null(df_pw)){
  #This loop will generate the lines and signficances
  for(i in 1:dim(df_pw)[1]){
    p<-p+geom_path(inherit.aes=F,aes(x,y),data = data.frame(x = c(which(levels(df[,grouping_column])==as.character(df_pw[i,"from"])),which(levels(df[,grouping_column])==as.character(df_pw[i,"to"]))), y = c(as.numeric(as.character(df_pw[i,"y"])),as.numeric(as.character(df_pw[i,"y"]))), measure=c(as.character(df_pw[i,"measure"]),as.character(df_pw[i,"measure"]))), color="black",lineend = "butt",arrow = arrow(angle = 90, ends = "both", length = unit(0.1, "inches")))
    p<-p+geom_text(inherit.aes=F,aes(x=x,y=y,label=label),data=data.frame(x=(which(levels(df[,grouping_column])==as.character(df_pw[i,"from"]))+which(levels(df[,grouping_column])==as.character(df_pw[i,"to"])))/2,y=as.numeric(as.character(df_pw[i,"y"])),measure=as.character(df_pw[i,"measure"]),label=as.character(cut(as.numeric(as.character(df_pw[i,"p"])),breaks=c(-Inf, 0.001, 0.01, 0.05, Inf),label=c("***", "**", "*", "")))),size=pairwise_text_size)
    if(exclude_pvalues_text_from_drawing){
      p<-p+geom_text(inherit.aes=F,aes(x=x,y=y,label=label),data=data.frame(x=(which(levels(df[,grouping_column])==as.character(df_pw[i,"from"]))+which(levels(df[,grouping_column])==as.character(df_pw[i,"to"])))/2,y=as.numeric(as.character(df_pw[i,"y"])),measure=as.character(df_pw[i,"measure"]),label=paste("p=",as.character(as.numeric(as.character(df_pw[i,"p"])))),sep=""),size=pairwise_text_size,vjust=-1)
    }
  }
}
if(use_provided_colors){
  p<-p+scale_color_manual(grouping_column,values=colours)
}
#if crashes at panel.margin change it to panel.spacing, if crashes at panel.spacing, change it to panel.margin
p<-p+theme(strip.background = element_rect(fill="white"))+theme(panel.spacing = unit(2, "lines"),
                                                                strip.text = element_text(size=strip_text_size),
                                                                legend.text=element_text(size=legend_text_size),
                                                                text = element_text(size=text_size),
                                                                axis.text=element_text(size=axis_text_size),
                                                                axis.title=element_text(size=axis_title_size),
                                                                axis.text.x = element_text(angle = 90, hjust = 1))
if(legends_position_bottom){
  p<-p+theme(legend.key = element_blank(),  #removes the box around each legend item
             legend.position = "bottom", #legend at the bottom
             legend.direction = "horizontal",
             legend.box = "horizontal",
             legend.box.just = "centre")
}
if(exclude_legends){
  p<-p+guides(colour="none")
}

p <- p + theme(legend.position = "none")

#pdf(paste("ANOVA_diversity_",which_level,"_",label,".pdf",sep=""),height=height_image,width=width_image)
print(p)
#dev.off()

