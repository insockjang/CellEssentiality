# This should be working in your local machine otherwise you have to open cytoscape in your server, then it will be okay to run
load("~/Synthetic_Lethality_Prediction/OV/NetProp_CGH_MUT.Rdata")
# 
load("~/NetworkAnalysis_RandomWalk/STRING/AdjMat_90percent.Rdata")

CUTOFF = 0.0005
gene.apair.CGH<-names(which(Test[,1]>=CUTOFF))
gene.anotherpair.CGH.gain<-names(which(Test[,2]>=CUTOFF))
gene.anotherpair.CGH.loss<-names(which(Test[,3]>=CUTOFF))
gene.apair.MUT<-names(which(Test[,4]>=CUTOFF))
gene.anotherpair.MUT<-names(which(Test[,5]>=CUTOFF))

seed.apair.CGH = rownames(f.0)[which(f.0[,1]==1)]
seed.anotherpair.CGH.gain = rownames(f.0)[which(f.0[,2]==1)]
seed.anotherpair.CGH.loss = rownames(f.0)[which(f.0[,3]==1)]
seed.apair.MUT = rownames(f.0)[which(f.0[,4]==1)]
seed.anotherpair.MUT = rownames(f.0)[which(f.0[,5]==1)]

# check how many neighbors are identified by network propagation
rbind(c(length(gene.apair.CGH),
        length(gene.anotherpair.CGH.gain),
        length(gene.anotherpair.CGH.loss),
        length(gene.apair.MUT),
        length(gene.anotherpair.MUT)),
      c(length(which(f.0[,1]==1)),
        length(which(f.0[,2]==1)),
        length(which(f.0[,3]==1)),
        length(which(f.0[,4]==1)),
        length(which(f.0[,5]==1))))

# check if seeds are included
setdiff(rownames(f.0)[which(f.0[,1]==1)],gene.apair.CGH)
setdiff(rownames(f.0)[which(f.0[,2]==1)],gene.anotherpair.CGH.gain)
setdiff(rownames(f.0)[which(f.0[,3]==1)],gene.anotherpair.CGH.loss)
setdiff(rownames(f.0)[which(f.0[,4]==1)],gene.apair.MUT)
setdiff(rownames(f.0)[which(f.0[,5]==1)],gene.anotherpair.MUT)

hit.name<-union(gene.apair.CGH,union(gene.anotherpair.CGH.gain,union(gene.anotherpair.CGH.loss,union(gene.apair.MUT,gene.anotherpair.MUT))))
seed.name<-union(seed.apair.CGH,union(seed.anotherpair.CGH.gain,union(seed.anotherpair.CGH.loss,union(seed.apair.MUT,seed.anotherpair.MUT))))

# Venn Diagram for 4 sets : candidate(apair), MUT(NP), CGH.gain(NP), CGH.loss(NP)
library(Vennerable)
VD1<-list(Candidate.CGH.NP = gene.apair.CGH,
          CGH.gain.NP = gene.anotherpair.CGH.gain,
          CGH.loss.NP = gene.anotherpair.CGH.loss)
VD2<-list(Candidate.MUT.NP = gene.apair.MUT,
          MUT.NP = gene.anotherpair.MUT)

Vstem1 <- Venn(VD1)
Vstem2 <- Venn(VD2)
png("~/Synthetic_Lethality_Prediction/OV/VD_NP.CGH.png",width = 640, height = 640)
par(cex =2)
plot(Vstem1,  show = list(Universe = FALSE))
# grid.text("Common Hits", vp = viewport(x=0.5, y=0.95, w=unit(1, "npc"), h=unit(1, "npc")),gp = gpar(cex=2))
dev.off()

png("~/Synthetic_Lethality_Prediction/OV/VD_NP.MUT.png",width = 640, height = 640)
par(cex =2)
plot(Vstem2,  show = list(Universe = FALSE))
dev.off()

VD.seed1<-list(Candidate.CGH = seed.apair.CGH,
               CGH.gain = seed.anotherpair.CGH.gain,
               CGH.loss = seed.anotherpair.CGH.loss)

VD.seed2<-list(Candidate.MUT = seed.apair.MUT,
               MUT = seed.anotherpair.MUT)

Vstem.seed1 <- Venn(VD.seed1)
Vstem.seed2 <- Venn(VD.seed2)

png("~/Synthetic_Lethality_Prediction/OV/VD_seed.CGH.png",width = 640, height = 640)
par(cex =2)
plot(Vstem.seed1,  show = list(Universe = FALSE))
dev.off()

png("~/Synthetic_Lethality_Prediction/OV/VD_seed.MUT.png",width = 640, height = 640)
par(cex =2)
plot(Vstem.seed2,  show = list(Universe = FALSE))
dev.off()


########################## This part should be carefully run in belltown, otherwise you might run this part in your  local machine(WARNING!!! Cytoscape should be running backend in order to generate cytoscape plot)
require(graph)
myMat<-MAT
# define hit.name
hit.name<-union(gene.apair.MUT,seed.anotherpair.MUT)
length(intersect(rownames(myMat),hit.name))
a<-match(hit.name,rownames(myMat))
b<-a[!is.na(a)]
hid<-myMat[b,b]
g<-as(hid,"graphNEL")


Gold<-gene.apair.MUT
Orange<-seed.apair.MUT
Red<-seed.anotherpair.MUT

color.red <- intersect(g@nodes,Red)
color.gold <- intersect(g@nodes,Gold)
color.orange <- intersect(g@nodes,Orange)



library(Rgraphviz)
library(RCytoscape)
g <- initEdgeAttribute (g, "weight", "numeric", 3.33)
g <- initNodeAttribute (g, "label", "char", "undefined")



cy = CytoscapeConnection()
pluginVersion (cy)
cy = CytoscapeConnection ()

window.title = 'OV_MUT_Synthetic Lethality Pairs'

if (window.title %in% as.character (getWindowList (cy)))
  deleteWindow (cy, window.title)
cw = new.CytoscapeWindow (window.title, g)

displayGraph(cw)

redraw(cw) 


for(k in 1:length(color.gold)){
  setNodeColorDirect(cw,color.gold[k],"FFFF11")
}
for(k in 1:length(color.orange)){
  setNodeColorDirect(cw,color.orange[k],"FFA500")
}

for(k in 1:length(color.red)){
  setNodeColorDirect(cw,color.red[k],"FF0000")
}



#orangered FF4500
#purple 800080
redraw(cw) 


for(k in 1:length(color.red)){
  setNodeSizeDirect(cw,color.red[k],80)
}

redraw(cw)

setDefaultEdgeColor(cw, "#D3D3D3")
setDefaultBackgroundColor (cw, '#FFFFFF')
redraw(cw)

setDefaultEdgeLineWidth(cw,0)
# edge color change
edge.names = as.character(cy2.edge.names(cw@graph))
edgenames = as.character(names(cy2.edge.names(cw@graph)))
K<-strsplit(edgenames,"~")
KK<-list()
for(k in 1:length(K)){
  k1<-match(K[[k]][1],rownames(hid))
  k2<-match(K[[k]][2],rownames(hid))
  if(hid[k1,k2]>0){
    setEdgeLineWidthDirect(cw,edge.names[k],hid[k1,k2]/100)
  }  
}

redraw(cw)


################# shortest path
require(igraph)

hid1<-hid
# make edge weight to be identical, otherwise shortest path considers the weights and looks weird.
hid1[which(hid1!=0)]<-1
g<-as(hid1,"graphNEL")
gg<-igraph.from.graphNEL(g) 
shortest.paths(gg)

### do further to find P2's synthetic lethality candidate pairs' subset 