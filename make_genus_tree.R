make_genus_tree <- function (phytree,placement = placement,species = species) { 
  phytree$tip.label <- gsub("[.]","_",phytree$tip.label)
  species <- gsub("[.]","_",species)
  species <- gsub(" ","_",species)
  
  genus.tree <- unique(unlist(lapply(strsplit(phytree$tip.label,"_"),function(x) x[[1]])))
  
  species.representative <- NULL
  
  for (i in 1:length(genus.tree)) {
    genus.tip <- phytree$tip.label[grep(genus.tree[[i]],phytree$tip.label)]
    temp <- genus.tip[sample(1:length(genus.tip),1)]
    species.representative <- cbind(temp,species.representative)
  }
  
  genustree <- keep.tip(phytree,species.representative)
  genustree$tip.label <- word(genustree$tip.label,sep="_",1)
  phytree1 <- genustree
  phytree1$tip.label <- paste0(phytree1$tip.label,"_")
  
  for (k in 1:nrow(placement)) {
    message(k)
    temp_row <- placement[k,]
    Add_tip <- paste0(temp_row[["Genus"]],"_")
    Sister <- c(temp_row[["Sister1"]],temp_row[["Sister2"]])
    Sister <- Sister[Sister != ""]
    
    ii <- pos <- list()
    
    for (l in 1:length(Sister)) {
      ii[[l]] <- grep(paste0(Sister[[l]],"_"), phytree1$tip.label)
    }
    
    pos <- ifelse(length(unlist(ii)) > 1,
                  findMRCA(phytree1,phytree1$tip.label[unlist(ii)]),
                  ii[[l]])
    
    nn <- pos
    
    if (length(Sister) > 1) {
      message("Sister >1")
      tt <- splitTree(phytree1, list(node = nn, bp =  0))
      wherenode <- tt[[2]]$edge[1,1]
      tt[[2]] <- bind.tip(tt[[2]],Add_tip,where=wherenode,position=0)
      phytree1 <- paste.tree(tt[[1]], tt[[2]])
    } else {
      message("ii =1")
      phytree1 <- bind.tip(phytree1, Add_tip, where = nn, position = 0)
    }
  }
  
  for (sp in unique(species)) {
    message(sp)
    phytree1 <- add.species.to.genus(phytree1,sp,where="root")
  }
  
  phytree1 <- drop.tip(phytree1,tip=phytree1$tip.label[!phytree1$tip.label %in% unique(species)])
  return(phytree1)
}
