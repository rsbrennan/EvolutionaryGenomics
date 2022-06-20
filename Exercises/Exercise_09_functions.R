plot_simulation <- function(initial.population.size,nloci,n.gens, migration.rate, selection,  plot){

  initial.population.size=400 #Initial population size
  if (selection == "none"){
    message("simulating with no selection")
    ########################
    ### Don't edit ###
    ########################
    e.v=0.001                    #Stochastic environmental variant
    d.v=1                       #Stochastic demographic variant
    bvs = t(array( seq(0,1, length = n.alleles.per.locus) ,c(n.alleles.per.locus,     add.loci))) #Breeding value of additive loci
    sex.ratio <- 0.5
    set.seed(2)
    start.1 <- start.2 <- initial.struct(initial.population.size,n.loci,n.alleles.per.locus   )
    start <- list(start.1, start.1) #Initial populations
    
    
    #b0 Maximum generated offspring
    #b1 Phenotypic optima
    #b2 Variance of the fitness curve
    #b3 Density-dependent demographic variant
    ##Parameters of the fitness function
    param.w1 <- list(b0 = 6,b1 = 0.5, b2 = 0.5, b3 = 0.01, d.v = d.v, add.loci = add.loci)
    param.w2 <- list(b0 = 6,b1 = 0.5, b2 = 0.5, b3 = 0.01, d.v = d.v, add.loci = add.loci)
    
    ##Parameters of the phenotype function
    param.z1 <- list(sex.ratio, fitness.pos, bvs, add.loci, e.v)
    param.z2 <- param.z1
    
  }
  
  
  if (selection == "weak"| selection == "strong"){
    ########################
    ### Don't edit ###
    ########################
    e.v=0.001                    #Stochastic environmental variant
    d.v=1                       #Stochastic demographic variant
    bvs = t(array( seq(0,1, length = n.alleles.per.locus) ,c(n.alleles.per.locus,     add.loci))) #Breeding value of additive loci
    sex.ratio <- 0.5
    set.seed(2)
    start.1 <- start.2 <- initial.struct(initial.population.size,n.loci,n.alleles.per.locus   )
    start <- list(start.1, start.1) #Initial populations
    
    
    #b0 Maximum generated offspring
    #b1 Phenotypic optima
    #b2 Variance of the fitness curve
    #b3 Density-dependent demographic variant
    ##Parameters of the fitness function
    if(selection == "strong"){
      message("simulating with strong selection")
      
    param.w1 <- list(b0 = 8,b1 = 0.25, b2 = 0.5, b3 = 0.01, d.v = d.v, add.loci = add.loci)
    param.w2 <- list(b0 = 8,b1 = 0.75, b2 = 0.5, b3 = 0.01, d.v = d.v, add.loci = add.loci)
    }
    if(selection == "weak"){
      message("simulating with weak selection")
      param.w1 <- list(b0 = 8,b1 = 0.4, b2 = 0.5, b3 = 0.01, d.v = d.v, add.loci = add.loci)
      param.w2 <- list(b0 = 8,b1 = 0.6, b2 = 0.5, b3 = 0.01, d.v = d.v, add.loci = add.loci)
    }
    ##Parameters of the phenotype function
    param.z1 <- list(sex.ratio, fitness.pos, bvs, add.loci, e.v)
    param.z2 <- param.z1
    
  }
  
  
  
  
  ###################
  ### Simulations ###
  ###################
  example <- evolve(
    x = start, 
    time = n.gens, 
    type = "additive", 
    recombination = "map", 
    recom.rate = recom.map, 
    migration.rate = migration.rate, 
    param.z = list(param.z1, param.z2), 
    param.w = list(param.w1, param.w2)
  )

  
  ###################
  ###     Fst     ###
  ###################
  ##Commputation of Fst values
  #We should first convert the output from class 'struct' to class 'loci'
  example.loci <- struct2pegas(example)
  
  #Estimating Fst with library 'pegas'
  FST <- pegas::Fst(example.loci)
  mean(FST[,2], na.rm=T)
  ####################
  #Figure
  if(plot ==TRUE){
  plot(loci.pos, FST[,"Fst"], pch=16, bty="l", xlab="Loci position", ylab="Fst", ylim=c(0,1), type="l", xaxs="i", yaxs="i", las=1, lwd=2  )
  # mark position of fitness alleles
  if(selection == "strong" | selection == "weak" ){
    points(x=fitness.pos, y=rep(0.02, length(fitness.pos)), col="seagreen2", pch=19)
  }

  }
  
  return(mean(FST[,2], na.rm=T))
  #return(head(FST))
  
}