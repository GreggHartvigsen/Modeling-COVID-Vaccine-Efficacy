###########################################################################
# Code for the paper: 
#   Hartvigsen, G. and Y. Dimitroff.
#   Modeling the SARS-CoV-2 Epidemic and the Efficacy of
#   Different Vaccines Across Different Network Structures
#   PLoS One 2025
# Gregg Hartvigsen
# hartvig@geneseo.edu
# SUNY Geneseo
#-----------------------------------

make.network = function(net.type,N,k,rewire.p){
  if (net.type == "Lattice") {
    g = make_lattice(dimvector = c(sqrt(N),sqrt(N)), circular = T)
  }
  if (net.type == "SW") {
    g = sample_smallworld(dim = 1, size = N, 
                          nei = round(k/2,0),p = rewire.p)
  }
  if (net.type == "Random") {
    g = sample_gnm(n = N,m = N*k/2) # random graph
  }
  if (net.type == "PA") {
    g = sample_pa(N, m = k/2, directed = F) # Barabasi - preferential attachment
  }
  return(g)
}

#--------------------------------------------------------
make.inds = function(N) {
  # inds$states
  state = rep("S", N)
  dayS = dayE = dayI = dayR = rep(0, N)
  ind.states = data.frame(state, dayS, dayE, dayI, dayR, stringsAsFactors = FALSE)
  #---------------------
  # inds$strains
  dayV1 = rep(0, N) # day received 1st vaccination
  dayV2 = rep(0, N) # day received 2nd vaccination
  ind.vacc.days = data.frame(dayV1, dayV2, stringsAsFactors = FALSE)
  #---------------------
  # inds$strains
  S1 = rep(-1,N)
  ind.inf.strains = data.frame(S1, stringsAsFactors = FALSE)
  #---------------------
  # inds$Vstrains
  VS1 = VS2 = rep(-1,N) # inf strains
  ind.vacc.strains = data.frame(VS1, VS2, stringsAsFactors = FALSE)
  #---------------------
  # Combine all dataframes into one list
  inds = list(ind.states, ind.inf.strains, ind.vacc.days, ind.vacc.strains)
  names(inds) = c("states","strains","daysV","Vstrains")
  return (inds)
}

#--------------------------------------------------------
innoculate.network = function(inds, num.E.T1, num.I.T1, Day) {
  if (num.E.T1 > 0) {
    innoculate = sample(which(inds$states$state == "S"), num.E.T1) # get num.infected.T1 "S" inds to be "I"
    inds$states$state[innoculate] = "E" # make exposed on day 1
    inds$states$dayE[innoculate] = Day # infected on day 1
    inds$strains$S1[innoculate] = 0 # give them strain 0
  }
  if (num.I.T1 > 0) {
    innoculate = sample(which(inds$states$state == "S"), num.I.T1) # get num.infected.T1 "S" inds to be "I"
    inds$states$state[innoculate] = "I" # make exposed on day 1
    inds$states$dayI[innoculate] = Day # infected on day 1
    inds$strains$S1[innoculate] = 0 # give them strain 0
  }
  return(inds)
}

#--------------------------------------------------------
get.states = function(inds) {
  NS = length(which(inds$states$state == "S"))
  NE = length(which(inds$states$state == "E"))
  NI = length(which(inds$states$state == "I"))
  NR = length(which(inds$states$state == "R"))
  return(c(NS,NE,NI,NR))
}

#---------------------------------------------------------
get.tot.num.strains = function(inds) {
  strains = as.numeric(unlist(inds$strains))
  strains = strains[-which(strains == -1)]
  return(length(unique(strains)))
}

get.num.cur.strains = function(inds) {
  EI.inds = which(inds$states$state == "E" | inds$states$state == "I")
  all.strains.now = numeric(length(EI.inds))
  if (length(EI.inds) > 0) {
    for (i in 1:length(EI.inds)) {
      strains = inds$strains[EI.inds[i],] # gets all strains for current ind
      strains = strains[length(strains)] # gets the last one, which should be current
      all.strains.now[i] = strains
    }
  }
  num.cur.strains = length(unique(all.strains.now))
  return(num.cur.strains)
}

#--------------------------------------------------------
make.time.step.data = function(inds) {
  temp.states = get.states(inds)
  tot.num.strains = get.tot.num.strains(inds)
  num.cur.strains = get.num.cur.strains(inds)
  dat = matrix(c(1, temp.states,tot.num.strains,num.cur.strains),ncol = 7)
  dat = as.data.frame(dat)
  names(dat) = c("Day","S","E","I","R","Num Curr Strains","Tot Num Strains")
  return(dat)
}

#--------------------------------------------------------
update.inds = function(inds, Day, num.days.E, num.days.I, num.days.R) {
  # Find all the infectious inds and try to recover
  exposed = which(inds$states$state == "E" & (Day - inds$states$dayE) >= num.days.E)
  if (length(exposed) > 0) {
    inds$states$state[exposed] = "I" # the become infectious
    inds$states$dayI[exposed] = Day
  }
  infectious = which(inds$states$state == "I" & (Day - inds$states$dayI) >= num.days.I)
  if (length(infectious) > 0) {
    inds$states$state[infectious] = "R"
    inds$states$dayR[infectious] = Day
  }
  recovered = which(inds$states$state == "R" & (Day - inds$states$dayR) >= num.days.R)
  if (length(recovered) > 0) {
    inds$states$state[recovered] = "S"
    inds$states$dayS[recovered] = Day
  }
  return (inds)
}

#-----------------------------------------------
n2bin = function(n, n.loci) {
  n = rev(as.numeric(intToBits(n)))
  return(n[-(1:(length(n) - n.loci))])
}
bin2n = function(n, n.loci) {
  n = paste(n, collapse = "")
  return (strtoi(n, base = 2))
}

#-----------------------------------------------
prop.possible.match = function(s1, s2, n.loci, n.loci.for.mismatch) {
  s1 = n2bin(s1, n.loci)
  s2 = n2bin(s2, n.loci)
  n.mismatches = sum(abs(s1-s2))
  partial.immunity = n.mismatches/n.loci.for.mismatch
  if (partial.immunity > 1) partial.immunity = 1
  # 1 means no protection, 0 = full protection
  return (1-partial.immunity)
}

#-------------------------------------------------------
# needs to include vaccine efficacy
vaccine.eff.protection = function(Day, I.strain, inds, ind.num, V1.days, 
                              vacc.efficacy, n.loci, n.loci.for.mismatch) {
  vacc.protection = 0 # not protected
  # vacc.strain = -1
  if (inds$daysV$dayV1[ind.num] > 0) { # received 1st vaccine
    vacc.strain = inds$Vstrains$VS1[ind.num]
    if (vacc.strain >= 0) { # has at least one vaccine
      if (inds$daysV$dayV2[ind.num] > 0) { # received 2nd vaccine
        vacc.strain = c(vacc.strain, inds$Vstrains$VS2[ind.num]) # add 2nd strain
        VE.row = Day - inds$daysV$dayV2[ind.num] 
        if (VE.row > V1.days) { # V2 is at max efficacy
          VE.row = V1.days
        }
        if (VE.row == 0) { # got vaccinated today!
          VE.row = 1
        }
        vacc.protection = vacc.efficacy[VE.row,2]
        
      } else { # only received 1st vaccine
        VE.row = Day - inds$daysV$dayV1[ind.num] 
        if (VE.row > V1.days) { # V2 is at max efficacy
          VE.row = V1.days
        }
        if (VE.row == 0) { # got vaccinated today!
          VE.row = 1
        }
        vacc.protection = vacc.efficacy[VE.row,1]
      }
    }
  }
  # Adjust vacc.protection based on how well vacc. matches
  v = prop.possible.match(I.strain, vacc.strain, n.loci, n.loci.for.mismatch)
  # 0 = not protected; 1 is perfect match at max efficacy (never?)
  vacc.protection = v*vacc.protection
  return (vacc.protection)
}

#------------------------------------------------
max.partial.immunity = function(I.strain, nbr.strains, n.loci, n.loci.for.mismatch) {
  nbr.strains = unique(as.numeric(nbr.strains)) # either vaccines or seens strains
  if (-1 %in% nbr.strains) {
    nbr.strains = nbr.strains[-which(nbr.strains == -1)]
  }
  if (length(nbr.strains) > 0) {
    matches = numeric(length(nbr.strains))
    for (i in 1:length(nbr.strains)) {
      matches[i] = prop.possible.match(I.strain, nbr.strains, n.loci, n.loci.for.mismatch)
    }
    return (max(matches)) # the maximum match (1 -> perfect match, 0 -> no match)
  } else {
    return (0) # nbr has no matches, so is completely susceptible
  }
}

#--------------------------------------------------------
spread = function(inds, Trans.prob, Day, doses.prop.per.day, V1.days, vacc.efficacy, 
                  mut.rate, n.loci, n.loci.for.mismatch) {
  infectious = which(inds$states$state == "I")
  if (length(infectious) > 0) {
    if (length(infectious) > 1) infectious = sample(infectious) # randomize them
    for (i in 1:length(infectious)) {
      I.ind = infectious[i] # ID of current infectious ind
      # get current inf. ind.'s strain
      I.strain = inds$strains[I.ind,]
      I.strain = I.strain[length(I.strain)] # get last strain
      # get all if this inds neighbors to challenge
      nbrs = unlist(ego(g, nodes = I.ind, mindist = 1))
      nbrs.to.infect = nbrs[which(runif(length(nbrs)) < Trans.prob)]
      # do only if some are getting challenged
      if (length(nbrs.to.infect) > 0) {
        nbrs.S = nbrs.to.infect[which(inds$states$state[nbrs.to.infect] == "S")] # nbrs that are S
        if (length(nbrs.S) > 0) {
          for (j in 1:length(nbrs.S)) {
            # Check for vaccinations/partial immunity
            # get closest strain match for each nbrs.R. "strain.immunity" holds prob. of NOT getting sick 
            #   1 = complete match, totally protected; 0 = no protection -> susceptible
            strain.immunity  = max.partial.immunity(I.strain, inds$strains[nbrs.S[j],],
                                                     n.loci, n.loci.for.mismatch)
            strain.immunity = 1 - strain.immunity # now strain.immunity: 0 for protected
            #   1 is full protection and 0 is no protection
            if (doses.prop.per.day > 0) {
              vaccine.immunity  = max.partial.immunity(I.strain, inds$Vstrains[nbrs.S[j],],
                                                       n.loci, n.loci.for.mismatch)
              vaccine.immunity = 1 - vaccine.immunity
              vacc.eff.protection = vaccine.eff.protection(Day, I.strain, inds, nbrs.S[j], V1.days,
                                                   vacc.efficacy, n.loci, n.loci.for.mismatch)
              vacc.eff.protection = 1 - vacc.eff.protection #  0 means protected
            } else {
              vacc.eff.protection = 1
            }
            if (runif(1) < strain.immunity * vacc.eff.protection) { # then NOT protected by having prior strains or useful vaccines
              # infect this nbrs.S[j] with strain 
              n.strains.columns = length(inds$strains[1,])
              nbr.strains = inds$strains[nbrs.S[j],]
              nbr.open.sites = which(nbr.strains == -1)
              if (length(nbr.open.sites) > 1) {
                nbr.open.site = nbr.open.sites[1] # just grab first location
              } else if (length(nbr.open.sites) == 0) { # if none the create new column
                temp = data.frame(matrix(rep(-1, length(inds$strains[,1])), ncol = 1))
                names(temp) = paste("S", n.strains.columns + 1, sep = "")
                inds$strains = cbind(inds$strains, temp)
                nbr.open.site = n.strains.columns + 1
              } else { # there's just one and use it
                nbr.open.site = nbr.open.sites
              }
              if (mut.rate > 0) {
                strain.temp = mutation(I.strain, mut.rate, n.loci)
              } else {
                strain.temp = I.strain
              }
              inds$strains[nbrs.S[j], nbr.open.site] = strain.temp # assign strain
              inds$states$state[nbrs.S[j]] = "E"
              inds$states$dayE[nbrs.S[j]] = Day
            }
          }
        }
      } # nbrs to infect
    } # for each infectious ind
  } # if there are infectious
  return (inds)
} # spread()

#---------------------------------------------------
get.vaccine.recipients = function(unvacc, num.to.vacc, vacc.strat, g) {
  if (vacc.strat == "random") {
    to.vaccinate = sample(unvacc, num.to.vacc)
    return (to.vaccinate)
  } else if (vacc.strat == "hubs") {
    # get the degrees in decreasing order
    verts = unvacc[order(degree(g)[unvacc], decreasing = T)]
    to.vaccinate = verts[1:num.to.vacc]
  return (to.vaccinate)
  }
  return(0)
}

#-----------------------------------------------------------
vaccinate = function (inds, Day, doses.prop.per.day, V1.days, 
                      vaccine.doses.per.person, vacc.strat, vacc.strain, vacc.only.S, g) {
  # if we're in here then we do at least one dose. A second dose is optional
  num.to.vacc = round(doses.prop.per.day * length(inds$states$state),0)
  if (num.to.vacc > 0) {
    if (vacc.only.S == TRUE) {
      unvacc = which(inds$daysV$dayV1 == 0 & inds$states$state == "S")
    } else {
      unvacc = which(inds$daysV$dayV1 == 0)
    }
    if (length(unvacc) > 0) {
      if (num.to.vacc > length(unvacc)) {
        num.to.vacc = length(unvacc)
      }
      to.vacc = get.vaccine.recipients(unvacc, num.to.vacc, vacc.strat, g) # should not already be vaccinated
      if (length(to.vacc) > 0) {
        for (i in 1:length(to.vacc)) { # for each ind to vaccinate add the strain and mark the day
          # get last strain for this ind
          inds$Vstrains[to.vacc[i],1] = vacc.strain # vacc strain to this ind
          inds$daysV$dayV1[to.vacc[i]] = Day # day the got this vaccine
        }
      }
    }
  }
  # Second dose - find all those needing second dose on this Day
  if (vaccine.doses.per.person == 2) {
    need.dose2 = which(inds$daysV$dayV1 > 0 & # have V1, and....
                         inds$daysV$dayV2 == 0 &  # haven't been V2'd, and...
                         (Day - inds$daysV$dayV1) == V1.days) # have waited long enough for V2
    if (length(need.dose2) > 0) {
      inds$daysV$dayV2[need.dose2] = Day
      inds$Vstrains$VS2[need.dose2] = vacc.strain
    }
  }
  return(inds)
}

# This function returns the probability that the vaccination is effective
#   for the person who received it. It is a sigmoid function

make.vacc.efficacy = function(V1.days, V2.days, dose1.max, dose2.max) {
  vacc.efficacy = matrix(0, nrow = V1.days, ncol = 2)
  vacc.efficacy[,1] = dose1.max/(1 + exp(-((1:(V1.days) - V1.days/2)/1.75)))
  vacc.efficacy[,2] = dose1.max + (dose2.max - dose1.max)/(1 + exp(-((1:(V2.days) - V2.days/2)/1.75)))
  return (vacc.efficacy)
}

#-----------------------------------------------
get.layout = function(net.type) {
  if (net.type == "PA") {
    my.layout = layout.fruchterman.reingold(g)
  } else if (net.type == "Lattice") {
    my.layout = layout.grid(g)
  } else {
    my.layout = layout_in_circle(g)
  }
  return (my.layout)
}

#-----------------------------------------------
plot.network = function(g, inds, Day, net.type, my.layout) {
  par(mfrow = c(1,1))
  my.cols = numeric(length(inds$states[,1]))
  my.cols[which(inds$states$state=="S")] = "green"
  my.cols[which(inds$states$state=="E")] = "orange"
  my.cols[which(inds$states$state=="I")] = "red"
  my.cols[which(inds$states$state=="R")] = "blue"
  plot(g, vertex.size = 5, vertex.color = my.cols,
       edge.color = "black", 
       main = paste("Day = ", Day),
       layout = my.layout, vertex.label = "")
  legend("topleft",legend = c("S","E","I","R"), 
         fill = c("green","orange","red","blue"))
}

#--------------------------------------------------------
make.ts.plot = function(time.step.data) {
  par(mar = c(5.1,4.1,4.1,2.1))
  matplot(time.step.data[,2:5], type = "l", lty = 1,lwd = 3, las = 1,
          ylim = c(0,N*1.3), ylab = "Abundance", xlab = "Day", cex.lab = 1.5,
          col = c("green","orange","red","blue"))
  leg.txt = c("S","E","I","R")	# this is the text for the legend
  legend("top",leg.txt, col = c("green","orange","red","blue"),
         lwd=3, horiz = T)
}

#-----------------------------------
summary.result = function(time.step.data){
  max.sick = max(rowSums(time.step.data[,3:4]))
  cat("Maximum number E + I at one time = ", max.sick,"\n")
  max.E
  max.I
  tot.sick = max(time.step.data$R)
  cat("The number of individuals that got sick = ",NR," out of ",N,"\n")
  cat("The epidemic lasted",length(time.step.data[,1]),"days\n")
}

#------------------------------------------------
mutation = function(strain, mut.rate, n.loci) {
  strain = n2bin(strain, n.loci)
  mut.sites = which(runif(n.loci) < mut.rate)
  if (length(mut.sites) > 0) {
    strain[mut.sites] = (strain[mut.sites] + 1) %% 2
  }
  strain = bin2n(strain,n.loci)
  return (strain)
}

time.step.graph = function(time.step.data) {
  matplot(time.step.data[,2:5], xlab = "Day", ylab = "N", las = 1, #log = "y",
          type = "l", lwd = 3, lty = 1, col = c("black","red","green","blue"),
          ylim = c(0, N*1.2))
  legend("top", horiz = T, legend = c("S","E","I","R"), 
         col = c("black","red","green","blue"), lty = 1, lwd = 5)
}

