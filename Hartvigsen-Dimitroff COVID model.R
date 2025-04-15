rand.seed = 1100 
###########################################################################
#   Hartvigsen, G. and Y. Dimitroff.
#   Modeling the SARS-CoV-2 Epidemic and the Efficacy of
#   Different Vaccines Across Different Network Structures
#   PLoS One 2025
# Gregg Hartvigsen
# hartvig@geneseo.edu
# SUNY Geneseo
#-----------------------------
#	This program builds and runs a discrete, individual-based
# version of an SEIR model. (S)usceptible individuals become
# (E)xposed when in contact with an (I)nfectious individual 
# with a probability (Trans.prob), derived from 
# a global parameter for Ro (R.naught). Infectious 
# individuals become (R)ecovered after num.days.I days.
#
# Vaccination assumptions include that vacc. begin when population
# reaches a threshold number of infectious inds (vacc.crit.pop so
# far has been set to zero.). Within an individual the vaccine
# efficacy increases over time according to a sigmoidal function. 
# This assumes it takes 21 days  for the first dose to become
# max.VE1 effective. Then another 21 days to become max.VE2
# effective. Vaccine effectiveness does not decrease with time.

# Currently vaccinating people regardless of their state. 

# This version keeps track of individuals in a list with 4
# dataframes. The structure is:
#
# inds:
#   states:
#     state (current state of ind)
#     dayS, dayE, dayI, dayR (days became these states)
#   strains:
#     S1 (this is the infected strain. New, evolved strains are added dynamically)
#   daysV: (when ind got each of the 2 vaccines. No further boosters are allowed now)
#     dayV1 (day ind got V1)
#     dayV2 (day ind got V2)
#   strainsV (vaccination strains recieved)
#     V1
#     V2

# The list "inds" is sent around and updated by various functions.

# Individuals are protected from reinfections for num.days.R days.
#   After this they return to the S class and may be infected (or
#   protected from their vaccination and previous strains). This
#   means that the number of R at any time is NOT the number that
#   have gotten the illness.
        
#------------------------------------------------------------
# This installs igraph (only if you don't have it) and loads it.
if (!require(igraph)) {
  install.packages("igraph")
  library(igraph)
} 
if (!require(plotrix)) {
  install.packages("plotrix")
  library(plotrix)
}
#------------------------------------------------------------
set.seed(rand.seed)
ver.num = 1.2
source("Vacc-Evol-fun1.2.R")
animate.network = F # do ONLY if doing 1 single run
Net.vis.max = 100 # don't visualize a networks larger than this (prob <= 250)
save.time.step.data = F # do ONLY if doing 1 rep
time.step.data.name = "time step - fewest sick.csv"
save.summary.data = T
save.summary.data.file.name = "sum1.csv"
show.time.step.graph = F # do ONLY if doing 1 rep
show.n.strains.graph = F # do ONLY if doing 1 rep
#VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
# Global parameters
N = 25000 # number of vertices
R.naught = c(2.5, 5.0) # R.naught (required range: 0 <= R.naught <= k)
num.E.T1 = 3
num.I.T1 = 3 # number  I at beginning
num.days.E = 3
num.days.I = 10 # number of days individuals are infectious
num.days.R = 10000 # c(30, 60, 100000) # An array. After this many days return to S state (very large # means don't)
net.type = c("SW") # "Lattice", Small-world ("SW"), "Random", and "PA" supported
rewire.p =  c(0.004, 0.01,  0.03, 0.1, 0.5) # used for SW network only, can hold multiple values (don't if not SW!)
k = 10 # mean num neighbors - must be even integer
nreps = 1
kill.on.this.day = 10000 # if it's going long. Make huge if needed
#------------------------------------------
# Vaccination parameters
doses.prop.per.day = c(0, 0.003, 0.006, 0.009, 0.012) # proportion of inds getting V1. V2 given after V1.days
vacc.strat = c("random","hubs")
vaccinate.only.S = c(FALSE) # if F then vaccinate anyone including R
vaccine.doses.per.person = c(1,2) # number vaccine doses tested
vacc.strain = 0 # strain to add like a strain inds recovered from
vacc.crit.pop = 0 # when pop of E+I+R exceeds this begin vaccinating
max.VE1 = seq(0.4, 0.75, by =0.05) # max efficacy of dose 1
max.VE2 = seq(0.75,0.9, by = 0.05) # max efficacy of dose 2
V1.days = 21 # num days for immune system to reach max effectiveness
V2.days = 21 # num days for immune system to reach max effectiveness
#------------------------------------------
# Evolution (achieved if mut.rate > 0)
mut.rate = 0 #seq(0, 0.001, by = 0.00025) #c(1e-3, 5e-4, 1e-4, 0)# prob. that a locus on strain is flipped when giving to host (should be small!)
n.loci = 12 # length of virus strain in binary (up to 31 alloweds)
n.loci.for.mismatch = 5 # more than this there is no partial immunity
# below is max num uniq seq w/o gradient of partial immunity
# max.num.uniq.seq = 2^(n.loci - n.loci.for.mismatch) - 1
#---------------------------------------------
tot.sims = nreps*length(rewire.p)*length(doses.prop.per.day)*
  length(vaccine.doses.per.person)*length(vacc.strat)*length(mut.rate)*
  length(max.VE1)*length(max.VE2)*length(num.days.R)*length(vaccinate.only.S) *
  length(R.naught)
sim.count = 0
start.time = Sys.time()
for (Ro in 1:length(R.naught)) {
  for (vs in 1:length(vacc.strat)) {
    for (nvd in 1:length(vaccine.doses.per.person)) {
      for (vacc.only.S in vaccinate.only.S) {
        for (rp in 1:length(rewire.p)) {
          for (ndoses in 1:length(doses.prop.per.day)) {
            for(mr in 1:length(mut.rate)){
              for (maxVE1 in 1:length(max.VE1)) {
                for (maxVE2 in 1:length(max.VE2)) {
                  for (dR in 1:length(num.days.R)) {
                    vacc.efficacy = make.vacc.efficacy(V1.days, V2.days, max.VE1[maxVE1], max.VE2[maxVE2])
                    for (the.rep in 1:nreps) {
                      sim.count = sim.count + 1
                      cat("Running sim",sim.count,"of",tot.sims,"\n")
                      Day = 1 # set Day variable
                      inds = make.inds(N)
                      inds = innoculate.network(inds, num.E.T1, num.I.T1, Day)
                      time.step.data = make.time.step.data(inds) # assumes Day = 1
                      g = make.network(net.type,N,k,rewire.p[rp])
                      Trans.prob = 1 - (1 - R.naught[Ro]/k)^(1/num.days.I) # transmission prob.
                      #-------------------------
                      if (animate.network == T) {
                        my.layout = get.layout(net.type)
                        if (N <= Net.vis.max) {
                          plot.network(g, inds, Day, net.type, my.layout)
                        }
                      }
                      #-------------------------
                      if (R.naught[Ro] > k) {
                        cat("Ro > k. Luke says: You ask the impossible!\n")
                        break
                      }
                      #--------------------------
                      #  Initial Vaccination
                      if (doses.prop.per.day[ndoses] > 0 & # vaccinating inds
                          sum(get.states(inds)[2:3]) > vacc.crit.pop & # we have more sickened than vacc.crit.pop 
                          vaccine.doses.per.person[nvd] != 0) { # vaccine doses are more than none
                        inds = vaccinate(inds, Day, doses.prop.per.day[ndoses], 
                                         V1.days, vaccine.doses.per.person[nvd], vacc.strat[vs], 
                                         vacc.strain, vacc.only.S, g)
                      }
                      #---------------------------------------------------------
                      # 2. Run simulation
                      while(1) {
                        Day = Day + 1
                        # 3. Check whether E->I and I->R
                        inds = update.inds(inds, Day, num.days.E, num.days.I, num.days.R[dR])
                        #----------------------------------
                        # 4. Check if there are still infected inds. If not, record stats and break
                        states = get.states(inds)
                        if (sum(states[2:3]) == 0) {
                          num.strains = get.tot.num.strains(inds)
                          current.strains = get.num.cur.strains(inds)
                          time.step.data = rbind(time.step.data,c(Day,states,num.strains,current.strains))
                          break
                        }
                        #--------------------------
                        # 5. Vaccinate
                        if (doses.prop.per.day[ndoses] > 0 & # vaccinating inds
                            sum(get.states(inds)[2:4]) > vacc.crit.pop & # we have more sickened than vacc.crit.pop 
                            vaccine.doses.per.person[nvd] != 0) { # vaccine doses are more than none
                          inds = vaccinate(inds, Day, doses.prop.per.day[ndoses],
                                           V1.days, vaccine.doses.per.person[nvd], vacc.strat[vs], 
                                           vacc.strain, vacc.only.S, g)
                        }
                        #--------------------------
                        # 6. Try to infect neighbors of infectious inds
                        inds = spread(inds, Trans.prob, Day, doses.prop.per.day[ndoses],
                                      V1.days, vacc.efficacy, mut.rate[mr], n.loci, n.loci.for.mismatch)
                        #--------------------------
                        # 7. Record stats
                        states = get.states(inds)
                        num.strains = get.tot.num.strains(inds)
                        current.strains = get.num.cur.strains(inds)
                        time.step.data = rbind(time.step.data,c(Day,states,num.strains,current.strains))
                        if (animate.network == T & N <= Net.vis.max) {
                          plot.network(g, inds, Day, net.type, my.layout)
                        }
                        if (Day == kill.on.this.day) break # used to break long runs
                      } # while(1)
                      states = get.states(inds)
                      num.sick = length(which(inds$states$dayR > 0))
                      max.num.I = max(time.step.data[,4])
                      max.I.day = which(time.step.data[,4] == max.num.I)[1] # just the first day
                      num.days = length(time.step.data[,1])
                      num.vacc.1 = length(which(inds$daysV$dayV1 > 0))
                      num.vacc.2 = length(which(inds$daysV$dayV2 > 0))
                      num.vacc.doses.per.day = round(doses.prop.per.day[ndoses]*N,0)
                      n.doses.per.person = vaccine.doses.per.person[nvd]
                      vac.strat = vacc.strat[vs]
                      num.strains = get.tot.num.strains(inds)
                      n.d.R = num.days.R[dR]
                      mutrate = mut.rate[mr]
                      MVE1 = max.VE1[maxVE1]
                      MVE2 = max.VE2[maxVE2]
                      R.knot = R.naught[Ro]
                      if (sim.count == 1) {
                        summary.data = data.frame(ver.num, N, R.knot, k, n.d.R, MVE1, MVE2, 
                                                  num.vacc.doses.per.day,n.doses.per.person, 
                                                  n.loci, n.loci.for.mismatch,
                                                  rewire.p[rp], num.vacc.1, num.vacc.2,
                                                  vac.strat, vacc.only.S, num.strains, mutrate,
                                                  num.sick, max.num.I, max.I.day, num.days)
                      } else {
                        summary.data = rbind(summary.data, 
                                             c(ver.num, N, R.knot, k, n.d.R, MVE1, MVE2, 
                                               num.vacc.doses.per.day,n.doses.per.person, 
                                               n.loci, n.loci.for.mismatch,
                                               rewire.p[rp], num.vacc.1, num.vacc.2,
                                               vac.strat, vacc.only.S, num.strains, mutrate,
                                               num.sick, max.num.I, max.I.day, num.days))
                      } # summary data
                      if (show.time.step.graph) {
                        time.step.graph(time.step.data)
                      }
                    } # the rep
                  } # num.days.R before returning to S
                } # max.VE2
              } # max.VE1
            } # mutation rate
          } # doses of vacc. given as proportion of pop per day
        } # rewire.p values
      } # vacc.only.S
    } # vacc.doses.per.person
  } # vacc.strat
} # Ro loop
#-----------------------------------------------------

time.step.data = as.data.frame(time.step.data)
names(time.step.data) = c("Day","S","E","I","R","Tot.Num.Strains","Num.DailyStrains")
if (save.time.step.data) {
  write.csv(time.step.data, time.step.data.name, row.names = FALSE)
}

if (save.summary.data) {
  write.csv(summary.data, save.summary.data.file.name, row.names = F)
}

#-----------------------------------------------------
# quick and dirty graphs


if (show.n.strains.graph) {
  par(mfrow = c(1,2))
  plot(time.step.data$Tot.Num.Strains, type = "l", las = 1, cex.lab = 1.5,
       xlab = "Time Step", ylab = "Total Number of Unique Strains")
  plot(time.step.data$Num.DailyStrains, type = "l", las = 1, cex.lab = 1.5,
       xlab = "Time Step", ylab = "Current Number Strains Circulating")
  par(mfrow = c(1,1))
}
################################################################

cat("The model ran for",format(Sys.time() - start.time),"\n")
summary.data
