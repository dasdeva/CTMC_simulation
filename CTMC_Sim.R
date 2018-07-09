Rcpp::sourceCpp('FastSim.cpp')

Fast_Sim = function() {
    SimRes =  UniformSimulation()
    States = SimRes$States
    Jump_times = SimRes$JumpTimes
    N = length(Jump_times) 
    change_indi = which(diff(States) != 0)
    return (list(N = N, ini = States[1],
                 states = c(States[change_indi +1]),
                 times = Jump_times[change_indi]))
}
