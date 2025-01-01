module constants

#geometry parameters
export dr, h
const dr = 3.9e-3 		     #average particle distance (decrease to make finer simulation)
const h = 2.5*dr		     #size of kernel support

#physical parameters
export rho0, vol, m
const rho0 = 1.0
const vol = dr*dr            #particle volume
const m = rho0*vol           #particle mass

#tail parameters
export L, W, pull_time
const L = 0.25               #length of the tail
const W = 0.01               #width of the tail
const pull_time = 0.5

end