#= 
helpful commands:
upload file to cluster:
download file from cluster:
=#

module turek_hron

using SmoothedParticles
include("tail.jl")

const folder_name = "results/turek_hron"

#=
Declare constants
=#

#geometry parameters
const dr = 3.9e-3 		     #average particle distance (decrease to make finer simulation)
const h = 2.5*dr		     #size of kernel support
const chan_l = 0.8 #2.2      #length of the channel
const chan_w = 0.41          #width of the channel
const cyl1 = dr*round(0.2/dr)  #x coordinate of the cylinder
const cyl2 = dr*round(0.2/dr)  #y coordinate of the cylinder
const cyl_r = 0.05           #radius of the cylinder
const wall_w = 2.5*dr        #width of the wall
const inflow_l = 3.0*dr      #width of inflow layer

#tail parameters
const L = 0.25
const W = 0.01

#physical parameters
const rho0 = 1.0
const vol = dr*dr            #particle volume
const m = rho0*vol           #particle mass


#temporal parameters
const dt = 0.1
const t_end = 1.0

#particle types
const FLUID = 0.
const WALL = 1.
const INFLOW = 1.
const OBSTACLE = 2.
const TAIL = 3.

#=
Define variables to be stored in a Particle 
=#

mutable struct Particle <: AbstractParticle
    x::RealVector #position    
    v::RealVector #velocity
    a::RealVector #acceleration
    rho::Float64 #density
    Drho::Float64 #rate of density
    P::Float64 #pressure
    #tail parameters
    X::RealVector #Lag. position    
    A::RealMatrix  #distortion
    H::RealMatrix  #correction matrix
    B::RealMatrix  #derivative of energy wrt A
    e::Float64     #fronorm squared of eta
    type::Float64 #particle type
    Particle(x, type) = begin
        return new(x, VEC0, VEC0, rho0, 0., 0., x, MAT0, MAT0, MAT0, 0., type)
    end
end

#=
Geometry
=#

function make_system()
    domain = Rectangle(-inflow_l, -10*wall_w, chan_l, chan_w + 10*wall_w)
    sys = ParticleSystem(Particle, domain, h)
    grid = Grid(dr, :square)

    #define geometry
    obstacle = Circle(cyl1, cyl2, cyl_r)
    pipe = Rectangle(-inflow_l, 0., chan_l, chan_w)
    wall = BoundaryLayer(pipe, grid, wall_w)
    wall = Specification(wall, x -> (-inflow_l <= x[1] <= chan_l))      # Why is this needed? -> Removes the left and right wall
    tail = Rectangle(cyl1 + cyl_r, cyl2, cyl1 + cyl_r + L, cyl2 + W)
    #generate particles
    generate_particles!(sys, grid, wall, x -> Particle(x, WALL))
    generate_particles!(sys, grid, obstacle, x -> Particle(x, OBSTACLE))
    generate_particles!(sys, grid, tail, x -> Particle(x, TAIL))
    return sys
end

#=
Physics
=#

function update_v!(p::Particle)
    p.v += 0.5*p.a*dt
    # dirichlet BC - fixes the base of the tail
    if p.X[1] < h + cyl1 + cyl_r
        p.v = VEC0
    end    
end

function update_x!(p::Particle)
    p.x += p.v*dt
    if p.type == TAIL
        #reset vars
        p.H = MAT0
        p.A = MAT0
        p.f = VEC0
        p.e = 0.
    end
end

function find_a!(p::Particle)
    if p.type == TAIL
        tail.pull!(p, L)
    end
end


#=
# Modified Verlet scheme
    1. Calculate v(t+0.5*Δt) = v(t) + 0.5*a(t)*Δt
    2. Calculate x(t+0.5*Δt) = x(t) + v(t+0.5*Δt)*Δt
    3. Derive a(t+Δt) from the interaction potential using x(t+Δt)
    4. v(t+Δt) = 0.5*v(t+0.5*Δt) + 0.5*a(t+Δt)*Δt
=#

function main()
    sys = make_system()
    out = new_pvd_file(folder_name)
    #Verlet scheme
    for k = 0 : Int64(round(t_end/dt))
        t = k * dt
        # update v
        # update x
        # calculate a
        # update v
        save_frame!(out, sys, :type)
    end    
    save_pvd_file(out)
    println("i did something")
end

end