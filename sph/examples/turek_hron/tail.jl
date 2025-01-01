module tail

using SmoothedParticles

#=
Declare constants
=#

const pull_force = 1.5 #pulling force [N]
const pull_time = 3.0  #for how long we pull

const c_l = 20.0   #longitudinal sound speed
const c_s = 200.0  #shear sound speed
const c_0 = sqrt(c_l^2 + 4/3*c_s^2)  #total sound speed
const rho0 = 1.0   #density
const nu = 1.0e-4    #artificial viscosity (surpresses noise but is not neccessary)

#= 
Algebra
=#

@inbounds function outer(x::RealVector, y::RealVector)::RealMatrix
    return RealMatrix(
        x[1]*y[1], x[2]*y[1], 0., 
        x[1]*y[2], x[2]*y[2], 0.,
        0., 0., 0.
    )
end

@inbounds function det(A::RealMatrix)::Float64
    return A[1,1]*A[2,2] - A[1,2]*A[2,1]
end

@inbounds function inv(A::RealMatrix)::RealMatrix
    idet = 1.0/det(A)
    return RealMatrix(
        +idet*A[2,2], -idet*A[2,1], 0., 
        -idet*A[1,2], +idet*A[1,1], 0.,
        0., 0., 0.
    )
end

@inbounds function trans(A::RealMatrix)::RealMatrix
    return RealMatrix(
        A[1,1], A[1,2], 0., 
        A[2,1], A[2,2], 0.,
        0.,  0., 0.
    )
end

@inbounds function dev(G::RealMatrix)::RealMatrix
    λ = 1/3*(G[1,1] + G[2,2] + 1.0)
    return RealMatrix(
        G[1,1] - λ,G[2,1], 0.0,
        G[1,2], G[2,2] - λ, 0.0,
        0.0, 0.0, 1.0 - λ
    )
end

#=
Physics
=#

function find_A!(p::AbstractParticle, q::AbstractParticle, r::Float64)
    ker = wendland2(h,r)
    x_pq = p.x - q.x
    X_pq = p.X - q.X
    p.A += -ker*outer(X_pq, x_pq)
    p.H += -ker*outer(x_pq, x_pq)
end

function find_B!(p::AbstractParticle)
    Hi = inv(p.H)
    p.A = p.A*Hi
    At = trans(p.A)
    G = At*p.A
    P = c_l^2*(det(p.A)-1.0)
    p.B = m*(P*inv(At) + c_s^2*p.A*dev(G))*Hi
end

function find_a!(p::AbstractParticle, q::AbstractParticle, r::Float64)
    ker = m*wendland2(h,r)
    rDker = m*rDwendland2(h,r)
    x_pq = p.x - q.x
    X_pq = p.X - q.X
    #acceleration
    p.a += -ker*(trans(p.A)*(p.B*x_pq))
    p.a += -ker*(trans(q.A)*(q.B*x_pq))
    #"eta" correction (remove this -> energy will not be conserved!)
    k_pq = +trans(p.B)*(X_pq - p.A*x_pq)
    k_qp = -trans(q.B)*(X_pq - q.A*x_pq)
    p.a += rDker*dot(x_pq, k_pq)*x_pq + ker*k_pq
    p.a -= rDker*dot(x_pq, k_qp)*x_pq + ker*k_qp
    #artificial_viscosity
    p.a += 2*m*vol*rDker*nu*(p.v - q.v)
end

function pull!(p::AbstractParticle, L)
    if p.X[1] > L-h
        p.a += m*RealVector(0., (vol*pull_force)/(h*W), 0.)
    end
end

end