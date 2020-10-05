module reproSim
using Pkg
Pkg.add("Distributions")
Pkg.add("StatsBase")
Pkg.add("CSV")
Pkg.add("DataFrames")
using Distributions
using StatsBase
using CSV
using DataFrames
function getPop(n)
    """return array of length n. with index, x, y chomosome values"""
return hcat(1:n,rand(Normal(0,0.5),n),rand(Normal(0,0.5),n))
end

function fselect(pop,n)
    """"returns population with additional column for fecundity"""
        # get samples
    fitList = sample(pop[:,1], Weights(pop[:, 4]), n, replace=true)
    #add fecundity column to pop.
 pop = hcat(pop,[sum(i .== fitList) for i in pop[:,1]])
 #return pop. without non-reproducers
return pop[pop[:,5] .!= 0.6,:]
end

function sselect(pop,n)
    """" return half of pop. weighted by fitness"""


return pop[Int.(sample(Int.(1:n), Weights(pop[:, 5]), Int(n/2), replace=false)),:]
end

function ggselect(pop,n)
    """" return half of pop. weighted by fitness"""

return pop[Int.(sample(Int.(1:n/2), Weights(pop[:, 1]), Int(n/2), replace=true)),:]
end

function clone(pop,n,mut)
    x = []
    y = []

    for mother in 1:size(pop)[1] # for each mother


        for fecundity in 1:pop[mother,7] #according to fecundity
            push!(x,pop[mother,2])
            push!(y,pop[mother,3])



        end
    end


    return hcat(1:n,x .+ rand(Normal(0,mut),n),y .+ rand(Normal(0,mut),n))



end

function sex(pop,fathers,n,mut)
    x = []
    y = []



    for mother in 1:size(pop)[1] # for each mother
        father = fathers[mother,:] # pick a father

        for fecundity in 1:pop[mother,7] #according to fecundity
            push!(x,pop[mother,rand(2:3)])
            push!(y,father[rand(2:3)])
        end
    end

    return hcat(1:n,x .+ rand(Normal(0,mut),n), y .+ rand(Normal(0,mut),n))
end

function fitPop(pop,f,s,gg,theta)
    z = (pop[:,2] .+ pop[:,3]) ./2

    for c in [1-f,1-s,1-gg]
        pop= hcat(pop, c .+ ((1 .-c) .*exp.((-1/2) .*(z .- theta) .^2)))
    end

    return pop
end

function getW(pop,theta)
    z = pop[:,2] .+ pop[:,3] ./2


    return [mean(exp.((-1/2) .*(z .- theta) .^2)) var(z)]
end

function sim(vartheta,lambda,rMode,mut,n,generations,f,s,gg)
    theta = 0
    noise = vartheta*(1-(lambda^2))
    pop = getPop(n)
    pop = fitPop(pop,f,s,gg, theta)
    dat = Array{Float64}(undef, 0, 2)
    r = Int(n/2)



    for t in 1:generations


        pop = sselect(pop,n)


        if rMode == "clone"
            pop = fselect(pop,n)
            pop = clone(pop,n,mut)
        else
            fathers = ggselect(pop,n)
            pop = fselect(pop,n)
            pop = sex(pop,fathers,n,mut)
        end

        theta = theta*lambda +rand(Normal(0,noise))
        pop = fitPop(pop,f,s,gg, theta)
        dat = vcat(dat,getW(pop,theta))
    end

    wg = geomean(dat[:,1])
    wa = mean(dat[:,1])
    meanv = mean(dat[:,2])
    vmean = var(dat[:,1])

    return [wg wa meanv vmean vartheta lambda rMode mut n generations f s gg]
end

end # module
