module TSPSubgradient

using LightGraphs

"Calculates a so called one tree. This tree is a MST with one additional
edge closing the smallest circle from the root. distmx is a matrix of
weights as it is used by the Prim algorithm for MSTs from LightGraphs.
Returns an undirected graph of the tree."
function one_tree(g, distmx, root)
    mst = prim_mst(g, distmx)
    mstg = Graph(length(g.vertices))
    for e in mst
        add_edge!(mstg, e)
    end
    prefs = sortperm(distmx[root,:])
    ns = Set(neighbors(mstg, root))
    push!(ns, root)
    one_circle_node = filter(x->!in(x, ns), prefs)[1]
    #println(!any([i==one_circle_node for i in mstg.fadjlist[root]]))
    add_edge!(mstg, root, one_circle_node)
    return mstg
end

"This subgradient approximation for the traveling salesman problem
uses a relaxation of the degree conditions on each node. For the TSP
a node must have degree(node)==2 which yields lambda*(degree(node)-2)
in the Lagrangian. Once these constraints are relaxed, the objective
becomes the calculation of a one-tree, which can be efficiently done
by calculating an MST with Prim on the transformed edge weights:
twe[i,j] = we[i,j] + lambdas[i] + lambdas[j]. This yields a gradient
update for the lambdas. Tau is the initial learning rate and alpha is
its decay. The function returns the costs for each step, the one-tree
calculated in each step and the Lagrangian multipliers of each step."
function tsp_subgradient(g, distmx, m, root; tau=0.1, alpha=0.99999)
    # 1)
    t = zeros(m)
    N = length(g.vertices)

    # 2)
    t[1] = tau
    lambdas = zeros(m,N)
    d = zeros(m, N)
    ots = []
    costs = []

    # 3)
    for k in 2:m
        # calculate current distance
        dist = [if i!=j distmx[j,i] + lambdas[k-1,i] + lambdas[k-1,j]
                else 0.0 end
                for j in 1:N, i in 1:N]
        # 3.1)
        ot = one_tree(g, dist, root)
        push!(ots, ot)
        # 3.2)
        d[k,:] = [degree(ot, i) - 2 for i in 1:N]
        # 3.3)
        for i in 1:N
            # bundle (or momentum) averaging
            lambdas[k, i] = lambdas[k-1, i] + t[k-1]*(0.7*d[k,i] + 0.3*d[k-1,i])
            lambdas[k, root] = 0
        end
        # 3.4)
        t[k] = alpha*t[k-1]

        cost = sum([distmx[e[1], e[2]] for e in edges(ot)])
        push!(costs, cost)
    end
    return costs, ots, lambdas
end

export one_tree, tsp_subgradient

end
