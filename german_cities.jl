include("TSPSubgradient.jl")

using TSPSubgradient, LightGraphs, GraphLayout, Plots

# distance matrix between cities (symmetric, undirected)
distmx = Float64[
    0   91  80  259 70  121;
    91  0   77  175 27  84;
    80  77  0   232 47  29;
    259 175 232 0   189 236;
    70  27  47  189 0   55;
    121 84  29  236 55  0;
]


N = size(distmx)[1]

root = 6

g = Graph(N)

for i in 1:N
    for j in 1:N
        if i!=j
            add_edge!(g, i, j)
        end
    end
end

ot = one_tree(g, distmx, 6)

iters = 300

costs, ots, lambdas = tsp_subgradient(g, distmx, iters, 6, tau=0.5)


# plot one tree
for k in 10:10:iters-1
    am = full(adjacency_matrix(ots[k]))
    loc_x, loc_y = layout_spring_adj(am)
    draw_layout_adj(am, loc_x, loc_y, labels=Array(1:6), filename=string("ot",k,".svg"))
end


using Plots
plotly()

plot(costs, linewidth=2,title="Cost")

plot(lambdas,linewidth=2,title="Lagrange multipliers")




