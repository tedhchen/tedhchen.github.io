This document outlines how to use the `multilayer.ergm` package to fit
exponential random graph models (ERGM) in `R`. The `multilayer.ergm`
package extends functionalities of the `ergm` package to multilayer
network applications by providing a way to count local network
configurations, or network motifs, that span more than one network
layer. This document assumes you have a working knowledge of ERGMs and
multilayer networks as described in [Chen
(2019)](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3189835).

Installation
------------

First, install the package from GitHub. To do so, you will need the
`devtools` package. If you are on a Windows machine, you will also need
the standalone [Rtools](https://cran.r-project.org/bin/windows/Rtools/).
If you do not already have the `network`, `ergm`, and `statnet.common`
packages, installing `multilayer.ergm` will install them as well.

Install the package by entering the following into your R console. The
argument `build_opts` can be removed entirely if you do not require the
package tutorial.

    devtools::install_github("tedhchen/multilayer.ergm", 
                             build_opts = c("--no-resave-data", "--no-manual"))

Of course, we want to load the package, which will also load the
`network`, `ergm`, and `statnet.common` packages.

    library(multilayer.ergm)

ChemG Policy Network
--------------------

In this exercise, we will be using data from the policy network
surrounding legislation over control of the chemical industry in Germany
in 1980 (ChemG). This data was used by [Leifeld and Schneider
(2012)](http://dx.doi.org/10.1111/j.1540-5907.2011.00580.x) in their
analysis of policy communication. To see the data associated with this
policy network, we can refer to its help file.

    help(chemg_data)

From the help file, we can see that there are five network data files.
Each of these networks have the same number of nodes, representing the
30 actors in this policy system. We might study these networks
independently as monoplex networks, but to understand the entire policy
process, which is a system that includes interactions within and across
all of these networks, we need to move to the multilayer approach.

Creating a Multilayer Network
-----------------------------

The `multilayer.ergm` package includes utility functions that help
create multilayer networks and check `network` objects to see whether
they satisfy the criteria for multilayer networks.

A multilayer network has multiple layers, which is defined by the user.
Technically, a `network` object is compatible with the `multilayer.ergm`
package as long as it has a vertex attribute called `layer.mem`, which
stands for *layer membership*. This vertex attribute should be a numeric
vector that identifies the layer each node belongs to. The following
line of code will assign the vertex attribute to the `network` object.

    network%v%'layer.mem' <- layer.ids

### Policy Multiplex Network

For this exercise, let us construct a multiplex network with one layer
being a network of scientific communication among the actors (`sci`) and
the other layer being perception of influence among the actors
(`infrep`).

In the matrix representation of a multiplex network, the network layers
are on the main diagonal of the block matrix, and the off-diagonal
blocks are identity matrices that serve to tie actors across layers. The
`to.multiplex` function is a utility function that assembles networks as
layers in a multilayer network, assigns layer ids to the nodes, and sets
the off-diagonal matrices as identity matrices.

    mnet <- to.multiplex(sci, infrep, output = "network")
    mnet

    ##  Network attributes:
    ##   vertices = 60 
    ##   directed = TRUE 
    ##   hyper = FALSE 
    ##   loops = FALSE 
    ##   multiple = FALSE 
    ##   bipartite = FALSE 
    ##   total edges= 406 
    ##     missing edges= 0 
    ##     non-missing edges= 406 
    ## 
    ##  Vertex attribute names: 
    ##     layer.mem vertex.names 
    ## 
    ## No edge attributes

We can check whether the `network` object is compatible with multilayer
terms manually, or with the `check.multilayer` function.

    mnet%v%'layer.mem'

    ##  [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2
    ## [36] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2

    check.multilayer(mnet)

    ## Layer membership for all nodes are properly specified.

ERGMs for Multilayer Networks
-----------------------------

Now that we have constructed our multilayer network `mnet`, we can fit
an ERGM. Let us work our way up from a monoplex approach to the
scientific communication network.

### Scientific Communication Monoplex Model

First, we create the network.

    scinet <- network(sci)
    scinet

    ##  Network attributes:
    ##   vertices = 30 
    ##   directed = TRUE 
    ##   hyper = FALSE 
    ##   loops = FALSE 
    ##   multiple = FALSE 
    ##   bipartite = FALSE 
    ##   total edges= 63 
    ##     missing edges= 0 
    ##     non-missing edges= 63 
    ## 
    ##  Vertex attribute names: 
    ##     vertex.names 
    ## 
    ## No edge attributes

Then fit an ERGM using the `ergm` function from the `ergm` package. Let
us include the baseline `edge` term, a term for how similar a pair of
actors' preferences are (`edgecov(prefsim)`), and a term for the
tendency for actors to reciprocate scientific communication (`mutual`).

    mod.monoplex <- ergm(scinet ~ edges + edgecov(prefsim) + mutual, 
                         control = control.ergm(seed = 206424))

(In case anyone is interested `206424` is based on weather statistics in
State College, PA when I wrote this vignette: -2°C, 0% precipitation,
64% humidity, and a wind speed of 24km/h.)

Let us look at the output:

    summary(mod.monoplex)

    ## 
    ## ==========================
    ## Summary of model fit
    ## ==========================
    ## 
    ## Formula:   scinet ~ edges + edgecov(prefsim) + mutual
    ## 
    ## Iterations:  2 out of 20 
    ## 
    ## Monte Carlo MLE Results:
    ##                 Estimate Std. Error MCMC % z value Pr(>|z|)    
    ## edges           -2.88849    0.29500      0  -9.792   <1e-04 ***
    ## edgecov.prefsim -0.01555    0.10156      0  -0.153    0.878    
    ## mutual           2.30823    0.44119      0   5.232   <1e-04 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##      Null Deviance: 1206  on 870  degrees of freedom
    ##  Residual Deviance:  429  on 867  degrees of freedom
    ##  
    ## AIC: 435    BIC: 449.3    (Smaller is better.)

Next, we move to the multiplex with the scientific communication network
on one layer, and preceptions of influence on the other layer.

### Policy Multiplex, Within-layer Model

We already created this network earlier, but before we can fit the ERGM,
we need to construct a network of sampling constraints that gives the
multilayer network its multiplex structure. The `to.multiplex` function
can do this as well. Make sure to specify `offzeros = TRUE`.

    free <- to.multiplex(matrix(1, ncol = 30, nrow = 30), matrix(1, ncol = 30, nrow = 30), 
                         output = "network", offzeros = TRUE)

Now we can fit ERGMs of the multiplex network. Let us look at what
happens when we fit an ERGM to a multiplex network but do not include
any terms that capture dependence across layers.

This within-layer model has the same terms as the previous monoplex
model, for both network layers. The `multilayer.ergm` package includes a
set of within-layer terms for multilayer networks. They are usually
named after their `ergm` counterparts, but with `_layer` appended to the
end. They also usually require an argument `layer`, which indicates
which layer the term should be computed on. For a list of network terms
that are available for multilayer networks, see `?multilayer_terms`.

    mod.within <- ergm(mnet
                       ~ edges_layer(layer = 1) + edges_layer(layer = 2)
                       + edgecov_layer(prefsim, layer = 1) + edgecov_layer(prefsim, layer = 2)
                       + mutual("layer.mem", diff = T),
                       control = control.ergm(seed = 206424),
                       constraints = ~fixallbut(free))

This is the output:

    summary(mod.within)

    ## 
    ## ==========================
    ## Summary of model fit
    ## ==========================
    ## 
    ## Formula:   mnet ~ edges_layer(layer = 1) + edges_layer(layer = 2) + edgecov_layer(prefsim, 
    ##     layer = 1) + edgecov_layer(prefsim, layer = 2) + mutual("layer.mem", 
    ##     diff = T)
    ## 
    ## Iterations:  2 out of 20 
    ## 
    ## Monte Carlo MLE Results:
    ##                         Estimate Std. Error MCMC % z value Pr(>|z|)    
    ## edges_layer.1           -2.88179    0.28836      0  -9.994  < 1e-04 ***
    ## edges_layer.2           -1.27679    0.18001      0  -7.093  < 1e-04 ***
    ## edgecov.layer.1.prefsim -0.01759    0.09886      0  -0.178  0.85875    
    ## edgecov.layer.2.prefsim  0.17504    0.06211      0   2.818  0.00483 ** 
    ## mutual.same.layer.mem.1  2.30022    0.44410      0   5.179  < 1e-04 ***
    ## mutual.same.layer.mem.2  0.33416    0.21998      0   1.519  0.12875    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##      Null Deviance: 2412  on 1740  degrees of freedom
    ##  Residual Deviance: 1516  on 1734  degrees of freedom
    ##  
    ## AIC: 1528    BIC: 1560    (Smaller is better.)

Remember that layer 1 is `sci` and layer 2 is `infrep`. From the output
we see that the coefficients for the scientific layer is effectively the
same as those from the monoplex model. (In fact, the only difference is
due to randomness in the estimation procedure.) This is to be expected,
because while two layers are included, we did not model any dependence
between them, meaning that they do not influence each other at all.

### Policy Multiplex, Cross-layer Model

Now, let us fit a model that includes cross-layer dependence terms.
Specifically, let us include, using `duplexdyad`, a term for cross-layer
reinforcement and a term for cross-layer reciprocity. From the code
below, we see that the `duplexdyad` term takes an argument
`layers = list(1, 2)`. It means that the two layers included in the
dependence term are layers 1 and 2. This format is common to cross-layer
dependence terms in the `multilayer.ergm` package.

    mod.cross <- ergm(mnet 
                      ~ edges_layer(layer = 1) + edges_layer(layer = 2)
                      + edgecov_layer(prefsim, layer = 1) + edgecov_layer(prefsim, layer = 2)
                      + mutual("layer.mem", diff = T)
                      + duplexdyad(c("e", "f"), layers = list(1, 2)),
                      control = control.ergm(seed = 206424),
                      constraints = ~fixallbut(free))

And the output:

    summary(mod.cross)

    ## 
    ## ==========================
    ## Summary of model fit
    ## ==========================
    ## 
    ## Formula:   mnet ~ edges_layer(layer = 1) + edges_layer(layer = 2) + edgecov_layer(prefsim, 
    ##     layer = 1) + edgecov_layer(prefsim, layer = 2) + mutual("layer.mem", 
    ##     diff = T) + duplexdyad(c("e", "f"), layers = list(1, 2))
    ## 
    ## Iterations:  2 out of 20 
    ## 
    ## Monte Carlo MLE Results:
    ##                         Estimate Std. Error MCMC % z value Pr(>|z|)    
    ## edges_layer.1            -3.3173     0.3186      0 -10.410   <1e-04 ***
    ## edges_layer.2            -1.3762     0.1837      0  -7.493   <1e-04 ***
    ## edgecov.layer.1.prefsim  -0.0646     0.1023      0  -0.631   0.5278    
    ## edgecov.layer.2.prefsim   0.1798     0.0633      0   2.841   0.0045 ** 
    ## mutual.same.layer.mem.1   2.3703     0.4709      0   5.033   <1e-04 ***
    ## mutual.same.layer.mem.2   0.3224     0.2210      0   1.459   0.1446    
    ## duplexdyad.e              1.2958     0.2860      0   4.530   <1e-04 ***
    ## duplexdyad.f             -0.1568     0.3133      0  -0.501   0.6167    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ##      Null Deviance: 2412  on 1740  degrees of freedom
    ##  Residual Deviance: 1493  on 1732  degrees of freedom
    ##  
    ## AIC: 1509    BIC: 1553    (Smaller is better.)

The `duplexdyad.e` term in the output is the cross-layer reciprocity
effect and the `duplexdyad.f` term is the cross-layer reinforcement
effect. (For more detail on the different configurations in the duplex
dyad census, see Fig. 4 of [Chen
(2019)](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3189835).)
As these results show, an actor's tendency to communicate scientific
information to another actor is reinforced by whether they think that
other actor is influential and vice versa.
