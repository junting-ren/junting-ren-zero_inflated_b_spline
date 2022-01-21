mean(mapply(function(X,Y) cor(X, Y), X = MCMCfit1$b_delta, Y = MCMCfit1$b_eta))
