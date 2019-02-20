##########
# Mathematical functions used in SGD fit and simulations
# same notations as those used in the paper
##########



g <- function(t,b,d)
{
	denominator <- ( 1. - (d/b)*exp(-(b-d)*t) )
	result <- (b-d) / denominator
	return(result)
}

m <- function(t, b, d, nu)
{
	test_zero <- b-d-nu
	# if we are in the neighborhood of zero, we use the Taylor series. Otherwise, the true function
	if(test_zero < 1e-8 && test_zero > -1e-8)
	{
		numerator <- b - d*exp(-(b-d)*t)
		denominator <- (1 + b*t)*(b-d)
		result <- numerator/denominator
	}
	else
	{
		numerator <- ( (1-b/(nu+d))*(exp((-b+d)*t) - b/d) )
		denominator <- (1-b/d)*(exp((-b+d+nu)*t) - b/(d+nu))
		result <- numerator / denominator
	}
	return(result)
}

rho_1mute0 <- function(t, b, d, nu)
{
	return( nu*m(t,b,d,nu)/(1-m(t,b,d,nu)) )
}
	
rho_1to1 <- function(t, b, d, nu)
{
	return( g(t, b, d)*(1 - m(t,b,d,nu)) )
}

rho_0to1 <- function(t, b, d, nu)
{
	return( 2*g(t, b, d)*(1 - m(t,b,d,nu)) )
}

rho_0tofixe <- function(t, b, d, nu)
{
	return( g(t, b, d)*m(t,b,d,nu) )
}

derivative_leaf <- function(t, y, parms)
{
	# Compute the derivative at a given state y and time t, with parameters parms = (b,d,nu,f)
	b <- parms[1]
	d <- parms[2]
	nu <- parms[3]
	f <- parms[4]

	d_uf0 <- -(nu + g(t, b, d))*y[1] + 2*g(t,b,d)*y[1]*y[2] + (1-f)*g(t, b, d)*m(t, b, d, nu)**2
	d_uf1 <- nu*y[1] - g(t, b, d)*y[2]*(1-y[2])
	d_Lf0 <- -(nu + g(t, b, d))*y[3] + 2*g(t,b,d)*y[2]*y[3] + 2*g(t,b,d)*y[1]*y[4] + f*g(t, b, d)*m(t, b, d, nu)**2
	d_Lf1 <- nu*y[3] - g(t, b, d)*y[4] + 2*g(t,b,d)*y[2]*y[4]

	return( list(c(d_uf0, d_uf1, d_Lf0, d_Lf1)) )
}

derivative_edge <- function(t, y, parms)
{
	# Compute the derivative at a given state y and time t, with parameters parms = (b,d,nu,f)
	b <- parms[1]
	d <- parms[2]
	nu <- parms[3]
	f <- parms[4]

	d_uf0 <- -(nu + g(t, b, d))*y[1] + 2*g(t,b,d)*y[1]*y[2] + (1-f)*g(t, b, d)*m(t, b, d, nu)**2
	d_uf1 <- nu*y[1] - g(t, b, d)*y[2]*(1-y[2])
	d_Lt0 <- -(nu + g(t, b, d))*y[3] + 2*g(t,b,d)*y[2]*y[3] + 2*g(t,b,d)*y[1]*y[4]
	d_Lt1 <- nu*y[3] - g(t, b, d)*y[4] + 2*g(t,b,d)*y[2]*y[4]
	
	return( list(c(d_uf0, d_uf1, d_Lt0, d_Lt1)) )
}



##########
# Inference functions
##########

LikelihoodSGDFromLambert <- function(tree, b, d, nu, f)
{
	# Compute the likelihood of any ultrametric tree in lambert representation under the SGD model with parameters given.
	if (length(tree) == 1)
	{
		t1 <- 0
		t2 <- tree[1]
		yt1 <- c(1-f, 0, f, 0)
		yt2 <- ode(yt1, c(t1, t2), derivative_leaf, parms = c(b,d,nu,f), atol=1e-100 )[2,2:5]

	}
	else{
		imax <- which.max( tree[2:length(tree)] ) +1
		t1 <- tree[imax]
		t2 <- tree[1]
		if( imax == 2 )
		{
			R1 <- c(tree[imax])
		}
		else{
			R1 <- c(tree[imax], tree[2:(imax-1)])
		}
		R2 <- tree[imax:length(tree)]

		yR1 <- LikelihoodSGDFromLambert(R1, b, d, nu, f)
		yR2 <- LikelihoodSGDFromLambert(R2, b, d, nu, f)

		gt1 <- g(t1, b, d)
		L0 <- (yR1[3]*yR2[4] + yR1[4]*yR2[3])*gt1
		L1 <- yR1[4]*yR2[4]*gt1

		yt1 <- c( yR2[1], yR2[2], L0, L1 )
		yt2 <- ode(yt1, c(t1, t2), derivative_edge, parms = c(b,d,nu,f), atol=1e-150 )[2,2:5]
	}
	return(yt2)
}

##########
# Simulating SGD phylogenetic trees
##########

Inverse_FB <- function(u, tau, b, d)
{
	return(log( u*(exp((b-d)*tau) - d/b) + d/b ) / (b-d))
}

Inverse_FB2 <- function(u, tau, b, d)
{
	return(log( sqrt(u)*(exp((b-d)*tau) - d/b) + d/b ) / (b-d))
}

Simulate_SGD_Phylogeny <- function(tau, remainder, lineage_type, b, d, nu)
{
	if (lineage_type == 0)
	{	
		# On tire le prochain temps auquel il se passe une naissance de type 1 :
		TB1 <- 0
		T <- tau
		while ((TB1 == 0) && (T > 0))
		{
			u <- runif(1, 0, 1)
			T <- Inverse_FB2(u, T, b, d)
			u2 <- runif(1, 0, 1)
			if ( u2 <= (1 - m(T,b,d,nu)) )
			{
				TB1 <- T
			}
		}
		
		# On tire le prochain temps auquel il se passe une fixation :
		TBfixe <- 0
		T <- tau
		while ((TBfixe == 0) && (T > TB1))
		{
			u <- runif(1, 0, 1)
			T <- Inverse_FB(u, T, b, d)
			u2 <- runif(1, 0, 1)
			if ( u2 <= m(T,b,d,nu) )
			{
				TBfixe <- T
			}
		}
		
		if ((TB1 <= 0) && (TBfixe <= 0))
		{
			resultat <- c(tau + remainder)
		}
		else{
			# On a donc eu la mutation avant, donc on repasse en type 0 au temps s1 :
			if (TBfixe > TB1)
			{
				resultat <- c(tau + remainder)
			}
			else{
			# On a eu d'abord une naissance de type 1, qui vient se brancher.
				resultat <- c( Simulate_SGD_Phylogeny(TB1, remainder + (tau-TB1), 0, b, d, nu) , Simulate_SGD_Phylogeny(TB1, 0, 1, b, d, nu) )
			}
		}

	}
	else if (lineage_type == 1){
		# On tire le prochain temps auquel il se passerait une naissance de type 1 :
		TB <- 0
		T <- tau
		while ((TB == 0) && (T > 0))
		{
			u <- runif(1, 0, 1)
			T <- Inverse_FB(u, T, b, d)
			u2 <- runif(1, 0, 1)
			if ( u2 <= (1 - m(T,b,d,nu)) )
			{
				TB <- T
			}
		}

		# On tire le prochain temps auquel il se passerait une mutation avec passage en type 0 :
		ds <- (tau-TB)/100.0
		s <- tau
		u <- runif(1, 0, 1)
		Inf_Int <- log(u)
		Int <- 0
		while ((s > 0) && (s > TB) && (Int > Inf_Int))
		{
			Int <- Int - rho_1mute0(s,b,d,nu)*ds
			s <- s - ds
		}
		TM <- s
			
		# On regarde ce qu'il se passe au premier evenement qui survient :
		if ((TB <= 0) && (TM <= 0))
		{
			resultat <- c(tau + remainder)
		}
		else{
			# On a donc eu la mutation avant, on repasse en type 0 :
			if (TM > TB)
			{
				resultat <- Simulate_SGD_Phylogeny(TM, remainder + (tau-TM), 0, b, d, nu)
			}
			else{
			# On a eu d'abord une naissance de type 1, qui vient se brancher.
				resultat <- c( Simulate_SGD_Phylogeny(TB, remainder + (tau-TB), 1, b, d, nu) , Simulate_SGD_Phylogeny(TB, 0, 1, b, d, nu) )
			}		
		}	
	}
	else{
		u <- runif(1, 0, 1)
		# S'il y a survie du type clonal sachant survie totale
		if (u <= m(tau, b, d, nu))
		{
			resultat <- Simulate_SGD_Phylogeny(tau, remainder, 0, b, d, nu)
		}
		else{
			resultat <- Simulate_SGD_Phylogeny(tau, remainder, 1, b, d, nu)
		}
	}

	return(resultat)
}



##########
# Translating tree representations
##########

Phylo2Lambert <- function(phylo)
{
	# find the root node
	for (node in phylo$edge[,1])
	{
		if (node %in% phylo$edge[,2])
		{
			next
		}else{
			root_node <- node
			break
		}
	}

	# translate into the "lambert" format
	counted_nodes <- c()
	lambert <- c()
	for (i in 1:length(phylo$tip.label))
	{
		node <- i
		branch_length <- 0
		while (!(node %in% counted_nodes) && (node != root_node))
		{
			counted_nodes <- c(counted_nodes, node)
			vector_bool <- phylo$edge[,2] == node
			branch_length <- branch_length + phylo$edge.length[vector_bool]
			node <- phylo$edge[,1][vector_bool]
		}
		lambert <- c(lambert, branch_length)
	}
	return(lambert)	
}

Lambert2Newick <- function(lamb)
{
	# recursively transforms the lambert representation into a newick tree
	if (length(lamb) == 1)
	{
		newick <- paste(":", lamb[1], sep="")
	}else{
		imax <- which.max( lamb[2:length(lamb)] ) + 1
		remains <- lamb[1] - lamb[imax]
		if( imax == 2 )
		{
			R1 <- c(lamb[imax])
		}else{
			R1 <- c(lamb[imax], lamb[2:(imax-1)])
		}
		R2 <- lamb[imax:length(lamb)]
		
		if( remains == 0 )
		{
			newick <- paste("(", Lambert2Newick(R1), ",", Lambert2Newick(R2), ");", sep="")
		}else{
			newick <- paste("(", Lambert2Newick(R1), ",", Lambert2Newick(R2), "):", remains, sep="")
		}
	}
	return(newick)
}

