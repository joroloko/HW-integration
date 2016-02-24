
module HW_int


	# question 1 b)
	# here are the packages I used

 using ApproxFun
	using FastGaussQuadrature
	using Roots
	using Sobol
	using PyPlot
	using Distributions
 ##using Gadfly
	# here are some functions I defined for useage
	# in several sub questions

	# demand function

function demandfunction(p)
			 return  2p^(-0.5)
	end

x= linspace(0.5,4.5,100)
y= map(demandfunction,x)
plot(x,y,linewidth=2)
plot(ones(2),[1,2],color="red")
plot([1,4],[1,1],color="red")
title("The demand function")
xlabel("Price")
ylabel("Demand")

	# gauss-legendre adjustment factors for map change

	# eqm condition for question 2
	# this is the equilibrium condition: total demand = supply,


	# weighted sum for integration from the slides.



	function question_1b(n)
 S=Chebyshev([-1,1])
 nodes = gausschebyshev(n)  # generate n Chebyshev Nodes
 function myfunct(p) ##### The demand function when its argument is mapped [1,4]-> [-1,1]
		return demandfunction((3p/2) +2.5)
	end
 f=Fun(ApproxFun.transform(S,map(myfunct,nodes[1])),S) ### We approximate the function using at data its values at the Chebbyshev nodes
 ApproxFun.plot(f) ###### We now plot the original function versus the interpolated function.
 x= linspace(-1,1,100)
 y= map(myfunct,x)
 plot(x,y,linewidth=2,color="red") #### THe true function is in red
 g= cumsum(f) ### The integral function
 myintegral= 1.5*(g(1) - g(-1)) #### We multiply by 1.5 because we shrank the domain when mapping [1,4]-> [-1,1].
 return myintegral-3 #### Finally, we take away 3 because we are interested in the change in consumer surplus.
	end


	function question_1c(n)

#### We start by generating random points in [1,4] x[1,2], which is [domain]x[range] of our function
 myvectorx=  1+3*rand(n)
 myvectory = 1+ rand(n)
 myvaluesy= map(demandfunction, myvectorx)
### Now we simply take the expectancy multiplied by the area of [domain]x[range]
myexpectation= mean(myvaluesy .> myvectory)*3

x= linspace(0.5,4.5,100)
y= map(demandfunction,x)
###### Plot the function as well as the sampled area
plot(x,y,linewidth=2)
plot(ones(2),[1,2],color="red")
plot([1,4],[1,1],color="red")
scatter(myvectorx,myvectory,s=1)
	 return myexpectation
	end


	function question_1d(n)
		mypoints=[]
		t = SobolSeq(2) 
		for(i =1:n)
   mypoints = [mypoints, next(t)] #### We get the list of Sobol points in R2
		end
		mypoints=reshape(mypoints, 2,n) ######### We make the point matrix a 2xn matrix
		##### We  now apply the same scaling procedure as before
		myvectorx = 1+3*mypoints[1,:]
		myvectory = 1+ mypoints[2,:]
		myvaluesy= map(demandfunction, myvectorx)
		myexpectation= mean(myvaluesy .> myvectory)*3
		x= linspace(0.5,4.5,100)
		y= map(demandfunction,x)
		###### Plot the function as well as the sampled area
		plot(x,y,linewidth=2)
		plot(ones(2),[1,2],color="red")
		plot([1,4],[1,1],color="red")
		scatter(myvectorx,myvectory,s=1)
		return myexpectation

	end



function question_2a(np)
	np = 10
	rules = gausshermite(np)
	### We build the nodes matrx as all combinations for the 100 nodes
	col1= repeat(rules[1],inner=[1],outer=[np])
	col2 = repeat(rules[1],inner=[np],outer=[1])
	nodes= [col1 col2]
##### Since the nodes are correlated, we need to compute the cholesky inverse of the variance matrix and premultiply the nodes by it
sigma= [[.002 .001],[ .001 .001]]
sigma= cholfact(sigma)[:U]
nodes=nodes*sigma
	#### The corresponding weights will be the product of the original weights
	weights1 = repeat(rules[2],inner=[1],outer=[np])
	weights2 = repeat(rules[2],inner=[np],outer=[1])
	weights= weights1.*weights2
	weights=weights/sum(weights)  #### Finally, we normalize the weights to sum to 1.
	function getroot(a) #### The function to apply to each row of nodes
		function minifunct(x) #### The function ap^(-1) +b p^(-1/2) -2 =0
			 return exp(a[1])*x +exp(a[2])*sqrt(x)-2  ###### Since the shock distributions are given in logs, we need to take the exponential to get the price equations
			end
		try return fzero(minifunct,0,5)^(-2) #### Try to get the root, and if there is no solution return 0; Since the solution is sqrt(1/p), we square it to get p(-1), and tke its inverse to get p
		catch return 0							             #### return 0 if no solution
	 end
	end
 values= zeros(np*np)
	for(i in 1:(np*np))
		values[i]= getroot(nodes[i,:])
	end
 meanprice= sum(weights.*values)
	secondmoment= weights.*(values.*values)
 varprice= sum(secondmoment)-(meanprice^2)
	return Dict("mean"=>meanprice,"var"=>varprice)
	end


function question_2b(n)
 echan= MvNormal([[.002 .001],[ .001 .001]]) ##### The distributions of log theta1 and log theta2
	echan =exp(rand(echan,n))		##### Now that we have the distributions of Theta_1 and Theta_2,
################### We now apply the same procedure as in the last exercise
 function getroot(a) #### The function to apply to each row of nodes
		function minifunct(x) #### The function ap^(-1) +b p^(-1/2) -2 =0
			 return a[1]*x +a[2]*sqrt(x)-2  ###### Since the shock distributions are given in logs, we need to take the exponential to get the price equations
			end
		try return fzero(minifunct,0,5)^(-2) #### Try to get the root, and if there is no solution return 0; Since the solution is sqrt(1/p), we square it to get p, and tke its inverse to get p
		catch return 0							             #### return 0 if no solution
	 end
	end
 values= zeros(n)                 ### Sampled prices
	for(i=1:n)
		values[i]=getroot(echan[:,i])
	end
	return Dict("mean"=>mean(values),"var"=>var(values))
	end


	# function to run all questions
	function runall(n=10)
		println("running all questions of HW-integration:")
		println("results of question 1:")
		question_1b(n)	# make sure your function prints some kind of result!
		question_1c(n)
		question_1d(n)
		println("")
		println("results of question 2:")
		q2 = question_2a(n)
		println(q2)
		q2b = question_2b(n)
		println(q2b)
		println("end of HW-integration")
	end

end
