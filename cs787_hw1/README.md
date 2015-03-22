hw Problem description  
--  
In class, we learned how dynamic programming had its origins in Fermat’s principle of least  
time. This problem deals with another connection between computer science and physics. It  
is known (*) that when a particle of mass m moves in a force field, it “chooses” its path so as  
to minimize the so-called action integral  

	A = Integral_t1_t2(T - V)dt

where T and V denote the kinetic and potential energies at time t. This means we can apply  
shortest path algorithms to solve first-year physics problems.  
As an example, consider a particle of mass m in a uniform gravitational field, say a ball dropped  
from a high place. Let the height coordinate be y, and the gravitational strength g.  
a) Suppose at time t, the particle is at y, and moves to position y −∆y at time t+∆t. Write  
down expressions for T and V that are accurate to first order in ∆t.  
b) Imagine a grid with time in the horizontal direction and height in the vertical dimension.  
Explain how you can sweep through the times, and find (for every grid point), an approximation  
to the action-minimizing path to that point. If the grid is m ×n, how many steps  
does your algorithm need?  
c) Using a computer and your answer to b), make a plot of the approximate path. If you  
didn’t get a parabolic shape, something is wrong and you have a bug in your program.  
(How you implement this is up to you.)  
You are free to substitute another problem for free fall if you wish. (Some ideas: harmonic  
motion with non-linear restoring force, planetary motion with funny force fields, ...) Plotting  
will be easiest, however, if the configuration space for your problem is one-dimensional.  
(*) There is some fine print on this one. Technically, the path only has to be “first order  
stationary” (similar to the derivative vanishing in calculus), so it might not be a minimizer.  
For this problem, and many others, though, you can just imagine it is. See Lecture 19 of R.  
P. Feynman, Lectures on Physics, v. 2, if you want to learn more.  
