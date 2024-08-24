These are original codes for some of the figures in Xie&Xiao PDRK paper (2014arXiv or 2016PST). Note, these codes may be different from the final version. Some figures (doppler/low hybrid) have been removed due to that the original data are too large.

I provide them here is only to hope that they can help the readers understand how PDRK can be used, i.e., how to sort the solution, how to remove artificial solutions, how to make surface plot.

Defaultly, the code uses ‘eig()’ to solve all the solutions. You can also use it as the other root finding solver with initial guess using ‘eigs()’ in Matlab to obtain one or several solutions around the initial guess.

Hua-sheng XIE (huashengxie@gmail.com)
08:39 2018/2/22