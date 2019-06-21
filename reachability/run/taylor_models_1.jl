using TaylorModels

f(x) = sum([n*sin(x^2 - x + 1) for n in 1:15])
a =  -1.0 .. 1.0 # Domain 
x0 = mid(a)     # Expansion point
# Independent variables for Taylor models, order 6,7,8
tm6 = TaylorModel1(6, interval(x0), a) 
tm7 = TaylorModel1(7, interval(x0), a)
tm8 = TaylorModel1(8, interval(x0), a)
# Taylor models corresponding to f(x) of order 6,7,8
ftm6 = f(tm6)
ftm7 = f(tm7)
ftm8 = f(tm8)

# Now the plot
using Plots; pyplot()
plot(range(inf(a), stop=sup(a), length=1000), x->f(x), label="f(x)", lw=2, xaxis="x", yaxis="f(x)")
plot!(ftm6, label="6th order")
plot!(ftm7, label="7th order")
plot!(ftm8, label="8th order")
savefig("taylor_model_plot_1.png")
print("done\n")
