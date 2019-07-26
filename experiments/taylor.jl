using TaylorModels

wh(x) = print("type is: $(typeof(x))\n")
disp(msg) = print("$(msg)\n")

f(x) = sum([n*sin(x^2 - x + 1) for n in 1:15])
A =  -1.0 .. 1.0 # Domain 
x0 = mid(A)     # Expansion point
tm = TaylorModel1(7, A) 
ftm = f(tm)

wh(tm)
wh(ftm)
disp( ftm(1) )

# Now the plot
using Plots; pyplot()
# /bin/bash: q: command not found
plot!(ftm, label="taylor model")
savefig("taylor_plot.png")
print("done\n")
