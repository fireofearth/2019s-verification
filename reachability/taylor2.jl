using TaylorModels

wh(x) = print("type is: $(typeof(x))\n")
disp(msg) = print("$(msg)\n")
header(msg) = print("\n$(msg)\n")

header("Initialization")
X = -0.5 .. 0.5
Y = 1.0 .. 2.0
t1 = TaylorModel1(3, X)
t2 = TaylorModel1(3, Y)
disp(t1)
showfull(t1)
disp(t2)
showfull(t2)

header("Addition")
disp(t1 * t2)

header("Multiplication")
disp(t1 * t2)


