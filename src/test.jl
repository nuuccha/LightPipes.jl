using FFTW
using Plots
@time U = LPBegin(0.1, 1e-6, 2000)
#U1 = LPBegin(0.1, 1e-6, 512)
#U1 = LPCopy(U)
U = LPCircAperture(0.03, 0.0, 0.0, U)
U = LPCircAperture(0.03, 0.005, 0.0, U)
U = LPCircScreen(0.005, 0.02, 0.02, U)
U = LPCircScreen(0.005, -0.02, 0.02, U)
U = LPShift(0.0,0.002, U)

U = LPRotate(15,U)
U= LPLens(20.0,U)
U= LPForvard(20,U)
U= LPLens(5.0,U)
U= LPForvard((5*20)/(20 - 5),U)
U = LPRotate(180,U)

U= LPInterpol(0.1, 500, 0.0, -0.0, U)
intens = (abs.(U.field)) .* (abs.(U.field))
display(heatmap(intens))
