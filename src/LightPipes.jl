#module LightPipes
using FFTW
using Plots

mutable struct LPField
    field::Matrix{ComplexF64}
    size::Float64
    lambda::Float64
    dim::Int32
    curv::Float64
end


function LPBegin(size,lambda,dim)
    field = Matrix{ComplexF64}(undef,dim,dim)
    field .= 1.0
    lambda = lambda
    size = size
    curv = 0.0
    return LPField(field,size,lambda,dim,curv) 
end



function LPCircAperture(rad, xc, yc, U::LPField)
    dx = U.size /(U.dim - 1.0)
    for j in 1:U.dim
        y = (j - U.dim/2.0)*dx - yc
        for i in 1:U.dim
            x = (i - U.dim/2.0)*dx - xc
            if x^2+y^2 > rad^2
                U.field[i,j]=0.0
            end
        end
    end
    return U
end
"""  Bilinear interpolation on 4 points in square grid """
function Interpol4(coordX::Float64, coordY::Float64, v00::ComplexF64, v01::ComplexF64, v10::ComplexF64, v11::ComplexF64 )
    cx= 1.0 - coordX
    cy = 1.0 - coordY
    f = cx * cy * v00 + coordX * cy * v10 + cx * coordY * v01 + coordX * coordY * v11
end   
""" Interpolation of the field into a new grid """
function LPInterpol(new_size, new_dim, x_shift , y_shift, U::LPField)
    U1 = Matrix{ComplexF64}(undef, new_dim, new_dim)
    dx = U.size /(U.dim - 1.0)
    dx1 = new_size/(new_dim - 1.0)
    xc = x_shift
    yc = y_shift
    for j in 1:new_dim
        y1 = (j - new_dim / 2.0) * dx1 + yc
        j1 = floor(y1/dx) 
        yy=  y1/dx - j1  
        j1 +=  U.dim / 2 
        j2 = j1 + 1
        j1 =Int(round(j1))
        j2=j1 + 1
        for i in 1:new_dim
            x1 = (i - new_dim / 2.0) * dx1 + xc
            i1 =floor(x1/dx)
            xx = x1/dx - i1  
            i1  +=  U.dim / 2
            i2 = i1 + 1
            i1 =Int(round(i1))
            i2=i1 + 1
            if 1 <= i1  <= (U.dim - 1) && 1 <= j1  <= (U.dim - 1) 
                #println(i1," ",j1)
                U1[i,j] = Interpol4(xx, yy, U.field[i1,j1], U.field[i1,j2], U.field[i2,j1], U.field[i2,j2])
            else
                U1[i,j] = 0.0 + 0.0im
            end
        end
    end
    return LPField(U1, new_size, U.lambda, new_dim, U.curv)
end

function LPFFT(U::LPField)
    U.field = fftshift(U.field)
    U.field = fft(U.field)
    U.field = fftshift(U.field)
    return U
end


function LPIFFT(U::LPField)
    U.field = fftshift(U.field)
    U.field = ifft(U.field)
    U.field = fftshift(U.field)
    return U
end

function LPForvard(z::Float64,U::LPField)
    nc = U.dim / 2 + 1
    wave_num = 2.0 * pi / U.lambda
    U = LPFFT(U)
    for j in 1:U.dim
        jc = j-nc
        for i in 1:U.dim
            ic = i - nc
            serv = (wave_num * ic * U.lambda / U.size)^2
            serv = serv + (wave_num * jc * U.lambda / U.size)^2
            serv = z * sqrt(wave_num^2 - serv)
            U.field[i,j]  = U.field[i,j] * exp((0.0 + 1.0im) * serv)
        end
    end
    U = LPIFFT(U)
    return U
end
""" Lens function """
function LPLens(ff::Float64,U::LPField)
    wave_num =2.0 * pi / U.lambda
    dx = U.size /(U.dim - 1.0)
    nc = U.dim/2+1
    for j in 1:U.dim
        y = (j - nc)*dx
        for i in 1:U.dim
            x = (i - nc)*dx
            serv = exp((0.0 - 1.0im)* (x*x + y*y) *  wave_num/2.0/ff)
            U.field[i,j] = U.field[i,j] * serv
        end
    end
    return U
end

function LPMix(U::LPField , U1::LPField)
    if U.size == U1.size && U.dim == U1.dim 
        U.field = U.field .+ U1.field
        return U
    else
        println("LPMix error: fields with different parameters")
    end
end

""" Introducing tilts in the wavefront 
 tiltx means tilt in the direction x,  around Y axis  """
function LPTilt(tiltx , tilty, U::LPField)
    wave_num =2.0 * pi / U.lambda
    dx = U.size /(U.dim - 1.0)
    for j in 1:U.dim
        y = (j - U.dim/2.0)*dx
        for i in 1:U.dim
            x = (i - U.dim/2.0)*dx
            serv = exp((0.0 - 1.0im)* (tiltx * x + tilty * y) *  wave_num)
            U.field[i,j] = U.field[i,j] * serv
        end
    end
    return U
end

function LPCopy(U::LPField)
    U1 = deepcopy(U)
    return U1
end

U = LPBegin(0.05, 1e-6, 1024)
#U1 = LPBegin(0.1, 1e-6, 512)
#U1 = LPCopy(U)
U = LPCircAperture(0.004, 0.00, 0.00, U)
#U1 = LPCircAperture(0.05, 0.0, 0.00, U1)
U = LPLens(5.0, U)
#U=LPMix(U , U1)
 U = LPForvard(5.0, U) 
 #U1= LPInterpol(0.4, 2100, 0.1, -0.0, U)
 #U = LPForvard(10.0 , U1) 
 #U1= LPInterpol(0.3, 2000, -0.1, 0.02, U)
 #=
 for ii in 1:30
    global U1= LPInterpol(0.3,299, 0., -0.0, U)
    global U= LPInterpol(2.3,2001, 0.0, -0U1= LPInterpol(0.3,400, 0., -0.0, U).0, U1)
 end
 =#
 
display(heatmap((abs.(U.field)).^2))


#@time U = LPFFT(U)
#display(heatmap(abs.(U.field)))
#=
for i in 1:U.dim
    println(real.(U.field[i,:]))
end
=#



#end # module

#=
""" example of the same but with immutable struct
in general looks less efficient, due to more memory and garbage collection
"""
struct LPFieldd
    field::Matrix{ComplexF64}
    size::Float64
    lambda::Float64
    dim::Int32
    curv::Float64
end

function LPBegind(size,lambda,dim)
    
    field = Matrix{ComplexF64}(undef,dim,dim)
    field .= 1.0
    lambda = lambda
    size = size
    curv = 0.0
    return LPFieldd(field,size,lambda,dim,curv) 
end

function LPFFT(U::LPFieldd)
    U1= fftshift(U.field)
    U1 = fft(U1)
    U1 = fftshift(U1)
    U.field .= (U1) # this is the trick to update a mutable vector !!!
    return U
end
=#