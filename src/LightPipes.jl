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


function LPMix(U::LPField , U1::LPField)
    if U.size == U1.size && U.dim == U1.dim 
        U.field = U.field .+ U1.field
        return U
    else
        println("LPMix error: fields with different parameters")
    end
end

function LPTilt(tiltx , tilty, U::LPField)
    wave_num =2.0 * pi / U.lambda
    dx = U.size /(U.dim - 1.0)
    for j in 1:U.dim
        y = (j - U.dim/2.0)*dx
        for i in 1:U.dim
            x = (i - U.dim/2.0)*dx
            serv = exp((0.0 + 1.0im)* (tiltx * x + tilty * y) *  wave_num)
            U.field[i,j] = U.field[i,j] * serv
        end
    end
    return U
end

function LPCopy(U::LPField)
    U1 = deepcopy(U)
    return U1
end

U = LPBegin(0.2, 1e-6, 512)
#U1 = LPBegin(0.1, 1e-6, 512)
U1 = LPCopy(U)
U = LPCircAperture(0.04, 0.00, 0.00, U)
U1 = LPCircAperture(0.05, 0.0, 0.00, U1)
U1 = LPTilt(5e-5, 5e-5 , U1)
U=LPMix(U , U1)
 U = LPForvard(0.9*0.02^2/1e-6 , U) 
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