"""
the structure LPField contains the field data. 
Square matrix with complex intensity, size is the side of the grid, in [m], [um], etc.
lambda - wavelength, same units as size, 
dim - grid dimension,  for instance dim =1024 for a 1024x1024 grid, 
curv - geometrical curvature of the beam. Default curv = 0.0 
"""
mutable struct LPField
    field::Matrix{ComplexF64}
    size::Float64
    lambda::Float64
    dim::Int32
    curv::Float64
end



"""
Initializing LPField with the user's data:
the structure LPField contains the field data. 
field: Square complex matrix, 
size:  the side of the grid, in [m], [um], etc.
lambda: - wavelength, same units as size, 
dim: - grid dimension,  for instance dim =1024 for a 1024x1024 grid, 
curv: - geometrical curvature of the beam WF. Default curv = 0.0 

The coordinate systems, as shown on screen: 
x - horizontal, left to the right
y - vertical, top to  bottom
angle: positive when rotating CLOCKWISE from x to y
angle defined in degrees. 360 degrees = 2*pi rad. 
"""
function LPBegin(size,lambda, dim)
    if iseven(dim)  # the array dimension should be even
        field = Matrix{ComplexF64}(undef,dim,dim)
        field .= 1.0
        lambda = lambda
        size = size
        curv = 0.0
        return LPField(field,size,lambda,dim,curv) 
    else
            println("LPBegin: the array dimension must be even, like 1002, 888, etc.")
    end
end


""" 
Introducing circular aperture in the existing field 
rad - aperture radius
xc,yc - coordinates of the aperture center
"""
function LPCircAperture(rad, xc, yc, U::LPField)
    yc = -yc
    dx = U.size /(U.dim - 1.0)
    dy = dx
    for j in 1:U.dim
        #jj = U.dim -j +1
        y = (j - 1 - U.dim/2.0)*dy - yc
        for i in 1:U.dim
            x = (i - 1 - U.dim/2.0)*dx - xc
            if x^2+y^2 > rad^2
                U.field[j,i]=0.0
            end
        end
    end
    return U
end



"""
Same as aperture, but screen
"""
function LPCircScreen(rad, xc, yc, U::LPField)
    yc = -yc
    dx = U.size /(U.dim - 1.0)
    dy = dx
    for j in 1:U.dim
        #jj = U.dim -j +1
        y = (j - 1 - U.dim/2.0)*dy - yc
        for i in 1:U.dim
            x = (i - 1 - U.dim/2.0)*dx - xc
            if x^2+y^2 < rad^2
                U.field[j,i]=0.0
            end
        end
    end
    return U
end

"""  Bilinear interpolation on 4 points in square grid """
function Interpol4(coordX::Float64, coordY::Float64, v00::ComplexF64, v01::ComplexF64, v10::ComplexF64, v11::ComplexF64 )
    cx= 1.0 - coordX
    cy = 1.0 - coordY
    #if cx < 0 || cy < 0  println(cx,"   ",cy)  end
    f = cx * cy * v00 + coordX * cy * v10 + cx * coordY * v01 + coordX * coordY * v11
end   
""" Interpolation of the field into a scaled grid """
function LPInterpol(new_size, new_dim, x_shift , y_shift, U::LPField)
    y_shift = -y_shift
    U1 = Matrix{ComplexF64}(undef, new_dim, new_dim)
    dx = U.size /(U.dim - 1.0)
    dx1 = new_size/(new_dim - 1.0)
    dy = dx
    dy1 = dx1
    xc = x_shift
    yc = y_shift
    for j in 1:new_dim
        #jj = new_dim - j + 1  # inversion of the "new" y axis
        y1 = (j  -1 - new_dim / 2.0) * dy1 + yc
        j1 = floor(y1/dy) 
        yy=  y1/dy - j1  
        j1 +=  U.dim / 2 
        j1 =Int(round(j1)) +1
        #j1 = U.dim - j1 + 1 # matching inversion of the "old" y axis
        j2 = j1 + 1
        #j1, j2 = j2, j1

        for i in 1:new_dim
            x1 = (i  -1   - new_dim / 2.0) * dx1 + xc
            i1 =floor(x1/dx)
            xx = x1/dx - i1  
            i1  +=  U.dim / 2
            i1 =Int(round(i1)) + 1
            i2=i1 + 1
            if 1 <= i1  <= (U.dim - 1) && 1 <= j1  <= (U.dim - 1) 
                #println(i1," ",j1)
                U1[j,i] = Interpol4( yy, xx, U.field[j1,i1], U.field[j1,i2], U.field[j2,i1], U.field[j2,i2])
            else
                U1[j,i] = 0.0 
            end
        end
    end
    return LPField(U1, new_size, U.lambda, new_dim, U.curv)
end

""" Shifts the coordinate system origin to the point dx,dy """
function LPShift(dx::Number,dy::Number,U::LPField)
    U1 = LPInterpol(U.size, U.dim, dx, dy, U)
end

""" Rotates the Field CLOCKWISE by angle [degrees] 
use negative angle to rotate counterclockwise """
function LPRotate(angle::Number, U::LPField)
   
    U1 = Matrix{ComplexF64}(undef, U.dim, U.dim)
    dx = U.size /(U.dim - 1.0)
    dy = dx
    ca = cos(angle * 2.0 *pi /360.)
    sa = sin(angle * 2.0 *pi /360.)
    for j in 1:U.dim
        y = (j - 1 - U.dim/2.0)*dy 
        for i in 1:U.dim
            x = (i - 1 - U.dim/2.0)*dx
            x1 = x*ca - y*sa
            y1 = x*sa + y*ca
            
            i1 =floor(x1/dx)
            xx = x1/dx - i1  
            i1  +=  U.dim / 2
            i1 =Int(round(i1)) + 1
            i2=i1 + 1

            j1 = floor(y1/dy) 
            yy=  y1/dy - j1  
            j1 +=  U.dim / 2 
            j1 =Int(round(j1)) +1
            #j1 = U.dim - j1 + 1 # matching inversion of the "old" y axis
            j2 = j1 + 1
            #j1, j2 = j2, j1
            if 1 <= i1  <= (U.dim - 1) && 1 <= j1  <= (U.dim - 1) 
                #println(i1," ",j1)
                U1[j,i] = Interpol4( yy, xx, U.field[j1,i1], U.field[j1,i2], U.field[j2,i1], U.field[j2,i2])
            else
                U1[j,i] = 0.0 
            end
        end
    end
    return LPField(U1, U.size, U.lambda, U.dim, U.curv)
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

function LPForvard(z::Number,U::LPField)
    nc = U.dim / 2 + 1
    wave_num = 2.0 * pi / U.lambda
    U = LPFFT(U)
    for j in 1:U.dim
        jj = U.dim - j +1
        jc = jj-nc
        for i in 1:U.dim
            ic = i - nc
            serv = (wave_num * ic * U.lambda / U.size)^2
            serv = serv + (wave_num * jc * U.lambda / U.size)^2
            serv = z * sqrt(wave_num^2 - serv)
            U.field[j,i]  = U.field[j,i] * exp(1.0im * serv)
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

""" Introducing tilts in the wavefront 
 tiltx means tilt in the direction x = (rotation around Y axis)  """
function LPTilt(tiltx , tilty, U::LPField)
    tilty = -tilty
    wave_num =2.0 * pi / U.lambda
    dx = U.size /(U.dim - 1.0)
    dy = dx
    for j in 1:U.dim
        #jj = U.dim -j +1
        y = (j - 1 - U.dim/2.0)*dy
        for i in 1:U.dim
            x = (i - 1 - U.dim/2.0)*dx
            serv = exp((0.0 + 1.0im)* (tiltx * x + tilty * y) *  wave_num)
            U.field[j,i] = U.field[j,i] * serv
        end
    end
    return U
end

function LPLens(F, U::LPField)
    wave_num =2.0 * pi / U.lambda
    dx = U.size /(U.dim - 1.0)
    dy = dx
    ud2 = U.dim/2.
    for j in 1:U.dim
        #jj = U.dim -j +1
        
        y = (j - 1.0 - ud2)*dy
        y2 = y*y
        for i in 1:U.dim
            x = (i - 1.0 - ud2)*dx
            rho2 = x^2 +y2
            F2=F^2
            serv = exp((1.0im)* sign(F) *sqrt(F2 - rho2) * wave_num) # spherical WF 
            #serv = exp((-1.0im)* rho2/2.0/F * wave_num) # parabolic WF 
            U.field[j,i] = U.field[j,i] * serv
        end
    end
    return U
end

function LPCopy(U::LPField)
    U1 = deepcopy(U)
    return U1
end
