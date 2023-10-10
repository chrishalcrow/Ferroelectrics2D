
function getP(P,i,j)
    return SVector{3,Float64}(
        P.field[i,j,1],
        P.field[i,j,2],
        P.field[i,j,3]
         )
end

function getDP(P,i,j)

    return SMatrix{2,3,Float64,6}(

    dx(P,i,j,1),
    dy(P,i,j,1),

    dx(P,i,j,2),
    dy(P,i,j,2),

    dx(P,i,j,3),
    dy(P,i,j,3)

    )

end

function getDDP(P,i,j)

    return SMatrix{3,3,Float64, 9}(

        d2x(P,i,j,1),  # [1,1]
        d2y(P,i,j,1), #[2,1]
        dxdy(P,i,j,1),

        d2x(P,i,j,2), #[1,2]
        d2y(P,i,j,2),
        dxdy(P,i,j,2),

        d2x(P,i,j,3),
        d2y(P,i,j,3),
        dxdy(P,i,j,3)

    )
    
end
function dx(P, i, j, a)
    @inbounds (-P.field[i+2,j,a] + 8.0*P.field[i+1,j,a] - 8.0*P.field[i-1,j,a] + P.field[i-2,j,a] )/(12.0*P.ls[1])
end

function dy(P, i, j, a)
    @inbounds (-P.field[i,j+2,a] + 8.0*P.field[i,j+1,a] - 8.0*P.field[i,j-1,a] + P.field[i,j-2,a] )/(12.0*P.ls[2])
end


function d2x(P, i, j, a)
    @inbounds (-P.field[i+2,j,a] + 16.0*P.field[i+1,j,a] - 30.0*P.field[i,j,a] + 16.0*P.field[i-1,j,a] - P.field[i-2,j,a])/(12.0*P.ls[1]^2)
end
function d2y(P, i, j, a)
    @inbounds (-P.field[i,j+2,a] + 16.0*P.field[i,j+1,a] - 30.0*P.field[i,j,a] + 16.0*P.field[i,j-1,a] - P.field[i,j-2,a])/(12.0*P.ls[2]^2)
end

function dxdy(P, i, j, a)
    @inbounds 0.5*( dxdydiff(P, i, j, a) - d2x(P, i, j, a) - d2y(P, i, j, a) )
end
function dxdydiff(P, i, j, a)
    @inbounds (-P.field[i+2,j+2,a] + 16.0*P.field[i+1,j+1,a] - 30.0*P.field[i,j,a] + 16.0*P.field[i-1,j-1,a] - P.field[i-2,j-2,a])/(12.0*P.ls[1]*P.ls[2])
end


#=

function d2xD(phi, a,i)
    @inbounds (-phi.field[a,i+2] + 16.0*phi.field[a,i+1] - 30.0*phi.field[a,i] + 16.0*phi.field[a,i-1] - phi.field[a,i-2] )/(12.0*phi.ls^2)
end





function dxD2D(phi, a, i, j, lsx)
    @inbounds (-phi[a,i+2,j] + 8.0*phi[a,i+1,j] - 8.0*phi[a,i-1,j] + phi[a,i-2,j])/(12.0*lsx)
end
function dyD2D(phi, a, i, j, lsy)
    @inbounds (-phi[a,i,j+2] + 8.0*phi[a,i,j+1] - 8.0*phi[a,i,j-1] + phi[a,i,j-2])/(12.0*lsy)
end

function d2xD2D(phi, a, i, j, lsx)
    @inbounds (-phi[a,i+2,j] + 16.0*phi[a,i+1,j] - 30.0*phi[a,i,j] + 16.0*phi[a,i-1,j] - phi[a,i-2,j])/(12.0*lsx^2)
end
function d2yD2D(phi, a, i, j, lsx)
    @inbounds (-phi[a,i,j+2] + 16.0*phi[a,i,j+1] - 30.0*phi[a,i,j] + 16.0*phi[a,i,j-1,] - phi[a,i,j-2])/(12.0*lsx^2)
end

function dxdydiffD2D(phi, a, i, j, lsx, lsy)
    @inbounds (-phi[a,i+2,j+2] + 16.0*phi[a,i+1,j+1] - 30.0*phi[a,i,j] + 16.0*phi[a,i-1,j-1] - phi[a,i-2,j-2])/(12.0*lsx*lsy)
end


function dxD2Du(phi, a, b, i, j, lsx)
    @inbounds (-phi[a][b,i+2,j] + 8.0*phi[a][b,i+1,j] - 8.0*phi[a][b,i-1,j] + phi[a][b,i-2,j])/(12.0*lsx)
end
function dyD2Du(phi, a, b, i, j, lsy)
    @inbounds (-phi[a][b,i,j+2] + 8.0*phi[a][b,i,j+1] - 8.0*phi[a][b,i,j-1] + phi[a][b,i,j-2])/(12.0*lsy)
end

function d2xD2Du(phi, a, b, i, j, lsx)
    @inbounds (-phi[a][b,i+2,j] + 16.0*phi[a][b,i+1,j] - 30.0*phi[a][b,i,j] + 16.0*phi[a][b,i-1,j] - phi[a][b,i-2,j])/(12.0*lsx^2)
end
function d2yD2Du(phi, a,b, i, j, lsx)
    @inbounds (-phi[a][b,i,j+2] + 16.0*phi[a][b,i,j+1] - 30.0*phi[a][b,i,j] + 16.0*phi[a][b,i,j-1,] - phi[a][b,i,j-2])/(12.0*lsx^2)
end

function dxdydiffD2Du(phi, a,b, i, j, lsx, lsy)
    @inbounds (-phi[a][b,i+2,j+2] + 16.0*phi[a][b,i+1,j+1] - 30.0*phi[a][b,i,j] + 16.0*phi[a][b,i-1,j-1] - phi[a][b,i-2,j-2])/(12.0*lsx*lsy)
end



function getP2D(P,i,j)
    return [ P[1,i,j], P[2,i,j], P[3,i,j] ]
end


function getDP2D!(DPpt,P,i,j,dx,dy)
    
    for a in 1:3
        DPpt[a,1] = dxD2D(P,a,i,j,dx)
        DPpt[a,2] = dyD2D(P,a,i,j,dy)
    end
    
end

function getDDP!(DPpt,P,i,j,dx,dy)
    
    for a in 1:3

        DPpt[a,1,1] = d2xD2D(P,a,i,j,dx)
        DPpt[a,2,2] = d2yD2D(P,a,i,j,dy)
        DPpt[a,1,2] = (dxdydiffD2D(P, a, i, j, dx, dy) - DPpt[a,1,1] - DPpt[a,2,2])/2.0
        DPpt[a,2,1] = DPpt[a,1,2]

    end
    
end


=#