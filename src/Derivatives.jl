
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
    @inbounds (-P.field[i+2,j,a] + 8.0*P.field[i+1,j,a] - 8.0*P.field[i-1,j,a] + P.field[i-2,j,a] )/(12.0*P.grid.ls[1])
end

function dy(P, i, j, a)
    @inbounds (-P.field[i,j+2,a] + 8.0*P.field[i,j+1,a] - 8.0*P.field[i,j-1,a] + P.field[i,j-2,a] )/(12.0*P.grid.ls[2])
end


function d2x(P, i, j, a)
    @inbounds (-P.field[i+2,j,a] + 16.0*P.field[i+1,j,a] - 30.0*P.field[i,j,a] + 16.0*P.field[i-1,j,a] - P.field[i-2,j,a])/(12.0*P.grid.ls[1]^2)
end
function d2y(P, i, j, a)
    @inbounds (-P.field[i,j+2,a] + 16.0*P.field[i,j+1,a] - 30.0*P.field[i,j,a] + 16.0*P.field[i,j-1,a] - P.field[i,j-2,a])/(12.0*P.grid.ls[2]^2)
end

function dxdy(P, i, j, a)
    @inbounds 0.5*( dxdydiff(P, i, j, a) - d2x(P, i, j, a) - d2y(P, i, j, a) )
end
function dxdydiff(P, i, j, a)
    @inbounds (-P.field[i+2,j+2,a] + 16.0*P.field[i+1,j+1,a] - 30.0*P.field[i,j,a] + 16.0*P.field[i-1,j-1,a] - P.field[i-2,j-2,a])/(12.0*P.grid.ls[1]*P.grid.ls[2])
end

