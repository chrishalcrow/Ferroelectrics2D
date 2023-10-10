


function gradient_flow!(P; shift=2, steps = 1, dt = P.ls[1]/500.0, tolerance = 0.0, checks = max(100,steps), print_stuff = true, dEdp = zeros(P.lp[1],P.lp[2],3) )


    if tolerance == 0 && checks > steps
        checks = steps
    end
    
    if print_stuff == true
        println("initial: energy: ", energy(P, shift=shift) )
    end

    counter = 0
    prev_error = 1.0e9
    
    while counter < steps
        
        gradient_flow_for_n_steps!(P,dEdp,checks,dt)
        
        err = max_abs_err(dEdp)
        if err > 3*prev_error
            error("Suspected numerical blowup. Please use a smaller dt. Currently, dt = ", dt)
        end
        prev_error = err

        counter += checks
        
        if print_stuff == true
            println("after ", counter, " steps, error = ", round(err, sigdigits=4))
        end

        if tolerance != 0.0    # => we are in tol mode    
            if err < tolerance
                counter = steps + 1    # => end the while loop
            else
                steps += checks    # => continue the while loop
            end
        end

		if isnan(err)
			println("error = NaN, choose smaller dt")
			break
		end

    end

    if print_stuff == true
        println("final energy: ", energy(P, shift=shift) )
    end

    return

end


function gradient_flow_for_n_steps!(ϕ,dEdp,n,dt)
    for _ in 1:n
        gradient_flow_1_step!(ϕ,dEdp,dt)
    end
end

function gradient_flow_1_step!(P, dEdp, dt)

	getdEdP!(P,dEdp)
    P.field .-= dt.*dEdp;
   
end


function getdEdP!(P,dedp)

	G4 = P.parameters.G4
	A2s = P.parameters.A2s
	A4s = P.parameters.A4s

	Threads.@threads for i in 3:P.lp[1]-2
        for j in 3:P.lp[2]-2

            @inbounds for a in 1:3
                dedp[i,j,a] = 0.0
            end

		Ppt = getP(P,i,j)
		DDPpt = getDDP(P,i,j)
					
		@inbounds for a in 1:3, b in 1:3
		
			
			dedp[i,j,a] += A2s[a,b]*Ppt[b]
            
            dedp[i,j,a] += -(G4[1,b,1,a] + G4[1,a,1,b])*DDPpt[1,b]
            dedp[i,j,a] += -(G4[2,b,2,a] + G4[2,a,2,b])*DDPpt[2,b]
            dedp[i,j,a] += -(G4[1,b,2,a] + G4[1,a,2,b] + G4[2,b,1,a] + G4[2,a,1,b])*DDPpt[3,b]
			
			for c in 1:3, d in 1:3
				dedp[i,j,a] += (A4s[a,b,c,d])*Ppt[b]*Ppt[c]*Ppt[d]
			end
			
		end
	
        end

    end
        
end



function max_abs_err(A)
    return maximum(abs, A) 
end






"""
    arrested_newton_flow!(skyrmion; skyrmion_dot, steps = n, tolerance = tol, dt=ls^2/80.0, frequency_of_checking_tolerance = freq, print_stuff = true)
    
Applies an arrested Newton flow to `skyrmion` whose initial time derivative field is skyrmion_dot with timestep `dt`, either for `n` steps or until the error falls below `tol`. The error is checked every `checks` steps.

See also [`gradient_flow!`, `newton_flow!`]
"""
function arrested_newton_flow!(ϕ; ϕd=zeros(3, ϕ.lp), dt=ϕ.ls/50.0, steps=1, tolerance = 0.0, checks = max(100,steps), print_stuff = true)

    if tolerance == 0 && checks > steps
        checks = steps
    end

    energy_density = zeros(ϕ.lp)
    old_field = deepcopy(ϕ.field);

    dEdp = zeros(3, ϕ.lp)
    dEdp2 = zeros(3, ϕ.lp)
    dEdp3 = zeros(3, ϕ.lp)
    dEdp4 = zeros(3, ϕ.lp)
    ϕ2 = deepcopy(ϕ)


    counter = 0
    while counter < steps

        arrested_newton_flow_for_n_steps!(ϕ,ϕ2,ϕd,old_field,dEdp,dEdp2,dEdp3,dEdp4,dt,energy_density,checks, EnergyANF(ϕ,energy_density))
        error = max_abs_err(dEdp)
        counter += checks

        if print_stuff == true 
            println("after ", counter, " steps, error = ", round(error, sigdigits=4), " energy = ", round(sum(energy_density)*ϕ.ls, sigdigits=8) )
            #println( round(error, sigdigits=8), "," )
        end

        if tolerance != 0.0    # => we are in tol mode
            if error < tolerance 
                counter = steps + 1    # => end the while loop
            else
                steps += checks
            end
        end

    end

    return

end
 
function arrested_newton_flow_for_n_steps!(ϕ,ϕ2,ϕd,old_field,dEdp1,dEdp2,dEdp3,dEdp4,dt,energy_density,n, initial_energy)

    new_energy = initial_energy
    
    for _ in 1:n

        old_energy = new_energy
        old_field .= ϕ.field

        newton_flow_for_1_step!(ϕ,ϕ2,ϕd,dEdp1,dEdp2,dEdp3,dEdp4,dt)
        new_energy = EnergyANF(ϕ,energy_density)

		#println("Old, ", old_energy, " and new, ", new_energy)

        if new_energy > old_energy

			println("ARREST")

            fill!(ϕd, 0.0);
            ϕ.field .= old_field;

            if new_energy > 1.2*old_energy
                error("Suspected numerical blow-up. Please use smaller dt. Currently, dt = ", dt)
            end
    
        end

    end

end

function newton_flow_for_1_step!(sk, sk2, skd ,dEdp1, dEdp2, dEdp3, dEdp4, dt)

    getdEdP!(sk, dEdp1)
    sk2.field .= sk.field .+ (0.5*dt).*skd

    #sk.pion_field .+= (0.5*dt).*skd
    #getdEdp!(sk2, dEdp2, sk, dEdp1, dt)
    getdEdP!(sk2, dEdp2)
	sk.field .= sk2.field - (0.5*dt)^2 .*dEdp1

    ##sk.pion_field .+= (0.5*dt)^2 .*dEdp1
    #getdEdp3!(sk, dEdp3, sk2, skd, dEdp1, dEdp2, dt)
    getdEdP!(sk, dEdp3)
	sk2.field .= sk.field .+ (0.5*dt).*skd .- (0.5*dt)^2 .*(4.0 .* dEdp2 .- dEdp1)

    ##sk.pion_field .+= (0.5*dt).*skd .+ (0.5*dt)^2 .*(4.0 .* dEdp2 .- dEdp1)
    #getdEdp4!(sk2, dEdp4, sk, dEdp1, dEdp2, dEdp3, skd, dt)
    getdEdP!(sk2, dEdp4)
	sk.field .= sk2.field + dt.*(  (5/6*dt).*dEdp2  .- (dt/6).*( dEdp1 .+ dEdp3 ) )  

	skd .-= (dt/6.0).*(dEdp1 .+ 2.0.*dEdp2 .+ 2.0.*dEdp3 .+ dEdp4)

    #RESET field to original value: sk.pion_field .-= (0.5*dt).*skd + (0.5*dt)^2 *(4.0 .* dEdp2 - dEdp1) + (0.5*dt)^2 .*dEdp1 + (0.5*dt).*skd
    #Then update: sk.pion_field .+= dt.*(skd + dt/6.0 .*( dEdp1 .+ dEdp2 .+ dEdp3 ) ), combined into:
    ###sk.pion_field .-=  dt.*(  (5/6*dt).*dEdp2  .- (dt/6).*( dEdp1 .+ dEdp3 ) )
    
    #
   
    #orthog_skd_and_sk_and_normer!(skd,sk)
    #normer!(sk)

    ##orthog_skd_and_sk_and_normer!(skd,sk)
   
end
