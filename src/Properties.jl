

function energy(P; shift=2)
	
	ED = zeros(P.lp[1], P.lp[2])
	
	Threads.@threads for j in (1+shift):(P.lp[2]-shift)
		for i in (1+shift):(P.lp[1]-shift)

			Ppt = getP(P,i,j)
			DPpt = getDP(P,i,j)

			ED[i,j] = engpt(Ppt,DPpt, P.parameters)

		end
	end

	total_energy = sum(ED)*P.ls[1]*P.ls[2]

	#if P.is_electrostatic
	#	total_energy += DipoleEnergy(P, sums, type, shift)
	#end
	
	return total_energy

end

function engpt(Ppt,DPpt,pars)
	
	return PotEng(Ppt, pars.A2t, pars.A4t, pars.V0) + DerEng(DPpt,pars.G4)
	
end


function DerEng(DPpt,G4)
	
	derengpt = 0.0

	for i in 1:2, j in 1:2, a in 1:3, b in 1:3
		derengpt += G4[i,a,j,b]*DPpt[i,a]*DPpt[j,b]
	end
		
	return derengpt
	
end


function PotEng(Ppt, A2t, A4t, V0)

	engpt = 0.0

	for a in 1:3, b in 1:3
		
		engpt += A2t[a,b]*Ppt[a]*Ppt[b]
	
		for c in 1:3, d in 1:3
			engpt += A4t[a,b,c,d]*Ppt[a]*Ppt[b]*Ppt[c]*Ppt[d]
		end
	end

	engpt -= V0
	
	return engpt
	
end

