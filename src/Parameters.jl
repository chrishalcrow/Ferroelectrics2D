function δ2(i,j)
    if i == j
        return 1.0
    else
        return 0.0
    end
end

function δ4(i,j,k,l)
    if i == j && j ==k && k==l
        return 1.0
    else
        return 0.0
    end
end

function δ6(i,j,k,l,m,n)
    if i == j && j ==k && k==l && l==m && m==n
        return 1.0
    else
        return 0.0
    end
end

function getA4s(A4t)
	A4s = zeros(3,3,3,3)
	for a in 1:3, b in 1:3, c in 1:3, d in 1:3
		A4s[a,b,c,d] = A4t[a,b,c,d] + A4t[b,a,c,d] + A4t[b,c,a,d] + A4t[b,c,d,a]
	end
	return A4s
end

function getA2s(A2t)
	A2s = zeros(3,3)
	for a in 1:3, b in 1:3
		A2s[a,b] = A2t[a,b] + A2t[b,a]
	end
	return A2s
end

function getLS(T)
	return sqrt( abs( 5.1e-10/(334000*(-381 + T)) ))*10^9
end


function set_parameters_lithium!(R2, material; is_electrostatic=false)
	

	TT = zeros(6,3,3)

	TT[1,1,1] = 1.0
	TT[2,2,2] = 1.0
	TT[3,3,3] = 1.0
#in the folliwing terms 0.5 due to assuming symmetry of initial matrix/tensor
	TT[4,2,3] = 0.5
	TT[4,3,2] = 0.5

	TT[5,1,3] = 0.5
	TT[5,3,1] = 0.5

	TT[6,2,1] = 0.5
	TT[6,1,2] = 0.5

    A2 = zeros(3,3)
	A2o = zeros(3,3)
	A2t = zeros(3,3)

	A4 = zeros(3,3,3,3)
	A4o = zeros(3,3,3,3)
	A4ot = zeros(3,3,3,3)
	A4t = zeros(3,3,3,3)
	
	Q4 = zeros(3,3,3,3)
	#Q4o = zeros(3,3,3,3)
	Q3 = zeros(6,3,3)
	Q3o = zeros(6,3,3)
	
	q4 = zeros(3,3,3,3)
	q4o = zeros(3,3,3,3)
    q3o = zeros(6,3,3)
	q3 = zeros(6,3,3)
	
	G2 = zeros(3,3)
	G4o = zeros(3,3,3,3)
	G4 = zeros(3,3,3,3)

	C4 = zeros(3,3,3,3)
	C4o = zeros(3,3,3,3)
	C2 = zeros(6,6)
    
	if (material == "CLN") || (material == "SLN")
		(C11,C12,C13,C14,C33,C44) = (1.9886e11, 0.5467e11, 0.6726e11, 0.0783e11, 2.3370e11, 0.5985e11)  #LiNbO3
		(e15, e22, e31, e33) = (3.658, 2.406, 0.373, 1.923)  #LiNbO3
		(eps11,eps22,eps33) = (84.3, 84.3, 28.9)  #LiNbO3
		Pvac = 0.75  #LiNbO3 C/m^2
		G11 = 0.0442 * 3.98e-11 #LiNbO3 N*m^4/C^2  it is g1 in Gopalan's paper
	elseif (material == "CLT") || (material == "SLT")
    	(C11,C12,C13,C14,C33,C44) = (2.3305e11, 0.4644e11, 0.8346e11, -0.1075e11, 2.7522e11, 0.9526e11)  #LiTaO3
        (e15, e22, e31, e33) = (2.628, 1.833, -0.188, 1.941)  #LiTaO3
        (eps11,eps22,eps33) = (52.7, 52.7, 44.0)  #LiTaO3
        Pvac = 0.525    #LiTaO3 C/m^2
        G11 = ((material == "CLT") ? 0.06 : 0.05) * 2.53e-11 #LiTaO3 N*m^4/C^2  it is g1 in Gopalan's paper
	end

    eps0 = 8.8541878128e-12
    C66 = (C11-C12)/2.

	  for a in 1:3, b in 1:3
	    for c in 1:3, d in 1:3
			G4o[a,b,c,d] =  2.0*G11*(δ2(a,b)*δ2(c,d) + δ2(a,c)*δ2(b,d) + δ2(a,d)*δ2(b,c) )/3.0
	    end
	end

	
      
      C4o[1,1,1,1] = C11;  C4o[2,2,2,2] = C11
      C4o[3,3,3,3] = C33
      C4o[2,3,2,3] = C44;  C4o[2,3,3,2] = C44;    C4o[3,2,3,2] = C44;    C4o[3,2,2,3] = C44
      C4o[1,3,1,3] = C44;  C4o[1,3,3,1] = C44;    C4o[3,1,3,1] = C44;    C4o[3,1,1,3] = C44
      C4o[1,2,1,2] = C66;  C4o[1,2,2,1] = C66;    C4o[2,1,2,1] = C66;    C4o[2,1,1,2] = C66

      C4o[1,1,2,2] = C12;  C4o[2,2,1,1] = C12
      C4o[1,1,3,3] = C13;  C4o[3,3,1,1] = C13
      C4o[1,1,2,3] = C14;  C4o[1,1,3,2] = C14;    C4o[3,2,1,1] = C14;    C4o[2,3,1,1] = C14

      C4o[2,2,3,3] = C13;  C4o[3,3,2,2] = C13
      C4o[2,2,2,3] =-C14;  C4o[2,2,3,2] =-C14;    C4o[3,2,2,2] =-C14;    C4o[2,3,2,2] =-C14

      C4o[1,2,1,3] = C14;  C4o[1,2,3,1] = C14;    C4o[2,1,1,3] = C14;    C4o[2,1,3,1] = C14
      C4o[1,3,1,2] = C14;  C4o[1,3,2,1] = C14;    C4o[3,1,1,2] = C14;    C4o[3,1,2,1] = C14

      C2o = [ C11 C12  C13 C14  0    0 ; # defining C_αβ
              C12 C11  C13 -C14 0    0 ;
              C13 C13  C33 0    0    0 ;
              C14 -C14 0   C44  0    0 ;
              0   0    0   0    C44 C14 ;
              0   0    0   0    C14 C66 ]
      
      C2oi = inv(C2o)
                        
      
      e = zeros(3,6)
      e = [ 0 0 0 0 e15 -e22;
            -e22 e22 0 e15 0 0;
            e31 e31 e33 0 0 0 ]
      
      d = zeros(3,6)
      for i in 1:3, j in 1:6, m in 1:6
            d[i,j] += e[i,m]*C2oi[m,j]
      end
      
      d15 = (d[1,5] + d[2,4])/2.
      d22 = (-d[1,6] -d[2,1] + d[2,2])/4.
      d31 = (d[3,1] + d[3,2])/2.
      d33 = d[3,3]
      
      Q13 = d31/(Pvac*eps0*(eps33-1.))
      Q24 = d22/(Pvac*eps0*(eps22-1.))
      Q33 = d33/(2. *Pvac*eps0*(eps33-1.))
      Q44 = d15/(2. *Pvac*eps0*(eps22-1.))
      
      g1 = C11*Q13 + C12*Q13 + C13*Q33
      g2 = 2. *C13*Q13 + C33*Q33
      g3 = -C11*Q24 + C12*Q24 + C14*Q44
      g4 = -2. *C14*Q24 + C44*Q44
      
      
      q4o[1,1,3,3] = g1
      q4o[2,2,3,3] = g1
      q4o[3,3,3,3] = g2
      
      q4o[1,1,2,3] = g3/2.;	q4o[1,1,3,2] = g3/2.
      q4o[2,2,2,3] =-g3/2.;	q4o[2,2,3,2] =-g3/2.
      q4o[2,3,2,3] = g4/2.;	q4o[2,3,3,2] = g4/2.;	q4o[3,2,2,3] = g4/2.;	q4o[3,2,3,2] = g4/2.	
      q4o[1,3,1,3] = g4/2.;	q4o[1,3,3,1] = g4/2.;	q4o[3,1,1,3] = g4/2.;	q4o[3,1,3,1] = g4/2.
      q4o[1,2,1,3] = g3/2.;	q4o[1,2,3,1] = g3/2.;	q4o[2,1,1,3] = g3/2.;	q4o[2,1,3,1] = g3/2.;

      
	  alpha21 = 1. /(2. * eps0* (eps33 -1.))
      alpha41 = alpha21/Pvac^2+ 2. *(2. *(C11 + C12)*Q13^2 + 4. *C13*Q13*Q33 + C33*Q33^2) #Gopalan's result, whcih coincides with ours alpha[1]/Pvac^2 - 4. *(beta1*psi2^2 + 4. * beta2*psi1^2 + 2. *beta4*psi1*psi2 + 2. *g1*psi1 + g2*psi2)
	  if material == "CLN"
		alpha42 = 1.7302e11 + 4. *(2. *C11*Q24^2 - 2. *C12*Q24^2 + Q44*(-4. *C14*Q24 + C44*Q44)) #addition is 1.11e9
      	alpha44 = 8300.0e9
	  elseif material == "SLN"
		alpha42 = 9.8771e11 + 4. *(2. *C11*Q24^2 - 2. *C12*Q24^2 + Q44*(-4. *C14*Q24 + C44*Q44))
	  	alpha44 = 271.0e12
	  elseif material == "CLT" 
		alpha42 = 2.4469578654e11 + 4. *(2. *C11*Q24^2 - 2. *C12*Q24^2 + Q44*(-4. *C14*Q24 + C44*Q44)) #addition is 3.2e9
      	alpha44 = 14000.0e9
	  elseif material == "SLT"
		alpha42 = 3.2564828067598423e12  + 4. *(2. *C11*Q24^2 - 2. *C12*Q24^2 + Q44*(-4. *C14*Q24 + C44*Q44))
	  	alpha44 = 2445.0e12
	  end

	  alpha22 = 1. /(eps0* (eps11 -1.)) - alpha21*alpha42/alpha41

      A2o[1,1] = alpha22/2.
      A2o[2,2] = alpha22/2. 
      A2o[3,3] = -alpha21/2.
      A4o[3,3,3,3] = alpha41/4.
      A4o[1,1,3,3] = (alpha42/2.)/6.;   A4o[1,3,1,3] = (alpha42/2.)/6.;   A4o[1,3,3,1] = (alpha42/2.)/6.;   A4o[3,1,1,3] = (alpha42/2.)/6.;   A4o[3,1,3,1] = (alpha42/2.)/6.;  A4o[3,3,1,1] = (alpha42/2.)/6.
      A4o[2,2,3,3] = (alpha42/2.)/6.;   A4o[2,3,2,3] = (alpha42/2.)/6.;   A4o[2,3,3,2] = (alpha42/2.)/6.;   A4o[3,2,2,3] = (alpha42/2.)/6.;   A4o[3,2,3,2] = (alpha42/2.)/6.;  A4o[3,3,2,2] = (alpha42/2.)/6.
      A4o[1,1,1,1] = alpha44/4.
      A4o[2,2,2,2] = alpha44/4.
      A4o[1,1,2,2] = 2. * (alpha44/4.)/6.;  A4o[1,2,1,2] = 2. * (alpha44/4.)/6.;    A4o[1,2,2,1] = 2. * (alpha44/4.)/6.;    A4o[2,1,1,2] = 2. * (alpha44/4.)/6.;    A4o[2,1,2,1] = 2. * (alpha44/4.)/6.;    A4o[2,2,1,1] = 2. * (alpha44/4.)/6.
            
	  g3 = 2. *G11 #factor 2 comes because in Chris paper in free energy 1/2
      g5 = 0. *g3
      g8 = 1. * g3
	  g10 = 1. * g3

	  G4o[1,3,1,3] = g3 
      G4o[2,3,2,3] = g3
	  G4o[1, 1, 2, 3] = g5/2; G4o[2, 3, 1, 1] = g5/2; G4o[2, 2, 2, 3] = -g5/2; G4o[2, 3, 2, 2] = -g5/2; G4o[1, 2, 1, 3] = g5/2; G4o[1, 3, 1, 2] = g5/2; G4o[2, 1, 1, 3] = g5/2; G4o[1, 3, 2, 1] = g5/2;

	  G4o[1, 1, 1, 1] = g8; G4o[2, 2, 2, 2] = g8; G4o[1, 1, 2, 2] = g8; G4o[2, 2, 1, 1] = g8;
	  G4o[2, 1, 2, 1] = g10; G4o[1, 2, 1, 2] = g10; G4o[1, 2, 2, 1] = -g10; G4o[2, 1, 1, 2] = -g10;      

	for a in 1:6, b in 1:3, c in 1:3, d in 1:3, e in 1:3
	    q3o[a,b,c] += TT[a,d,e]*q4o[d,e,b,c]	# defining q_αbc using Voigt contraction Eq. (4)
	end

	for a in 1:6, b in 1:3, c in 1:3, d in 1:6
	    Q3o[a,b,c] += C2oi[a,d]*q3o[d,b,c]		# defining Q_αbc Eq. (6)
	end

	for e in 1:3, f in 1:3, g in 1:3, h in 1:3
	    A4ot[e,f,g,h] = A4o[e,f,g,h]
	    for α in 1:6, β in 1:6 
	        A4ot[e,f,g,h] -= 0.5*C2o[α,β]*Q3o[α,e,f]*Q3o[β,g,h]	# A(tilda)_ijkl with corrections due to Eq. (7)
	    end
	end

	

    #  f(p) = Vo([0,0,p],A2o,A4ot)
    #  D(f) = p -> ForwardDiff.derivative(f,float(p))
	#p0 = find_zero( D(f), (0.1,1.5) ) #it coincides with Pvac
	p0 = Pvac
    PV = zeros(3)
    PV .= [0,0,p0]
	
	println("The ground state is PV = ", PV)

	#renormalizing (rescaling) parameters
	A2o ./= abs(g1) #g1 = q[1,1,3,3]
	A4o ./= abs(g1)/p0^2
	q4o ./= abs(g1)
	G4o  /= G11
	C4o ./= abs(g1)*p0^2	
	PV ./= p0 # it will be [0,0,1]
	length_scale = sqrt(G11/abs(g1))


	A2 = zeros(3,3)
	
	PV = R2*PV
	 
	for a in 1:3, b in 1:3, c in 1:3, d in 1:3
	    A2[a,b] += R2[a,c]*R2[b,d]*A2o[c,d]
	    for e in 1:3, f in 1:3, g in 1:3, h in 1:3
	        A4[a,b,c,d] += R2[a,e]*R2[b,f]*R2[c,g]*R2[d,h]*A4o[e,f,g,h]
	        C4[a,b,c,d] += R2[a,e]*R2[b,f]*R2[c,g]*R2[d,h]*C4o[e,f,g,h]
			G4[a,b,c,d] += R2[a,e]*R2[b,f]*R2[c,g]*R2[d,h]*G4o[e,f,g,h]
	        q4[a,b,c,d] += R2[a,e]*R2[b,f]*R2[c,g]*R2[d,h]*q4o[e,f,g,h]
	    end
	end
	
	
	

	for a in 1:6, b in 1:3, c in 1:3, d in 1:3, e in 1:3
	    q3[a,b,c] += TT[a,d,e]*q4[d,e,b,c]
	    Q3[a,b,c] += TT[a,d,e]*Q4[d,e,b,c]
	end

	for a in 1:6, b in 1:6, c in 1:3, d in 1:3, e in 1:3, f in 1:3
	    C2[a,b] += TT[a,c,d]*TT[b,e,f]*C4[c,d,e,f]
	end



	for α in 1:6, β in 1:6, b in 1:3, c in 1:3
	    Q3[α,b,c] += inv(C2)[α,β]*q3[β,b,c]
	end

	qll = zeros(6,3,3)
	qperp = zeros(6,3,3)
	
	Qll = zeros(6,3,3)
	Qperp = zeros(6,3,3)
	
	Cll = zeros(6,6)
	Cperp = zeros(6,6)
	Cm = zeros(6,6)

	qperp[2,:,:] .= q3[2,:,:]
	qperp[3,:,:] .= q3[3,:,:]
	qperp[4,:,:] .= q3[4,:,:]

	qll[1,:,:] .= q3[1,:,:]
	qll[5,:,:] .= q3[5,:,:]
	qll[6,:,:] .= q3[6,:,:]

	Qperp[2,:,:] .= Q3[2,:,:]
	Qperp[3,:,:] .= Q3[3,:,:]
	Qperp[4,:,:] .= Q3[4,:,:]

	Qll[1,:,:] .= Q3[1,:,:]
	Qll[5,:,:] .= Q3[5,:,:]
	Qll[6,:,:] .= Q3[6,:,:]

	for a in 1:6, b in 1:6
    
	    if a in [1,5,6] && b in [1,5,6]
	        Cll[a,b] = C2[a,b]
	    elseif a in [2,3,4] && b in [2,3,4]
	        Cperp[a,b] = C2[a,b]
	    else
	        Cm[a,b] = C2[a,b]
	    end
		
	end

	Cll33 = [ Cll[1,1] Cll[1,5] Cll[1,6] ;
		      Cll[5,1] Cll[5,5] Cll[5,6] ;
		      Cll[6,1] Cll[6,5] Cll[6,6] ]

	Cll33inv = inv(Cll33)

	Cll66inv=[ Cll33inv[1,1] 0 0 0 Cll33inv[1,2] Cll33inv[1,3] ;
				0 0 0 0 0 0 ;
				0 0 0 0 0 0 ;
				0 0 0 0 0 0 ;
			    Cll33inv[2,1] 0 0 0 Cll33inv[2,2] Cll33inv[2,3] ;
			    Cll33inv[3,1]  0 0 0 Cll33inv[3,2] Cll33inv[3,3] ]

	ϵV = zeros(6)

	for α in 1:6, b in 1:3, c in 1:3
	    ϵV[α] += Q3[α,b,c]*PV[b]*PV[c]
	end

	H = zeros(6,3,3)
	K = zeros(6)
	
	

	for α in [1,5,6], b in 1:3, c in 1:3
	    H[α,b,c] = Qll[α,b,c]
	    for β in 1:6, γ in 1:6
	        H[α,b,c] += Cll66inv[α,β]*Cm[β,γ]*Qperp[γ,b,c]
	    end
	end


	for α in [1,5,6], β in 1:6, γ in 1:6
	    K[α] += -Cll66inv[α,β]*Cm[β,γ]*ϵV[γ]
	end
	
	for α in [2,3,4]
	    K[α] = ϵV[α]
	end
		
	for a in 1:3, b in 1:3, c in 1:3, d in 1:3
		A4t[a,b,c,d] = A4[a,b,c,d]
		for alp in 1:6
			A4t[a,b,c,d] -= 1.0*q3[alp,a,b]*H[alp,c,d]
			for bet in 1:6
				A4t[a,b,c,d] += 0.5*H[bet,a,b]*C2[bet,alp]*H[alp,c,d]
			end
		end
	end
	

	
	for a in 1:3, b in 1:3
		A2t[a,b] = A2[a,b]
		for alp in 1:6
			A2t[a,b] -= q3[alp,a,b]*K[alp]
			for bet in 1:6
				A2t[a,b] += 0.5*(C2[alp,bet]+C2[bet,alp])*H[bet,a,b]*K[alp]
			end
		end
	end
	
	A0t = 0
	for alp in 1:6, bet in 1:6
		A0t += 0.5*K[alp]*C2[alp,bet]*K[bet]
	end
	
	
	V0 = 0.0

	if is_electrostatic
		println("ELECTROSTATIC ENERGY IS ON!!!! ")
	end

	
	return Parameter(PV, V0, A2, A4, G4, A2t, A4t, getA2s(A2t), getA4s(A4t), abs(g1)#=*p0^2=#, length_scale, is_electrostatic)
	
end









