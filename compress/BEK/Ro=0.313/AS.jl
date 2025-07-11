##depend packages
include("CRD_STA.jl")
include("Fun.jl")
using NonlinearEigenproblems
using DelimitedFiles
function Caculate_AS(Mr)
    #parameters
    N_cheb = 49
    gamma = 1.4
    sigma = 0.72
    Ro = 0.313
    Co = 2 - Ro - Ro^2
    Tw = 1
        #baseflow
        global u0,v0,w0,f,q,D,D2,x = baseflow_var(N_cheb,Ro);
        global H,T = T_ca(Mr,f,q,w0,gamma,Tw);
        global F,G,H,T,rho,z = interp(u0,v0,H,T,x,N_cheb,"sim");
        global lam = - (2/3) * T;
        global kappa = (1/sigma) * T;
        global al0 = 0.1
        global data_all = [0 0 0 0]

        #outloop of R
        for R = 120 : 0.5 : 140
            global mode = 2
            #innerloop of beta
            # for be = 0.1 : -0.005 : -0.05
            for be = 0.25 : -0.005 : -0.05
                if mode == 1
                    data_temp = [R be -1 -1]
                    global data_all = [data_all;data_temp]
                    writedlm("Dataall_$Tw _$Mr.dat",data_all[2:end,:])
                    continue
                end

                Ma = Mr/R
                PinPoint= []
                total = []
                writedlm("AS.dat",total)
                B0,B1 = Timemode(F,G,H,rho,lam,kappa,T,sigma,gamma,R,Ma,al0,be,N_cheb,Ro,Co)
                global C = eigen(B0,B1)
                global val = filter(x->-0.01<imag(x)<0.01&&abs(real(x))<0.2,C.values)

                if val == ComplexF64[]3
                    data_temp = [R be -1 -1]
                    global data_all = [data_all;data_temp]
                    writedlm("Dataall_$Tw _$Mr.dat",data_all[2:end,:])
                    continue
                end 

                for i = 1 : min(2,length(val))
                    indi = []
                    val_temp = val[i]
                    for al = 0.1 : 0.006 : 0.7
                        vec = eigvector(val[i],C.values,C.vectors)
                        B0,B1 = Timemode(F,G,H,rho,lam,kappa,T,sigma,gamma,R,Ma,al,be,N_cheb,Ro,Co)
                        val0,vec0 = RQI(B0,B1,val_temp,q0=vec)
                        indi = [indi;val0]
                        val_temp = val0
                        vec = vec0 
                        open("AS.dat", "a") do io
                        write(io,"al = $al,eig = $val0\n")
                        end
                    end
                    if i == 1

                        total = indi

                    else

                        total = [total indi]

                    end
                end

                for i = 1 : length(axes(total,2))

                    d2 = diff1(real(total[:,i]),0.006)

                    for j = 1 : length(d2)-1

                        if d2[j] * d2[j+1] < 0 && abs(d2[j+1])<0.005

                            PinPoint = [PinPoint;total[j,i]]

                        end

                    end

                end

                if PinPoint == []
                data_temp = [R be -1 -1]
                else
                data_temp = [R be real(PinPoint[findmax(imag(PinPoint))[2],1]) imag(PinPoint[findmax(imag(PinPoint))[2],1])]
                end
                
                global data_all = [data_all;data_temp]
                if  length(axes(data_all,1)) > 5 && data_all[end,4] < 0 && data_all[end-2,4] < 0 && data_all[end-3,4] < 0 && data_all[end-4,4] > 0 
                    global mode = 1
                end
                #io
                writedlm("Dataall_$Tw _$Mr.dat",data_all[2:end,:])

            end
            
        end
    end
Caculate_AS(0.3)