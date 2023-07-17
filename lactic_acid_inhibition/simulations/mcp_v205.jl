using JuMP
using LinearAlgebra
using Ipopt
using DelimitedFiles
using TickTock
using JLD2
import DataFrames
import CSV
import ForwardDiff
import Dierckx

include("/home/users/mgotsmy/julia/dFBA/custom/postprocessing_02.jl")

function fed_batch_dFBA(N0,V0,nCI,S,t_max,t_min,V_max,c_G,glu,atp,so4,xxx,dna,qmin,qmax,Kup,vlb,vub,d,verbose=4)
    
# number of reactions
    nR = size(S,2)
    # number of metabolites
    nM = size(S,1)
    # number of tracked metabolites 
    nN = length(N0)
    # number collocation points
    nCP = 3
    
    #--------------------------
    # SIMULATION HYPERPARAMETERS
    
    w=1e-20         # weigth for flux sum minimization on the OF
    phi_scale = 1e0 # 
    phi1=phi_scale
    phi2=phi_scale
    phi3=phi_scale

    #--------------------------
    # COLLOCATION AND RADAU PARAMETERS
    colmat = [0.19681547722366  -0.06553542585020 0.02377097434822;
              0.39442431473909  0.29207341166523 -0.04154875212600;
              0.37640306270047  0.51248582618842 0.11111111111111]
    radau  = [0.15505 0.64495 1.00000]

    #--------------------------
    # JuMP MODEL SETUP
    m = Model(optimizer_with_attributes(Ipopt.Optimizer, 
            "warm_start_init_point" => "yes", 
            "print_level" => verbose, 
            "linear_solver" => "ma27", 
            "tol" => 1e-4, 
            "acceptable_iter" => 5, 
            "acceptable_tol" => 1e-2))

    ATPM_base_value = 3.15
    kATPM = 196.15209181
    pH = 7
    pKa = 3.86
    
    #--------------------------
    # VARIABLE SET UP
    @variables(m, begin
        # Differential Equation Variables
        N[   1:nN, 1:nCI, 1:nCP]  # metabolite amounts [mmol]
        Ndot[1:nN, 1:nCI, 1:nCP]  # dN/dt [mmol/h]
        V[         1:nCI, 1:nCP]  # reactor volume [L]
        Vdot[      1:nCI, 1:nCP]  # dV/dt [L/h]
        q[   1:nR, 1:nCI       ]  # FBA fluxes [mmol/(g h)]

        # KKT Variables
        lambda[ 1:nM, 1:nCI]
        alpha_U[1:nR, 1:nCI]
        alpha_L[1:nR, 1:nCI]
        FO_U[1:nR, 1:nCI]
        FO_L[1:nR, 1:nCI]

        # Process Variabels
        tend        # h
        # length of control intervals
        rtCI[1:nCI] # unitless
        # feed
        F[  1:nCI]  # [L/h]
        # initial sulfate amount
        S0          # [mmol]
    end)

    println("Setting start values.")
    @JLD2.load "/home/users/mgotsmy/julia/230502_UoT_stuff/simulations/mcp_v200/variables.jld2" N_ Ndot_ V_ Vdot_ q_ lambda_ alpha_U_ alpha_L_ FO_U_ FO_L_ tCI_ F_ S0_ rtCI_ tend_
    
    for _N in 1:nN
        for _CI in 1:nCI
            for _CP in 1:nCP
                set_start_value(N[_N,_CI,_CP],N_[_N,_CI,_CP])
                set_start_value(Ndot[_N,_CI,_CP],Ndot_[_N,_CI,_CP])
                set_start_value(V[_CI,_CP],V_[_CI,_CP])
                set_start_value(Vdot[_CI,_CP],Vdot_[_CI,_CP])
            end
        end
    end
    
    for _CI in 1:nCI
        for _R in 1:nR
            set_start_value(q[_R,_CI],q_[_R,_CI])
            set_start_value(alpha_L[_R,_CI],alpha_L_[_R,_CI])
            set_start_value(alpha_U[_R,_CI],alpha_U_[_R,_CI])
            set_start_value(FO_L[_R,_CI],FO_L_[_R,_CI])
            set_start_value(FO_U[_R,_CI],FO_U_[_R,_CI])
        end
        set_start_value(rtCI[_CI],rtCI_[_CI])
        set_start_value(F[_CI],F_[_CI])
    end

    for _CI in 1:nCI
        for _M in 1:nM
            set_start_value(lambda[_M,_CI],lambda_[_M,_CI])
        end
    end
    
    set_start_value(S0,S0_)
    set_start_value(tend,tend_)
    println("Finished setting start values.")

    #--------------------------
    # SET UP OBJECTIVE FUNCTION
    
    # Uptake Rates as NLexpressions
    @NLexpressions(m, begin
        # glucose
        q_G[i=1:nCI], F[i]*c_G/N[2,i,2] # L/h * mmol/L / g = mmol/(g h)
        # absolute values for tCI
        tCI[i=1:nCI], tend/nCI*rtCI[i]
        # ATPM dependence on acid toxicity
        q_M[i=1:nCI], ATPM_base_value + kATPM * (N[4,i,2]/V[i,2])/(1 + 10^(pH-pKa))
    end)
    
    @NLobjective(m, Max, +
        N[4,end,end]/sum(tCI[i] for i in 1:nCI) * 1e1 -
        sum(
            sum(
                - phi1*FO_L[mc,i] 
                - phi3*FO_U[mc,i] for mc in 1:nR) 
            for i in 1:nCI)/nCI)
    
    @NLconstraints(m, begin
        NLc_q_G[ i=1:nCI], -q[glu,i] - q_G[i] == 0 # |q_G| >= |q[glu]| [mmol/(g h)]
        NLC_q_G2[i=1:nCI], qmin + (qmax - qmin)*(-1 + 2/(1 + exp(- Kup * q[xxx,i]))) - q_G[i] == 0
            
        FO1[mc=1:nR,i=1:nCI], FO_L[mc,i]  == (q[mc,i] -vlb[mc])*alpha_L[mc,i]
        FO2[mc=1:nR,i=1:nCI], FO_U[mc,i]  == (q[mc,i] -vub[mc])*alpha_U[mc,i]
            
        NLc_q_M[i=1:nCI], q_M[i] - q[atp,i] == 0
            
        # INTEGRATION BY COLLOCATION
        # set up collocation equations - 2nd-to-nth point
        coll_N[l=1:nN, i=2:nCI, j=1:nCP], N[l,i,j] == N[l,i-1,nCP]+tCI[i]*sum(colmat[j,k]*Ndot[l,i,k] for k in 1:nCP)
        coll_V[        i=2:nCI, j=1:nCP], V[i,j]   == V[  i-1,nCP]+tCI[i]*sum(colmat[j,k]*Vdot[  i,k] for k in 1:nCP)
        # set up collocation equations - 1st point
        coll_N0[l in [1,2,4], i=[1], j=1:nCP], N[l,i,j] == N0[l] + tCI[i]*sum(colmat[j,k]*Ndot[l,i,k] for k in 1:nCP)
        coll_S0[l in [3],     i=[1], j=1:nCP], N[l,i,j] == S0    + tCI[i]*sum(colmat[j,k]*Ndot[l,i,k] for k in 1:nCP)
        coll_V0[              i=[1], j=1:nCP], V[  i,j] == V0    + tCI[i]*sum(colmat[j,k]*Vdot[  i,k] for k in 1:nCP)
    end)

    #------------------------#
    # SET UP BOUNDS
    for mc in 1:nR
        for i in 1:nCI
            set_lower_bound(q[mc,i],vlb[mc])
            set_upper_bound(q[mc,i],vub[mc])
            set_lower_bound(alpha_U[mc,i],0)
            set_upper_bound(alpha_L[mc,i],0)
            for j in 1:nCP
                for n in 1:nN
                    set_lower_bound(N[n,i,j],0)
                end
                set_lower_bound(V[  i,j],V0)
                set_upper_bound(V[  i,j],V_max)
            end
        end
    end
    set_upper_bound(S0,30)
    
    for i in 1:nCI
        set_lower_bound(rtCI[i],0.8)
        set_upper_bound(rtCI[i],1.2)
    end
    set_lower_bound(tend,t_min)
    set_upper_bound(tend,t_max)
    
    #------------------------#
    # SET UP OTHER CONSTRAINTS
    
    @constraints(m, begin
        # DIFFERENTIAL EQUATIONS
        # glucose
        m1[i=1:nCI, j=1:nCP], Ndot[1,i,j] == F[i]*c_G + q[glu,i]*N[2,i,j] # mmol/h = L/h * mmol/L + mmol/(g h) * g
        # biomass
        m2[i=1:nCI, j=1:nCP], Ndot[2,i,j] == q[xxx,i]*N[2,i,j] # g/h = g/(g h) * g
        # sulfate
        m3[i=1:nCI, j=1:nCP], Ndot[3,i,j] == q[so4,i]*N[2,i,j]
        # product
        m4[i=1:nCI, j=1:nCP], Ndot[4,i,j] == q[dna,i]*N[2,i,j]
        # volume
        v1[i=1:nCI, j=1:nCP], Vdot[i,j]   == F[i]

        # SYSTEM CONSTRAINTS
        c_S[ mc=1:nM,i=1:nCI], sum(S[mc,k] * q[k,i] for k in 1:nR) == 0
        #c_UB[mc=1:nR,i=1:nCI],  q[mc,i]  - vub[mc] <= 0 # implemented as set_lower_bound
        #c_LB[mc=1:nR,i=1:nCI], -q[mc,i]  + vlb[mc] <= 0 # implemented as set_lower_bound

        c_UB_S0                         , S0          == 30

        #--------------------------
        # KKT
        Lagr[mc=1:nR,i=1:nCI],  d[mc] + 
                                w*q[mc,i]  + 
                                alpha_L[mc,i] + 
                                alpha_U[mc,i] + 
                                sum(S[k,mc]*lambda[k,i] for k in 1:nM) == 0
                                

        alpha1_LB[mc=1:nR,i=1:nCI],   alpha_L[mc,i]    <= 0
        alpha1_UB[mc=1:nR,i=1:nCI],   alpha_U[mc,i]    >= 0    
            
        # here, for nfe = 6, the 2nd & 3rd fluxes are fixed to the first, and the 5th & 6th fluxes are fixed to the 4th
        #c_raj1[i=2:Int(nCI/2)],         q[xxx,i] - q[xxx,i-1] == 0
        #c_raj2[i=Int((nCI/2)+2):nCI],   q[xxx,i] - q[xxx,i-1] == 0
    end)

    #--------------------------
    # MODEL OPTIMIZATION
    println("Model Preprocessing Finished.")
    tick()
    solveNLP = JuMP.optimize!
    status = solveNLP(m)
    t = tok()
    #tock()
    println("Model Finished Optimization")

   # Get values for plotting
    N_ = JuMP.value.(N[:,:,:])
    Ndot_ = JuMP.value.(Ndot[:,:,:])
    V_ = JuMP.value.(V[:,:])
    Vdot_ = JuMP.value.(Vdot[:,:])
    q_ = JuMP.value.(q[:,:])
    lambda_ = JuMP.value.(lambda[:,:])
    alpha_U_ = JuMP.value.(alpha_U[:,:])
    alpha_L_ = JuMP.value.(alpha_L[:,:])
    FO_U_ = JuMP.value.(FO_U[:,:])
    FO_L_ = JuMP.value.(FO_L[:,:])
    tCI_ = JuMP.value.(tCI[:])
    F_ = JuMP.value.(F[:])
    S0_ = JuMP.value.(S0)
    rtCI_ = JuMP.value.(rtCI[:])
    tend_ = JuMP.value.(tend)
    
    @JLD2.save dirname*"/variables.jld2" N_ Ndot_ V_ Vdot_ q_ lambda_ alpha_U_ alpha_L_ FO_U_ FO_L_ tCI_ F_ S0_ rtCI_ tend_

return  N_, Ndot_, V_, Vdot_, q_, tCI_, FO_L_, FO_U_, m, t
end

function run_kkt_simulation()

    # load stoichiometric matrix & flux bounds
    S = readdlm("../Si.csv", ',')
    vlb2 = readdlm("../lbi.csv", ',')
    vlb=vlb2[:,1]
    vub2 = readdlm("../ubi.csv", ',')
    vub=vub2[:,1]
    
    # reaction indices
    glu = 12
    atp = 716
    so4 = 260
    xxx = 19
    dna = 92 # actually EX_lac__D_e

    # open d-lactic acid exchange
    vlb[dna] = 0
    vub[dna] = 1000

    # open glucose
    vlb[glu] = -10
    vub[glu] = 0

    # open biomass
    vlb[xxx] = 0
    vub[xxx] = 1000
    
    # initial values
    V0 = .5   # L
    G0 = 0.   # mmol
    X0 = 1 # g
    S0 = 30 # 8.9  # mmol
    c_G = 450/0.18015588 # g/L /(g/mmol) = mmol/L # basically max solubility

    N0 = [G0,X0,S0,0]

    t_max = 20 # h
    t_min = 1  # h
    V_max = 1  # L
    nCI   = t_max ÷ 1

    # KKT objective function: optimize product
    nR = size(S,2)
    d = 0.0*Vector{Float64}(undef,nR) 
    d[dna] = -1
    
    # parameters for uptake - growth relationship:
    qmin = 0.5
    qmax = 10.0
    Kup  = 5

    N_, Ndot_, V_, Vdot_, q_, tCI_, FO_L_, FO_U_, m_, t_ = fed_batch_dFBA(N0,V0,nCI,S,t_max,t_min,V_max,c_G,glu,atp,so4,xxx,dna,qmin,qmax,Kup,vlb,vub,d,5)
    println("DONE")
    println(termination_status(m_))
      
    #-------------------------
    # SAVING ALL INTERESTING FILES
    
    summary = "#----" * 
    "Termination Status" * "\n" * 
    string(termination_status(m_)) * "\n" * 
    "\nObjective Value" * "\n" * 
    string(JuMP.objective_value(m_)) * "\n" * 
    "\nNormalized Integrated Biomass" * "\n" * 
    string(sum(N_[2,:,:]/(6*3))) * "\n" * 
    "\nSum of growth rates" * "\n" * 
    string(sum(q_[xxx,:])) * "\n" * 
    "\nFinal pDNA amount" * "\n" * 
    string(N_[4,end,end]) * "\n" * 
    "\nNormalized Complementary Slackness" * "\n" * 
    string(-(sum(-FO_L_)+sum(-FO_U_))/nCI) * "\n" * 
    "\nComplementary Slackness per CTRL interval" * "\n" * 
    string([round(sum(FO_L_[:,i])+sum(FO_U_[:,i]),digits=2) for i in 1:nCI]) * "\n" * 
    "\nS0" * "\n" * 
    string(value(JuMP.variable_by_name(m_,"S0"))) * "\n" * 
    ("\nControl Intervals") * "\n" * 
    string([round(i,digits=2) for i in tCI_]) * "\n" * 
    ("\nProcess End") * "\n" * 
    string(round(sum(tCI_),digits=2)) * "\n" * 
    ("\nFeed Rates") * "\n" * 
    string([round(i,digits=3) for i in Vdot_[:,1]]) * "\n" *
    "\nTime Elapsed" * "\n" * string(t_)
    println(summary)
    
    # POST PROCESSING

    df = DataFrames.DataFrame(hcat(
            get_time_points(tCI_),
            get_points(V_,tCI_,V0),
            get_points(N_,tCI_,[N0[1],N0[2],value(JuMP.variable_by_name(m_,"S0")),N0[4]]),
            get_points(Vdot_,tCI_),
            get_points(Ndot_,tCI_),
            get_fluxes(q_,[glu,xxx,so4,dna,atp],tCI_)),
        ["t","V","G","X","S","P","r_V","r_G","r_X","r_S","r_P","q_G","q_X","q_S","q_P","q_M"]);
    
    write(   dirname*"/summary.txt",summary)
    writedlm(dirname*"/q.csv",  q_, ',')
    writedlm(dirname*"/vlb.csv",vlb,',')
    writedlm(dirname*"/vub.csv",vub,',')
    CSV.write(dirname*"/df.csv",df)
    
    
    
end

println("Start Script")
scriptname = PROGRAM_FILE
dirname = scriptname[begin:end-3]
println("Creating results directory: ",dirname)

mkdir(dirname)

debug = false
if debug == true
    run_kkt_simulation()
else
    global worked = false
    global nTRY = 0
    while worked == false
        try
            run_kkt_simulation()
            global worked = true
        catch e
            # print(e)
            if isa(e,ErrorException)
                global nTRY += 1
                println("Error During Initialization. Retrying ... (",nTRY,")")
            else
                rethrow(e)
            end
        end
    end
end

println("Script Ended")
