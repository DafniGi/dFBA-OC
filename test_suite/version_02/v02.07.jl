using JuMP
using LinearAlgebra
using Ipopt
using DelimitedFiles
using JLD2
import DataFrames
import CSV
import ForwardDiff
import Dierckx

include("/home/users/mgotsmy/julia/230502_UoT_stuff/230601_CP_FBA/postprocessing_03.jl")

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
    # PROCESS PARAMETERS
    
    ATPM_base_value = 3.15
    kATPM = 196.15209181
    pH = 7
    pKa = 3.86
    
    #--------------------------
    # SIMULATION HYPERPARAMETERS
    
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
            "acceptable_iter" => 15, 
            "acceptable_tol" => 1e-2))
    
    #--------------------------
    # SET UP VARIABLES
    @variables(m, begin
        # Differential Equation Variables
        N[   1:nN, 1:nCI, 1:nCP]  # metabolite amounts [mmol]
        Ndot[1:nN, 1:nCI, 1:nCP]  # dN/dt [mmol/h]
        V[         1:nCI, 1:nCP]  # reactor volume [L]
        Vdot[      1:nCI, 1:nCP]  # dV/dt [L/h]
        q[   1:nR, 1:nCI, 1:nCP]  # FBA fluxes [mmol/(g h)]

        # KKT Variables
        lambda[ 1:nM, 1:nCI, 1:nCP]
        alpha_U[1:nR, 1:nCI, 1:nCP]
        alpha_L[1:nR, 1:nCI, 1:nCP]
        FO_U[   1:nR, 1:nCI, 1:nCP]
        FO_L[   1:nR, 1:nCI, 1:nCP]

        # Process Variabels
        tend        # h
        # length of control intervals
        rtCI[1:nCI] # unitless
        # glucose uptake rate
        q_G[1:nCI]  # [L/h]
        # initial sulfate amount
        S0          # [mmol]
    end)

    #--------------------------
    # SET UP STARTING VALUES
    println("Setting start values.")
    @JLD2.load "/home/users/mgotsmy/julia/230502_UoT_stuff/230601_CP_FBA/v026/variables.jld2" N_ Ndot_ V_ Vdot_ q_ lambda_ alpha_U_ alpha_L_ FO_U_ FO_L_ tCI_ F_ S0_ tCI_
    nN_,nCI_,nCP_ = size(N_)
    
    for _N in 1:nN
        for _CI in 1:nCI
            for _CP in 1:nCP
                set_start_value(N[_N,_CI,_CP],N0[_N])
                set_start_value(V[_CI,_CP],V0)
            end
        end
    end
    set_start_value(S0,S0_)
    set_start_value(tend,t_max)
    println("Finished setting start values.")
    
    #--------------------------
    # SET UP TIME POINTS
    @NLexpressions(m,begin
        tCI[i=1:nCI], tend/nCI*rtCI[i]
        end)
    
    #--------------------------
    # SET UP OBJECTIVE FUNCTION
    @NLobjective(m, Max, +
        N[4,end,end]/sum(tCI[i] for i in 1:nCI) *  1.0  -
        sum(
            sum(
                sum(
                        - phi1*FO_L[mc,i,j] 
                        - phi3*FO_U[mc,i,j] 
                    for mc in 1:nR) 
                for j in 1:nCP)
            for i in 1:nCI)/(nCI*nCP))
    
    #--------------------------
    # SET UP PROCESS MODEL
    @expressions(m, begin
        # feed rate
        F[  i=1:nCI,j=1:nCP], q_G[i]/c_G*N[2,i,j]
        end)
    
    @NLexpressions(m, begin
        # change of ATPM due to toxicity
        q_M[i=1:nCI,j=1:nCP], ATPM_base_value + kATPM * (N[4,i,j]/V[i,j])/(1 + 10^(pH-pKa))
        # growth as a function of glucose uptake
        q_X[i=1:nCI,j=1:nCP], log((-q_G[i]-qmax+2qmin)/(q_G[i]-qmax))*1/Kup
    end)
    
    @NLconstraints(m, begin
        # fix ATPM in FBA
        NLc_q_M1[i=1:nCI,j=1:nCP],   q[atp,i,j] - q_M[i,j] == 0
        # fix q_G in FBA
        NLc_q_G1[i=1:nCI,j=1:nCP], - q[glu,i,j] - q_G[i] == 0 # |q_G| >= |q[glu]| [mmol/(g h)]
        # fix biomass as function of q_G in FBA
        NLc_q_G2[i=1:nCI,j=1:nCP],   q[xxx,i,j] - q_X[i] == 0
        end)
    @constraints(m, begin
        # initial sulfate concentration
        c_UB_S0, S0 == 30
        end)
    
    #--------------------------
    # SET UP DIFFERENTIAL EQUATIONS
    @constraints(m, begin
        # glucose
        m1[i=1:nCI, j=1:nCP], Ndot[1,i,j] == F[i,j]*c_G + q[glu,i,j]*N[2,i,j]
        # biomass
        m2[i=1:nCI, j=1:nCP], Ndot[2,i,j] == q[xxx,i,j]*N[2,i,j]
        # sulfate
        m3[i=1:nCI, j=1:nCP], Ndot[3,i,j] == q[so4,i,j]*N[2,i,j]
        # product
        m4[i=1:nCI, j=1:nCP], Ndot[4,i,j] == q[dna,i,j]*N[2,i,j]
        # volume
        v1[i=1:nCI, j=1:nCP], Vdot[i,j]   == F[i,j]
        end)
    
    
    #--------------------------
    # SET UP KKT
    @constraints(m, begin
        # complementary slackness
        FO1[mc=1:nR,i=1:nCI,j=1:nCP], FO_L[mc,i,j]  == (q[mc,i,j] -vlb[mc])*alpha_L[mc,i,j]
        FO2[mc=1:nR,i=1:nCI,j=1:nCP], FO_U[mc,i,j]  == (q[mc,i,j] -vub[mc])*alpha_U[mc,i,j]
        # FBA constraints 
        c_S[ mc=1:nM,i=1:nCI,j=1:nCP], sum(S[mc,k] * q[k,i,j] for k in 1:nR) == 0
        Lagr[mc=1:nR,i=1:nCI,j=1:nCP],  d[mc] + 
                                        alpha_L[mc,i,j] + 
                                        alpha_U[mc,i,j] + 
                                        sum(S[k,mc]*lambda[k,i,j] for k in 1:nM) == 0
        # KKT condition
        alpha1_LB[mc=1:nR,i=1:nCI,j=1:nCP],   alpha_L[mc,i,j]    <= 0
        alpha1_UB[mc=1:nR,i=1:nCI,j=1:nCP],   alpha_U[mc,i,j]    >= 0    
        end)
    
    
    #--------------------------
    # SET UP COLLOCATION
    @NLconstraints(m, begin
        # set up collocation equations - 2nd-to-nth point
        coll_N[l=1:nN, i=2:nCI, j=1:nCP], N[l,i,j] == N[l,i-1,nCP]+tCI[i]*sum(colmat[j,k]*Ndot[l,i,k] for k in 1:nCP)
        coll_V[        i=2:nCI, j=1:nCP], V[i,j]   == V[  i-1,nCP]+tCI[i]*sum(colmat[j,k]*Vdot[  i,k] for k in 1:nCP)
        # set up collocation equations - 1st point
        coll_N0[l in [1,2,4], i=[1], j=1:nCP], N[l,i,j] == N0[l] + tCI[i]*sum(colmat[j,k]*Ndot[l,i,k] for k in 1:nCP)
        coll_S0[l in [3],     i=[1], j=1:nCP], N[l,i,j] == S0    + tCI[i]*sum(colmat[j,k]*Ndot[l,i,k] for k in 1:nCP)
        coll_V0[              i=[1], j=1:nCP], V[  i,j] == V0    +       tCI[i]*sum(colmat[j,k]*Vdot[  i,k] for k in 1:nCP)
    end)

    #------------------------#
    # SET UP BOUNDS
    for mc in 1:nR
        for i in 1:nCI
            for j in 1:nCP
                set_lower_bound(q[mc,i,j],vlb[mc])
                set_upper_bound(q[mc,i,j],vub[mc])
                set_lower_bound(alpha_U[mc,i,j],0)
                set_upper_bound(alpha_L[mc,i,j],0)
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
        #c_raj1[i=2:Int(nCI/2)],         q[xxx,i] - q[xxx,i-1] == 0
        #c_raj2[i=Int((nCI/2)+2):nCI],   q[xxx,i] - q[xxx,i-1] == 0
    end)

    #--------------------------
    # MODEL OPTIMIZATION
    println("Model Preprocessing Finished.")
    solveNLP = JuMP.optimize!
    status = solveNLP(m)
    println("Model Finished Optimization")

   # Get values for plotting
    N_ = JuMP.value.(N[:,:,:])
    Ndot_ = JuMP.value.(Ndot[:,:,:])
    V_ = JuMP.value.(V[:,:])
    Vdot_ = JuMP.value.(Vdot[:,:])
    q_ = JuMP.value.(q[:,:,:])
    lambda_ = JuMP.value.(lambda[:,:,:])
    alpha_U_ = JuMP.value.(alpha_U[:,:,:])
    alpha_L_ = JuMP.value.(alpha_L[:,:,:])
    FO_U_ = JuMP.value.(FO_U[:,:,:])
    FO_L_ = JuMP.value.(FO_L[:,:,:])
    tCI_ = JuMP.value.(tCI[:])
    rtCI_ = JuMP.value.(rtCI[:])
    q_G_ = JuMP.value.(q_G[:])
    F_ = JuMP.value.(F[:,:])
    S0_ = JuMP.value.(S0)
    tend_ = JuMP.value.(tend)
    
    @JLD2.save dirname*"/variables.jld2" tCI_ tend_ N_ Ndot_ V_ Vdot_ q_ lambda_ alpha_U_ alpha_L_ FO_U_ FO_L_ rtCI_ F_ S0_ q_G_

return  N_, Ndot_, V_, Vdot_, q_, tCI_, FO_L_, FO_U_, m
end

function run_kkt_simulation()
    S = readdlm("/home/users/mgotsmy/julia/dFBA/DC_dFBA/Ecoli_iJO1366/Si.csv", ',');
    vlb2 = readdlm("/home/users/mgotsmy/julia/dFBA/DC_dFBA/Ecoli_iJO1366/lbi.csv", ',');
    vlb=vlb2[:,1]
    vub2 = readdlm("/home/users/mgotsmy/julia/dFBA/DC_dFBA/Ecoli_iJO1366/ubi.csv", ',');
    vub=vub2[:,1]
    
    
    # reaction indices
    glu = 12
    atp = 716 # lb is 6.86
    so4 = 260
    xxx = 19
    dna = 92 # actually EX_lac__D_e

    # set kappa to 200%
    vlb[dna] = 0
    vub[dna] = 1000

    # open glucose
    vlb[glu] = -10
    vub[glu] = 0

    # open biomass
    vlb[xxx] = 0
    vub[xxx] = 1000
    ;
    
    # initial values
    V0 = .5   # L
    G0 = 0.   # mmol
    X0 = 1 # g
    S0 = 30 # 8.9  # mmol
    c_G = 450/0.18015588 # g/L /(g/mmol) = mmol/L # basically max solubility

    N0 = [G0,X0,S0,0]

    t_max =  10 
    t_min = 1  # h
    V_max = 1  # L
    nCI   =  15 
    nCP   = 3

    nR = size(S,2)
    # these are some KKT vectors
    d=0.0*Vector{Float64}(undef,nR) 
    d[dna] = -1
    ;
    
    # relationship:
    qmin = 0.5
    qmax = 10.0
    Kup = 5

    t_ = @elapsed  N_, Ndot_, V_, Vdot_, q_, tCI_, FO_L_, FO_U_, m_ = fed_batch_dFBA(N0,V0,nCI,S,t_max,t_min,V_max,c_G,glu,atp,so4,xxx,dna,qmin,qmax,Kup,vlb,vub,d,5)
    println("DONE")
    println(termination_status(m_))
      
    #-------------------------
    # SAVING ALL INTERESTING FILES

    summary = "#----" * 
    "Termination Status" * "\n" * 
    string(termination_status(m_)) * "\n" * 
    "\nObjective Value" * "\n" * 
    string(JuMP.objective_value(m_)) * "\n" * 
    "\nSum of growth rates" * "\n" * 
    string(sum(q_[xxx,:,:])) * "\n" * 
    "\nFinal pDNA amount" * "\n" * 
    string(N_[4,end,end]) * "\n" * 
    "\nNormalized Complementary Slackness" * "\n" * 
    string(-(sum(-FO_L_)+sum(-FO_U_))/(nCI*nCP)) * "\n" * 
    "\nComplementary Slackness per CTRL interval" * "\n" * 
    string([[round(sum(FO_L_[:,i,j])+sum(FO_U_[:,i,j]),digits=2) for j in 1:nCP] for i in 1:nCI]) * "\n" * 
    "\nS0" * "\n" * 
    string(value(JuMP.variable_by_name(m_,"S0"))) * "\n" * 
    ("\nControl Intervals") * "\n" * 
    string([round(i,digits=2) for i in tCI_]) * "\n" * 
    ("\nProcess End") * "\n" * 
    string(round(sum(tCI_),digits=2)) * "\n" * 
    ("\nFeed Rates") * "\n" * 
    string([[round(Vdot_[i,j],digits=3) for j in 1:nCP] for i in 1:nCI]) * "\n" *
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

mkpath(dirname)

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
