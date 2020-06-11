import MPSPyLib as mps
import numpy as np
import sys
import os.path
import matplotlib.pyplot as plt
import time

def main(PostProcess=False):
    """

    """
    epsilon_e= 0.5
    epsilon_nuc = 0.5

    exchange = 0

    eigenstate = 0
    totalstate = eigenstate + 1

    specflag = 'S'

    t = 1.0
    L = 12
    N = 6
    eefactor = 1
    nfactor = 2

    testflag = 'Finite'
    max_sweep = 10
    int_range = L
    lambda_up = 20
    lambda_lo = 0
    itr = 21

    # Build operators
    Operators = mps.BuildFermiOperators()

    # Hopping terms
    H = mps.MPO(Operators)
    H.AddMPOTerm('bond', ['fdagger','f'], hparam='t', weight=1.0, Phase=True)


    
    # Adding the electron-electron interaction
    
    
    eeinvr = []

    for eerange in range(int_range):
        eeinvr.append( 1.0 /(eerange + 1.0 + epsilon_e))

    H.AddMPOTerm('FiniteFunction', ['nftotal', 'nftotal'], f=eeinvr, hparam='lambda_e', weight=1.0)  


    nucinvr = []

    for nucrange in range(int_range):
        nucinvr.append( 1.0 /(nucrange + 1.0 + epsilon_nuc )) 


    H.AddMPOTerm('FiniteFunction', ['nftotal', 'I'], f=nucinvr, hparam='lambda_nuc', weight= -1.0)  



    '''



    # Adding the nucleus-electron interaction
    
    # infinite function for the nuclear interaction

    eeinvr = lambda x : 1.0 /(x+epsilon_e)
    H.AddMPOTerm('InfiniteFunction', ['nftotal','nftotal'], func=eeinvr, hparam='lambda_e', weight=1.0, L=1000, tol=1e-12, maxnterms=100) 



    nucinvr = lambda x : 1.0 /(x+epsilon_nuc)
    H.AddMPOTerm('InfiniteFunction', ['nftotal','I'], func=nucinvr, hparam='lambda_nuc', weight= -1.0, L=1000, tol=1e-12, maxnterms=100) 
    
    # Alternative way of defining the term using Finite functions

    
    '''

    H.printMPO()

    # Observables
    myObservables = mps.Observables(Operators)
    # Correlation functions
    #myObservables.AddObservable('corr', ['fdagger','f'], 'bspdm')
    #myObservables.AddObservable('corr', ['fdagger','f'], 'spdm', Phase=True)

    # Convergence data
    myConv = mps.MPSConvParam(max_bond_dimension=500, max_num_sweeps=max_sweep, max_num_lanczos_iter=100, variance_tol=1e-12, method=specflag, warmup_bond_dimension=500)

    #lambda_list = np.logspace(lambda_lo, lambda_up, itr)
    #lambda_list=np.linspace(lambda_lo,lambda_up,itr)


    lambda_list = np.array([ 0.0, 0.03, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.25, 1.5, 1.75, 2, 2.5, 3, 4, 5, 6, 7,
     8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 30, 40, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000])


    if not eefactor:
        lambda_list = np.array([0, 0.03, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.25, 1.5, 1.75, 2, 10, 20, 30, 40, 50, 100, 200, 500, 1000 ])

    if not nfactor:
        lambda_list = np.array([0, 0.03, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.25, 1.5, 1.75, 2, 2.5, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 30, 40, 50, 100])


    parameters= []

    start_time = time.time()

    for lambda_e in lambda_list: 


    # Define statics
        lambda_nuc = lambda_e

        ID = 'var' + testflag + 'L_' + str(L) + 'N' + str(N) + 'int' + str(int_range) +  'lambda' + str(lambda_nuc) + 'exc' + str(exchange) + 'eig' + \
         str(eigenstate) + 'factor' + str(eefactor) + '_' + str(nfactor)  + specflag


        dirID = 'var' + testflag + 'L_' + str(L) + 'N' + str(N) + 'int' + str(int_range)  + 'exc' + str(exchange) + 'eig' + \
         str(eigenstate) + 'factor' + str(eefactor) + '_' + str(nfactor) + specflag +  '/'

        parameters.append({ 
        # Directories
            'simtype'                   : testflag,
            'job_ID'                    : 'GWB_ET',
            'unique_ID'                 : ID,
            'Write_Directory'           : dirID,
            'Output_Directory'          : 'OUTPUTS' + dirID,
        # System size and Hamiltonian parameters
            'L'                         : L,
            't'                         : t, 
            'lambda_e'                  : eefactor * ( 1 - exchange ) * lambda_e, 
            'lambda_nuc'                : nfactor * lambda_nuc, 
            'n_excited_states'          : eigenstate,
            'verbose'                   : 2,
            'logfile'                   : True,
        # Specification of symmetries and good quantum numbers
            'Abelian_generators'        : ['nftotal'],
            'Abelian_quantum_numbers'   : [N],
            'MPSObservables'            : myObservables,
            'MPSConvergenceParameters'  : myConv
        })

    # Write Fortran-readable main files
    MainFiles = mps.WriteFiles(parameters, Operators, H,
                               PostProcess=PostProcess)

    # Run the simulations
    if(not PostProcess):
        if os.path.isfile('./Execute_MPSMain'):
            RunDir = './'
        else:
            RunDir = None
        mps.runMPS(MainFiles, RunDir=RunDir)
        return

    # Postprocessing
    # --------------
    #Outputs = mps.Result(parameters, 'eMPS', 1)
    Outputs = mps.ReadStaticObservables(parameters)
    energies = np.zeros(lambda_list.shape[0]*totalstate)

    end_time = time.time()

    elapsed_time = end_time - start_time

    print(elapsed_time)
    #print(len(Outputs))

    itr_count = min(100, lambda_list.shape[0])
    #print(itr_count)
    for lambda_i in range(itr_count*totalstate):

        '''
        spdm = Outputs[lambda_i]['spdm']
        spdmeigs, U = np.linalg.eigh(spdm)
        bspdm = Outputs[lambda_i]['bspdm']
        bspdmeigs, U = np.linalg.eigh(bspdm)
        '''

        energies[lambda_i] = Outputs[lambda_i]['energy']

        #print('Eigenvalues of <f^{\dagger}_i f_j> with Fermi phases', spdmeigs)
        #$print('Eigenvalues of <f^{\dagger}_i f_j> without Fermi phases', bspdmeigs)
        #print(Outputs[lambda_i]['energy'])
        #print(Outputs[lambda_i]['state'])


    gsenergy = np.zeros(itr_count)
    newenergy = np.zeros((totalstate, itr_count))
    for state in range(totalstate):
        for lambda_i in range(itr_count):
            newenergy[state][lambda_i] = energies[lambda_i*totalstate + state] - energies[lambda_i*totalstate]
            gsenergy[lambda_i] = energies[lambda_i*totalstate]
    print(lambda_list)
    print(gsenergy)


    file_name = testflag + str(L) + '_' + str(N) + 'res_ex_' + str(exchange)  + '_eig' + str(eigenstate) +'_range' +  str(int_range) +'_' + str(max_sweep) + '.dat'
    gs_name = 'gsvar' +  testflag + str(L) + '_' + str(N) + 'int' + str(int_range) + '_ex_' + str(exchange)  + '_eescale_' + str(eefactor) + '_nscale_' + str(nfactor) + specflag + '.dat'
    plot_name = testflag + str(L) + '_' + str(N) + 'energy_ex_' + str(exchange)   + '_eig' +  str(eigenstate) + '_range' +  str(int_range) +'_' + str(max_sweep) + '.png'

    np.savetxt(file_name, newenergy, fmt='%10.5f')
    np.savetxt(gs_name, gsenergy, fmt='%8.5f')
    #eplot(lambda_list, newenergy, plot_name, itr_count)


    return


def eplot(lambda_list, energies, plot_name, itr_count):
    colors = ['red', 'red', 'blue', 'green']

    fig = plt.figure()
    ax = plt.gca()
    for eig in range(1, energies.shape[0]):
        ax.scatter(lambda_list[:itr_count] ,energies[eig],  alpha=1, color=colors[eig], edgecolors='none')
    #ax.set_yscale('log')
    #ax.set_xscale('log')
    plt.ylim([0, 2.5])
    plt.xlim([0, 23])
    plt.savefig(plot_name)
    #plt.loglog(lambda_list, energies)
    #plt.ylim(1e-2, 1e8)
    #plt.show()

if(__name__ == '__main__'):
    # Check for command line arguments
    Post = False
    for arg in sys.argv[1:]:
        key, val = arg.split('=')
        if(key == '--PostProcess'): Post = (val == 'T') or (val == 'True')

    # Run main function
    main(PostProcess=Post)
