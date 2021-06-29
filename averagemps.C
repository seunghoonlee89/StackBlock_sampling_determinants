/*
 * =====================================================================================
 *
 *       Filename:  readrotatationmatrix.C
 *
 *    Description: Read rotation matrix from external DMRG code to build spinblocks.
 *
 *        Version:  1.0
 *        Created:  02/13/2017 12:29:47 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Sheng Guo 
 *
 * =====================================================================================
 */


#include "MatrixBLAS.h"
#include "Stackspinblock.h"
#include "initblocks.h"
#include "input.h"
#include "rotationmat.h"
#include "global.h"
#ifndef SERIAL
#include "mpi.h"
#include <boost/mpi.hpp>
#include "readrotationmatrix.h"
#include "sweep.h"
#include "Stackdensity.h"
#include "pario.h"
#include "Stackwavefunction.h"
#include "stackguess_wavefunction.h"
#include "operatorfunctions.h"
#endif

void ReadInput(char* conf);
void restart(double sweep_tol, bool reset_iter);

using namespace SpinAdapted;
using namespace SpinAdapted::Sweep;

void getComplementarySites(std::vector<int> &sites, std::vector<int> &complementarySites) ;

int test(int argc, char* argv[])
{
  std::vector<SpinQuantum> quanta;
  std::vector<int> quantaStates;
  std::vector<int> newquantaStates;
  std::vector<double> oldmatrix;
  for(int dotsite = 0; dotsite < 10; dotsite++){
    cout << "Site : " << dotsite <<endl;
    ReadQuanta(quanta, quantaStates, newquantaStates, dotsite+1);
    cout << "Quanta " <<endl;
    int states = 0;
    for(int i=0; i< quanta.size(); ++i)
    {
      states += quantaStates[i];

      cout << quanta[i] <<" " <<quantaStates[i] << " " <<newquantaStates[i]<<endl;
    }

    std::vector<Matrix> RotationM;
    ReadRotationMatrix(RotationM, quantaStates, newquantaStates, dotsite);
    cout << " rotation matrix " <<endl;
    for(int i=0; i< RotationM.size(); ++i)
    {

      cout << RotationM[i]<<endl;
    }
  }
  return 0;

}

int main(int argc, char* argv[])
{
  //test(argc,argv);
#ifndef SERIAL
  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world;
  if (world.rank() == 0) {
    pout << "Runing with " << world.size() << " processors" << endl;
  }
#endif

  ReadInput(argv[1]);
  dmrginp.matmultFlops.resize(numthrds, 0.);

  double* stackmemory = new double[dmrginp.getMemory()];
  Stackmem.resize(numthrds);
  Stackmem[0].data = stackmemory;
  Stackmem[0].size = dmrginp.getMemory();

  dmrginp.initCumulTimer();


  //BuildBlock(-1);
  BuildFromRotationMatrix(0);
  BuildFromRotationMatrix(-1);
  SweepParams sweepParams;
  sweepParams.set_sweep_iter() = 0;
  sweepParams.current_root() = -1;
	Sweep::CanonicalizeWavefunction(sweepParams, false, 0);
	Sweep::CanonicalizeWavefunction(sweepParams, false, -1);
	Sweep::CanonicalizeWavefunction(sweepParams, true, 0);
	Sweep::CanonicalizeWavefunction(sweepParams, true, -1);
	Sweep::CanonicalizeWavefunction(sweepParams, false, 0);
	Sweep::CanonicalizeWavefunction(sweepParams, false, -1);
  sweepParams.savestate(true,1);
  //do_one();
	//restart(1e-7, false);
  return 0;
}


double SpinAdapted::AverageMPS::do_one(SweepParams &sweepParams, const bool &warmUp, const bool &forward, const bool &restart, const int &restartSize, int targetState=-1)
{
  double finalError = 0.0;
  int integralIndex = 0;
  StackSpinBlock system;
  const int nroots = dmrginp.nroots(sweepParams.get_sweep_iter());

  sweepParams.set_sweep_parameters();
  // a new renormalisation sweep routine
  pout << endl;
  if (forward)
    { pout << "\t\t\t Starting sweep "<< sweepParams.set_sweep_iter()<<" in forwards direction"<<endl; }
  else
    { pout << "\t\t\t Starting sweep "<< sweepParams.set_sweep_iter()<<" in backwards direction" << endl; }
  pout << "\t\t\t ============================================================================ " << endl;

  //TODO
  restart = false;
  InitBlocks::InitStartingBlock (system,forward, -1, 0, sweepParams.get_forward_starting_size(), sweepParams.get_backward_starting_size(), restartSize, restart, warmUp, integralIndex);
  if(!restart)
    sweepParams.set_block_iter() = 0;

 
  p2out << "\t\t\t Starting block is :: " << endl << system << endl;

  for (int i=0; i < nroots; i++)
    StackSpinBlock::store (forward, system.get_sites(), system, -1, i); // if restart, just restoring an existing block --
  sweepParams.savestate(forward, system.get_sites().size());
  bool dot_with_sys = true;
  vector<int> syssites = system.get_sites();

  if (dmrginp.outputlevel() > 0)
    mcheck("at the very start of sweep");  // just timer

  for (; sweepParams.get_block_iter() < sweepParams.get_n_iters(); ) // get_n_iters() returns the number of blocking iterations needed in one sweep
    {
      pout << "\n\t\t\t Block Iteration :: " << sweepParams.get_block_iter() << endl;
      pout << "\t\t\t ----------------------------" << endl;

      if (sweepParams.get_block_iter() != 0) 
	      sweepParams.set_guesstype() = TRANSFORM;
      else
        sweepParams.set_guesstype() = TRANSPOSE;


      
      p1out << "\t\t\t Blocking and Decimating " << endl;
	  
      StackSpinBlock newSystem; // new system after blocking and decimating

      //Need to substitute by:
      if (warmUp )
      {
        for (int i=0; i < nroots; i++)
	        Startup(sweepParams, syssites, newSystem, 0, i);
      }
      else {
	      BlockDecimateAndCompress (sweepParams, false, dot_with_sys, syssites);
      }
      

      //Need to substitute by?

      if (!warmUp ){
	      finalError = max(sweepParams.get_lowest_error(),finalError);
      }
      
      dot_with_sys = true;
      p1out << "\t\t\t saving state " << syssites.size() << endl;
      ++sweepParams.set_block_iter();
      
#ifndef SERIAL
      mpi::communicator world;
      calc.barrier();
#endif
      sweepParams.savestate(forward, syssites.size());
      if (dmrginp.outputlevel() > 0)
         mcheck("at the end of sweep iteration");
    }


  //when we are doing twodot, we still need to do the last sweep to make sure that the
  //correctionVector and base wavefunction are propogated correctly across sweeps
  //especially when we switch from twodot to onedot algorithm
  if (!sweepParams.get_onedot() && !warmUp) {
      pout << "\n\t\t\t Block Iteration :: " << sweepParams.get_block_iter() << endl;
      pout << "\t\t\t ----------------------------" << endl;
      if (forward)
	{ p1out << "\t\t\t Current direction is :: Forwards " << endl; }
      else
	{ p1out << "\t\t\t Current direction is :: Backwards " << endl; }
    sweepParams.set_onedot() = true;
    sweepParams.set_env_add() = 0;
    bool dot_with_sys = true;
    WavefunctionCanonicalize(sweepParams, system, warmUp, dot_with_sys, targetState, baseState);
    sweepParams.set_onedot() = false;
    sweepParams.set_env_add() = 1;
  }

  //system.deallocate();

  pout << "\t\t\t Largest Error for Sweep with " << sweepParams.get_keep_states() << " states is " << finalError << endl;
  sweepParams.set_largest_dw() = finalError;
  

  pout << "\t\t\t ============================================================================ " << endl;

  // update the static number of iterations

  ++sweepParams.set_sweep_iter();

  return finalError;
}

void SpinAdapted::AverageMPS::BlockDecimateAndCompress (SweepParams &sweepParams, const bool &useSlater, const bool& dot_with_sys, std::vector<int>& systemsites)
{
  int sweepiter = sweepParams.get_sweep_iter();
  const int nroots = dmrginp.nroots(sweepParams.get_sweep_iter());

  p2out << "\t\t\t dot with system "<<dot_with_sys<<endl;
  p1out <<endl<< "\t\t\t Performing Blocking"<<endl;

  dmrginp.guessgenT -> start();
  bool forward = (system.get_sites() [0] == 0);
  StackSpinBlock systemDot;
  StackSpinBlock environment, environmentDot, newEnvironment;
  StackSpinBlock envDot, big;
  int systemDotStart, systemDotEnd;
  int environmentDotStart, environmentDotEnd, environmentStart, environmentEnd;
  int systemDotSize = sweepParams.get_sys_add() - 1;
  int environmentDotSize = sweepParams.get_env_add() -1;


  StackDensityMatrix tracedMatrix;

  if (forward)
  {
    systemDotStart = dmrginp.spinAdapted() ? *system.get_sites().rbegin () + 1 : (*system.get_sites().rbegin ())/2 + 1 ;
    systemDotEnd = systemDotStart + systemDotSize;
    environmentDotStart = systemDotEnd + 1;
    environmentDotEnd = environmentDotStart + environmentDotSize;
  }
  else
  {
    systemDotStart = dmrginp.spinAdapted() ? system.get_sites()[0] - 1 : (system.get_sites()[0])/2 - 1 ;
    systemDotEnd = systemDotStart - systemDotSize;
    environmentDotStart = systemDotEnd - 1;
    environmentDotEnd = environmentDotStart - environmentDotSize;
  }
  systemDot = StackSpinBlock(systemDotStart, systemDotEnd, system.get_integralIndex(), false);//singleSiteBlocks[system.get_integralIndex()][systemDotStart];
  environmentDot = StackSpinBlock(environmentDotStart, environmentDotEnd, system.get_integralIndex(), false);//singleSiteBlocks[system.get_integralIndex()][environmentDotStart];


  for(int i=0;i<nroots;i++)
  {
    Sweep::makeSystemEnvironmentBigOverlapBlocks(systemsites, systemDot, environmentDot, system, 
        newSystem, environment, newEnvironment, big, sweepParams, dot_with_sys, useSlater, 0, -1, i);

    StackSpinBlock newbig;
    //TODO
    if (sweepParams.get_onedot() && !dot_with_sys)
    {
      newSystem.set_integralIndex() = perturbationsystem.get_integralIndex();
      newSystem.initialise_op_array(OVERLAP, false);
      newSystem.setstoragetype(DISTRIBUTED_STORAGE);
      newSystem.BuildSumBlock(NO_PARTICLE_SPIN_NUMBER_CONSTRAINT, system, systemdot);
      InitBlocks::InitBigBlock(newSystem, environment, newbig); 

      perturbationBig.get_rightBlock()->clear();
      perturbationBig.clear();
    }
    else
      newbig = big;

    StackWavefunction solution, projection; 
    solution.initialise(dmrginp.effective_molecule_quantum_vec(), newbig.get_leftBlock()->get_ketStateInfo(), newbig.get_rightBlock()->get_ketStateInfo(), sweepParams.get_onedot());
    solution.set_onedot(sweepParams.get_onedot());
    solution.Clear();

    projection.initialise(dmrginp.effective_molecule_quantum_vec(), newbig.get_leftBlock()->get_braStateInfo(), newbig.get_rightBlock()->get_brastateInfo(), sweepParams.get_onedot()); 
    projection.solution.set_onedot(sweepParams.get_onedot());
    projection.Clear();


    GuessWave::guess_wavefunctions(solution, e, big, sweepParams.set_guesstype(), sweepParams.get_onedot(), dot_with_sys, 1, 0.0, baseState); 
    if (mpigetrank() == 0)
      newbig.multiplyOverlap(solution, iwave[i], MAX_THRD);

    projection.SaveWavefunctionInfo (newbig.get_ketStateInfo(), newbig.get_leftBlock()->get_sites(), i);
    if (i==0) 
      tracedMatrix.allocate(perturbationnewbig.get_leftBlock()->get_braStateInfo());
    operatorfunctions::MultiplyWithOwnTranspose (projection, tracedMatrix, 1.0);  

    projection.deallocate();
    solution.deallocate();
    newEnvironment.deallocate();
    environment.deallocate();
    newSystem.deallocate();
    system.deallocate();
  }

  std::vector<Matrix> brarotatematrix;
  int largeNumber = 1000000;
	double error = makeRotateMatrix(tracedMatrix, brarotatematrix, largeNumber, sweepParams.get_keep_qstates());
  

  //Transform and store operators.
  for(int i=0;i<nroots;i++)
  {

    Sweep::makeSystemEnvironmentBigOverlapBlocks(systemsites, systemDot, overlapenvironmentdot, 
						 overlapsystem, overlapnewsystem, overlapenvironment, 
						 overlapnewenvironment, overlapBig, sweepParams, dot_with_sys, useSlater,
						 0, targetState, projectors[l]);

    Sweep::makeSystemEnvironmentBigOverlapBlocks(system, systemDot, newSystem, environment, environmentDot, newEnvironment, big, sweepParams, dot_with_sys, useSlater, system.get_integralIndex(), -1, i);
    StackSpinBlock newbig;
    //TODO
    if (sweepParams.get_onedot() && !dot_with_sys)
    {
      newSystem.set_integralIndex() = perturbationsystem.get_integralIndex();
      newSystem.initialise_op_array(OVERLAP, false);
      newSystem.setstoragetype(DISTRIBUTED_STORAGE);
      newSystem.BuildSumBlock(NO_PARTICLE_SPIN_NUMBER_CONSTRAINT, system, systemdot);
      InitBlocks::InitBigBlock(newSystem, environment, newbig); 

      perturbationBig.get_rightBlock()->clear();
      perturbationBig.clear();
    }
    else
      newbig = big;

    if (i==0)
      SaveRotationMatrix (newbig.get_leftBlock()->get_sites(), brarotateMatrix, -1);
    if (i==nroot-1)
     systemsites = newbig.get_leftBlock()->get_sites()

    SpinAdapted::LoadRotationMatrix (newSystem.get_sites(), rotatematrix, i);
    newSystem.transform_operators(brarotatematrix, rotatematrix, false, false);

    newEnvironment.deallocate();
    environment.deallocate();
    newSystem.deallocate();
    system.deallocate();

  }

  kettracedMatrix.deallocate();
  bratracedMatrix.deallocate();

#ifndef SERIAL
  broadcast(calc, brarotateMatrix, 0);
#endif




  if (dmrginp.outputlevel() > 0)
    mcheck("after rotation and transformation of block");

  p2out << *dmrginp.guessgenT<<" "<<*dmrginp.multiplierT<<" "<<*dmrginp.operrotT<< "  "<<globaltimer.totalwalltime()<<" timer "<<endl;
  p2out << *dmrginp.makeopsT<<" makeops "<<endl;
  p2out << *dmrginp.datatransfer<<" datatransfer "<<endl;
  p2out <<"oneindexopmult   twoindexopmult   Hc  couplingcoeff"<<endl;  
  p2out << *dmrginp.oneelecT<<" "<<*dmrginp.twoelecT<<" "<<*dmrginp.hmultiply<<" "<<*dmrginp.couplingcoeff<<" hmult"<<endl;
  p2out << *dmrginp.buildsumblock<<" "<<*dmrginp.buildblockops<<" build block"<<endl;

}



void SpinAdapted::AverageMPS::Startup (SweepParams &sweepParams, std::vector<int> syssites, int stateA, int stateB)
{

  StackSpinBlock system;
  if (dmrginp.outputlevel() > 0) 
    mcheck("at the start of block and decimate");
  p1out << "\t\t\t Performing Blocking"<<endl;
  dmrginp.guessgenT -> start();
  // figure out if we are going forward or backwards  
  bool forward = (syssites[0] == 0);

  StackSpinBlock::restore(forward, syssites, system, -1, stateB);


  StackSpinBlock systemDot;
  int systemDotStart, systemDotEnd;
  int systemDotSize = sweepParams.get_sys_add() - 1;
  if (forward)
  {
    systemDotStart = dmrginp.spinAdapted() ? *system.get_sites().rbegin () + 1 : (*system.get_sites().rbegin ())/2 + 1 ;
    systemDotEnd = systemDotStart + systemDotSize;
  }
  else
  {
    systemDotStart = dmrginp.spinAdapted() ? system.get_sites()[0] - 1 : (system.get_sites()[0])/2 - 1 ;
    systemDotEnd = systemDotStart - systemDotSize;
  }
  vector<int> spindotsites(2); 
  spindotsites[0] = systemDotStart;
  spindotsites[1] = systemDotEnd;

  systemDot = StackSpinBlock(systemDotStart, systemDotEnd, system.get_integralIndex(), true);//singleSiteBlocks[system.get_integralIndex()][systemDotStart];

  const int nexact = forward ? sweepParams.get_forward_starting_size() : sweepParams.get_backward_starting_size();

  dmrginp.guessgenT -> stop();

  newSystem.set_integralIndex() = system.get_integralIndex();
  newSystem.initialise_op_array(OVERLAP, false);
  newSystem.setstoragetype(DISTRIBUTED_STORAGE);
  newSystem.BuildSumBlock (NO_PARTICLE_SPIN_NUMBER_CONSTRAINT, system, systemDot);


  pout << "\t\t\t System  Block"<<newSystem;
  newSystem.printOperatorSummary();

  std::vector<Matrix> leftrotateMatrix, rightrotateMatrix;

  LoadRotationMatrix (newSystem.get_sites(), leftrotateMatrix, stateA);
  SaveRotationMatrix (newSystem.get_sites(), leftrotateMatrix, -1);
  LoadRotationMatrix (newSystem.get_sites(), rightrotateMatrix, stateB);

    
#ifndef SERIAL
  mpi::communicator world;
  broadcast(calc, leftrotateMatrix, 0);
  broadcast(calc, rightrotateMatrix, 0);
#endif

  p1out <<"\t\t\t Performing Renormalization "<<endl<<endl;

  dmrginp.operrotT->start();
  newSystem.transform_operators(leftrotateMatrix, rightrotateMatrix);
  dmrginp.operrotT->stop();

  StackSpinBlock::store(forward, newSystem.get_sites(), newSystem, -1, stateB);
  newSystem.deallocate();
  system.deallocate();

}

