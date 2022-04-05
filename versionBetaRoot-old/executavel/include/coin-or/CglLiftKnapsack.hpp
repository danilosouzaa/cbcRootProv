// Copyright (C) 2000, International Business Machines
// Corporation and others.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#ifndef CglLiftKnapsack_H
#define CglLiftKnapsack_H

#include <string>
#include <numeric>

#include "CglCutGenerator.hpp"
#include "CglTreeInfo.hpp"
#include "CoinCutPool.hpp"
#include "OsiClpSolverInterface.hpp"



/** Knapsack Cover Cut Generator Class */
class CGLLIB_EXPORT CglLiftKnapsack : public CglCutGenerator {
   friend CGLLIB_EXPORT void CglLiftKnapsackUnitTest(const OsiSolverInterface * siP, const std::string mpdDir );

public:
   /** A method to set which rows should be tested for knapsack covers */

   /**@name Generate Cuts */
  //@{
  /** Generate knapsack cover cuts for the model of the solver interface, si. 
      Insert the generated cuts into OsiCut, cs.
  */
  virtual void generateCuts(const OsiSolverInterface & si, OsiCuts & cs, const CglTreeInfo info = CglTreeInfo());

  //@}

  /**@name Constructors and destructors */
  //@{
  /// Default constructor 
  CglLiftKnapsack ();
 
  /// Copy constructor 
  CglLiftKnapsack (const CglLiftKnapsack &);

  /// Clone
  virtual CglCutGenerator * clone() const;

  /// Assignment operator 
  CglLiftKnapsack & operator=(const CglLiftKnapsack& rhs);
  
  /// Destructor 
  virtual
    ~CglLiftKnapsack ();

  /// This can be used to refresh any information
  virtual void refreshSolver(OsiSolverInterface * solver);
  //@}


  /**@name Sets and gets */
  //@{
  /// Set limit on number in knapsack
  inline void setMaxInKnapsack(int value)
           { if (value>0) maxInKnapsack_ = value;}
  /// get limit on number in knapsack
  inline int getMaxInKnapsack() const
           {return maxInKnapsack_;}
  /// Switch off expensive cuts
  inline void switchOffExpensive()
  { expensiveCuts_=false;}
  /// Switch on expensive cuts
  inline void switchOnExpensive()
  { expensiveCuts_=true;}
  //verify if value is integer
  inline bool isValueInteger(double v){
    return ceil(v)==floor(v);
  }


  inline size_t getCutsAgg() const
    {return numCutsAgg_;}
  inline size_t getCutsAgr() const
    {return numCutsAgr_;}
  
private:
  
 // Private member methods

  int deriveAKnapsack(
    const OsiSolverInterface & si, 
    OsiCuts & cs,
    CoinPackedVector & krow,
    bool treatAsLRow,
    double & b,
    int *  complement,
    double *  xstar,
    int rowIndex,
    int numberElements,
    const int * index,
    const double * element);

  int deriveAKnapsack(
    const OsiSolverInterface & si, 
    OsiCuts & cs,
    CoinPackedVector & krow,
    double & b,
    int *  complement,
    double *  xstar,
    int rowIndex,
    const CoinPackedVectorBase & matrixRow);


  std::string generateCpp(FILE *fp);


  /**@name Private methods */
  //@{

  /** deriveAKnapsack 
                 returns 1 if it is able to derive
                 a (canonical) knapsack inequality
                in binary variables of the form ax<=b 
                 from the rowIndex-th  row in the model, 
                returns 0 otherwise.
  */

//-------------------------------------------------
//Method Danilo Souza
//-------------------------------------------------
  bool valueIsInteger(double v);
  
  


  void quicksortCof(double *values, int *idc, int began, int end);

  void shuffleVectorInt(int *vec, int sz);


void quicksortInts(int *values, double *idc, int began,
                                   int end);

 
  bool tryConvertCoefFracInInt(double *temp,int sz);



  void backTracking(int *v, double *p, int idCandidato, double excesso, int *solutions, int *numberSolutions, int szSolution);

  void quicksortDoubleA(double *values, std::vector<int> &idc, int began, int end) ;

  void quicksortIncludeComplement(double *values, int *idx, int *complement, int began,
                                   int end);

  int* algoritm_4_1(int *w, int c, int sz, double *p);

  int* algoritm_4_2(int *w, int c, int sz, double *p);

  int *algoritm_4_2_double(double *w, double c, int sz, double *p);

  int* algoritm_4_3(int *a, int a_0, int sz, std::vector<int> C1, std::vector<int> C2, std::vector<int> F,std::vector<int> R);
  

  //No applied lifting method
  int *returnCoverMinimal_4_5(double *a_temp, double a_0_temp, double *xAsterisc, int sz);
  /*
  double *algoritm_4_3_double(double *a, double a_0, int sz, std::vector<int> C1,
                                   std::vector<int> C2, std::vector<int> F,
                                   std::vector<int> R);*/
  
  double *algoritm_4_4(double *a, double a_0, int sz, std::vector<int> C, std::vector<int> N_C);

  void quicksortIntA(int *values, std::vector<int> &idc, int began, int end);

  void quicksortVector(double *values, std::vector<int> &idc, double *a, int began, int end);

  void quicksortVectorTwoCriterio(double *values, std::vector<int> &C, int *a);

  void quicksortVectorTwoCriterioDouble(double *values,std::vector<int> &C, double *a);

  int *algoritm_4_5(int *a, int a_0, double *xAsterisc, int sz);

  double run4_5_LCI_Letchford_and_Souli(double *a, double a_0, double *xAsterisc, int sz, std::vector<std::vector<double>> &poolSolution);

  int *algoritm_4_6(int *a, int a_0, double *xAsterisc, int sz);

  double *algoritm_4_7(double *a, double a_0, double *xAsterisc, int sz);

  void SAPerConstraints(double *a, double b, double *xAsterisc, int szConst, double timeLeft,std::vector<std::vector<double>> &poolSolution);

  double LCIAdamPerConstraints(int *coverSolution, double *a, double a_0, double *xAsterisc, int szConstraints, std::vector<double> &sol);

  int *newNeighborhoodOne(int *solutionCurrent, int sz);

  int *newNeighborhoodTwo(int *solutionCurrent, int sz);
  int *newNeighborhoodThree(int *solutionCurrent, int sz);
  int *newNeighborhoodFour(int *solutionCurrent, int sz);
  int *newNeighborhoodFive(int *solutionCurrent, int sz);

  void BackTrackingPerConstraints(
    double *a, double a_0, double *xAsterisc, int szConst, double timeLeft,
    std::vector<std::vector<double>> &poolSolution);


  double mixLiftCover(int *coverSolution, double *a, double a_0,
                                     double *xAsterisc, int szConstraints,
                                     std::vector<double> &sol);

  int isMinimalDouble(double *a, double a_0, int sz, std::vector<int> C);


int agroupBKColCliqueEqualCoef(double *a, int *idxCol, int sz, 
                                  int *complement, long int *idxClique,
                                  std::vector<double> &xStarClique,
                                  const OsiSolverInterface &si);

void separateCuts(double *w, double w_0, int sz, int *idxCol, 
                                    int *complement, const OsiSolverInterface &si
                                    , double timeLeft, OsiCuts &cs, const CglTreeInfo &info);

void separateCutsLSCDC(double *w, double w_0, int sz, int *idxCol,
                                   int *complement, const OsiSolverInterface &si,
                                   double timeLeft, OsiCuts &cs, const CglTreeInfo &info);

void fixVar(int idxCol, int complement);

int *greedyProposedCover(double *w, double rhs, int sz, double *p);

bool randomClq(std::set<size_t> &clqCur, std::map<size_t, std::vector<size_t>> dephElements, size_t maxClq, size_t nElements, double *fitnessBest, double *minCoefBest,double *maxCoefBest, std::map<size_t, std::vector<double>> p, const CoinConflictGraph *ppcg, double &dltMax, double pMax);

size_t heuristicCliqueFound(long int *idxClique, size_t *idOrig, size_t *posIndVtx, size_t nCols, const CoinConflictGraph *ppcg, double *xstar, double *w, std::vector<size_t> &itemRep, std::vector<double> &xstarClq, double dltMaxIni);

std::set<size_t> neighborhoodOneClq(std::set<size_t> oriClq, std::map<size_t, std::vector<size_t>> dephElements, size_t maxClq, size_t nElements, double *fitnessBest, double *minCoefBest, double *maxCoefBest, std::map<size_t, std::vector<double>> p, const CoinConflictGraph *ppcg, double &dltMax, double pMax);
std::set<size_t> neighborhoodTwoClq(std::set<size_t> oriClq, std::map<size_t, std::vector<size_t>> dephElements, size_t maxClq, size_t nElements, double *fitnessBest, double *minCoefBest, double *maxCoefBest, std::map<size_t, std::vector<double>> p, const CoinConflictGraph *ppcg, double &dltMax, double pMax);
std::set<size_t> neighborhoodThreeClq(std::set<size_t> oriClq, std::map<size_t, std::vector<size_t>> dephElements, size_t maxClq, size_t nElements, double *fitnessBest, double *minCoefBest, double *maxCoefBest, std::map<size_t, std::vector<double>> p, double &dltMax);
std::set<size_t> neighborhoodFourClq(std::set<size_t> oriClq, std::map<size_t, std::vector<size_t>> dephElements, size_t maxClq, size_t nElements, double *fitnessBest, double *minCoefBest, double *maxCoefBest, std::map<size_t, std::vector<double>> p, const CoinConflictGraph *ppcg, double &dltMax, double pMax);
std::set<size_t> neighborhoodFiveClq(std::set<size_t> oriClq, std::map<size_t, std::vector<size_t>> dephElements, size_t maxClq, size_t nElements, double *fitnessBest, double *minCoefBest, double *maxCoefBest, std::map<size_t, std::vector<double>> p, const CoinConflictGraph *ppcg, double &dltMax, double pMax);

bool insertClqPool(std::map<size_t, std::vector<size_t>> &elClqPool,std::map<size_t, std::vector<size_t>> &clqPool, std::set<size_t> clq, size_t &szPoolClq,double *w, std::vector<double> &sumCofClq);

std::vector<size_t> updatePoolDisjointClq(std::map<size_t, std::vector<size_t>> &clqPool, size_t &szPoolClq, std::vector<double> listClqWeight, size_t nCols, double dltMax,  std::vector<double> sumCofClq);

int agroupColDisjointCliques(double *a, double b, int *idxCol, int sz,
                                              int *complement, long int *idxClique, std::vector<size_t> &itemRep,
                                              std::vector<double> &xStarClique, 
                                              const OsiSolverInterface &si/*, double deltaCoef*/);

//----------------------------------------------------


public:
  /** Creates cliques for use by probing.
      Only cliques >= minimumSize and < maximumSize created
      Can also try and extend cliques as a result of probing (root node).
      Returns number of cliques found.
  */




    /**
   * Number of cuts separated.
   **/
  static size_t sepCuts_;
    /**
   * Execution time spent for the LCI cut separator.
   **/
  static double sepTime_;


private:


  // Private member data

  //Danilo ;
  size_t maxIteSrcClqHeu_;
  size_t maxExactClq_;
  double timeRemaining_;
  /**@name Private member data */
  //@{
  /// epsilon
  double epsilon_;  
  /// Tolerance to use for violation - bigger than epsilon_
  double epsilon2_;
  /// 1-epsilon
  double onetol_;  
  /// Maximum in knapsack
  int maxInKnapsack_;
  /** which rows to look at. If specified, only these rows will be considered
      for generating knapsack covers. Otherwise all rows will be tried */
  int numRowsToCheck_;
  int* rowsToCheck_;
  /// exactKnapsack can be expensive - this switches off some
  bool expensiveCuts_;
  /// Cliques
  /// **** TEMP so can reference from listing
  const OsiSolverInterface * solver_;
  int whichRow_;
  int * complement_;
  double * elements_;
  /// Number of cliques
  /// Clique type

  /// Number of columns
  int numberColumns_;
size_t numCutsAgg_;
size_t numCutsAgr_;


  /**
  * Auxiliary structure used to temporary
  * store a cut.
  **/
  OsiRowCut osrc_;
  /** For each column with nonzero in row copy this gives a clique "number".
      So first clique mentioned in row is always 0.  If no entries for row
      then no cliques.  If sequence > numberColumns then not in clique.
  */
  //CliqueEntry * cliqueRow_;
  /// cliqueRow_ starts for each row
  //int * cliqueRowStart_;
  //@}
};

//#############################################################################
/** A function that tests the methods in the CglLiftKnapsack class. The
    only reason for it not to be a member method is that this way it doesn't
    have to be compiled into the library. And that's a gain, because the
    library should be compiled with optimization on, but this method should be
    compiled with debugging. */
CGLLIB_EXPORT
void CglLiftKnapsackUnitTest(const OsiSolverInterface * siP,
			      const std::string mpdDir );
  
#endif
