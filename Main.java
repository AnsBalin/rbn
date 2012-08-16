
import java.util.*;
import java.io.*;

public class Main {
  
  private ArrayList<Network> library = new ArrayList<Network>();
  private ArrayList<Molecule> bucket = new ArrayList<Molecule>();
  private ArrayList<Reaction> reactions = new ArrayList<Reaction>();
  private int numReactions;
  private static double temperature = 1;
  private static boolean reuse = false;
  private static String path = "/data/library/";
  private static String suffix = "";
  
  private int product;
  
  // Constructor adds null network to library, with index 0
  public Main(int numSpecies, int numNodes, int numBondingSites, int numMolecules){
    
    System.out.printf("\n\tPopulating library...\n");
    populateLibrary(numSpecies, numNodes, numBondingSites);
    System.out.printf("\tLibrary successfully populated with:\n\t\t%d Species\n\t\t%d nodes each\n\t\t%d bonding sites each", numSpecies, numNodes, numBondingSites);
    
    System.out.printf("\n\n\tPopulating world...\n");
    populateWorld(numMolecules, numSpecies);
    System.out.printf("\tWorld successfully populated with %d atoms.\n", numMolecules);
    
    numReactions=0;
    stats();
    
  }
  
  public Main(int numSpecies, int A, int B, int C, int numMolecules, double concentration){
    
    ArrayList<Network> temp = new ArrayList<Network>();
    
    for(int i=0; i<numSpecies; i++){
      
      DataOutput in = new DataOutput(path+i+".txt");
      try{
        temp.add(in.readFile());
        temp.get(0).setPop(0);
      }catch(IOException e){
        e.printStackTrace();
      }catch(ClassNotFoundException e){
        e.printStackTrace();
      }
      
    }
    library.add(temp.get(0));
    library.add(temp.get(A));
    library.add(temp.get(B));
    library.add(temp.get(C));
    
    Molecule nullMolecule = new Molecule();
    bucket.add(nullMolecule);

    for(int i =1; i < (int) (numMolecules*(1 - concentration)); i+=2){
      
      bucket.add(new Molecule(1, nullMolecule, nullMolecule));
      library.get(1).incrPop();
      bucket.add(new Molecule(2, nullMolecule, nullMolecule));
      library.get(2).incrPop();
      
      bucket.get(i).calculateMoleculeID(0,0, new int[0]);
      bucket.get(i+1).calculateMoleculeID(0,0, new int[0]);
      
      product = 0;
    
    }
  
    for(int i = (int)(numMolecules*(1 - concentration)); i<numMolecules; i++){
    
      bucket.add(new Molecule(3, nullMolecule, nullMolecule)); 
      library.get(3).incrPop();
      bucket.get(i).calculateMoleculeID(0,0, new int[0]);
      
    }
    
  }
  
  public int getCycleLength(int index){return library.get( bucket.get(index).getID() ).getCycleLength();}
  public double getActivity(int index){return library.get( bucket.get(index).getID() ).getActivity();}
  public int getPopulation(){return bucket.size();}
  public int getNumSpecies(){return library.size();}
  public int getNumReactions(){return numReactions;}
  
 public ArrayList<Reaction> getReactions(){return reactions;}
  public ArrayList<Network> getLibrary(){return library;}
  public ArrayList<Molecule> getBucket(){return bucket;}
  
  public void populationDistribution(){
    
    int[][] species = new int[library.size()][2];
    for(int i=0; i<bucket.size(); i++){
      
      Molecule M = bucket.get(i);
      species[M.getID()][0] = M.getID();
      species[M.getID()][1]++;
    }
    int k;
    for(int j=0; j<library.size(); j++){
      
      int a = species[j][1];
      k=j;
      while(k>0 && species[k-1][1] > a){
        species[k][1] = species[k-1][1];
        k--;
      }
      species[k][1]=a;
    }
    for(int j=0; j<library.size(); j++){
      
      System.out.printf("%d\n", species[j][1]);
    }
  }
  
  // If generateNew = 0     Generate new species
  // If generateNew != 0    Use species from previous run
  public void populateLibrary(int numSpecies, int numNodes, int numBondingSites){
    
    if(!reuse){
      library.add(new Network(0));
    
      for(int i=1; i<numSpecies; i++){
        library.add(new Network(numNodes, i, numBondingSites));
      
      }
    
      for(int i=0; i<numSpecies; i++){
        DataOutput out = new DataOutput("/data/library/"+i+".txt");
      
        try{
          out.writeToFile(library.get(i));
        }catch(IOException e){System.out.println("Could not write file.");}
      }
    }
    else{
      
      for(int i=0; i<numSpecies; i++){
          
        DataOutput in = new DataOutput(path+i+".txt");
        try{
          library.add(in.readFile());
          library.get(0).setPop(0);
        }catch(IOException e){
          e.printStackTrace();
        }catch(ClassNotFoundException e){
          e.printStackTrace();
        }
        
      }
      
    }
    
    
  }
  
  public void populateWorld(int numMolecules, int numSpecies){
    
    int rand;
    
    Molecule nullMolecule = new Molecule();
    bucket.add(nullMolecule);
    
    for(int i=1; i<numMolecules; i++){
      rand = (int) (Math.random()*(double)(numSpecies-1)) + 1;
      bucket.add(new Molecule(rand , nullMolecule, nullMolecule));
      bucket.get(i).calculateMoleculeID(0,0, new int[0]);
      bucket.get(i).setState(library.get(rand).getState());
      //String str = bucket.get(i).toString();
      //System.out.printf(str);
    }
  }
  
  public void stats(){
    
    int[] cycleLengths = new int[bucket.size()];
    double[] activities = new double[bucket.size()];
    
    double cycleLengthMean, activityMean;
    int cycleLengthMedian, cycleLengthUQ, cycleLengthLQ;
    double activityMedian, activityUQ, activityLQ, activityVar, activityStdDev, cycleLengthVar, cycleLengthStdDev;
    System.out.printf("\n\t==========================================");
    System.out.printf("\n\t                STATISTICS   ");
    System.out.printf("\n\t==========================================\n");
    System.out.printf("\n\tCalculating averages...\n");
    cycleLengthMean=0;
    activityMean=0;
    for(int i=0; i<bucket.size(); i++){
      
      cycleLengths[i] = getCycleLength(i);
      cycleLengthMean += (double)cycleLengths[i];
      
      activities[i] = getActivity(i);
      activityMean += activities[i];
    }
    
    cycleLengthMean = cycleLengthMean/bucket.size();
    activityMean = activityMean/bucket.size();
    
    System.out.printf("\tSorting arrays...\n");
    Arrays.sort(cycleLengths);
    Arrays.sort(activities);
    
    cycleLengthLQ = cycleLengths[(int) bucket.size()/4];
    cycleLengthMedian = cycleLengths[(int) bucket.size()/2];
    cycleLengthUQ = cycleLengths[(int) 3*bucket.size()/4];
    
    activityLQ = activities[(int) bucket.size()/4];
    activityMedian = activities[(int) bucket.size()/2];
    activityUQ = activities[(int) 3*bucket.size()/4];
    
    
    System.out.printf("\tCalculating variance...");
    
    activityVar=0;
    cycleLengthVar=0;
    cycleLengthStdDev=0;
    activityStdDev=0;
    for(int i=0; i<bucket.size(); i++){
      
      activityVar += (getActivity(i) - activityMean)*(getActivity(i) - activityMean);
      cycleLengthVar += ((double)getCycleLength(i) - cycleLengthMean)*((double)getCycleLength(i) - cycleLengthMean);
    }
    
    activityVar = activityVar/(bucket.size()-1);
    cycleLengthVar = cycleLengthVar/(bucket.size()-1);
    
    activityStdDev = Math.sqrt(activityVar);
    cycleLengthStdDev = Math.sqrt(cycleLengthVar);
    
    
    
    System.out.printf("\n\n\t\t\tCyclelength\tActivity\n");
    System.out.printf("\t==========================================\n");
    System.out.printf("\tMean\t\t%.2f\t\t%.2f\n", cycleLengthMean, activityMean);
    System.out.printf("\tMedian\t\t%d\t\t%.2f\n", cycleLengthMedian, activityMedian);
    System.out.printf("\tUpper Q\t\t%d\t\t%.2f\n", cycleLengthUQ, activityUQ);
    System.out.printf("\tLower Q\t\t%d\t\t%.2f\n", cycleLengthLQ, activityLQ);
    System.out.printf("\tVariance\t%.2f\t\t%.2f\n", cycleLengthVar, activityVar);
    System.out.printf("\tStd. Dev.\t%.2f\t\t%.2f\n\n", cycleLengthStdDev, activityStdDev);
    
    
  }
  
  public int[] molecularSizeDistrb(){
    
    int maxSize=0;
    
    for (int i=0; i<bucket.size(); i++){
      maxSize = (bucket.get(i).getSize() > maxSize)? bucket.get(i).getSize() : maxSize;
    }
    
    int[] sizeDistrb = new int[maxSize+1];
  
    for (int i=1; i<bucket.size(); i++){
    
      sizeDistrb[bucket.get(i).getSize()] += bucket.get(i).getSize();
      
    }
    
    return sizeDistrb;
    
  }
  
  public int[] cycleLengthDistrb(){
    
    int maxCL=0;
    Network A;
    
    for(int i=0; i<bucket.size(); i++){
      A = library.get( bucket.get(i).getID() );
      maxCL = (A.getCycleLength() > maxCL)? A.getCycleLength() : maxCL;
    }
    
    int[] CLDistrb = new int[maxCL+1];
    
    for(int i=0; i<bucket.size(); i++){
      A = library.get( bucket.get(i).getID() );
      
      CLDistrb[A.getCycleLength()]++;
      
    }
    
    return CLDistrb;
    
  }
  
  public void reactionAnalysis(){
    //System.out.printf("%d\n", reactions.size());
    for(int i=0; i<reactions.size(); i++){
      if(reactions.get(i).getCount() > 1){
        System.out.printf("%d\t", reactions.get(i).getCount());
        System.out.printf(""+reactionToString(reactions.get(i))+"\n");
      }
    }  
  }
  
  // selects n random molecules
  public ArrayList<Molecule> selectRandMols(int n){
    ArrayList<Molecule> randMols = new ArrayList<Molecule>();
    int max = bucket.size();
    int prevRand = 0;
    //System.out.printf("Max: %d\n", max);
    int rand;
    for(int i=0; i<n; i++){
      rand =  (int) (Math.random()*(double)max);
      while(rand==0 || rand==prevRand){rand =  (int) (Math.random()*(double)max);}
      randMols.add(bucket.get(rand));
      //System.out.printf("%d\t%d\n", rand, randMols.get(i).getID());
      prevRand = rand;
    }
    
    return randMols;
    
  }
  
  public ArrayList<Network> moleculeArrToNetworkArr(ArrayList<Molecule> molecules){
    
    ArrayList<Network> networks = new ArrayList<Network>();
    
    for(int i=0; i<molecules.size(); i++){
      
      networks.add(retrieve(molecules.get(i)));
    }
    return networks;
  }
  
  
  //Checks to see if Network A is in the library and if it isn't, add's it
  //returns index of where Network A is in the library if/when added
  public int libraryUpdate(Network A){
    int index=-1;
    int i = 0;
    while (index==-1 && i<library.size()){
      
      index = (A.equals(library.get(i))) ? i : -1;
      i++;
    }
    
    if(index == -1){
      library.add(A);
      index = library.size()-1;
      
    }
    
    return index;
  }
  
  public int reactionsUpdate(Reaction A){
  
    int index=-1;
    int i=0;
    
    while( index==-1 && i<reactions.size() ){
    
      index = ( A.getID().equals( reactions.get(i).getID() ) )? i:-1;
      i++;
      
    }
    
    if(index==-1){
      reactions.add(A);
      index = reactions.size()-1;
    }
    
    A.setIntID(index);
    
    return index;
    
    
  }
  
  
  public Network retrieve(int NetworkID){
    
    return library.get(NetworkID);
  }
  
  public Network retrieve(Molecule A){
    //System.out.printf("...%d\n", A.getID());
    return retrieve( A.getID() );
  }
  
  public Network bond(Molecule Asub, Molecule Bsub, int bondingSite1, int bondingSite2){
    
    Network netAsub, netBsub, netResult;
    netAsub = retrieve(Asub);
    netBsub = retrieve(Bsub);
    
    netResult = netAsub.bond(netBsub, bondingSite1, bondingSite2);
    return netResult;
    
  }
  
  public String reactionToString(Reaction R){
    
    String str = "";
    
    for(int i=0; i<R.getReactants().length; i++){
      str = str+R.getMolReactants().get(i).toStringf()+" ";
      if(i!=R.getReactants().length-1){
        str = str+"+ ";
      }
    }
    str=str+"---> ";
    for(int i=0; i<R.getProducts().length; i++){
      str = str+R.getMolProducts().get(i).toStringf()+"  ";
      if(i!=R.getProducts().length-1){
        str = str+"+ ";
      }
    }
    return str;
  
  }
  
  
  public boolean formBond(Molecule A, Molecule B, Molecule Asub, Molecule Bsub, int bondingSite1, int bondingSite2){
    
    molsA = A.split(Asub.getMolID());
    A.unFlagAll();
    molsB = B.split(Bsub.getMolID());
    B.unFlagAll();
    
    
    after.addAll( moleculeArrToNetworkArr(molsA) );
    after.addAll( moleculeArrToNetworkArr(molsB) );
    after.add( bond(Asub, Bsub, bondingSite1, bondingSite2) );

    
    numReactions++;
    ArrayList<Molecule> reactants = new ArrayList<Molecule>();
    ArrayList<Molecule> products = new ArrayList<Molecule>();
    
    reactants.add(A);
    reactants.add(B);
    
    //System.out.printf("\n\t"+A.toStringf()+" + "+B.toStringf()+" --> ");
    bucket.remove(A);
    bucket.remove(B);
    library.get(A.getID()).decrPop();
    library.get(B.getID()).decrPop();
    
    
    
    for(int i=0; i<molsA.size(); i++){
      
      bucket.add(molsA.get(i));
      products.add(molsA.get(i));
      library.get(molsA.get(i).getID()).incrPop();
      //System.out.printf(""+molsA.get(i).toStringf()+" + ");
    }
    
    for(int i=0; i<molsB.size(); i++){
        
        bucket.add(molsB.get(i));
        products.add(molsB.get(i));
        library.get(molsB.get(i).getID()).incrPop();
        //System.out.printf(""+molsB.get(i).toStringf()+" + ");
      }
    
    
    
    int index = libraryUpdate(after.get(after.size()-1));
    bucket.add(new Molecule(index, Asub, Bsub));
    library.get(index).incrPop();
    bucket.get(bucket.size()-1).setSize(Asub.getSize()+Bsub.getSize());
    
    products.add(bucket.get(bucket.size()-1));
    
    checkCatalysis(reactants, products);
    
    Reaction R = new Reaction(reactants, products);
    int reactionIndex = reactionsUpdate(R);
    reactions.get(reactionIndex).incrCount();
    if(reactions.get(reactionIndex).getCount()==1){
        String str1="";
        str1 += reactions.get(reactionIndex).getMolReactants().get(0).toStringf()+" + "+
                reactions.get(reactionIndex).getMolReactants().get(1).toStringf()+" --> ";
        
        str1 += reactions.get(reactionIndex).getMolProducts().get(0).toStringf();
        for(int i=1; i<reactions.get(reactionIndex).getMolProducts().size(); i++){
          str1 += " + "+reactions.get(reactionIndex).getMolProducts().get(i).toStringf();
        }
        
        System.out.printf(""+str1+"\n");
      }
    //System.out.printf(""+products.get(0).toStringf()+"\t%d\n", products.get(0).getID());
    //System.out.printf(""+bucket.get(bucket.size()-1).toStringf()+"\n");
    //System.out.printf("index: %d\nreaction #: %d\n", index, products.get(products.size()-).getID());
      //System.out.printf("\t"+reactions.get(reactionIndex).getID()+"\n");

    
  }
  
  public void checkCatalysis( ArrayList<Molecule> reactants, ArrayList<Molecule> products ){
    boolean containsMon, containsDi;
    
    containsMon=false;
    containsDi=false;
    int ABIndex, CIndex, ACIndex;
    double A1, A2, A3;
    ABIndex=0;
    CIndex=0;
    ACIndex=0;
    
    // Check reaction for form: AC + B --> AB + C
    if( (products.size() == 2) && (products.get(0).getID() != reactants.get(0).getID()) ){
      
      for(int i=0; i<products.size(); i++){
        if(products.get(i).getSize()==1){containsMon=true; CIndex=i;}
        if(products.get(i).getSize()==2){containsDi=true; ABIndex = i;}
        if(reactants.get(i).getSize()==2){ACIndex = i;}
      }
      if(containsMon && containsDi){
        
        Molecule AB = new Molecule(products.get(ABIndex));
        Molecule A = new Molecule(AB.getChildren(0));
        Molecule B = new Molecule(AB.getChildren(1));
        Molecule C = new Molecule(products.get(CIndex));
        Molecule AC = new Molecule(reactants.get(ACIndex));
        
        // Get activation energy level A1 for reaction A + B -> AB
        double cAB = (double) retrieve(AB).getCycleLength();
        double aAB = retrieve(AB).getActivity();
        double cA = (double) retrieve(A).getCycleLength();
        double aA = retrieve(A).getActivity();
        double cB = (double) retrieve(B).getCycleLength();
        double aB = retrieve(B).getActivity();
        
        A1 = (cAB + cA + cB)/2 + (aAB + aA + aB)/3 - (cA + cB);
        
        
        // Get compound activation energy level A2+A3 for reactions A + C -> AC and AC + B -> AB + C
        double cC = (double) retrieve(C).getCycleLength();
        double aC = retrieve(C).getActivity();
        double cAC = (double) retrieve(AC).getCycleLength();
        double aAC = retrieve(AC).getActivity();
        
        A2 = (cAC + cA + cC)/2 + (aAC + aA + aC)/3 - (cA + cC);
        A3 = (cAC + cB + cAB + cC)/2 + (aAC + aB + aAB + aC)/4 - (cAC + cB);
        
        if((A2+A3 < A1) && (aC !=0)){
          
          System.out.printf("\n\tMolecule %d catalyses the reaction %d + %d --> %d\n", C.getID(), A.getID(), B.getID(), AB.getID());
          System.out.printf("\t\tCyLen.\tActiv.\n");
          System.out.printf("\tA\t%.2f\t%.2f\n", cA, aA);
          System.out.printf("\tB\t%.2f\t%.2f\n", cB, aB);
          System.out.printf("\tC\t%.2f\t%.2f\n", cC, aC);
          System.out.printf("\tAB\t%.2f\t%.2f\n", cAB, aAB);
          System.out.printf("\tAC\t%.2f\t%.2f\n", cAC, aAC);
          System.out.printf("\tA1: %f\n", A1);
          System.out.printf("\tA2+A3: %f\n", A2+A3);
          
        
        }
        
        
      }
      
    
    
    
    }
  
  
  }
  
  public void collide(Molecule A, Molecule B, boolean bspar){
  
    
    ArrayList<Molecule> arrayA = new ArrayList<Molecule>();
    ArrayList<Molecule> arrayB = new ArrayList<Molecule>();
    ArrayList<Network> before = new ArrayList<Network>();
    ArrayList<Network> after = new ArrayList<Network>();
    ArrayList<Molecule> molsA = new ArrayList<Molecule>();
    ArrayList<Molecule> molsB = new ArrayList<Molecule>();
    
    Molecule Asub, Bsub;
    Network net1, net2;
    
    int[][] P;
    int[] bondingSites1, bondingSites2;
    int bondingSite1, bondingSite2;
    
    A.sortIntoArray(arrayA);
    B.sortIntoArray(arrayB);
    
    P = new int[arrayA.size()][arrayB.szie()];
    
    before.add(retrieve(A));
    before.add(retrieve(B));
    
    if(bspar){}
    
    else{
      
      for(int i=0; i<arrayA.size(); i++){
      
        for(j=0; j<arrayB.size(); j++){
        
          Asub = new Molecule(arrayA.get(i));
          Bsub = new Molecule(arrayB.get(j));
          
          net1 = retrieve(Asub);
          net2 = retrieve(Bsub);
          
          net1.setState(Asub.getState());
          net2.setState(Bsub.getState());
          
          bondingSites1 = net1.getBondingSites();
          bondingSites2 = net2.getBondingSites();
          
          //
          // Ordered choosing of bonding site.
          //
          
          for(int i = bondingSites1.length-1; i>=0; i--){
            bondingSite1 = net1.getNodes().get(bondingSites1[i]).getFilled()? bondingSite1 : bondingSites1[i];
          }
          
          for(int i = bondingSites2.length-1; i>=0; i--){
            bondingSite2 = net2.getNodes().get(bondingSites2[i]).getFilled()? bondingSite2 : bondingSites2[i];
            
          }
          
          
          molsA = A.split(arrayA.get(i).getMolID());
          A.unFlagAll();
          molsB = B.split(arrayB.get(j).getMolID());
          B.unFlagAll();
          
          
          after.addAll( moleculeArrToNetworkArr(molsA) );
          after.addAll( moleculeArrToNetworkArr(molsB) );
          after.add( bond(Asub, Bsub, bondingSite1, bondingSite2) );
          
          avActivity = 0;
          
          for(int i=0; i< before.size(); i++){
            
            avCycleLengthBefore += (double) before.get(i).getCycleLength();
            avActivity += before.get(i).getActivity();
            
          }
          
          for(int i=0; i<after.size(); i++){
            
            avCycleLengthAfter += (double) after.get(i).getCycleLength();
            avActivity += after.get(i).getActivity();
            
          }
          
          //avCycleLengthAfter = avCycleLengthAfter/after.size();
          //avCycleLengthBefore = avCycleLengthBefore/before.size();
          avActivity = avActivity/(before.size()+after.size());
          
          double exponent = ((avCycleLengthAfter + avCycleLengthBefore)/2 + avActivity)
          - avCycleLengthBefore;
          
          //System.out.printf("%f\t%f\t%f\n", avCycleLengthBefore, exponent+avCycleLengthBefore, avCycleLengthAfter);
          
          Z += Math.exp(-exponent/temperature);
          P[i][j] = Z;
          
        
        }
        
        
      }
  
      rand = Math.rand();
      
      
      for(int i=0; (i<arrayA.size())&&(!found); i++){
        for(int j=0; (i<arrayB.size())&&(!found); i++){
          if(P[i][j]/(Z+1)>rand){
            
            formBond(A, B, arrayA.get(i), arrayB.get(j), bondingSite1, bondingSite2);
            
            found = true;
          
          }
        }
      }
      
    }
  
  
  }
  
  
  public void breakUp(Molecule A){
    
    Molecule Asub;
    Network networkA;
    ArrayList<Molecule> fragmentsMol = new ArrayList<Molecule>();
    ArrayList<Network> fragments = new ArrayList<Network>();
    
    int[] bondingSites;
    int bondingSite;
    
    int[] AsubMolID;
    
    bondingSite=0;
    
    AsubMolID = A.selectRandom();
    
    Asub = new Molecule(A.getFromMoleculeID(AsubMolID, AsubMolID.length, 0));
    
    fragmentsMol = A.split(AsubMolID);
    fragmentsMol.add(Asub);
    fragments = moleculeArrToNetworkArr(fragmentsMol);
    
    A.unFlagAll();
    
    networkA = library.get(A.getID());
    networkA.setState(A.getState());
    
    double avActivity = networkA.getActivity();
    double avCycleLengthBefore = networkA.getCycleLength();
    double avCycleLengthAfter = 0;
    double P = 0;
    
    for(int i=0; i<fragments.size(); i++){
      
      avActivity += (double) fragments.get(i).getActivity();
      avCycleLengthAfter += (double) fragments.get(i).getCycleLength();
      
    }
    
    avActivity = avActivity/(fragments.size()+1);
    //avCycleLengthAfter = avCycleLengthAfter/fragments.size();
    
    
    double exponent = ((avCycleLengthAfter + avCycleLengthBefore)/2 + avActivity)
    - avCycleLengthBefore;
    
    P = Math.exp(-exponent/temperature);
    
    /*if(avActivityAfter < ActivityBefore){
      
      P = 1;
      
    }
    else{
      
      //P=0;
      //P = exp[(int)avCycleLengthAfter-(int)avCycleLengthBefore];
      P = Math.exp( ((ActivityBefore-avActivityAfter)/ (1*temperature) ) );
      //P=1;
    }*/
    
    boolean bondBroken = (Math.random()<P)? true:false;
    
    if(bondBroken){
      //System.out.printf("\t"+A.toStringf()+" ---> %d\n", fragments.size());
      numReactions++;
      ArrayList<Molecule> reactants = new ArrayList<Molecule>();
      ArrayList<Molecule> products = new ArrayList<Molecule>();
      
      reactants.add(A);
      bucket.remove(A);
      library.get(A.getID()).decrPop();
      
      for(int i=0; i<fragmentsMol.size(); i++){
        
        bucket.add(fragmentsMol.get(i));
        products.add(fragmentsMol.get(i));
        library.get(fragmentsMol.get(i).getID()).incrPop();
        
      }
      
      //products.add(bucket.get(bucket.size()-1));
      //library.get(bucket.get(bucket.size()-1).getID()).incrPop();
      
      
      Reaction R = new Reaction(reactants, products);
      int reactionIndex = reactionsUpdate(R);
      reactions.get(reactionIndex).incrCount();
      
    
    }
  
  }
  
  public static void main(String args[]){
    
    ArrayList<Molecule> reactants = new ArrayList<Molecule>();
    ArrayList<Molecule> unstable = new ArrayList<Molecule>();
    ArrayList<DataOutput> data = new ArrayList<DataOutput>();
    
    int[] arr = new int[2];
    int[] CLDistrb10 = new int[0];
    int[] CLDistrb100 = new int[0];
    
    ////////////
    //DEFAULTS//
    ////////////
    int initPop = 5000;
    int numSpecies = 50;
    reuse = false;      // Will generate new RBNs on each run
    path = "/data/library/";   //Default path to load RBNs from
    
    //
    // Command-line arguments
    //
    System.out.printf("\n\n\tInterpreting input arguments...\n\n");
    for(int i=0; i<args.length; i++){
      String str = args[i];
      if(str.startsWith("temp=")){
        temperature = Double.parseDouble(str.substring(5));
        System.out.printf("\t\tTemperature = %.2f\n",temperature);
      }
      
      else if(str.startsWith("initPop=")){
        initPop = Integer.parseInt(str.substring(8));
        System.out.printf("\t\tInitial population = %d\n",initPop);
      }
      
      else if(str.startsWith("numSpecies=")){
        numSpecies = Integer.parseInt(str.substring(11));
        System.out.printf("\t\tInitial # of species = %d\n", numSpecies );
      }
      else if (str.startsWith("fileSuffix=")){
        suffix=str.substring(11);
        System.out.printf("\t\tSaving files with suffix"+suffix);
      }
      
      else if(str.startsWith("reuse")){
      
        if(str.length()>5){
          path = str.substring(5);
          reuse = true;
          System.out.printf("\t\tRe-using previous RBNs saved in path "+path+"\n");
          
        }
        
        else{
          reuse = true;
          System.out.printf("\t\tRe-using previous RBNs saved in path "+path+"\n");
        }
      }
    }
    
    data.add( new DataOutput("/data/numReactionsPop"+suffix+".txt", true) );
    data.add( new DataOutput("/data/molecularSizeDistrb"+suffix+".txt", true) );
    data.add( new DataOutput("/data/librarySize"+suffix+".txt", true) );
    data.add( new DataOutput("/data/tempSizeDistrb_"+suffix+".txt", true) );
    data.add( new DataOutput("/data/cyclelengths"+suffix+".txt", true) );
    data.add( new DataOutput("/data/catalysis_"+suffix+".txt", true));
    
    for(int i=0; i<data.size(); i++){
      try{
        data.get(i).clearFile();
      }catch(IOException e){
        System.out.println("Could not clear file "+data.get(i).getPath()+"\n");
      }  
    }
    
    Main m = new Main(numSpecies, 10, 2, initPop);
    //Main m = new Main(numSpecies, 39, 18, 32, 100000, Double.parseDouble(suffix)/100);
    
    
    //m.populationDistribution();
    int temp=0;
    
    
    for(int i=0; i<200*initPop; i++){
      //System.out.printf("bucket size %d\n", m.getNumSpecies());
      
      if(i%50==0){
        int[] sizeDistrb = m.molecularSizeDistrb();
        double[] tempSizeDistrb = new double[sizeDistrb.length+1];
        tempSizeDistrb[0]=temperature;
        for(int k=0;k<sizeDistrb.length; k++){
          
          tempSizeDistrb[k+1]= (double)sizeDistrb[k];
        
        }
        //System.out.printf("%d\t%d\n", i, m.getLibrary().size());
        
        double tempCL=0;
        for(int k=0; k<m.getPopulation(); k++){
        
          Molecule M = m.getBucket().get(k);
          tempCL += m.retrieve(M).getCycleLength()*M.getSize();
          
        }
        
        double[] CLarr = new double[2];
        CLarr[0] = i;
        CLarr[1] = tempCL;
        
        String catStr = "";
        
        for(int k=1; k<m.getLibrary().size(); k++){
          catStr += ""+((double)m.getLibrary().get(k).getPop()/(double)m.getPopulation())+"\t";
        }
        
        String librarySize = ""+i+"\t"+m.getLibrary().size();
        
        arr[0] = m.getNumReactions()-temp;
        arr[1] = m.getPopulation();
        //System.out.printf("%d\t%d\t%d\n",i, m.getNumReactions()-temp, m.getPopulation());
        temp=m.getNumReactions();
        try{
          data.get(5).writeToFile( catStr );
          data.get(4).writeToFile( CLarr );
          data.get(3).writeToFile( tempSizeDistrb );
          data.get(2).writeToFile( librarySize );
          data.get(1).writeToFile( sizeDistrb );
          data.get(0).writeToFile( arr );
        }
        catch(IOException e){
          System.out.println("Could not write to file "+data.get(0).getPath());
        }
        
        /*temperature -= 0.001;
        if(temperature <0.001){
          temperature = 0.001;
        }*/
        
        
      }
      
      reactants = m.selectRandMols(2);
      
      Molecule A, B;
      
      A = reactants.get(0);
      B = reactants.get(1);
      
      if(! (A.getID()==0 || B.getID() == 0)){
        m.collide(A, B);
      }
      
      
      
      unstable = m.selectRandMols(2);
      
      A = unstable.get(0);
      B = unstable.get(1);
      
      //System.out.printf("\n\n|||||\t%d\n", m.getPopulation());
      if( (A.getID()!=0) && (A.getSize()!=1) ){
        m.breakUp(A);
      }
      if( (B.getID()!=0) && (B.getSize()!=1) ){
        m.breakUp(B);
      }
      
      //System.out.printf("|||||\t%d\n", m.getPopulation());
      
      
      
    }
    
  
    //m.reactionAnalysis();
    
    
    /*CLDistrb100 = m.cycleLengthDistrb();
    
    for(int i=0; (i<CLDistrb10.length || i<CLDistrb100.length); i++){
    
      if(i>=CLDistrb10.length){
        System.out.printf("0\t");
      }
      else{
        System.out.printf("%d\t", CLDistrb10[i]);
        
      }
      if(i>=CLDistrb100.length){
        System.out.printf("0\n");
      }
      else{
        System.out.printf("%d\n", CLDistrb100[i]);
        
      }
      
      
    }*/
    
    //m.populationDistribution();
    //m.molecularSizeDistrb();
    
  }
  
}












