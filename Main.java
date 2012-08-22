
import java.util.*;
import java.io.*;
import java.lang.*;

public class Main {
  
  private ArrayList<Network> library = new ArrayList<Network>();
  private ArrayList<Molecule> bucket = new ArrayList<Molecule>();
  private ArrayList<Reaction> reactions = new ArrayList<Reaction>();
  private int numReactions;
  private static int catalysisProductIndex=-1;
  private static boolean foundCatalysisProduct  = false;
  private static double temperature = 1;
  private static boolean reuse = false;
  private static String path = "/data/library/";
  private static String suffix = "";
  
  
  // Matrix 1 - numOutputs
  // Matrix 2 - Function #
  private ArrayList<ArrayList<ArrayList<Double>>> matrix1 = new ArrayList<ArrayList<ArrayList<Double>>>();
  private ArrayList<ArrayList<ArrayList<Double>>> matrix2 = new ArrayList<ArrayList<ArrayList<Double>>>();
  
  private int product;
  
  // Constructor adds null network to library, with index 0
  public Main(int numSpecies, int numNodes, int numBondingSites, int numMolecules){
    
    System.out.printf("\n\tPopulating library...\n");
    populateLibrary(numSpecies, numNodes, numBondingSites);
    System.out.printf("\tLibrary successfully populated with:\n\t\t%d Species\n\t\t%d nodes each\n\t\t%d bonding sites each", numSpecies, numNodes, numBondingSites);
    
    System.out.printf("\n\n\tPopulating world...\n");
    populateWorld(numMolecules, numSpecies);
    System.out.printf("\tWorld successfully populated with %d atoms.\n\n", numMolecules);
    
    numReactions=0;
    //stats();
    
    for(int i=0; i<11; i++){
      matrix1.add(new ArrayList<ArrayList<Double>>());
      for(int j=0; j<11; j++){
        matrix1.get(i).add(new ArrayList<Double>());
      }
    }
    
    for(int i=0; i<16; i++){
      matrix2.add(new ArrayList<ArrayList<Double>>());
      for(int j=0; j<16; j++){
        matrix2.get(i).add(new ArrayList<Double>());
      }
    }
    
    
    
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
  public ArrayList getMatrix1(){return matrix1;}
  public ArrayList getMatrix2(){return matrix2;}
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
  
  public void progress(String process, int i, int N){
    
    int k=0;
    
    System.out.printf("\r");
    String outputString;
    
    outputString = "\t";
    for(k=0; k<=i; k+=N/50){
      
      outputString+="█";
    
    }
    while(k<=N){
    
      outputString+=" ";
      k+=N/50;
    }
    System.out.printf("%20s "+outputString+"%d%%", process, (int)(100*((double)i/(double)N)));
    
  
  }
  
  // If generateNew = 0     Generate new species
  // If generateNew != 0    Use species from previous run
  public void populateLibrary(int numSpecies, int numNodes, int numBondingSites){
    
    if(!reuse){
      library.add(new Network(0));
    
      for(int i=1; i<numSpecies; i++){
        library.add(new Network(numNodes, i, numBondingSites));
        if(i % (numSpecies/50)==0){
          progress("Adding...", i, numSpecies);
        }
      
      }
      progress("Adding...", numSpecies, numSpecies);
      System.out.printf("\n");
      for(int i=0; i<numSpecies; i++){
        DataOutput out = new DataOutput("/data/library/"+i+".txt");
        if(i % (numSpecies/50)==0){
          progress("Saving atoms", i, numSpecies);
        }
        try{
          out.writeToFile(library.get(i));
        }catch(IOException e){System.out.println("Could not write file.");}
      }
      progress("Adding...", numSpecies, numSpecies);
      System.out.printf("\n");
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
      if(i % (numMolecules/50) == 0){progress("Adding...", i, numMolecules);}
      rand = (int) (Math.random()*(double)(numSpecies-1)) + 1;
      bucket.add(new Molecule(rand , nullMolecule, nullMolecule));
      bucket.get(i).calculateMoleculeID(0,0, new int[0]);
      bucket.get(i).setState(library.get(rand).getState());
      //String str = bucket.get(i).toString();
      //System.out.printf(str);
    }
    progress("Adding...", numMolecules, numMolecules);
    
    System.out.printf("\n");
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

  
  public void collide(Molecule A, Molecule B, boolean bspar){
  
    ArrayList<Network> before = new ArrayList<Network>();
    ArrayList<Molecule> arrayA = new ArrayList<Molecule>();
    ArrayList<Molecule> arrayB = new ArrayList<Molecule>();
    boolean found=false;
    Molecule Asub, Bsub;
    Network net1, net2;
    double[][] P;
    int[] bondingSites1, bondingSites2;
    int bondingSite1, bondingSite2;
    double avActivity, avCycleLengthBefore, avCycleLengthAfter, avActivityAfter, Z;
    bondingSite1=0;
    bondingSite2=0;
    Z=0;
    
    A.sortIntoArray(arrayA);
    B.sortIntoArray(arrayB);
    
    P = new double[arrayA.size()][arrayB.size()];
    
    before.add(retrieve(A));
    before.add(retrieve(B));
      
    for(int i=0; i<arrayA.size(); i++){
    
      for(int j=0; j<arrayB.size(); j++){
        
        
        ArrayList<Molecule> molsA = new ArrayList<Molecule>();
        ArrayList<Molecule> molsB = new ArrayList<Molecule>();
        
        ArrayList<Network> after = new ArrayList<Network>();
        
        avCycleLengthAfter=0;
        avCycleLengthBefore=0;
        avActivityAfter=0;
        
      
        //Asub = new Molecule(arrayA.get(i));
        //Bsub = new Molecule(arrayB.get(j));
        Asub = arrayA.get(i);
        Bsub = arrayB.get(j);
        
        net1 = retrieve(Asub);
        net2 = retrieve(Bsub);
        
        net1.setState(Asub.getState());
        net2.setState(Bsub.getState());
        
        bondingSites1 = net1.getBondingSites();
        bondingSites2 = net2.getBondingSites();
        
    
        
        //
        // Ordered choosing of bonding site.
        //
        
        //System.out.printf("%d\t%d\n", bondingSites1.length, bondingSites2.length);
        
        for(int k = bondingSites1.length-1; k>=0; k--){
          //System.out.printf("\t\t%d\n", bondingSite1);
          bondingSite1 = net1.getNodes().get(bondingSites1[k]).getFilled()? bondingSite1 : bondingSites1[k];
        }
        
        for(int k = bondingSites2.length-1; k>=0; k--){
          // System.out.printf("\t\t%d\n", k);
          bondingSite2 = net2.getNodes().get(bondingSites2[k]).getFilled()? bondingSite2 : bondingSites2[k];
          
        }
        
        
        molsA = A.split(arrayA.get(i).getMolID());
        A.unFlagAll();
        molsB = B.split(arrayB.get(j).getMolID());
        B.unFlagAll();
        
        
        after.addAll( moleculeArrToNetworkArr(molsA) );
        after.addAll( moleculeArrToNetworkArr(molsB) );
        after.add( bond(Asub, Bsub, bondingSite1, bondingSite2) );
        
        avActivity = 0;
        
        for(int k=0; k< before.size(); k++){
          
          avCycleLengthBefore += (double) before.get(k).getCycleLength();
          avActivity += before.get(k).getActivity();
          
        }
        
        for(int k=0; k<after.size(); k++){
          
          avCycleLengthAfter += (double) after.get(k).getCycleLength();
          avActivity += after.get(k).getActivity();
          
        }
        
        //avCycleLengthAfter = avCycleLengthAfter/after.size();
        //avCycleLengthBefore = avCycleLengthBefore/before.size();
        avActivity = avActivity/(before.size()+after.size());
        
        double exponent = ((avCycleLengthAfter + avCycleLengthBefore)/2 + avActivity)
        - avCycleLengthBefore;
        
        //System.out.printf("%f\t%f\t%f\n", avCycleLengthBefore, exponent+avCycleLengthBefore, avCycleLengthAfter);
        
        
        P[i][j] = Math.exp(-exponent/temperature);
        Z += P[i][j];
        //System.out.printf("\tbondingSite1:   %d\n", bondingSite1);
        //System.out.printf("\tbondingSite2:   %d\n", bondingSite2);
        
        
        matrix2.get(net1.getNodes().get(bondingSite1).getFunction()).get(net2.getNodes().get(bondingSite2).getFunction()).add(new Double(P[i][j]/(P[i][j]+1)));
        matrix1.get(net1.getNodes().get(bondingSite1).getNumOutputs()).get(net2.getNodes().get(bondingSite2).getNumOutputs()).add(new Double(P[i][j]/(P[i][j]+1)));
        
      }
      
                        
    }
  
  
  }
  
  public void save(ArrayList<DataOutput> data, ArrayList<ArrayList<ArrayList<Double>>> matrix, boolean range, int fileNo){
  
    int size=matrix.size();
    
    for(int i=0; i<size; i++){
      
      double[] matrixZ = new double[size];
      
      for(int j=0; j<size; j++){
        //System.out.printf("[%d, %d]\t", i, j);
        Object[] arrList = matrix.get(i).get(j).toArray();
        double dbl=0;
        
        if(arrList.length>0){
          Arrays.sort(arrList);
          if(range){
            dbl = new Double(arrList[(int)(0.75*(double)arrList.length)].toString()) - new Double (arrList[(int)(0.25*(double)arrList.length)].toString());
          }
          else{
            dbl = new Double(arrList[(int)(arrList.length/2)].toString());
          }
          
        }
        else{}
        matrixZ[j] = dbl;
        
      }
      
      try{
        
        data.get(fileNo).writeToFile(matrixZ);
        
      }catch(Exception e){}
    }
  }
    
  public void cleanUp(){
    
    //System.out.printf("\tcleanUp() initiated...\n");
    int trash=0;
    
    for(int i=0; i<11; i++){
      for(int j=0; j<11; j++){
        
        int matrixSize = matrix1.get(i).get(j).size();
        if(matrixSize>10000){
          matrix1.get(i).get(j).subList(10000, matrixSize-1).clear();
          trash += matrixSize - 10000;
          //System.out.printf("\n\r\tcleanUp() collected %d items.", trash);
        }
      }
    
    }
  
    for(int i=0; i<16; i++){
      for(int j=0; j<16; j++){
        
        int matrixSize = matrix2.get(i).get(j).size();
        if(matrixSize > 10000){
          matrix2.get(i).get(j).subList(10000, matrixSize-1).clear();
          trash += matrixSize - 10000;
          //System.out.printf("\r\tcleanUp() collected %d items.", trash);
      
        }
      }
    }
    
    System.out.printf("\tcleanUp() collected %d items.", trash);
    
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
    
    
    data.add( new DataOutput("/data/numOutputsRange"+suffix+".dat", true));
    data.add( new DataOutput("/data/numOutputsMean"+suffix+".dat", true));
    data.add( new DataOutput("/data/functionRange"+suffix+".dat", true));
    data.add( new DataOutput("/data/functionMean"+suffix+".dat", true));
    
    for(int i=0; i<data.size(); i++){
      try{
        data.get(i).clearFile();
      }catch(IOException e){
        System.out.println("Could not clear file "+data.get(i).getPath()+"\n");
      }  
    }
    
    Main m = new Main(numSpecies, 10, 2, initPop);
    
    ArrayList<ArrayList<Double>> matrixTemp = m.getMatrix1();
    int ii=0;
    while(ii<100*initPop){
      
      if(ii%10000==0){m.progress("\tColliding...", ii, 100*initPop);m.cleanUp();}
      
      reactants = m.selectRandMols(2);
      
      Molecule A, B;
      
      A = reactants.get(0);
      B = reactants.get(1);
      
      if(! (A.getID()==0 || B.getID() == 0)){
        m.collide(A, B, false);
      }
      ii++;
    }
    System.out.println("\n\tSaving files...");
    m.save(data, m.getMatrix1(), true, 0);
    m.save(data, m.getMatrix1(), false, 1);
    m.save(data, m.getMatrix2(), true, 2);
    m.save(data, m.getMatrix2(), false, 3);
    
    System.out.printf("\tCompleted successfully\n\n");
    
    /*for(int i=0; i<16; i++){
      
      int foo = matrix.get(i).size();
      double[] probability = new double[min];
      if(foo<min){break;}
      for(int j=0; j<min; j++){
        probability[j] = matrix.get(i).get(j).doubleValue();
        System.out.printf("(%d, %d)\t%.4f\n", i, j, probability[j]);
        
      }
      
      try{
        
        data.get(0).writeToFile(probability);
      }catch(IOException e){
        System.out.printf("Could not write to file.\n");
        break;
      }
      
      
    }*/
    
    
    /*System.out.println("size 1: "+matrix.get(0).size());
    System.out.println("size 2: "+matrix.get(1).size());
    System.out.println("size 3: "+matrix.get(2).size());
    System.out.println("size 4: "+matrix.get(3).size());
    System.out.println("size 5: "+matrix.get(4).size());
    System.out.println("size 6: "+matrix.get(5).size());
    System.out.println("size 7: "+matrix.get(6).size());
    System.out.println("size 8: "+matrix.get(7).size());*/
    
    
    

  }
  
}












