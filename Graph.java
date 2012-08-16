import java.util.*;

public class Graph{
  
  public class GraphNode{
    
    private Object contents;
    private Class type;
    private int ID;
    private String StringID;
    private ArrayList children;
    private ArrayList parents;
    private int colour; // -1: black, 0: white, 1+:grey
    
    public GraphNode(){
      contents = null;
      ID=-2;
      StringID="Null";
      children = null;
      parents = null;
    }
    
    public GraphNode(Molecule o){
      
      contents = new Molecule(o);
      type = o.getClass();
      children = new ArrayList<GraphNode>();
      parents = new ArrayList<GraphNode>();
      colour = 0;
      ID = o.getID();
      StringID = o.getStringID();
      
    }
    public GraphNode(Reaction o){
      
      contents = new Reaction(o);
      type = o.getClass();
      children = new ArrayList<GraphNode>();
      parents = new ArrayList<GraphNode>();
      colour = 0;
      ID = o.getIntID();
      StringID = o.getID();
      
    }
    
    public int getColour(){
      return colour;
    }
    
    public int getID(){
      return ID;
    }
    
    public Class getType(){
      return type;
    }
    
    public void print(){
    
      System.out.println(""+StringID);
    }
    
    public ArrayList<GraphNode> getChildren(){
      return children;
    }
    
    public void colour(int i){
    
      colour=i;
      
    }
    
    public void addChild(GraphNode GNChild){
      
      children.add(GNChild);
      
    }
    
    public void addParent(GraphNode GNParent){
    
      parents.add(GNParent);
    }
    
    //Forms directed edge from this to GN
    public void connect(GraphNode GN){
      
      this.addChild(GN);
      GN.addParent(this);
      
    }
  } 
  
  private class Traverse{
    
    
    
  }
  
  private ArrayList<GraphNode> graphNodes;
  private int numCycles;
  
  public Graph(){
    
    graphNodes = new ArrayList<GraphNode>();
    numCycles = 0;
    
  }
  
  public int checkID(int IDCheck, Object o){
    int index=-1;
    for(int i=0; i<graphNodes.size(); i++){
    
      if( (graphNodes.get(i).getID()==IDCheck) && (graphNodes.get(i).getType()==o.getClass())){
        index=1;
        break;
      }
    
    }
    return index;
    
  }
  
  public int getNumCycles(){
    return numCycles;
  }
  
  public void add(Reaction R, Main m){
    
    Molecule M = new Molecule();
    
    int IDCheck;
    int index = checkID(R.getIntID(), R);
    graphNodes.add(new GraphNode());
    
    if(index==-1){
      index=graphNodes.size();
      graphNodes.add(new GraphNode(new Reaction(R)));
    }
    
    
        
    for(int i=0; i<R.getReactants().length; i++){
      
      IDCheck = R.getReactants()[i];
      int reactantIndex = checkID(IDCheck, M);
      
      if(reactantIndex==-1){
        
        for(int k=0; k<m.getLibrary().size(); k++){
          
          if(k==R.getReactants()[i]){
            reactantIndex = graphNodes.size();
            graphNodes.add(new GraphNode(new Molecule(m.getLibrary().get(k).getMol())));
            break;
          
          }
        }
      }
      
      graphNodes.get(reactantIndex).connect(graphNodes.get(index));
      
    }
    for(int i=0; i<R.getProducts().length; i++){
      
      IDCheck = R.getProducts()[i];
      int productIndex = checkID(IDCheck, M);
      
      if(productIndex==-1){
        
        for(int k=1; k<m.getLibrary().size(); k++){
          
          if(k==R.getProducts()[i]){
            productIndex = graphNodes.size();
            graphNodes.add(new GraphNode(new Molecule(m.getLibrary().get(k).getMol())));
            break;
            
          }
        }
      }
      
      graphNodes.get(index).connect(graphNodes.get(productIndex));
      
    }
  
  }
  
  // Performs a depth-first-search of this graph
  public void dfs(){
  
    for(int i=1; i<graphNodes.size(); i++){
      if((graphNodes.get(i).getColour()==0)&&(graphNodes.get(i).getChildren()!=null)){
        dfsVisit(graphNodes.get(i), 1);
        
      }
    }
  
  }
  
  public void dfsVisit(GraphNode GN, int n){
  
    GN.colour(n);
    //System.out.printf("Printing %d\n", n);
    //GN.print();
    ArrayList<GraphNode> children = GN.getChildren();
    
    for(int i=0; i<children.size(); i++){
    
      if(children.get(i).getColour()>0){
        numCycles++;
        System.out.printf("%d\t%d\n", n, n-children.get(i).getColour()+1);
      }
      if((children.get(i).getColour()==0)&&(children.get(i).getChildren()!=null)){
        dfsVisit(children.get(i), n+1);
      }
      else{}
    
    }
    
    GN.colour(-1);
    
  }
  
  public static void main(String args[]){
  
    Main m = new Main(3, 10, 2, 10);
    
    for(int i=0; i<100; i++){
      if(i%100==0){
        //System.out.printf("%d\n", i);
      }
      for(int j=0; j<m.getBucket().size(); j++){
        m.getLibrary().get( m.getBucket().get(j).getID() ).setMol(  m.getBucket().get(j) );
      }
      if(m.getBucket().size()==2){
        break;
      }

      ArrayList<Molecule> reactants;
      reactants = m.selectRandMols(2);
      
      Molecule A, B;
      
      A = reactants.get(0);
      B = reactants.get(1);
      
      if(! (A.getID()==0 || B.getID() == 0)){
        m.collide(A, B);
      }
    }
    
    Graph g = new Graph();
    System.out.printf("Building graph...\n");
    for(int i=0; i<m.getReactions().size(); i++){
      //if(i%1==0){System.out.printf("\t%d\n", i );}
      if(m.getReactions().get(i).getCount() > 0){
        g.add(m.getReactions().get(i), m);
      }
      
    }
    System.out.printf("Performing DFS...\n");
    g.dfs();
    System.out.printf("Number of cycles: %d\n", g.getNumCycles());
    
  
  
  }



}








