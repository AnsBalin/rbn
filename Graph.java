import java.util.*;

public class Graph{
    
  public static final int _REACTION_ = Integer.MIN_VALUE;
  public static final int _MOLECULE_ = Integer.MIN_VALUE + 1;
  public static final int _NULL_ = Integer.MIN_VALUE + 2;
  public static final int WHITE = 0;
  public static final int BLACK = -1;
  
  
  public class GraphNode{
    
    private Object contents;
    private int type;
    private int ID;
    private String StringID;
    private ArrayList children;
    private ArrayList parents;
    private int colour; // -1: black, 0: white, 1+:gray
    
    public GraphNode(){
      contents = null;
      ID=-2;
      StringID="Null";
      children = null;
      parents = null;
      type = _NULL_;
    }
    
    public GraphNode(Molecule o){
      
      contents = new Molecule(o);
      type = _MOLECULE_;
      children = new ArrayList<GraphNode>();
      parents = new ArrayList<GraphNode>();
      colour = 0;
      ID = o.getID();
      StringID = o.getStringID();
      
    }
    public GraphNode(Reaction o){
      
      contents = new Reaction(o);
      type = _REACTION_;
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
    
    public int getType(){
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
      
      addChild(GN);
      GN.addParent(this);
      
    }
  } 
  
  private ArrayList<GraphNode> graphNodes;
  private int numCycles;
  private int visited=0;
  
  public Graph(){
    
    graphNodes = new ArrayList<GraphNode>();
    graphNodes.add(new GraphNode());
    numCycles = 0;
    
  }
  
  public int checkReactionID(int IDCheck){
    int index=-1;
    for(int i=0; i<graphNodes.size(); i++){
      if( (graphNodes.get(i).getID()==IDCheck) && (graphNodes.get(i).getType() == _REACTION_) ){
        index=i;
      
      }
    }
    return index;
  }
  
  public int checkMoleculeID(int IDCheck){
    int index=-1;
    for(int i=0; i<graphNodes.size(); i++){
      if( (graphNodes.get(i).getID()==IDCheck) && (graphNodes.get(i).getType() == _MOLECULE_) ){
        index=i;
        
      }
    }
    return index;
  }
  
  public int getNumCycles(){
    return numCycles;
  }
  
  public ArrayList<GraphNode> getGraphNodes(){return graphNodes;}
  
  public void add(Reaction R, Main m){
    
    int IDCheck;
    int index = checkReactionID(R.getIntID());  // index = index of R in Graph. checks if R is in Graph, if it is, returns index.
   
    if(index==-1){                               // ie. if R isnt in Graph
      index=graphNodes.size();
      graphNodes.add(new GraphNode(new Reaction(R)));  // NEW INSTANCE OF R IS ADDED
    }
    
    for(int i=0; i<R.getReactants().length; i++){
      
      IDCheck = R.getReactants()[i];
      int reactantIndex = checkMoleculeID(IDCheck);   // Checks to see if reactant is already in graph.
      
      if(reactantIndex==-1){    // ie. if reactant isnt in graph
        reactantIndex = graphNodes.size();
        graphNodes.add(new GraphNode(new Molecule(m.getLibrary().get(IDCheck).getMol())));   // NEW INSTANCE OF MOLECULE IS ADDED
      }
      
      //System.out.printf("(%d) is a parent of (%d)\n", reactantIndex, index);
      graphNodes.get(reactantIndex).connect(graphNodes.get(index));      //Makes directed edge
      
    }
    
    for(int i=0; i<R.getProducts().length; i++){
      
      IDCheck = R.getProducts()[i];
      int productIndex = checkMoleculeID(IDCheck);
      
      if(productIndex==-1){
        productIndex = graphNodes.size();
        graphNodes.add(new GraphNode(new Molecule(m.getLibrary().get(IDCheck).getMol())));
      }
      //if(index==1){System.out.printf("\n(1) has another child` ie %d child!!!!\n", R.getProducts().length);}
      //if(index==0){System.out.printf("!!!!!!!!\n");}
      //System.out.printf("(%d) is a parent of (%d)\n", index, productIndex);
      graphNodes.get(index).connect(graphNodes.get(productIndex));
      
    }
  
  }
  
  // Performs a depth-first-search of this graph
  public void dfs(){
    int depth;
    for(int i=0; i<graphNodes.size(); i++){
      //if(graphNodes.get(i).getChildren()!=null){System.out.printf("// (%d) has %d children\n", i, graphNodes.get(i).getChildren().size());}
      if((graphNodes.get(i).getColour()==0)&&(graphNodes.get(i).getChildren()!=null)){
        depth = 1;
        dfsVisit(graphNodes.get(i), depth);
        
      }
    }
  
  }
  
  public void dfsVisit(GraphNode GN, int depth){
    GN.colour(depth);
    //System.out.printf("\n");
    //for(int z=0; z<depth; z++){System.out.printf("\t");}
    //System.out.printf("(%d) - depth = %d\n", graphNodes.indexOf(GN), depth);
    
    ArrayList<GraphNode> children = GN.getChildren();
    
    for(int i=0; i<children.size(); i++){
      
      if((children.get(i).getColour()==0)&&(children.get(i).getChildren()!=null)){
        //for(int z=0; z<depth; z++){System.out.printf("\t");}
        //System.out.printf("Visting child %d of (%d)\n", i, graphNodes.indexOf(GN));
       
        dfsVisit(children.get(i), depth+1);
      }
      //System.out.printf("!!%d\n", children.get(i).getColour());
      if(children.get(i).getColour()>0){
        numCycles++;
        //for(int z=0; z<depth; z++){System.out.printf("\t");}
        //System.out.printf("Child %d of (%d) is (%d)\n", i, graphNodes.indexOf(GN), graphNodes.indexOf(children.get(i)));
        //for(int z=0; z<depth; z++){System.out.printf("\t");}
        System.out.printf("%d\n", GN.getColour()- children.get(i).getColour()+1);
        
      }
      
      else{} 
    
    }
    
    GN.colour(-1);
    visited++;
    //Main.progress("Visiting nodes", visited, graphNodes.size());
    
  }  
  
  public boolean filter(Reaction R, int option, int arg ){
    
    boolean result = false;
    
    switch(option){
    
      case 0:
        if(R.getCount() > arg){
          result = true;
        }
        break;
    
    
    
    }
    
    return result;
  
  }
 
  
  public static void main(String args[]){
    
    int pop = 1000;
    Main m = new Main(5, 10, 2, pop);
    m.getReactions().add(new Reaction());
    m.temperature = 1000;
    
    for(int i=0; i<100*pop; i++){
      if(i%100==0){
        m.progress("Running...", i, 100*pop);
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
        m.collide(A, B, false);
      }
      
      for(int j=0; j<m.getBucket().size(); j++){
        m.getLibrary().get( m.getBucket().get(j).getID() ).setMol(new Molecule(  m.getBucket().get(j) ) );
      }
    }
    
    m.progress("Running...", 10000,10000);
    
    System.out.printf("\n");
    Graph g = new Graph();
    System.out.printf("# Reactions = %d", m.getReactions().size());
    for(int i=1; i<m.getReactions().size(); i++){
      //if(i%10==0){m.progress("Graph build", i, m.getReactions().size());}
      Reaction X = new Reaction();
     
      
      
      if(g.filter(m.getReactions().get(i), 0, 10)){
        g.add(m.getReactions().get(i), m);
      }
      
    }
    //m.progress("Graph build", m.getReactions().size(), m.getReactions().size());
    System.out.printf("\n");
    System.out.printf("\tPerforming DFS...\n");
    g.dfs();
    System.out.printf("\tNumber of cycles: %d\n", g.getNumCycles());
    
  }

}








