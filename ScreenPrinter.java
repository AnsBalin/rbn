/*
 
 Printer for dynamic command line interface.
 
 See:
 http://www.termsys.demon.co.uk/vtansi.htm
 for full explanation of escape sequences and colour codes etc.
 
*/

import java.lang.*;
import java.util.*;


public class ScreenPrinter{
  
  public class Box{
    
    public boolean frozen;
    private int width;
    private int height;
    private int[] position = new int[2];
    
    public Box(int width, int height, int row, int column){
      
      this.width = width;
      this.height = height;
      position[0] = row;
      position[1] = column;
      
    }
    
    public int[] getPosition(){return position;}
    public int getHeight(){return height;}
    public int getWidth(){return width;}
    
    public void printLine(String line, int lineNumber){
      
      cursor(position[0]+lineNumber-1, position[1]);
      if(1==1){
        
        System.out.print(line);
        for(int i=0; i< width-line.length(); i++){
          System.out.print(" ");
        }
        
      }
      else{
        
        line = line.substring(0, width);
        System.out.print(line);
      }
    
    }
  
  
  }
  
  public static final int width = 80;     // #columns
  public static int height;               // #rows
  public static final char esc = 0x1B;    // ASCII escape character
  private int numMolecules;
  private int numSpecies;
  private double temperature;
  private int time;
  private int generations;
  private ArrayList<Box> boxes;
  public Box debug;
  
  public ScreenPrinter(){
    
    
    boxes = new ArrayList<Box>();
    
    
    // This can be edited by user to produce boxes of different sizes
    // and positions. 
    boxes.add(new Box(30, 4, 3, 10)); // "Populating library..."
    boxes.add(new Box(30, 2, 3, 40)); // "Populating bucket..."
    boxes.add(new Box(30, 3, 8, 10)); // temperature, library size
    boxes.add(new Box(30, 3, 8, 40)); // collisions, population, # reactions
    boxes.add(new Box(60, 4, 12, 10)); // progress bar
    
    height = boxes.get(boxes.size()-1).getPosition()[0] + boxes.get(boxes.size()-1).getHeight() +1;  
    
    debug = new Box(80, 10, height+1, 1);
    
    clear();
    makeTitle("RBN WORLD");
    makeBorders();
    
  }
  
  // Clears screen by reseting terminal
  public void clear(){
    System.out.print(String.format("%cc", esc));
  }
  
  // Set formatting attributes.
  public void format( boolean underscore, int foreground, int background ){
  
    if(underscore){
      System.out.print(String.format("%c[4;%d;%dm", esc, foreground, background));
    }
    
    else{
      System.out.print(String.format("%c[%d;%dm", esc, foreground, background));
    }
  }
  
  public String formatString( boolean underscore, int foreground, int background ){
    
    if(underscore){
      return String.format("%c[4;%d;%dm", esc, foreground, background);
    }
    
    else{
      return String.format("%c[%d;%dm", esc, foreground, background);
    }
  }

  
  // Reset all attributes
  public void reformat(){
    System.out.print(String.format("%c[0m", esc));
  }
  
  // Reset all attributes
  public String reformatString(){
    return String.format("%c[0m", esc);
  }

  
  // Positions cursor at given row/column
  public void cursor( int row, int column ){
    System.out.print(String.format("%c[%d;%df", esc, row, column));
  }
  
  public void makeTitle(String title){
    
    int padding = (width-title.length())/2;
    cursor(1, 1);
    reformat();
    for (int i=0; i<padding; i++){
      
      System.out.printf("═");
    
    }
    
    format(false, 31, 47);
    System.out.printf(title);
    reformat();
    for(int i=padding+title.length(); i<width; i++){
      System.out.printf("═");
    }
  
  }
  
  public void makeBorders(){
  
    int column = 1;
    for(int i=2; i<height; i++){
      cursor(i, column);
      System.out.print("║");
    
    }
    column = width;
    for(int i=2; i<height; i++){
      cursor(i, column);
      System.out.print("║");
      
    }
    for(int i=1; i<=width; i++){
      cursor(height, i);
      System.out.print("═");
    }
    cursor(1,1);
    System.out.print("╔");
    cursor(1,width);
    System.out.print("╗");
    cursor(height,1);
    System.out.print("╚");
    cursor(height,width);
    System.out.print("╝");
  }
  
  // String can begin at positions 1-8 (columns 0 to 80). This method pads 
  // a string starting at a given position and pads it on either side
  public void padString(String str, int position){
    String padString = "";
    for(int i=1; i<(position-1)*10; i++){
      padString += " ";
    }
    
    padString += str;
    
  }
  
  public void printLine(String str){
    
  }
  
  // Main printer function
  public void print(String line, int boxNumber, int lineNumber){
  
    boxes.get(boxNumber).printLine(line, lineNumber);
    
  }
  
  public String progress( int i, int N){
    
    double k=0;
    
    //System.out.printf("\r");
    String outputString;
    
    outputString = "";
    
    debug.printLine(String.format("%d", outputString.length()),7);
    int temp=0;
    for(k=0; k<=(double)i; k+=(double)N/50){
      
      outputString+="█";
      temp++;
    }
    
    while(k<=N){
      
      outputString+="-";
      k+=(double)N/50;
      temp++;
    }
    
    debug.printLine(String.format("temp = %d", temp), 6);

    debug.printLine(String.format("i = %d", i), 1);
    debug.printLine(String.format("N = %d", N), 2);
    debug.printLine(String.format("i/N = %d", (int)(100*((double)i/(double)N))), 3);
    debug.printLine(outputString, 4);
    debug.printLine(String.format("%d", outputString.length()), 5);
    
    String one = "-----";
    String two = "█████";
    
    debug.printLine(String.format("length of one = %d\nlength of two= %d", one.codePointCount(0, one.length()), two.codePointCount(0, two.length())), 8);
    return String.format(outputString+" %d%%", (int)(100*((double)i/(double)N)));
    
    
  }

  
  public static void main(String args[]){
  
    ScreenPrinter sp = new ScreenPrinter();
    
    sp.cursor(height, width);
    System.out.printf("\n");
    
  }
  
}







