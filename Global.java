public class Global{
  
  public static boolean[][] lookUp = {
  //
  //INPUTS AB can be encoded as ints
  //
  //LOGIC     FF      FT      TF      TT
  //BINARY    00      01      10      11
  //INDEX     0       1       2       3
  {false, false,  false,  false},
  {false, false,  false,  true},
  {false, false,  true,   false},
  {false, false,  true,   true},
  
  {false, true,   false,  false},
  {false, true,   false,  true},
  {false, true,   true,   false},
  {false, true,   true,   true},
  
  {true,  false,  false,  false},
  {true,  false,  false,  true},
  {true,  false,  true,   false},
  {true,  false,  true,   true},
  
  {true,  true,   false,  false},
  {true,  true,   false,  true},
  {true,  true,   true,   false},
  {true,  true,   true,   true},
  };
  
  public static String[] chars = {"0","A","B","C","D","E","F","G","H",
  "I","J","K","L","M","N","O","P","Q",
  "R","S","T","U","V","W","X","Y","Z",
  "a", "b", "c", "d", "e", "f", "g", "h",
  "i", "j", "k", "l", "m", "n", "p", "q",
  "r", "s", "t", "u", "v", "w", "x", "y", "z",
  "0", "1", "2", "3", "4", "5", "6", "7", "8",
  "9", 
  
  };
  
  public static void main(String args[]){
    
    for(int x=0;x<16;x++){
      for(int y=0;y<4; y++){
        System.out.printf("%b\t", lookUp[x][y]);
      }
      System.out.printf("\n");
      
    }
  }
  
}