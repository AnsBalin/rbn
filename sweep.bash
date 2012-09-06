#!/bin/bash/

#   This script performs repeated executions of Main with variable input arguments
#   It utilises /usr/bin/osascript script on OS X to open a new Terminal window
#+  for each execution.
#   To port this into linux, replace the AppleScript code with the command:
#+  gnome-terminal -e "insert script here"

echo Begin

n=0
temperature=1
while [ $n -le 50 ]
do
  echo "======================================================="
  var1=$(ls -la | awk '{print$9, "\t", $5}' | grep 'log.txt' | awk '{print $2}')
  #size of log file
  
  temp="$temperature"
  suffix="$temperature"
  echo Suffix and temperature are: $suffix, $temp


  ####### AppleScript Code here #######
  button=`/usr/bin/osascript << EOT
  tell app "Terminal" 
    do script "
      cd ~/Dropbox/rbn
      java Main reuse scriptTemp numSpecies=10 temp=$temp fileSuffix=$suffix
      exit 0
      "
  end tell
  EOT`
  ######## End of AppleScript ########


  var2=$(ls -la | awk '{print$9, "\t", $5}' | grep 'log.txt' | awk '{print $2}')  
  #size of log file
  
  while [ $var1 -eq $var2 ]
  do  #while the log file remains the same size, ie. while Main is still running
    sleep 2;
    var2=$(ls -la | awk '{print$9, "\t", $5}' | grep 'log.txt' | awk '{print $2}');
    
  done #when log file is updated, ie. Main has terminated

  echo "Run number $n"
  echo "temperature = temp"
  n=`expr $n + 1`
  temperature=`expr $temperature + 1`

done

echo End