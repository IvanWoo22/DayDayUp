#!/bin/bash 

function read_dir(){
  for file in `ls $1`
  do
    if [ -d $1"/"$file ]  ##Here blank is very important.
    then 
      read_dir $1"/"$file 
    else
      echo $1"/"$file
    fi
  done
}


read_dir $1
