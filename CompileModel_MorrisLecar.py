import os

print("1. Running ./model2diffit models/MorrisLecar.mdf")
os.system('./model2diffit ./models/MorrisLecar.mdf') 
print("2. Copying MorrisLecar.cc to model.cc")
os.system("rm ./model.cc") #have to do this first b/c of a permissions issue with my account (rickymorse)
os.system("cp ./MorrisLecar.cc ./model.cc")
print("3. Running make")
os.system("make clean; make all")
