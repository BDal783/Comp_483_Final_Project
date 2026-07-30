
from pathlib import Path 
import torch
import torch.nn as nn
import torch.optim as optim 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import time
import datetime
import threading
import os
import math
import random
from Bio import SeqIO
from Bio import Seq
import shutil

#This code makes arrays for refrenced later. I dont think I ended up using the seccond one tho
frames = ["",".","..","..."]
al=["a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z"]

#This code sets up some variables to be used later
done=0
ready=0
space="     "
past8=""
seq=0
a=[""]
out=[]
rtnd=None

#These next two functions create empty arrays of type None of type String
def initArrayStr(amt):
    ary=[]
    for x in range(amt):
        ary.append("")
    return ary

def initArrayObj(amt):
    ary=[]
    for x in range(amt):
        ary.append(None)
    return ary

#This runs the loading screen you see
def cuteLoadingScreen():
    count=0
    frames = ["",".","..","..."]
    while(ready==0):
        count+=1
        print(f"\rLoading{frames[count%4]:<4}",end="",flush=True)
        time.sleep(1/3)
        if count>=4:
            count=0
    print(f"\r",end="",flush=True)    

#This function is something I stole from stack overflow that finds the nth occurence of somthing in a string
def findNth(haystack, needle, n):
    parts= haystack.split(needle, n+1)
    if len(parts)<=n+1:
        return -1
    return len(haystack)-len(parts[-1])-len(needle)

#This returns a array of all indexes of ">" character
def count(fileLoc,chr=">"):
    fileName=fileLoc
    with open (fileName,"r") as text:
        file = text.read()
        amt = file.count(chr)
        rtn=[amt,[]]
        for x in range(amt):
            file2=file
            pos=findNth(file,">",x)
            rtn[1].append(pos)

        rtn[1].append(len(file))

        return rtn

#This code counts the nucleotides in part of a file
def countSeq(fileLoc,posIn,next,idx):

    global done, started,spli,out,seq,done,ready,past8
    fileName=fileLoc
    with open (fileName,"r") as f:   
        rtn=""
        seq=0
        f.read(posIn+3)
        next_char=f.read(0)
        past8=""

        pos=posIn+2

        while seq==0:
            pos+=1

            if len(past8)>=8:
                past8=past8[1:9]
            past8+=next_char
            if (past8[-1:]=="\n"):
                seq=1
            else:
                next_char=f.read(1)
        
        rtn=f.read(next-pos)
    rtn=rtn.replace(" ","")
    rtn=rtn.replace("\n","")
    out[idx]=rtn
    done+=1

#This runs countSeq on every part of the file in multiple threads
def readFile(fileLoc):

    global out,rtnd, all
  
    fileName=fileLoc
    rtnd=count(fileName)
    threads=initArrayObj(rtnd[0])
    out=initArrayStr(rtnd[0])
    done=0
    ready=0
    space="     "
    al=["a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z"]
    past8=""
    seq=0
    frames = ["",".","..","..."]
    a=[""]

#THIS RUNS PROCESSING OF FILE
    for x in range(rtnd[0]):
        #print(rtnd[1][x],end="     ")
        threads[x]=threading.Thread(target=countSeq,args=(fileName,rtnd[1][x],rtnd[1][x+1],x))
        threads[x].start()
        time.sleep(0.01)
    for x in threads:
        x.join()
    time.sleep(0.1)

#END PROCESSING



#ALL BELOW PRINTS OUTPUT
#NEVERMIND I REMOVED THE PART THAT PRINTS THE OUTPUT IT NOW JUST RETURNS DATA
    aa=out

    for xxx in aa:
        x=xxx.replace(" ","").lower()
        a[0]=a[0]+x

    b=[]
    for x in a:
        x=x.lower()
        b=([float(math.floor( (len(x)/10)+0.5)/100),x.count("a"),x.count("g"),x.count("c"),x.count("t")])
    ready=1
    return b

#This resets variables used in countSeq, count, and readFile
def resetReadFile():
    global all
    done=0
    ready=0
    space="     "
    past8=""
    seq=0
    a=[""]
    out=[]
    rtnd=None

#This creates an array of random numbers that acts as the seed for randTargets and randInputs so they have the same randomization
def randArray(amount=10):
    amount=min(amount,len(virusNames))
    rtn = [0]*amount
    for x in range(amount):
        ran=random.randint(1,len(virusNames))
        while ran in rtn:
            ran=random.randint(1,len(virusNames))
        rtn[rtn.index(0)]=ran

    rtn.sort()

    return rtn

#This randomizes the targets for training
def randTargets(randOut):
    randTarg=[]

    for x in randOut:
        randTarg.append( [( targets.tolist() )[x-1]])
    randTargTens=torch.tensor(randTarg)
    return randTargTens

#This randomises the outputs for training
def randInputs(randOut):
    randIn=[]

    for x in randOut:
        randIn.append((inputs.tolist())[x-1])

    randInTens=torch.tensor(randIn)
    return randInTens

#This returns the output of randInputs and randTargets with the same seed so that the inputs and targets match up
def randData(a=10):
    a=min(a,len(virusNames))
    pHold=randArray(a)
    return[[randInputs(pHold),randTargets(pHold)],pHold]

#This gets the data from training to be used for the output writing
def getData():
    rtn=[[],[],[],[],[],[],[],[],[]]
    count=0
    tar=targets.tolist()
    for x in tar:
        rtn[0].append(virusNames[len(rtn[0])])
        rtn[1].append(x)
        phold=model(inputs[count]).item()
        rtn[2].append(phold)
        rtn[3].append(x-phold)
        for y in range(5):
            rtn[4+y].append(inputs[len(rtn[4+y])][y])
        count+=1
    return rtn

#This writes to the output file (MORE COMMENTS INSIDE FUNCTION)
def writeOutput(inRTN):
    #This sets up variables
    file_name=default_dir_txt
    global all
    out1=[[],[],[]]
    #This opens the file
    with open(file_name,"a") as f:
        #This starts writing to the file
        time = datetime.datetime.now()
        ctime=time.strftime("%c")
        f.write("_____________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________________\n=============================================================================================================================================================================================================================================================\n")
        f.write(f"Date/Time: {ctime}\nLearning Rate: {learning_rate}\nEpochs: {runs}\nFinal Loss: {round(loss.item(),5)}\n\n")

        header = ["Name","Mutation Rate*","Predicted Rate*","Error","Genome Size (KB)*","A%**","G%**","C%**","T%**"]
        f.write(f"{header[0]:^50} | {header[1]:^25} | {header[2]:^25} | {header[3]:^25} | {header[4]:^25} | {header[5]:^25} | {header[6]:^25} | {header[7]:^25} | {header[8]:^25}")
        f.write("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|~~~~~~~~~~~~~~~~~~~~~~~~~~~|~~~~~~~~~~~~~~~~~~~~~~~~~~~|~~~~~~~~~~~~~~~~~~~~~~~~~~~|~~~~~~~~~~~~~~~~~~~~~~~~~~~|~~~~~~~~~~~~~~~~~~~~~~~~~~~|~~~~~~~~~~~~~~~~~~~~~~~~~~~|~~~~~~~~~~~~~~~~~~~~~~~~~~~|~~~~~~~~~~~~~~~~~~~~~")
        for x in range(len(inRTN[0])):
                  space=" "
                  #This adds the data to an output array for graphing in the graph
                  out1[0].append(inRTN[1][x])
                  out1[1].append(inRTN[2][x])
                  out1[2].append(inRTN[4][x])
                  if inRTN[3][x]<0:
                      space=""
                #This writes the data to the file
                  f.write(f"\n{inRTN[0][x]:^50} | {inRTN[1][x]:<25} | {inRTN[2][x]:<25} | {space+str(inRTN[3][x]):<25} | {inRTN[4][x]:<25} | {inRTN[5][x]:<25} | {inRTN[6][x]:<25} | {inRTN[7][x]:<25} | {inRTN[8][x]:<25}")

        #THIS PART ALLOWS YOU TO RUN THE AI ON OTHER DATA NOT FROM THE TRAINING DATA
        out2=[[],[],[]]
        if len(virusNamesTest) != 0:
            f.write("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|~~~~~~~~~~~~~~~~~~~~~~~~~~~|~~~~~~~~~~~~~~~~~~~~~~~~~~~|~~~~~~~~~~~~~~~~~~~~~~~~~~~|~~~~~~~~~~~~~~~~~~~~~~~~~~~|~~~~~~~~~~~~~~~~~~~~~~~~~~~|~~~~~~~~~~~~~~~~~~~~~~~~~~~|~~~~~~~~~~~~~~~~~~~~~~~~~~~|~~~~~~~~~~~~~~~~~~~~~")
        for X in range(len(virusNamesTest)):
            #This part gets input data from testing data
            nameO = str(virusNamesTest[X])
            #This saves the mutation rate diffrently if a rate is or is not given
            if isinstance(targetsTest[X],float):
                mutRO = math.log(targetsTest[X])
            else: 
                mutRO = targetsTest[X]

            #This part gets input data from testing data, again
            gSizeO = math.log(inputsTest[X][0])

            aCountO = float(inputsTest[X][1])
            gCountO = float(inputsTest[X][2])
            cCountO = float(inputsTest[X][3])
            tCountO = float(inputsTest[X][4])

            #This gets the total of the nuceltides and devides to make them into percentages, than devides the total by two to double the percentages since the model likes that
            total = (aCountO+gCountO+cCountO+tCountO)/2

            #This devides teh nuceltides by the total
            aCountO/=total
            gCountO/=total
            cCountO/=total
            tCountO/=total


            #This scales the inputs and then puts it into the model to get the prediction from the testing data
            reIn = torch.tensor([[gSizeO, aCountO, gCountO, cCountO, tCountO]])

            reIn = reIn.float()

            if not (0 in input_stds.tolist()):
                reIn = (reIn - input_means) / input_stds

            outO = model( reIn )


            #This only forwords data to be graphed if the mutation rate is given
            if isinstance(targetsTest[X],float):

                out2[0].append(mutRO)
                out2[1].append(outO.item())
                out2[2].append(gSizeO)
            
            #This writes testing data to the file
            f.write(f"\n{nameO:^50}")
            f.write(f" | {mutRO:<25}")
            f.write(f" | {outO.item():<25}")
            if isinstance(targetsTest[X],float):
                f.write(f" | {(mutRO-outO.item()):<25}")
            else:
                pholding="No Error Available"
                f.write(f" | {(pholding):<25}")                
            f.write(f" | {gSizeO:<25}")
            f.write(f" | {aCountO:<25}")
            f.write(f" | {gCountO:<25}")
            f.write(f" | {cCountO:<25}")
            f.write(f" | {tCountO:<25}")
            
        #This adds a footnote in the output file
        f.write("\n\n\n*The values listed is the natural log of the value indicated by the row header, for the rows containing an asterisk\n**Yup, these precentages sure are wrong. Nope, it wont be explained.\n\n\n")
        #This sends data to be graphed
        return [out1[0],out1[1],out1[2],out2[0],out2[1],out2[2]]
    
#This dates the output from writeOutput and inputs it into the graphing function
def getResults():
    vIn = writeOutput(getData())
    visuals(vIn[0],vIn[1],vIn[2],vIn[3],vIn[4],vIn[5])


#This graphs the output
def visuals(targ,guess,gsize,targ2,guess2,gsize2):
    matplotlib.use('Agg')
    global default_dir_img

    #This part sets up variables

    size=[]
    for x in range(len(targ)):
        size.append(1)
    for x in range(len(targ2)):
        targ.append(targ2[x])
        guess.append(guess2[x])
        gsize.append(gsize2[x])
        size.append(2)

    minGsize = min(gsize)

    maxGsize = max(gsize)

    lineP = [min(guess),max(guess)]
    linePX= [min(targ),max(targ)]


    lineN=min(lineP[0],linePX[0])
    lineX=max(lineP[1],linePX[1])

    #THERE IS REASON BEHIND "uhhWhat" HOWEVER IM TOO TIERD TO EXPLAIN IT.

    uhhWhat = 50 - 30 * ( 0.0* ( ( 1- ( ( 1/1.01 )**len(gsize) ) ) / ( 1- ( 1/0.01 ) ) ) )
  
    uhhList=[]
    for y in size:
        uhhList.append(uhhWhat*y)

    #This part does the actuall graphing

    plt.plot([lineN,lineX],[lineN,lineX],c="#000000",zorder=-100,lw=2,alpha=1/3)

    plt.scatter(targ,guess,c=gsize,cmap='viridis',vmin=minGsize,vmax=maxGsize,alpha=1,s=uhhList)
   
    plt.gcf().set_layout_engine('constrained')

    for x in range(len(targ)):
        av=(targ[x]/2+guess[x]/2)
        plt.plot([av,targ[x]],[av,guess[x]],lw=1.5,zorder=-10000,alpha=0.25,c="#000000")


    cbar = plt.colorbar(pad=0.2,label="Color scale based on genome size")

    plt.xlabel("Actual mutation rate of virus")

    plt.ylabel("AI predicted mutationrate of virus")

    plt.title("Graphical Representation of Mutation Rate AI")

    cbar.ax.yaxis.set_ticks_position("left")


    plt.savefig(default_dir_img)

#I forget what this part does, but it seems simple enough

def plotOut(Xpoints):
    visuals(getData(Xpoints))

#This sets up the lame old AI model

class TestNN(nn.Module):
    
    def __init__(self):
        super(TestNN, self).__init__()
        self.fc1 = nn.Linear(5,10)
        self.relu = nn.ReLU()
        self.dout = nn.Dropout(p=0)
        self.fc2 = nn.Linear(10,25)
        self.fc3 = nn.Linear(25,50)
        self.fc4 = nn.Linear(50,25)
        self.fc5 = nn.Linear(25,10)
        self.fc6 = nn.Linear(10,1)
         
    def go(self,x):
        x=self.relu(x)
        #x=self.dout(x)
        return x

    def forward(self,x):
        x = self.fc1(x)
        x = self.go(x) 
        x = self.fc2(x)
        x = self.go(x)
        x = self.fc3(x)
        x = self.go(x)      
        x =self.fc4(x)
        x = self.go(x)       
        x =self.fc5(x)
        x =self.fc6(x)
        return x

#This sets up the new awsome AI model

class smolNN(nn.Module):
    def __init__(self):
        super(smolNN, self).__init__()
        self.net = nn.Sequential(
            nn.Linear(5, 16),
            nn.ReLU(),
            nn.Linear(16, 24),
            nn.ReLU(),
            nn.Linear(24, 16),
            nn.ReLU(),
            nn.Linear(16, 5),
            #nn.ReLU(),
            nn.Linear(5, 1) 
        )
        
    def forward(self, x):
        return self.net(x)

def TestandTrainingSplit():
    # Read sequences
    records = list(SeqIO.parse("fasta/proteins.txt", "fasta"))

    if os.path.exists("Testing"):
        shutil.rmtree("Testing")

    # Create directories
    test_dir = Path("Testing")

    test_dir.mkdir(exist_ok=True)
    amino_to_codon = {
        'A': 'GCT', 'B': 'AAC', 'C': 'TGT', 'D': 'GAT', 'E': 'GAA', 'F': 'TTT',
        'G': 'GGT', 'H': 'CAT', 'I': 'ATT', 'K': 'AAA', 'L': 'CTT',
        'M': 'ATG', 'N': 'AAT', 'P': 'CCT', 'Q': 'CAA', 'R': 'CGT',
        'S': 'TCT', 'T': 'ACT', 'V': 'GTT', 'W': 'TGG', 'X': 'NNN', 'Y': 'TAT', 'Z':'CAA'
    }

    def protein_to_dna(protein_seq):
    # Convert sequence to uppercase and map each amino acid to its codon
        return "".join(amino_to_codon.get(aa, 'NNN') for aa in protein_seq.upper())

    #Write individual FASTA files into testing and training

    for record in records:
        name = record.id
        sequence = protein_to_dna(record.seq)
        with open(test_dir / f"{name}.fna", "w") as f:
            f.write(">")
            f.write(name)
            f.write(" ")
            f.write(record.description)
            f.write("\n")
            f.write(sequence)

#THIS PART RUNS THE LOADING SCREEN
ready=0
TestandTrainingSplit()

cuteLoadThreadt = threading.Thread(target=cuteLoadingScreen)
cuteLoadThreadt.start()

#THIS PART BELOW SETS UP THE AI
#IT TURNS OUT MAKING THE AI SMALL HELPS A LOT

model = smolNN()

learning_rate=0.003

criterion = nn.MSELoss()

optimizer = optim.SGD(model.parameters(),lr=learning_rate)

#This sets up file directories loctions

fileRel=os.path.realpath(__file__)[::-1]
fileRel=fileRel[fileRel.index("/"):][::-1]

directory_path = Path(f'{fileRel}/Training')
directory_path_test = Path(f'{fileRel}/Testing')

default_dir_img=Path(f"{fileRel}plot_mut_Size_ln.png")
default_dir_txt = Path(f"{fileRel}/outData_mut_size_DataSet1_Seq.txt")

file_count = sum(1 for item in directory_path_test.iterdir() if item.is_file())

#This sets up empty arrays

fileList=[]
filePath=[]
fileListTest=[]
filePathTest=[]
mutRate=[]
virusNames=[]
virusNamesTest=[]
inputs=[]
targets=[]
rateInNames=[]
rateInValue=[]
targetsTest = []
inputsTest  = []


#This part gets the file name of every training file

for path in directory_path.iterdir():
    if path.is_file():  
        fileList.append( str( path.name ) )
        filePath.append( str(    path   ) )  

#This part gets the file name of every testing file

for path in directory_path_test.iterdir():
    if path.is_file():  
        fileListTest.append( str( path.name ) )
        filePathTest.append( str(    path   ) )          

#This adds the file names to the virus list names after removing the .fna file heading

for x in fileList:
    virusNames.append(x[:x.index(".fna")])
for x in fileListTest:
    virusNamesTest.append(x[:x.index(".fna")])

#This gets the mutation rates for all data

with open (fileRel+"mutationRates.txt","r") as f:

    for lines in f:
        rateInNames.append(lines[:lines.index(":")-1])
        rateInValue.append(float(lines[lines.index(":")+1:].replace("\n","")))


#This part adds the input data and target data for training into their own variables

for x in range(len(virusNames)):
    targets.append(rateInValue[rateInNames.index(virusNames[x])])
    inputs.append (  readFile(filePath[x])) 

    resetReadFile()


#This part adds the input data and target data for testing into their own variables

for x in range(len(virusNamesTest)):
    if virusNamesTest[x] in rateInNames:
        targetsTest.append(rateInValue[rateInNames.index(virusNamesTest[x])])
    else:
        targetsTest.append("No Mutation Rate")
    inputsTest.append (  readFile(filePathTest[0])   ) 

    resetReadFile()

#This ends the loading screen D:

ready = 1
time.sleep(2/5)



#The following converts inputs into numbers usable by the ai, spesificly taking the natural log of the virus size and the mutation rate

dataLen = len(virusNames)

count=0
for x in inputs:
    inputs[count][0]=math.log(x[0])
    count+=1
count=0
for x in targets:
    targets[count]=math.log(x)
    count+=1

#The following covnerts input data for nucleotides total to percentage since otherwise the model explodes
#It also doubles the percentages for the funs

for x in range(len(inputs)):
    total=inputs[x][1]+inputs[x][2]+inputs[x][3]+inputs[x][4]
    for y in range(1,len(inputs[x])):
        inputs[x][y]=float(inputs[x][y])
        inputs[x][y]=2*inputs[x][y]/total

#This converts the targets and inputs into tensors

targets = torch.tensor(  targets )
inputs  =  torch.tensor(inputs)

#This ends the load screen. Again. for some reason. Yeah idk

ready=1
time.sleep(2/5)

#This standardises the input data

backup = inputs

inputs = inputs.float()

input_means = inputs.mean(dim=0)
input_stds = inputs.std(dim=0)
#This if statment accounts for division by zero
if not (0 in input_stds.tolist()):
    inputs = (inputs - input_means) / input_stds

#THIS PART BELOW RUNS TRAINING

#This prompts the user for loss quit values and scycles
default_loops=500
print(f"How many training scycles do you want? (Type zero for default of {default_loops})")
loops = int(input())

stopL = str(input("What loss value would you like to quit at? (Enter blank or -1 if you dont wish to quit)\n"))

if stopL == "":
    stopL=float(-1)
stopL=float(stopL)    

if loops ==0:
    loops = default_loops

#This sets up
go=1
runs = 0
avLossMtr=[0]*5

#This runs the training loop

for epoch in range(loops):
   #This if statment lets it only train if it has not yet hit the loss limit
    if go==1: 
        #This changes the data trained on every other epoch
        if epoch%2==0:
            #This picks 10 random data points from the training to be trained on for the next two epochs
            randot=randData(10)#//.tolist()
            randOt=randot[0]
        randIn=torch.tensor(randOt[0].tolist())
        randTr=torch.tensor(randOt[1].tolist())
        
        #This increases the runs variable. Idk what is does

        runs+=1

        #This is the actuall training part with th eback prop and stuff

        optimizer.zero_grad()

        outputs=model(randIn)
        #//print(outputs)
        #//print()
        #//print(randTr)
        #//print()
        loss = criterion(outputs,randTr)
        #//time.sleep(1000)
        #//print("test")
        loss.backward()
        optimizer.step()

        #This prints the epoch number, loss, and average loss

        lossInt=round(loss.item(),5)
        avLossMtr.pop(0)
        avLossMtr.append(lossInt)
        avLoss=round(sum(avLossMtr)/5,5)


        print(f'Epoch [{epoch+1}/{loops}]          Loss: {lossInt:<15} Average Loss: {avLoss:<15}')
        #//print(randIn)
        
        #This ends the loop if target loss is hit

        if stopL>=avLoss and go==1:
            go=0
            loopsEnd=loops
#This sets the loop variable to the loop you ended on for printing data
if go == 0:
    loops=loopsEnd
                
#THIS PART GIVES THE OUTPUT FROM THE AI

getResults()