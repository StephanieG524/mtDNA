# Stephanie Gardner
# sggardne
# BME160
# Final Project - graphicRep.py


import matplotlib.pyplot as plt
from sequenceAnalysis import FastAreader
from pandas import ExcelWriter
from pandas import ExcelFile
import pandas as pd
import math
import numpy as np


class GraphicRep:
    '''
    reads and prepare files to be implemented in a matplotlib graphical representation
    of mtDNA. Can eventially be used for different circular DNA
    '''

    myReader = FastAreader('HomoSapiensMitochondrion.fa')
    for head, seq in myReader.readFasta():
        head = head
        seq = seq

    excelFile = 'mtDNACoding.xlsx' #names workbook
    funcLoc = pd.read_excel(excelFile) #loads speadsheet

    excelFile1 = 'MitoMap.xlsx' #names workbook
    mitoMap = pd.read_excel(excelFile1)

    locDict = {}
    mutDict = {}


    for j in range(0,38): #parses through function localtoin spread sheet
        locus = funcLoc.iloc[j,0] #creates dictionary with key - locus and value - list of start, stop and name of locus
        locDict.setdefault(locus,[])
        locDict[locus].append(funcLoc.iloc[j])

    locDictValues = list(locDict.values()) #[locus][start][stop][desc] // [0][0][2] --> stop pos of locus 1

    for j in range(0,330): #parses through mutation spread sheet and saves to dictionary
        locus = mitoMap.iloc[j,1] #dictionary key - locus and value - start, stop, Disease,Allele,RNA,Homoplasmy,Heteroplasmy,pathogenicity
        mutDict.setdefault(locus,[])
        mutDict[locus].append(mitoMap.iloc[j])

    mutDictValues = list(mutDict.values()) #saves vales to a list to be used later on



    def optionA(self,locus):
        '''
        OPTION A:
        input desired locus for viewing/analysis
            base mutations in a certain Locus
                calculates % of mutations at that Locus
            pathogenictiy of certain locus
            saves all info to a text file
        '''
        if locus is None:
            pass
        else:
            with open(locus+'MutationInfo', 'w') as f: #prepare file to save

                try:
                    length = ((len(list(GraphicRep.mutDict[locus])))/330) * 100 #average of mutations at that locus over total mutations
                    f.write("Percent of total mutations found at this locus are: {0:0.2f} % \n".format(length))

                    paths = 0
                    i = 0
                    x=float('nan')

                    for item in GraphicRep.mutDict[locus]: #adds up all the pathogenicities to be averaged out
                        i += 1
                        if math.isnan(item[7]) == False:
                            path = item[7]
                            paths += float(path) #online

                    pathAvg = paths/i #average = sum of all pathogenicity levels over the count of pathogenicities taken

                    if pathAvg == 0.00:
                        f.write('No pathogenicity data is available for this locus')
                    else:
                        f.write("Average Pathogenicity percentage for this locus is: {0:0.2f}%\n".format(pathAvg))
                except KeyError:
                    f.write('No mutation data is available for this locus')


    def optionB(self,low,high,disease):
        '''
        OPTION B:
            input pathogenicity interval [0-100]
                prints loci that have mutations with that pathogenicity
                prints specific base mutations with that pathogenicity
        '''
        if (low or high) is None:
            pass
        else:
            with open('mtDNAPathogInfo', 'w') as f:

                if disease is False:
                    f.write("Sorted based on position\nPathogenicity ({},{})- Codes for - Base Mutation".format(low,high))
                    for i in range(0,25): #for each locus
                        for j in range(0,44): # for all mutations in each locus
                            try:
                                val = float(GraphicRep.mutDictValues[i][j][7])
                                if low <= val <= high: #if value is inbetween low and high interval save to file
                                    f.write("{}% - {} - {}\n".format(GraphicRep.mutDictValues[i][j][7],GraphicRep.mutDictValues[i][j][4],GraphicRep.mutDictValues[i][j][3])) #[locus][entry number][column head]
                            except IndexError:
                                pass
                if disease is True:
                    f.write("Sorted based on position\nPathogenicity ({},{})- Codes for - Base Mutation - Linked Disease\n".format(low,high))
                    for i in range(0,25): #for each locus
                        for j in range(0,44): # for all mutations in each locus
                            try:
                                val = float(GraphicRep.mutDictValues[i][j][7])
                                if low <= val <= high: #if value is inbetween low and high interval save to file
                                    f.write("{}% - {} - {} - {}\n".format(GraphicRep.mutDictValues[i][j][7],GraphicRep.mutDictValues[i][j][4],GraphicRep.mutDictValues[i][j][3],GraphicRep.mutDictValues[i][j][2])) #[locus][entry number][column head]
                            except IndexError:
                                pass
                    f.write('\n\nLHON-Leber Hereditary Optic Neuropathy\nMM-Mitochondrial Myopathy\nAD-Alzeimers Disease\nLIMM-Lethal Infantile Mitochondrial Myopathy\nADPD-Alzeimers Disease and Parkinsonss Disease\nMMC-Maternal Myopathy and Cardiomyopathy\nNARP-Neurogenic muscle weakness, Ataxia, and Retinitis Pigmentosa\nFICP-Fatal Infantile Cardiomyopathy Plus, a MELAS-associated cardiomyopathy\nMELAS-Mitochondrial Encephalomyopathy, Lactic Acidosis, and Stroke-like episodes\nLDYT-Lebers hereditary optic neuropathy and DYsTonia\nMERRF-Myoclonic Epilepsy and Ragged Red Muscle Fibers\nMHCM-Maternally inherited Hypertrophic CardioMyopathy\nCPEO-Chronic Progressive External Ophthalmoplegia\nKSS-Kearns Sayre Syndrome\nDM-Diabetes Mellitus\nDMDF-Diabetes Mellitus + DeaFness\nCIPO-Chronic Intestinal Pseudoobstruction with myopathy and Ophthalmoplegia\nDEAF-Maternally inherited DEAFness or aminoglycoside-induced DEAFness\nPEM-Progressive encephalopathy\nSNHL-SensoriNeural Hearing Loss\n')

    # def optionC(self,locus,plasmy):
    #     '''
    #     OPTION C:
    #         input hetero or homoplasmy
    #             prints base mutations that correspond with...
    #
    #     '''
    #     if plasmy == True:
    #         with open('Hetero/HomoplasmyInfo', 'w') as f:
    #
    #
    #             for i in range(0,25): #for each locus
    #                 for j in range(0,44): # for all mutations in each locus
    #                     try:
    #                         if locus == GraphicRep.mutDictValues[i][j][1]:
    #                             homoplasmy = GraphicRep.mutDictValues[i][j][5]
    #                             heteroplasmy = GraphicRep.mutDictValues[i][j][6]
    #                             if homoplasmy == '+' and homo == True:
    #                                 f.write('Homoplasmic Mutations: \n{} - {}\n').format(GraphicRep.mutDictValues[i][j][4],GraphicRep.mutDictValues[i][j][3]))
    #                             if heteroplasmy == '+' and hetero == True:
    #                                 f.write('Heteroplasmic Mutations: \n{}-{}\n').format(GraphicRep.mutDictValues[i][j][4],GraphicRep.mutDictValues[i][j][3]))
    #                     except IndexError:
    #                         pass
    #     else:
    #         pass

    def get_cmap(n, name='hsv'):
        '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct
        RGB color; the keyword argument name must be a standard mpl colormap name.
        https://stackoverflow.com/questions/14720331/how-to-generate-random-colors-in-matplotlib'''
        return plt.cm.get_cmap(name, n)


    def createMap(self):
        '''
        creates map of mitochondira DNA with color coded mutations and loci
        '''

        circumf = len(GraphicRep.seq) #circumference is equal to the length of the sequence
        baseLen = circumf/len(GraphicRep.seq) #each base is len of 1


        xs = circumf/len(GraphicRep.seq) #baseLen
        baseList = []
        baseList = np.arange(0,circumf,xs)  #populates base list with a list of base position
                                            #and their corresponding position on the circle figure


        cmap = GraphicRep.get_cmap(38) #creates a color map of 38 colors

        fig,ax = plt.subplots() #make figure
        ax = plt.subplot(polar=True) #make polar axis
        ax.set_aspect('equal') #set figure to print out with equal dimensions
        ax.set_xticklabels(['Origin','','','','']) #set origin to the top of the figure
        ax.grid() #clear theta grid from figure
        ax.set_rticks([]) #clear radius ticks from figure
        ax.set_theta_zero_location('N',offset = 0) #set the 0 location to the top of circle

        x=float('nan') #used in conditional (if x is not a number...)



        for i in range(0,25): #locus
            for j in range(0,44): #mutation (max number of mutations any of the loci have)
                try:
                    #makes each pathogenicity interval a special color
                    val = float(GraphicRep.mutDictValues[i][j][7])
                    if (0.0 <= val <= 20.0):
                        colors = 'purple'

                    elif (21.0 <= val <= 40.0):
                        colors = 'blue'

                    elif (41.0 <= val <= 60.0):
                        colors = 'green'

                    elif (61.0 <= val <= 80.0):
                        colors = 'yellow'

                    elif (81.0 <= val <= 99.0):
                        colors = 'orange'

                    elif (val == 100.0):
                        colors = 'red'

                    elif (math.isnan(val) == True):
                        colors = 'black'

                    mutPos = GraphicRep.mutDictValues[i][j][0]
                    theta=((np.pi*2)/len(baseList))*mutPos#calculate theta
                    #populate graph with mutation markers (|) that are roated the correspoding degree
                    ax.plot(theta,0,marker = (2,0,(theta*(180/np.pi))),color=colors,markersize=20,zorder = 6)
                except IndexError:
                    pass
                    #try except needed for locus that have fewer mutations than the max (44)

        for i in range(0,38): #for each locus
            locusName = GraphicRep.locDictValues[i][0][0]
            start = GraphicRep.locDictValues[i][0][1]
            stop = GraphicRep.locDictValues[i][0][2]
            thretaa = ((np.pi*2)/len(baseList))
            #plot gray markers where each locus starts
            ax.plot(thretaa*start,0,marker = (2,0,(thretaa*start)*(180/np.pi)),color='gray',markersize=60,zorder=6)
            for k in range(GraphicRep.locDictValues[i][0][1],GraphicRep.locDictValues[i][0][2]):
                theta = ((np.pi*2)/len(baseList))*k
                #plot all bases in a locus and switch colors at the end of the locus
                ax.plot(theta,0,marker = (2,0,(theta*(180/np.pi))),color=cmap(i),markersize=40,zorder=5)

        ax.annotate('D Loop',
        xy=(0, 0),  # theta, radius
        xytext=(.3, .9),    # fraction, fraction
        textcoords='figure fraction',
        zorder =8,#'offset points'   : Specify an offset (in points) from the xy value
        arrowprops=dict(width=.05,
        headwidth = 2,
        facecolor='black', shrink=.05),
        horizontalalignment='left',
        verticalalignment='bottom')

        y = .90
        for i in range(0,15): #for each locus
            locusNameL = GraphicRep.locDictValues[i][0][0]
            startL = GraphicRep.locDictValues[i][0][1]
            stopL = GraphicRep.locDictValues[i][0][2]
            thetaL = ((np.pi*2)/len(baseList))

            ax.annotate(locusNameL,
            xy=(thetaL*startL, .01),  # theta, radius
            xytext=(.1, y),    # fraction, fraction
            textcoords='figure fraction',
            zorder =8,#'offset points'   : Specify an offset (in points) from the xy value
            arrowprops=dict(width=.05,
            headwidth = 2,
            facecolor='black', shrink=.05),
            horizontalalignment='left',
            verticalalignment='bottom')
            y = y - .06

        x = .17
        for i in range(15,27):
            locusNameL = GraphicRep.locDictValues[i][0][0]
            startL = GraphicRep.locDictValues[i][0][1]
            stopL = GraphicRep.locDictValues[i][0][2]
            thetaL = ((np.pi*2)/len(baseList))

            ax.annotate(locusNameL,
            xy=(thetaL*startL, .01),  # theta, radius
            xytext=(x, .05),    # fraction, fraction
            textcoords='figure fraction',
            zorder =8,#'offset points'   : Specify an offset (in points) from the xy value
            arrowprops=dict(width=.05,
            headwidth = 2,
            facecolor='black', shrink=.05),
            horizontalalignment='left',
            verticalalignment='bottom')
            x = x + .06

        z = .14
        for i in range(27,38): #for each locus
            locusNameL = GraphicRep.locDictValues[i][0][0]
            startL = GraphicRep.locDictValues[i][0][1]
            stopL = GraphicRep.locDictValues[i][0][2]
            thetaL = ((np.pi*2)/len(baseList))

            ax.annotate(locusNameL,
            xy=(thetaL*startL, .01),  # theta, radius
            xytext=(.9, z),    # fraction, fraction
            textcoords='figure fraction',
            zorder =8,#'offset points'   : Specify an offset (in points) from the xy value
            arrowprops=dict(width=.05,
            headwidth = 2,
            facecolor='black', shrink=.05),
            horizontalalignment='left',
            verticalalignment='bottom')
            z = z + .07



        plt.savefig('/Users/stephaniegardner/Desktop/BME160/FinalProject/mtDNAMutations.png')





    def zoomPlot(self,locus):
        '''
        zooms in on a specified locus given in the command line and saves figure to file
        with appropriate locus name.
        First half of this method is the same as createMap.
        '''
        if locus is None:
            pass
        else:

            circumf = len(GraphicRep.seq)
            baseLen = circumf/len(GraphicRep.seq)

            xs = circumf/len(GraphicRep.seq) #baseLen
            baseList = []
            baseList = np.arange(0,circumf,xs)

            cmap = GraphicRep.get_cmap(38)

            figz = plt.figure()
            axz = figz.add_subplot(111, polar=True)
            axz.set_aspect('equal')
            axz.set_xticklabels(['Origin','','','',''])
            axz.grid()
            axz.set_rticks([])
            axz.set_theta_zero_location('N',offset = 0)


            x=float('nan')

            for i in range(0,25): #locus
                for j in range(0,44): #mutation
                    try:
                        val = float(GraphicRep.mutDictValues[i][j][7])
                        if (0.0 <= val <= 20.0):
                            colors = 'purple'
                        elif (21.0 <= val <= 40.0):
                            colors = 'blue'
                        elif (41.0 <= val <= 60.0):
                            colors = 'green'
                        elif (61.0 <= val <= 80.0):
                            colors = 'yellow'
                        elif (81.0 <= val <= 99.0):
                            colors = 'orange'
                        elif (val == 100.0):
                            colors = 'red'
                        elif (math.isnan(val) == True):
                            colors = 'black'
                        mutPos = GraphicRep.mutDictValues[i][j][0]
                        theta=((np.pi*2)/len(baseList))*mutPos
                        axz.plot(theta,0,marker = (2,0,(theta*(180/np.pi))),color=colors,markersize=200,zorder = 6,label=GraphicRep.mutDictValues[i][j][2])
                    except IndexError:
                        pass


            for i in range(0,38): #for each locus
                locusName = GraphicRep.locDictValues[i][0][0]
                start = GraphicRep.locDictValues[i][0][1]
                stop = GraphicRep.locDictValues[i][0][2]
                thretaa = ((np.pi*2)/len(baseList))
                axz.plot(thretaa*start,0,marker = (2,0,(thretaa*start)*(180/np.pi)),color='gray',markersize=150,zorder=6)
                for k in range(GraphicRep.locDictValues[i][0][1],GraphicRep.locDictValues[i][0][2]):
                    theta = ((np.pi*2)/len(baseList))*k
                    axz.plot(theta,0,marker = (2,0,(theta*(180/np.pi))),color=cmap(i),markersize=100,zorder=5)

            '''
            to zoom in on a certain locus, theta min and max must be reset to accomodate
            the locus for viewing.
            '''
            for p in range(0,38):
                if locus == GraphicRep.locDictValues[p][0][0]:
                    thetamin = (360/len(baseList))*GraphicRep.locDictValues[p][0][1]
                    thetamax = (360/len(baseList))*GraphicRep.locDictValues[p][0][2]
                    axz.set_thetamin(thetamin-5)
                    axz.set_thetamax(thetamax+5)

            plt.title(locus) #titles the figure with locus name

            plt.savefig('/Users/stephaniegardner/Desktop/BME160/FinalProject/'+locus+'zoomPlot.png')

class CommandLine() :

    def __init__(self, inOpts=None) :
        '''
        Implement a parser to interpret the command line string using argparse.
        '''
        import argparse
        self.parser = argparse.ArgumentParser(description = 'Program prolog - this program will print out mutation and disease related information for mtDNA',
                                                add_help = True, prefix_chars = '-',)
        self.parser.add_argument('-l', '--locus', type=str, dest = 'locus',
                                 choices= list(GraphicRep.locDict.keys()),
                                 action = 'store',
                                 help='locus')
        self.parser.add_argument('-low', '--lowPathogenicity', type=float, dest = 'low',
                                action = 'store',
                                 help='minimum pathogenicity percentage')
        self.parser.add_argument('-high', '--highPathogenicity', type=float, dest = 'high',
                                    action = 'store',
                                 help='maximum pathogenicity percentage')
        self.parser.add_argument('-disease', '--diseaseInfo', dest = 'disease',
                                    action = 'store', default = False,
                                 help='option for linked disease to be included in pathogenicitiy infor')
        # self.parser.add_argument('-plasmy', '--hetero/homoplasmy', dest = 'plasmy',
        #                             action = 'store', default = False,
        #                          help='option for viewing hetero/homoplasmy mutations at a certain locus')

        if inOpts is None :
            self.args = self.parser.parse_args()
        else :
            self.args = self.parser.parse_args(inOpts)


def main(inCL = None):
    '''
    Command line class that allows arguments given in the terminal to be
    parsed through.
    '''
    myCommandLine = CommandLine(inCL)


    rep = GraphicRep()

    rep.createMap()
    rep.zoomPlot(myCommandLine.args.locus)
    rep.optionA(myCommandLine.args.locus)
    rep.optionB(myCommandLine.args.low,myCommandLine.args.high,myCommandLine.args.disease)
    #rep.optionC(myCommandLine.args.locus,myCommandLine.args.plasmy)



main()




'''

citations:
https://matplotlib.org/gallery/shapes_and_collections/artist_reference.html#sphx-glr-gallery-shapes-and-collections-artist-reference-py
ax.set_aspect('equal') https://stackoverflow.com/questions/34902477/drawing-circles-on-image-with-matplotlib-and-numpy
https://stackoverflow.com/questions/14720331/how-to-generate-random-colors-in-matplotlib

'''
