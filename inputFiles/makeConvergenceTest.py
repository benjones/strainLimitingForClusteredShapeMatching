#!/usr/bin/python

topPart = '''{
	"particleFiles" : [ "inputFiles/blue_noise/thinBox-bn-1000.txt"],
"dt" : 0.0055555555555556,
	"neighborRadius" : 0.25,
	"alpha" : 1.0,
	"omega" : 0.1,
"mass" : 0.01,
	"numConstraintIters" : 20,
	"maxStretch" : 0.05,
'''
bottomPart = '''
	"springDamping": 0.2,
	"constraintPlanes" : [
		{ 
"x" : -0.4
		}


	],
	"gravity" : [0, -10, 0]
}'''

nClusters = 1000
while nClusters > 0:
    f = open("convergenceTest_%d.txt" % nClusters, "w+")
    f.write(topPart)
    f.write('"nClusters" : %d, \n' % nClusters)
    f.write(bottomPart)
    f.close()
    nClusters /= 2

