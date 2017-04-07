#!/usr/bin/python

topPart = '''{
	"particleFiles" : [ "inputFiles/blue_noise/thinBox-bn.txt"],

	"neighborRadius" : 0.25,
	"alpha" : 1.9,
	"omega" : 0.8,
"mass" : 0.01,
	"numConstraintIters" : 30,
	"maxStretch" : 0.03,
'''
bottomPart = '''
	"springDamping": 0.5,
	"constraintPlanes" : [
		{ 
"x" : -0.4
		}


	],
	"gravity" : [0, -10, 0]
}'''

nClusters = 5317
while nClusters > 0:
    f = open("convergenceTest_%d.txt" % nClusters, "w+")
    f.write(topPart)
    f.write('"nClusters" : %d, \n' % nClusters)
    f.write(bottomPart)
    f.close()
    nClusters /= 2

