
dims = [10,10,10]

for i in range(-dims[0],dims[0]+1) :
   for j in range(-dims[1],dims[1]+1) :
      for k in range(-dims[2],dims[2]+1) :
         x = float(i)/dims[0]
         y = float(j)/dims[1]
         z = float(k)/dims[2]
         print x, y, z
      #
   #
#
