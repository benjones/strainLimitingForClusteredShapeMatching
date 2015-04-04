
dims = [5,5,5]

for i in range(-dims[0],dims[0]+1) :
   excess = abs(i)-3
   for j in range(-dims[1]-excess,dims[1]+1+excess) :
      for k in range(-dims[2]-excess,dims[2]+1+excess) :
         x = float(i)/dims[0]/2
         y = float(j)/(dims[1]+excess)/2
         z = float(k)/(dims[2]+excess)/2


         print x, y, z
      #
   #
#
