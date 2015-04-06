
dims = [5,5,5]

for i in range(-dims[0],dims[0]+1) :
   for j in range(-dims[1],dims[1]+1) :
      for k in range(-dims[2],dims[2]+1) :
         x = pow(abs(float(i)/dims[0]),0.5)/2
         y = float(j)/dims[1]/2
         z = float(k)/dims[2]/2

         if i < 0 :
            x = -x
         #

         print x, y, z
      #
   #
#
